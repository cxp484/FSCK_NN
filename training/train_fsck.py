#!/usr/bin/env python3
"""Train a 7-input, 32-kappa PyTorch model directly from the binary FSCK database."""
import argparse
import copy
import itertools
import json
import math
from pathlib import Path
import random
import time

import onnx
import torch
from torch import nn

from fsck_backend import Backend, DEFAULT_DB, ROOT, build, files, source_digest, validate_raw

from pressure_selection import parse_pressure_selection

INPUTS = ['P_atm', 'xCO2', 'xH2O', 'xCO', 'fv', 'Tgas_K', 'Tref_K']

class KappaModel(nn.Module):
    """forward returns log10 kappa; predict_kappa returns native database units."""
    def __init__(self, width=128, depth=3, kappa_floor=1e-8):
        super().__init__()
        self.config = dict(width=width, depth=depth, kappa_floor=kappa_floor)
        self.register_buffer('kappa_floor', torch.tensor(kappa_floor))
        layers = [nn.Linear(7,width), nn.SiLU()]
        for _ in range(depth-1):
            layers.extend([nn.Linear(width,width), nn.SiLU()])
        layers.append(nn.Linear(width,32))
        self.network = nn.Sequential(*layers)

    def forward(self, inputs):
        x = inputs.clone()
        # Fixed physical scaling; no validation-data statistics are fitted.
        x[...,0] = (torch.log10(x[...,0])+1.0)/math.log10(800.0)
        x[...,4] = torch.log1p(x[...,4]/1e-9)/math.log1p(1e-5/1e-9)
        x[...,5:7] = (x[...,5:7]-300.0)/2700.0
        return self.network(x)

    def targets(self, kappas):
        return torch.log10(kappas.clamp_min(self.kappa_floor))

    def predict_kappa(self, inputs):
        return torch.pow(10.0,self(inputs))


def read_file(backend, entry):
    path,state = entry
    raw = validate_raw(path)
    # tensor copies the reusable ctypes buffer before the next Fortran call.
    k = torch.tensor(list(backend.evaluate(path,state)),dtype=torch.float64).reshape(784,32)
    reference = torch.tensor(raw,dtype=torch.float64).reshape(784,33)[:,:32]
    if not torch.allclose(k,reference,rtol=2e-10,atol=1e-12):
        raise ValueError(f'Fortran/database mismatch: {path}')
    if (k<0).any():
        raise ValueError(f'Negative kappa: {path}')
    x = torch.empty((784,7),dtype=torch.float32)
    x[:,:5] = torch.tensor(state,dtype=torch.float32)
    temperatures = torch.arange(300,3001,100,dtype=torch.float32)
    x[:,5] = temperatures.repeat_interleave(28)
    x[:,6] = temperatures.repeat(28)
    return x,k.float()


def split_files(entries, fraction, seed):
    if len(entries)<2:
        raise ValueError('Select at least two database files for train/validation splitting')
    ordered=list(entries)
    random.Random(seed).shuffle(ordered)
    count=max(1,min(len(ordered)-1,round(len(ordered)*fraction)))
    return ordered[count:],ordered[:count]


def batches(backend, entries, batch_size, chunk_files, seed, shuffle):
    ordered=list(entries)
    rng=random.Random(seed)
    generator=torch.Generator().manual_seed(seed)
    if shuffle: rng.shuffle(ordered)
    for start in range(0,len(ordered),chunk_files):
        tensors=[read_file(backend,e) for e in ordered[start:start+chunk_files]]
        x=torch.cat([t[0] for t in tensors])
        k=torch.cat([t[1] for t in tensors])
        if shuffle:
            order=torch.randperm(len(x),generator=generator)
            x,k=x[order],k[order]
        for offset in range(0,len(x),batch_size):
            yield x[offset:offset+batch_size],k[offset:offset+batch_size]


def epoch(model, backend, entries, args, epoch_number, optimizer=None):
    training=optimizer is not None
    model.train(training)
    squared=absolute=0.0
    elements=0
    last_report=time.monotonic()
    with torch.set_grad_enabled(training):
        for x,k in batches(backend,entries,args.batch_size,args.chunk_files,
                           args.seed+epoch_number,training):
            x,k=x.to(args.device),k.to(args.device)
            target=model.targets(k)
            prediction=model(x)
            loss=(prediction-target).square().mean()
            if not torch.isfinite(loss):
                raise RuntimeError('Nonfinite training/validation loss')
            if training:
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                nn.utils.clip_grad_norm_(model.parameters(),10.0)
                optimizer.step()
            squared+=loss.item()*k.numel()
            absolute+=(prediction-target).abs().sum().item()
            elements+=k.numel()
            if time.monotonic()-last_report>30:
                print(f'  {"train" if training else "validation"}: {elements//32} rows processed',flush=True)
                last_report=time.monotonic()
    return {'log10_mse':squared/elements,'log10_mae':absolute/elements,'rows':elements//32}


def save_checkpoint(path,payload):
    temporary=path.with_suffix('.partial')
    torch.save(payload,temporary)
    temporary.replace(path)


class KappaInference(nn.Module):
    """ONNX interface: physical inputs to kappa in native database units."""
    def __init__(self, model):
        super().__init__()
        self.model = model

    def forward(self, inputs):
        return self.model.predict_kappa(inputs)


def export_onnx(path, model):
    # Export a CPU copy without changing the live model's device or mode.
    inference = KappaInference(copy.deepcopy(model).cpu().eval())
    example = torch.tensor([[1.0, .1, .2, .05, 1e-6, 1000.0, 1500.0]])
    temporary = path.with_suffix('.onnx.partial')
    try:
        torch.onnx.export(
            inference, example, str(temporary), opset_version=17, dynamo=False,
            input_names=['physical_inputs'], output_names=['kappa'],
            dynamic_axes={'physical_inputs': {0: 'batch'}, 'kappa': {0: 'batch'}})
        graph = onnx.load(str(temporary))
        onnx.helper.set_model_props(graph, {
            'inputs': json.dumps(INPUTS),
            'outputs': '32 kappa values in native database units; not log10',
            'preprocessing': 'Included in graph',
        })
        onnx.checker.check_model(graph)
        onnx.save(graph, str(temporary))
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def train(args):
    torch.manual_seed(args.seed)
    torch.set_num_threads(args.threads)
    database=args.database.resolve()
    entries=list(itertools.islice(files(database,args.pressure),args.limit_files)) if args.limit_files else list(files(database,args.pressure))
    pressure_folders=sorted({e[0].parent.name for e in entries},key=float)
    print(f'Pressure folders used: {", ".join(pressure_folders)}',flush=True)
    training,validation=split_files(entries,args.validation_fraction,args.seed)
    if args.device=='auto':
        args.device='cuda' if torch.cuda.is_available() else 'cpu'
    output=args.output.resolve()
    if output==database or database in output.parents:
        raise ValueError('Training output must be outside the database')
    if output.exists() and any(output.iterdir()):
        raise ValueError('Output directory must be empty; choose a new directory for this run')
    backend=Backend(build(args.compiler,args.debug))
    output.mkdir(parents=True,exist_ok=True)
    model=KappaModel(args.width,args.depth,args.kappa_floor).to(args.device)
    optimizer=torch.optim.Adam(model.parameters(),lr=args.learning_rate)
    manifest={'inputs':INPUTS,'outputs':[f'kappa_{i:02}' for i in range(1,33)],
              'target_transform':'log10(max(kappa,kappa_floor))','model':model.config,
              'database':str(database),'fortran_sha256':source_digest(),
              'pressure_folders':pressure_folders,
              'train_files':[str(e[0].relative_to(database)) for e in training],
              'validation_files':[str(e[0].relative_to(database)) for e in validation],
              'settings':{key:str(value) if isinstance(value,Path) else value for key,value in vars(args).items()}}
    (output/'run.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'Training: {len(training)} files / {len(training)*784} rows; validation: {len(validation)} files / {len(validation)*784} rows; device: {args.device}',flush=True)
    best=float('inf')
    for n in range(1,args.epochs+1):
        tr=epoch(model,backend,training,args,n,optimizer)
        val=epoch(model,backend,validation,args,n)
        metrics={'epoch':n,'train':tr,'validation':val}
        with (output/'metrics.jsonl').open('a') as stream:
            stream.write(json.dumps(metrics)+'\n')
        payload={'model_state':model.state_dict(),'model_config':model.config,
                 'optimizer_state':optimizer.state_dict(),'epoch':n,'metrics':metrics,
                 'inputs':INPUTS,'g':list(backend.g),'w':list(backend.w),
                 'fortran_sha256':manifest['fortran_sha256']}
        save_checkpoint(output/'last.pt',payload)
        export_onnx(output/'last.onnx',model)
        if val['log10_mse']<best:
            best=val['log10_mse']
            save_checkpoint(output/'best.pt',payload)
            export_onnx(output/'best.onnx',model)
        print(f'Epoch {n}: train log10 MSE={tr["log10_mse"]:.6g}, validation log10 MSE={val["log10_mse"]:.6g}',flush=True)
    print(f'Saved model: {output/"best.pt"}',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--database',type=Path,default=DEFAULT_DB)
    p.add_argument('--output',type=Path,default=ROOT/'runs'/time.strftime('%Y%m%d_%H%M%S'))
    p.add_argument('--pressures','--pressure',dest='pressure',default='all',
                   help='Pressure folders in atm: all (default), 1, 1:5 (inclusive), or 0.5,1,3; ranges and values may be combined')
    p.add_argument('--limit-files',type=int,help='First N files in sorted order; useful for smoke tests only')
    p.add_argument('--epochs',type=int,default=20)
    p.add_argument('--batch-size',type=int,default=1024)
    p.add_argument('--chunk-files',type=int,default=8,help='Maximum files loaded together for row shuffling')
    p.add_argument('--width',type=int,default=128)
    p.add_argument('--depth',type=int,default=3)
    p.add_argument('--learning-rate',type=float,default=1e-3)
    p.add_argument('--kappa-floor',type=float,default=1e-8,help='Log-target floor in native database units')
    p.add_argument('--validation-fraction',type=float,default=0.1)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--device',default='auto',choices=['auto','cpu','cuda','mps'])
    p.add_argument('--compiler',default='gfortran')
    p.add_argument('--debug',action='store_true')
    args=p.parse_args()
    try:
        parse_pressure_selection(args.pressure)
    except ValueError as error:
        p.error(str(error))
    for key in ['epochs','batch_size','chunk_files','width','depth','threads']:
        if getattr(args,key)<1: p.error(f'{key} must be positive')
    if args.limit_files is not None and args.limit_files<2: p.error('--limit-files must be at least 2')
    if not 0<args.validation_fraction<1: p.error('--validation-fraction must be between 0 and 1')
    for key in ['learning_rate','kappa_floor']:
        if not math.isfinite(getattr(args,key)) or getattr(args,key)<=0: p.error(f'{key} must be finite and positive')
    train(args)

if __name__=='__main__': main()
