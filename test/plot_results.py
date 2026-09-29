#!/usr/bin/env python3
"""Plot existing Fortran comparison CSV files. No model loading or inference."""
import argparse
import json
import os
from pathlib import Path
import numpy as np

def metrics(reference, prediction):
    finite = np.isfinite(reference) & np.isfinite(prediction)
    r, p = reference[finite], prediction[finite]
    result = {'points': int(len(reference)), 'finite_points': int(finite.sum()),
              'nonfinite_points': int((~finite).sum())}
    if len(r):
        result.update(mae=float(np.mean(abs(p-r))), rmse=float(np.sqrt(np.mean((p-r)**2))))
        positive = (r>0) & (p>0)
        if positive.any():
            result['log10_rmse_positive'] = float(np.sqrt(np.mean((np.log10(p[positive])-np.log10(r[positive]))**2)))
    return result


def plot_results(data, output):
    os.environ.setdefault('MPLCONFIGDIR', str(output/'matplotlib-cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    summary = {}
    valid = data['nn_monotone'] == 1
    for name, label in [('k','Kappa (native database units)'), ('a','Stretching function a')]:
        reference, prediction = data[f'{name}_db'], data[f'{name}_nn']
        summary[name] = metrics(reference, prediction)
        fig, axes = plt.subplots(1,2,figsize=(11,5),layout='constrained')
        for ax, logarithmic in zip(axes, [False, True]):
            mask = np.isfinite(reference) & np.isfinite(prediction)
            if logarithmic:
                mask &= (reference>0) & (prediction>0)
            if name == 'a':
                for group, color, legend in [(valid,'tab:blue','Monotone NN curves'), (~valid,'tab:orange','Nonmonotone NN: diagnostic only')]:
                    pick=mask & group
                    if pick.any(): ax.scatter(reference[pick],prediction[pick],s=7,alpha=.35,color=color,label=legend,rasterized=True)
                if mask.any(): ax.legend(fontsize=8)
            else:
                ax.scatter(reference[mask],prediction[mask],s=7,alpha=.35,rasterized=True)
            if mask.any():
                lo=min(reference[mask].min(),prediction[mask].min())
                hi=max(reference[mask].max(),prediction[mask].max())
                if hi == lo:
                    lo,hi=(lo/2,hi*2) if logarithmic else (lo-1,hi+1)
                if logarithmic:
                    factor=(hi/lo)**.025
                    lo,hi=lo/factor,hi*factor
                else:
                    margin=(hi-lo)*.025
                    lo,hi=lo-margin,hi+margin
                ax.plot([lo,hi],[lo,hi],'k--',linewidth=1,label='1:1')
                if logarithmic: ax.set_xscale('log'); ax.set_yscale('log')
                ax.set_xlim(lo,hi); ax.set_ylim(lo,hi)
            else:
                ax.text(.5,.5,'No eligible points',ha='center',transform=ax.transAxes)
            ax.set_aspect('equal',adjustable='box')
            ax.set_xlabel('Direct database'); ax.set_ylabel('NN')
            ax.set_title(('Positive values, log scale' if logarithmic else 'Linear scale')+f' ({mask.sum()} points)')
            ax.grid(True,alpha=.2)
        fig.suptitle(label+' — 1:1 comparison')
        fig.savefig(output/f'{name}_parity.png',dpi=180)
        fig.savefig(output/f'{name}_parity.pdf')
        plt.close(fig)
    summary['a_monotone_nn_only'] = metrics(data['a_db'][valid],data['a_nn'][valid])
    per_state=data[::32]
    summary['states']=len(per_state)
    summary['nonmonotone_nn_states']=int((per_state['nn_monotone']==0).sum())
    summary['negative_a_nn_points']=int((data['a_nn']<0).sum())
    # Report actual discrete normalization, without renormalizing original afun output.
    for source in ('db','nn'):
        totals=np.sum((data['w']*data[f'a_{source}']).reshape(-1,32),axis=1)
        finite=totals[np.isfinite(totals)]
        summary[f'a_{source}_weighted_sum_range']= [float(finite.min()),float(finite.max())] if len(finite) else None
    return summary



def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('csv',type=Path)
    parser.add_argument('--output',type=Path,help='Plot directory; defaults to CSV directory')
    args=parser.parse_args()
    output=(args.output or args.csv.parent).resolve()
    output.mkdir(parents=True,exist_ok=True)
    data=np.atleast_1d(np.genfromtxt(args.csv,delimiter=',',names=True))
    required={'state_id','q','g','w','k_db','k_nn','a_db','a_nn','nn_monotone'}
    if not required.issubset(data.dtype.names or ()) or len(data)==0 or len(data)%32:
        raise ValueError('Expected nonempty comparison CSV with 32 rows per state')
    if not np.array_equal(data['q'],np.tile(np.arange(1,33),len(data)//32)):
        raise ValueError('CSV quadrature rows are out of order')
    summary=plot_results(data,output)
    summary['comparison_csv']=str(args.csv.resolve())
    (output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary,indent=2))
    print(f'Plots: {output}')

if __name__=='__main__': main()
