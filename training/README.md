# FSCK-NN: direct PyTorch training

Use `train_fsck.py` to train directly from the binary database. It calls the
Fortran kappa backend in memory and passes batches to PyTorch. No CSV export
or intermediate training dataset is required. `fsck_backend.py` supplies the shared binary-database reader and Fortran
compilation helpers.

Requirements: Python 3.9+, gfortran, PyTorch, and ONNX (install from `training/requirements.txt`).

## Train

From the `FSCK_NN` directory:

```bash
# Small pipeline check; these first files are NOT a representative training set.
python3 training/train_fsck.py --pressure 1 --limit-files 6 --epochs 2 --output training/runs/my_smoke

# Train on all files at 1 atm.
python3 training/train_fsck.py --pressure 1 --epochs 20 --output training/runs/one_atm

# Train on the full database at every pressure.
python3 training/train_fsck.py --epochs 20 --output training/runs/full_database
```

Select pressure folders numerically in atm (folder `01.0` is pressure `1`):

```bash
python3 training/train_fsck.py --pressures all
python3 training/train_fsck.py --pressures 1:5
python3 training/train_fsck.py --pressures 0.5,1,3
python3 training/train_fsck.py --pressures 0.5,2:5
```

`all` is the default. Ranges include both endpoints and select existing folders;
`1.5:3` includes folders at 2 and 3 atm if present. Values and ranges can be
combined, and overlapping selections are included only once. The original
`--pressure 1` option remains an alias. Each requested value/range must match
at least one folder; malformed, reversed, or unmatched selections raise errors.
Actual pressure folders used are printed and recorded in `run.json`.
`--limit-files` applies globally **after** pressure filtering, so a small limit
may cover only the lowest selected pressure.

Each output directory must be new or empty. A full run can be substantial:
every selected database record is visited once per epoch, including validation.
`--limit-files` selects the first N sorted files and is intended for smoke tests.

## Data and model

Inputs, in order: `P_atm, xCO2, xH2O, xCO, fv, Tgas_K, Tref_K`.
Targets: the 32 kappa values returned by the Fortran `get_k` / `k_interp` path.
Each file provides all 784 temperature pairs. Reference temperature remains an
input because it changes the Planck-weighted k distribution.

The default network is a fully connected MLP with three 128-unit hidden layers
and SiLU activation. Pressure uses logarithmic scaling; soot uses log1p scaling;
temperatures are scaled by the known 300–3000 K range. Mole fractions remain
unchanged. All preprocessing is inside the model; inference takes physical
input values in the units above.

The loss is MSE on `log10(max(kappa, 1e-8))`, in the database's native kappa
units. The floor handles zero values and is configurable with `--kappa-floor`.
Thus values below the floor are not resolved separately by the training loss.
`model.forward()` returns log10 kappa; `model.predict_kappa()` returns positive
kappa values in native database units. There are no a-function outputs. No
monotonicity constraint across the 32 outputs is imposed.

`--batch-size` (default 1024) controls optimizer batches. `--chunk-files`
(default 8) bounds the group of files loaded for in-memory shuffling. Files and
rows are shuffled each training epoch. Data is not loaded all at once. The
Fortran module holds mutable state, so calls are sequential in one process;
this program does not call it from concurrent DataLoader workers.

Validation holds out entire pressure/composition files, keeping all temperature
pairs of a file together. The split is reproducible using `--seed` and recorded
in `run.json`. Similar compositions and the same composition at another
pressure can occur across splits: this measures held-out-file performance, not
strict extrapolation to unseen composition families. No validation statistics
are used to fit input normalization.

Options include `--width`, `--depth`, `--learning-rate`, `--validation-fraction`,
`--threads`, `--device cpu|cuda|mps|auto`, `--compiler`, and `--debug`.
`auto` selects CUDA if available, otherwise CPU.

## Outputs and inference

- `best.pt`: lowest validation log10 MSE checkpoint.
- `last.pt`: final completed epoch checkpoint.
- `best.onnx` and `last.onnx`: inference exports matching the respective checkpoints.
- `metrics.jsonl`: training/validation log10 MSE, log10 MAE, and row counts.
- `run.json`: file split, model settings, input/output names, and Fortran hash.

Checkpoints include model configuration, preprocessing buffers, optimizer
state, epoch, and fixed quadrature metadata. The w values in metadata are not
trained targets. Automatic training resumption is not implemented.

ONNX exports are written automatically with each corresponding checkpoint.
They use opset 17, accept float32 `physical_inputs` of shape `(batch, 7)` in
the input order above, and return `kappa` of shape `(batch, 32)` in native
database units (not log10). The batch dimension is dynamic; preprocessing and
conversion from logarithms are included. Optimizer state and quadrature
metadata remain in the `.pt` files.

From the `FSCK_NN` directory, start Python and load a checkpoint:

```python
import sys
import torch

sys.path.insert(0, "training")
from train_fsck import KappaModel

checkpoint = torch.load("training/runs/my_smoke/best.pt", map_location="cpu", weights_only=True)
model = KappaModel(**checkpoint["model_config"])
model.load_state_dict(checkpoint["model_state"])
model.eval()
x = torch.tensor([[1.0, 0.1, 0.2, 0.05, 1e-6, 1000.0, 1500.0]])
with torch.no_grad():
    kappa = model.predict_kappa(x)  # shape (1, 32)
```

Checkpoints are generated locally and are not bundled. A two-epoch run on six
files verifies integration; it is not an accurate production surrogate.
Before deployment, validate errors across the operating range and radiative
quantities of interest, not just the aggregate logarithmic loss.

## Fortran backend

The backend reads 784 temperature-pair records per database file, each with
33 float32 values. It uses the first 32 kappas and omits the final maximum
kappa, matching the original V4 loader. Files are processed using singleton
pressure/composition axes and all 28 local/reference temperatures.

`fortran/row_bridge.inc` is inserted into a generated module copy under `training/build/` (relative to `FSCK_NN`).
The preserved numerical sources are unchanged. Each record calls `get_k` and
`k_interp`; no a-function is calculated. Quadrature uses
`quadgen2(.false.,g,w,32,2.d0)`, matching the V4Table runtime.

Original SRCS copyright/license notices are retained in the Fortran sources.
