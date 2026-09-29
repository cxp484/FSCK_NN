# FSCK-NN

Train a neural-network surrogate for the 32-point FSCK kappa database and
compare its ONNX predictions against direct database interpolation in Fortran.

- [`training/`](training/README.md): PyTorch training and ONNX export.
- [`test/`](test/README.md): Fortran ONNX inference, database comparison, and plots.
- [`fortran/`](fortran/): preserved numerical sources and database bridge.
- [`docs/fsck_nn_guide.tex`](docs/fsck_nn_guide.tex): theory and usage guide.

Run commands from this `FSCK_NN` directory:

```bash
python3 -m pip install -r training/requirements.txt
python3 training/train_fsck.py --database /path/to/FSKTableV4 \
  --pressure 1 --epochs 20 --output training/runs/one_atm
```

The radiation database, trained checkpoints, comparison results, compiled
binaries, and ONNX Runtime SDK are not included. Build the documentation with
`bash docs/compile.sh`. See the linked instructions for the Fortran compiler,
SDK setup and inference.

## Redistribution of bundled sources

The files in `fortran/` retain their original SRCS copyright/license notices.
Those notices specify non-commercial use and state that the files must not be
distributed without permission from the copyright holder. Obtain that permission
before uploading these bundled sources to GitHub or otherwise redistributing
them. This repository does not grant a new license to the third-party code.
