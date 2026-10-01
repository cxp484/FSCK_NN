# FSCK-NN documentation

`fsck_nn_guide.tex` contains five sections:

1. FSCK database filenames and binary format.
2. Direct PyTorch training with `train_fsck.py` and supporting files.
3. A simple program to read the FSCK database and export CSV.
4. Test the NN model using Fortran and ONNX Runtime.
5. FSCK reference-temperature theory, Planck weighting, and the a-function.

It includes TikZ diagrams and requires a LaTeX installation with TikZ, listings,
xurl, Latin Modern, and standard math/table packages. From this folder:

```bash
./compile.sh
```

Or from the `FSCK_NN` directory:

```bash
bash docs/compile.sh
```

The script runs LaTeX twice and writes `fsck_nn_guide.pdf` in this folder.
Temporary build files are cleaned up automatically. If compilation fails,
the previous PDF is preserved and diagnostics are saved to `compile-error.log`.
It finds `pdflatex` on PATH or at the standard MacTeX location. Set `PDFLATEX`
to an executable path to override it.

No external figures or shell escape are needed.

The compiled PDF is generated locally and excluded from Git.
