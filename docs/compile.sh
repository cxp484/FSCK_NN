#!/usr/bin/env bash
# Compile the guide from any working directory.
set -euo pipefail

docs_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
latex_bin="${PDFLATEX:-pdflatex}"

if ! command -v "$latex_bin" >/dev/null 2>&1; then
    if [[ "$latex_bin" == pdflatex && -x /Library/TeX/texbin/pdflatex ]]; then
        latex_bin=/Library/TeX/texbin/pdflatex
    else
        printf 'Error: pdflatex not found. Install TeX Live/MacTeX or set PDFLATEX to its executable path.\n' >&2
        exit 1
    fi
fi

build_dir="$(mktemp -d "${TMPDIR:-/tmp}/fsck-latex.XXXXXX")"
trap 'rm -rf -- "$build_dir"' EXIT

cd -- "$docs_dir"
for pass in 1 2; do
    printf 'Compiling fsck_nn_guide.tex (pass %s/2)...\n' "$pass"
    if ! "$latex_bin" -interaction=nonstopmode -halt-on-error -file-line-error \
        -output-directory="$build_dir" fsck_nn_guide.tex \
        > "$build_dir/compiler-output.txt" 2>&1; then
        cat "$build_dir/compiler-output.txt" >&2
        cp -- "$build_dir/compiler-output.txt" "$docs_dir/compile-error.log"
        printf 'Compilation failed. Details: %s/compile-error.log\n' "$docs_dir" >&2
        exit 1
    fi
done

cp -- "$build_dir/fsck_nn_guide.pdf" "$docs_dir/fsck_nn_guide.pdf"
printf 'Created %s/fsck_nn_guide.pdf\n' "$docs_dir"
