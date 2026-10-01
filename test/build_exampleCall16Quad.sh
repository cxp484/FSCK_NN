#!/usr/bin/env bash
# Build the direct-database 16-point runtime example; ONNX Runtime is not needed.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$script_dir/build/example16"
cd "$script_dir/build/example16"
sources=(commonData.f90 quadlib.f90 commonRoutines.f90 FSKDBPath.f90 libJCai.f90 fskTable.f90 fskTableV4.f90)
for i in "${!sources[@]}"; do sources[$i]="$script_dir/../fortran/${sources[$i]}"; done
"${FC:-gfortran}" -O2 -g -fcheck=all -fbacktrace -ffree-line-length-none \
    "${sources[@]}" "$script_dir/exampleCall16Quad.f90" -o ../exampleCall16Quad
printf 'Built %s/build/exampleCall16Quad\n' "$script_dir"
