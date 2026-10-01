#!/usr/bin/env bash
# Build the standalone CSV exporter from any working directory.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$script_dir/build"
"${FC:-gfortran}" -std=f2008 -O2 "$script_dir/FSCK_to_csv.f90" \
    -o "$script_dir/build/FSCK_to_csv"
printf 'Built %s/build/FSCK_to_csv\n' "$script_dir"
