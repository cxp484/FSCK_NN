#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ort_dir="${ONNXRUNTIME_ROOT:-}"
if [[ -z "$ort_dir" || ! -f "$ort_dir/include/onnxruntime_c_api.h" ]]; then
    echo 'Set ONNXRUNTIME_ROOT to an extracted ONNX Runtime C SDK (include/ and lib/).' >&2
    exit 1
fi
ort_dir="$(cd -- "$ort_dir" && pwd)"
mkdir -p "$script_dir/build"
cd "$script_dir/build"
"${CC:-cc}" -std=c11 -O2 -fPIC -I"$ort_dir/include" -c "$script_dir/onnx_bridge.c" -o onnx_bridge.o
sources=(commonData.f90 quadlib.f90 commonRoutines.f90 FSKDBPath.f90 libJCai.f90 fskTable.f90 fskTableV4.f90)
for i in "${!sources[@]}"; do sources[$i]="$script_dir/../fortran/${sources[$i]}"; done
flags=(-O2 -g -fcheck=all -fbacktrace -ffree-line-length-none)
"${FC:-gfortran}" "${flags[@]}" -c "${sources[@]}" "$script_dir/nn_model.f90" "$script_dir/comparison_grid.f90"
objects=(commonData.o quadlib.o commonRoutines.o FSKDBPath.o libJCai.o fskTable.o fskTableV4.o nn_model.o comparison_grid.o onnx_bridge.o)
"${FC:-gfortran}" "${flags[@]}" "$script_dir/compare_fsck.f90" "${objects[@]}" \
    -L"$ort_dir/lib" -Wl,-rpath,"$ort_dir/lib" -lonnxruntime -o compare_fsck
printf 'Built comparison executable in %s/build\n' "$script_dir"
