#!/usr/bin/env bash
# Configure, build and run the self-tests of a selection of demo driver
# directories against the installed library in ./install. Usage:
#
#     .devcontainer/run_selected_tests.sh demo_drivers/poisson/one_d_poisson [more dirs...]
#
# Each directory's build goes into <dir>/build; a summary is printed at the end.
set -uo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
jobs="${JOBS:-$(nproc)}"
summary=()

for dir in "$@"; do
  echo "================================================================"
  echo ">>> ${dir}"
  echo "================================================================"
  if cmake -G Ninja -S "${dir}" -B "${dir}/build" -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_C_COMPILER_LAUNCHER=ccache -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
        > "${dir}/build.configure.log" 2>&1 \
     && cmake --build "${dir}/build" -j "${jobs}" > "${dir}/build.compile.log" 2>&1 \
     && ctest --test-dir "${dir}/build" -j "${jobs}" --output-on-failure; then
    summary+=("PASS  ${dir}")
  else
    summary+=("FAIL  ${dir}  (see ${dir}/build.*.log and ${dir}/build/Validation)")
  fi
done

echo
echo "================================================================"
echo "Summary"
echo "================================================================"
printf '%s\n' "${summary[@]}"
