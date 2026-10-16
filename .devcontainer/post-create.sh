#!/usr/bin/env bash
# Runs once, inside the container, after it has been created (see
# "postCreateCommand" in devcontainer.json). Keep this light: it should set up
# the workspace, not build the whole library.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."

echo ">>> Toolchain"
for tool in cmake ninja g++ gfortran ccache python3; do
  printf '    %-10s %s\n' "$tool" "$(command -v "$tool" >/dev/null && "$tool" --version 2>/dev/null | head -n1 || echo 'MISSING')"
done
if command -v mpicxx >/dev/null; then
  printf '    %-10s %s\n' "mpicxx" "$(mpicxx --version | head -n1)"
fi

echo ">>> Prebuilt third-party libraries in ${OOMPH_TPL_INSTALL_DIR:-<unset>}"
if [ -f "${OOMPH_TPL_INSTALL_DIR:-}/cmake_flags_for_oomph_lib.json" ]; then
  python3 - <<'EOF'
import json, os
with open(os.path.join(os.environ["OOMPH_TPL_INSTALL_DIR"], "cmake_flags_for_oomph_lib.json")) as f:
    for k, v in sorted(json.load(f).items()):
        print(f"    {k} = {v}")
EOF
else
  echo "    (none found; run ./oomph_build.py without --no-build-tpl to build them in-tree)"
fi

# Make sure git doesn't complain about the bind-mounted checkout being owned by
# a different user than the one running inside the container.
git config --global --add safe.directory "$(pwd)" || true

# Keep the pre-commit hooks available if the user wants them
# (see .pre-commit-config.yaml).
if ! command -v pre-commit >/dev/null; then
  python3 -m pip install --quiet --user --break-system-packages pre-commit 2>/dev/null || true
fi

cat <<'EOF'

>>> Ready. To build oomph-lib against the prebuilt third-party libraries:

        oomph-build --oomph-OOMPH_INSTALL_HEADERS_AS_SYMLINKS=ON

    (this wraps ./oomph_build.py --no-build-tpl --enable-ccache
     --ext-CMAKE_INSTALL_PREFIX="$OOMPH_TPL_INSTALL_DIR"; any extra flags are
     forwarded, e.g. "-c Debug --oomph-OOMPH_ENABLE_PARANOID=ON --reuse-oomph").

    To build and run a demo driver's self-tests afterwards:

        cmake -G Ninja -S demo_drivers/poisson/one_d_poisson -B demo_drivers/poisson/one_d_poisson/build
        cmake --build demo_drivers/poisson/one_d_poisson/build
        ctest --test-dir demo_drivers/poisson/one_d_poisson/build

EOF
