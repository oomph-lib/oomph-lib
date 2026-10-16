# Local container setup for oomph-lib

This directory (together with the top-level [`Dockerfile`](../Dockerfile) and
[`.dockerignore`](../.dockerignore)) provides a reproducible Linux build
environment for oomph-lib. It is intended for local development; the GitHub
Actions workflows in `.github/workflows/` are unaffected.

## What's in the image

| Stage     | Contents                                                                                          |
| --------- | ------------------------------------------------------------------------------------------------- |
| `base`    | Ubuntu 24.04, gcc/g++/gfortran, CMake (>= 3.24), Ninja, ccache, Python 3, gdb, valgrind, clang-format, optional OpenMPI, non-root user `oomph` |
| `tpl`     | Third-party libraries from `external_distributions/` built and installed into `/opt/oomph-lib/tpl`. OpenBLAS + SuperLU by default; CGAL / MUMPS / Hypre / Trilinos via build args |
| `dev`     | **Default target.** `base` + the prebuilt libraries + an `oomph-build` helper. Expects the repository to be bind-mounted at `/workspaces/oomph-lib` |
| `library` | `dev` + a complete build/install of oomph-lib from the current checkout                           |

## VS Code / Dev Containers

Open the repository in VS Code and choose *Dev Containers: Reopen in
Container*. The image is built from the `dev` target with the build args in
[`devcontainer.json`](./devcontainer.json); `post-create.sh` prints the
toolchain and the prebuilt third-party libraries that were found.

Inside the container:

```bash
# Build and install oomph-lib (into ./install) against the prebuilt TPLs.
# Any extra arguments are forwarded to ./oomph_build.py.
oomph-build --oomph-OOMPH_INSTALL_HEADERS_AS_SYMLINKS=ON

# Debug build with PARANOID/RANGE_CHECKING, reusing existing build dirs
oomph-build -c Debug --oomph-OOMPH_ENABLE_PARANOID=ON --reuse-oomph

# Build + run a demo driver's self tests
cmake -G Ninja -S demo_drivers/poisson/one_d_poisson -B demo_drivers/poisson/one_d_poisson/build
cmake --build demo_drivers/poisson/one_d_poisson/build
ctest --test-dir demo_drivers/poisson/one_d_poisson/build

# ...or several at once (configure + build + ctest each, with a summary)
.devcontainer/run_selected_tests.sh \
    demo_drivers/navier_stokes/driven_cavity \
    demo_drivers/navier_stokes/rayleigh_channel
```

`oomph-build` is just

```bash
python3 oomph_build.py --no-build-tpl --enable-ccache \
    --OOMPH_ENABLE_MPI="$OOMPH_ENABLE_MPI" \
    --ext-CMAKE_INSTALL_PREFIX="$OOMPH_TPL_INSTALL_DIR" "$@"
```

so you can always call `./oomph_build.py` directly if you prefer (e.g. to
build a different set of third-party libraries in-tree).

The ccache directory is kept in a named Docker volume (`oomph-lib-ccache`) so
compilation results survive container rebuilds.

## Plain Docker

```bash
# Development image (default target)
docker build -t oomph-lib:dev .
docker run --rm -it -v "$PWD":/workspaces/oomph-lib oomph-lib:dev
oomph-build

# ...with MPI and MUMPS/Hypre prebuilt
docker build -t oomph-lib:dev-mpi \
    --build-arg ENABLE_MPI=ON --build-arg BUILD_MUMPS=ON --build-arg BUILD_HYPRE=ON .

# Fully built library from the current checkout
docker build --target library -t oomph-lib:library .
```

### Build arguments

| Arg                              | Default        | Purpose                                             |
| -------------------------------- | -------------- | --------------------------------------------------- |
| `BASE_IMAGE`                     | `ubuntu:24.04` | Base image (needs CMake >= 3.24 via apt, or adjust) |
| `ENABLE_MPI`                     | `OFF`          | Install OpenMPI; build SuperLU_DIST                 |
| `BUILD_CGAL`                     | `OFF`          | Prebuild CGAL (+ Boost)                             |
| `BUILD_MUMPS`                    | `OFF`          | Prebuild MUMPS                                      |
| `BUILD_HYPRE`                    | `OFF`          | Prebuild Hypre                                      |
| `BUILD_TRILINOS`                 | `OFF`          | Prebuild Trilinos (slow!)                           |
| `USERNAME`/`USER_UID`/`USER_GID` | `oomph`/1000   | Container user                                      |
| `TPL_JOBS`                       | `4`            | Parallel jobs for the third-party build             |
| `OOMPH_JOBS`                     | `4`            | Parallel jobs for the `library` target              |
