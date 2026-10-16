# syntax=docker/dockerfile:1
# ------------------------------------------------------------------------------
# oomph-lib development container
#
# Stages:
#
#   base    -- Ubuntu + the toolchain needed to build oomph-lib (gcc/g++/
#              gfortran, CMake >= 3.24, Ninja, ccache, Python 3, optional
#              OpenMPI) and a non-root user.
#
#   tpl     -- Builds the third-party libraries from external_distributions/
#              (OpenBLAS + SuperLU by default; CGAL/MUMPS/Hypre/Trilinos can be
#              switched on with build args) and installs them into
#              /opt/oomph-lib/tpl so that they don't have to be rebuilt every
#              time the source tree changes.
#
#   dev     -- (default target) base + prebuilt third-party libraries. Intended
#              to be used with the repository bind-mounted into
#              /workspaces/oomph-lib (this is what .devcontainer/ does).
#
#   library -- dev + a full build/install of oomph-lib itself from the current
#              checkout. Useful for a one-off, reproducible build:
#
#                  docker build --target library -t oomph-lib:library .
#
# Typical usage of the default (dev) image:
#
#     docker build -t oomph-lib:dev .
#     docker run --rm -it -v "$PWD":/workspaces/oomph-lib oomph-lib:dev
#     oomph-build                            # wrapper around oomph_build.py
#
# Build args (all optional):
#
#     ENABLE_MPI=ON|OFF           install OpenMPI and build MPI-enabled TPLs
#     BUILD_CGAL=ON|OFF           } which third-party libraries to prebuild
#     BUILD_MUMPS=ON|OFF          } (OpenBLAS + SuperLU are always built as
#     BUILD_HYPRE=ON|OFF          } they're required by oomph-lib)
#     BUILD_TRILINOS=ON|OFF       }
#     USERNAME/USER_UID/USER_GID  non-root user (match your host UID/GID so
#                                 that bind-mounted files keep sane ownership)
#     TPL_JOBS                    parallel jobs for the third-party build
# ------------------------------------------------------------------------------

ARG BASE_IMAGE=ubuntu:24.04

# ==============================================================================
FROM ${BASE_IMAGE} AS base
# ==============================================================================

ARG ENABLE_MPI=OFF
ARG USERNAME=oomph
ARG USER_UID=1000
ARG USER_GID=${USER_UID}

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# Toolchain. Ubuntu 24.04 ships CMake 3.28 which satisfies oomph-lib's
# requirement of CMake >= 3.24; if you change BASE_IMAGE to something older,
# use scripts/install_cmake.sh instead of the apt package.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        gfortran \
        cmake \
        ninja-build \
        ccache \
        git \
        python3 \
        python3-pip \
        pkg-config \
        ca-certificates \
        curl \
        wget \
        patch \
        unzip \
        less \
        vim-tiny \
        sudo \
        gdb \
        valgrind \
        clang-format \
    && if [ "${ENABLE_MPI}" = "ON" ]; then \
         apt-get install -y --no-install-recommends openmpi-bin libopenmpi-dev; \
       fi \
    && rm -rf /var/lib/apt/lists/*

# Non-root user. Ubuntu 24.04 already ships a user with UID 1000 ("ubuntu");
# remove it so that we can reuse the UID for our own user.
RUN if getent passwd "${USER_UID}" >/dev/null; then \
        userdel -r "$(getent passwd "${USER_UID}" | cut -d: -f1)" || true; \
    fi \
    && if ! getent group "${USER_GID}" >/dev/null; then \
        groupadd --gid "${USER_GID}" "${USERNAME}"; \
    fi \
    && useradd --uid "${USER_UID}" --gid "${USER_GID}" --create-home --shell /bin/bash "${USERNAME}" \
    && echo "${USERNAME} ALL=(ALL) NOPASSWD:ALL" > "/etc/sudoers.d/${USERNAME}" \
    && chmod 0440 "/etc/sudoers.d/${USERNAME}"

# Let MPI run as root / oversubscribe inside the container (harmless if MPI is
# not installed).
ENV OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_rmaps_base_oversubscribe=1

# Where the prebuilt third-party libraries live (consumed by the stages below
# and by the oomph-build helper).
ENV OOMPH_TPL_INSTALL_DIR=/opt/oomph-lib/tpl \
    OOMPH_ENABLE_MPI=${ENABLE_MPI}

# ==============================================================================
FROM base AS tpl
# ==============================================================================

ARG ENABLE_MPI=OFF
ARG BUILD_CGAL=OFF
ARG BUILD_MUMPS=OFF
ARG BUILD_HYPRE=OFF
ARG BUILD_TRILINOS=OFF
ARG TPL_JOBS=4

WORKDIR /tmp/oomph-lib

# Only what the third-party build needs; keeps this layer cacheable when the
# rest of the source tree changes.
COPY oomph_build.py ./
COPY external_distributions/CMakeLists.txt external_distributions/
COPY external_distributions/cmake external_distributions/cmake

RUN python3 oomph_build.py \
        --no-build-oomph \
        --OOMPH_ENABLE_MPI="${ENABLE_MPI}" \
        --ext-CMAKE_INSTALL_PREFIX="${OOMPH_TPL_INSTALL_DIR}" \
        --ext-OOMPH_BUILD_OPENBLAS=ON \
        --ext-OOMPH_BUILD_SUPERLU=ON \
        --ext-OOMPH_BUILD_SUPERLU_DIST="${ENABLE_MPI}" \
        --ext-OOMPH_BUILD_CGAL="${BUILD_CGAL}" \
        --ext-OOMPH_BUILD_MUMPS="${BUILD_MUMPS}" \
        --ext-OOMPH_BUILD_HYPRE="${BUILD_HYPRE}" \
        --ext-OOMPH_BUILD_TRILINOS="${BUILD_TRILINOS}" \
        -j "${TPL_JOBS}" \
        --verbose \
    && rm -rf /tmp/oomph-lib

# ==============================================================================
FROM base AS dev
# ==============================================================================

ARG USERNAME=oomph

COPY --from=tpl ${OOMPH_TPL_INSTALL_DIR} ${OOMPH_TPL_INSTALL_DIR}

# Small convenience wrapper: builds oomph-lib against the prebuilt third-party
# libraries. Any extra arguments are forwarded to oomph_build.py, e.g.
#
#     oomph-build --oomph-OOMPH_INSTALL_HEADERS_AS_SYMLINKS=ON
#     oomph-build -c Debug --oomph-OOMPH_ENABLE_PARANOID=ON --reuse-oomph
#
RUN printf '%s\n' \
        '#!/usr/bin/env bash' \
        'set -euo pipefail' \
        'if [ ! -f oomph_build.py ]; then' \
        '  echo "oomph-build: run this from the root of the oomph-lib checkout" >&2' \
        '  exit 1' \
        'fi' \
        'exec python3 oomph_build.py \' \
        '  --no-build-tpl \' \
        '  --enable-ccache \' \
        '  --OOMPH_ENABLE_MPI="${OOMPH_ENABLE_MPI:-OFF}" \' \
        '  --ext-CMAKE_INSTALL_PREFIX="${OOMPH_TPL_INSTALL_DIR}" \' \
        '  "$@"' \
        > /usr/local/bin/oomph-build \
    && chmod +x /usr/local/bin/oomph-build

ENV CCACHE_DIR=/home/${USERNAME}/.ccache \
    CCACHE_MAXSIZE=2G

RUN mkdir -p /workspaces/oomph-lib "${CCACHE_DIR}" \
    && chown -R "${USERNAME}:${USERNAME}" /workspaces "/home/${USERNAME}"

USER ${USERNAME}
WORKDIR /workspaces/oomph-lib

CMD ["/bin/bash"]

# ==============================================================================
FROM dev AS library
# ==============================================================================

ARG USERNAME=oomph
ARG OOMPH_JOBS=4

COPY --chown=${USERNAME}:${USERNAME} . /workspaces/oomph-lib

RUN oomph-build -j "${OOMPH_JOBS}" --oomph-OOMPH_INSTALL_HEADERS_AS_SYMLINKS=OFF --verbose
