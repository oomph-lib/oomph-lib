FROM ubuntu:20.04

# Docker containers work on a concept of "layers". Each time an extra step is
# run, this adds another "layer" to the container. If the build fails at step X
# and you try to rerun the build, docker will attempt to restart from the layer
# just before step X. Therefore, costly, and unlikely-to-break steps should be
# performed early in the installation so that the remaining steps are short in
# comparison. For this reason, we clone the oomph-lib repository as soon as we
# have a working copy of "git".

# Update package information from configured sources; get the timezone data
# early on to avoid the annoying interactive command questioning you about it
# you about it
RUN apt-get update
RUN apt-get install -y tzdata

# Install both general and oomph-lib-specific requirements
RUN apt-get install -y git wget

# Clone the repo early on. If the installation of our required packages fail later on,
# then when a cached container a
# cache of this container layer
RUN git clone https://github.com/puneetmatharu/oomph-lib.git
RUN cd oomph-lib && git checkout feature-convert-to-cmake-build-system
RUN cd oomph-lib && ls -l

# Install CMake v3.20+ in order to use CMake presets and add it to the PATH
RUN wget https://github.com/Kitware/CMake/releases/download/v3.20.5/cmake-3.20.5-Linux-x86_64.sh \
    -q -O /tmp/cmake-install.sh \
    && chmod u+x /tmp/cmake-install.sh \
    && mkdir /usr/bin/cmake \
    && /tmp/cmake-install.sh --skip-license --prefix=/usr/bin/cmake \
    && rm /tmp/cmake-install.sh
ENV PATH="/usr/bin/cmake/bin:${PATH}"

# Install general and oomph-lib-specific requirements
RUN apt-get install -y build-essential gfortran libopenblas-dev libopenmpi-dev ninja-build openmpi-bin

# Python requires sudo privileges
RUN apt-get install -y python3

# Run cmake
# RUN cd oomph-lib && cmake --preset optimised -B build
RUN cd oomph-lib && cmake -G Ninja \
    -D OOMPH_USE_OPENBLAS=ON \
    -D CMAKE_INSTALL_PREFIX="~/oomph-installation" \
    -B build

# Build and install project
RUN cd oomph-lib/build && ninja -j4 && ninja install