# Building oomph-lib in a docker container

Add description here...

## Installation

```bash
# docker build -t <name>:<version> -f <path-to-dockerfile> <location-to-build>
docker build -t oomphlib:v1 -f ../Dockerfile .
```

## Supported operating systems

By default, a Docker environment will build with Ubuntu. You can switch to macOS
if you wish.

### Instructions for switching to macOS

...Fill this in...

## Debugging

If the docker build fails at any step, you should inspect the failing build by restarting from the preceding step in an interactive session. By default, the docker build kit hides a large amount of output during the build, including the container ID created at each step. To make the docker build run "verbosely", disable the docker build kit by setting the `DOCKER_BUILDKIT` environment variable before the call to `docker build`:

```bash
DOCKER_BUILDKIT=0 docker build -t oomphlib:v1 -f ../Dockerfile .
```

### Example

```bash
...
Step 7/9 : RUN cd oomph-lib && ls -l
 ---> Using cache
 ---> df28dae4db64
Step 8/9 : RUN cd oomph-lib && cmake --preset optimised -B build
 ---> Running in e63a02b2d538
CMake Error: The source directory "/oomph-lib/optimised" does not exist.
Specify --help for usage, or press the help button on the CMake GUI.
```

Now restart from the container cached from the previous step

```bash
docker run -ti --rm df28dae4db64 sh
```

