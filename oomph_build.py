#!/usr/bin/env python3

import json
import subprocess
import sys
import time
from argparse import Namespace, ArgumentParser
from pathlib import Path


# ANSI escape sequences for bold green text
BOLD_GREEN = "\033[1;32m"
RESET_COLOR = "\033[0m"


def read_external_json_file(json_path):
    """
    Read a JSON file that contains the external CMake flags (e.g. cmake_flags_for_oomph_lib.json).
    """
    try:
        # data is expected to be a dict of { str -> str }, e.g.
        # {
        #   "CMAKE_BUILD_TYPE": "Release",
        #   "OOMPH_ENABLE_MPI": "ON",
        #   ...
        # }
        with open(json_path, "r") as f:
            data = json.load(f)
            return data
    except (json.JSONDecodeError, FileNotFoundError) as e:
        print(f"ERROR: Failed to read or parse {json_path}: {e}", file=sys.stderr)
        sys.exit(1)


def run_command(cmd: list[str], cwd: Path, verbose: bool):
    """
    Run a subprocess command given by list 'cmd'.
    If verbose is False, capture output; on failure, print the output and exit.
    """
    try:
        if not verbose:
            result = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        else:
            result = subprocess.run(cmd, cwd=cwd)
    except Exception as e:
        print(f"ERROR: failed to run command: {' '.join(cmd)}\n{e}", file=sys.stderr)
        sys.exit(1)

    if result.returncode != 0:
        print(f"ERROR: Command {' '.join(cmd)} failed with exit code {result.returncode}.", file=sys.stderr)
        if not verbose:
            sys.stderr.write(result.stdout)
        sys.exit(result.returncode)


def bold_green(text: str):
    return f"{BOLD_GREEN}{text}{RESET_COLOR}"


def print_progress(message, *args, pad_to: int | None = None, flush=True, **kwargs):
    """
    Print progress messages in bold green.
    """
    if pad_to is None:
        print(bold_green(message), *args, flush=flush, **kwargs)
    else:
        print(bold_green(f"{message:<{pad_to}}"), *args, flush=flush, **kwargs)


def print_time(time_elapsed: float, *args, verbose: bool = False, flush=True, **kwargs):
    if verbose:
        print(bold_green(f"Time taken [sec]: {time_elapsed:.2f}\n"), *args, flush=flush, **kwargs)
    else:
        print(bold_green(f"[{time_elapsed:.2f} s]"), *args, flush=flush, **kwargs)


def write_presets_file(presets_dict, json_path):
    """
    Write the given presets_dict to json_path as CMakeUserPresets.json
    """
    with open(json_path, "w") as f:
        json.dump(presets_dict, f, indent=2)


def generate_external_dist_preset(
    external_dist_build_dir: Path,
    ext_cache: dict,
    generator: str = "Ninja",
) -> dict:
    """
    Create a dictionary representing CMakeUserPresets.json for external_distributions project.
    ext_cache is a dict of {var_name: var_value} for external OOMPH_ flags.
    """
    preset_dict = {
        "version": 3,
        "cmakeMinimumRequired": {"major": 3, "minor": 24, "patch": 0},
        "configurePresets": [
            {
                "name": "tpl",
                "displayName": "Third-party libraries",
                "description": "Configure third-party libraries",
                "generator": generator,
                "binaryDir": str(external_dist_build_dir),
                "cacheVariables": ext_cache,
            }
        ],
        "buildPresets": [
            {
                "name": "tpl",
                "configurePreset": "tpl",
            }
        ],
    }
    return preset_dict


def generate_root_preset(
    root_build: Path,
    root_cache_vars: dict,
    generator = "Ninja",
) -> dict:
    """
    Create a dictionary representing CMakeUserPresets.json for the root oomph-lib project.
    root_cache_vars is a dict of {var_name: var_value}.
    """
    preset_dict = {
        "version": 3,
        "cmakeMinimumRequired": {"major": 3, "minor": 24, "patch": 0},
        "configurePresets": [
            {
                "name": "main",
                "displayName": "Main project",
                "description": "Configure main project",
                "generator": generator,
                "binaryDir": str(root_build),
                "cacheVariables": root_cache_vars,
            }
        ],
        "buildPresets": [
            {
                "name": "main",
                "configurePreset": "main",
            }
        ],
    }
    return preset_dict


def configure_build_and_install_external_libs(
    args: Namespace,
    external_dist_dir: Path,
    external_dist_build_dir: Path,
    verbose: bool,
):
    """
    Generate external project presets, run cmake configure and build+install, then
    return the path to cmake_flags_for_oomph_lib.json (if found).
    """
    # Convert ext- arguments to a dict of {OOMPH_*: value}
    ext_cache = {}
    for (opt, val) in vars(args).items():
        if (opt.startswith("ext_OOMPH")) and (val is not None):
            # Example: ext_OOMPH_ENABLE_MPI -> OOMPH_ENABLE_MPI
            flag_name = opt.replace("ext_OOMPH", "OOMPH")
            ext_cache[flag_name] = val

    # Generate CMakeUserPresets.json for external project
    presets_dict = generate_external_dist_preset(external_dist_build_dir, ext_cache, args.generator)
    preset_path = external_dist_dir / "CMakeUserPresets.json"
    write_presets_file(presets_dict, preset_path)

    # Construct configuration command
    print_progress(">>> Configuring third-party libraries", pad_to=60, end="\n" if verbose else "")
    config_cmd = ["cmake", "--preset", "tpl"]
    if args.ext_extra_flags:
        config_cmd += args.ext_extra_flags
    start_time = time.perf_counter()
    run_command(config_cmd, external_dist_dir, verbose)
    time_elapsed = time.perf_counter() - start_time
    print_time(time_elapsed, verbose=verbose)

    # Build and install external
    print_progress(">>> Building and installing", pad_to=60, end="\n" if verbose else "")
    build_cmd = ["cmake", "--build", "--preset", "tpl"]
    start_time = time.perf_counter()
    run_command(build_cmd, external_dist_dir, verbose)
    time_elapsed = time.perf_counter() - start_time
    print_time(time_elapsed, verbose=verbose)


def configure_build_and_install_root(
    args: Namespace,
    root_dir: Path,
    root_build: Path,
    ext_flags: dict,
    verbose: bool
):
    """
    Generate root project presets, run cmake configure and build+install.
    ext_flags is the dict of flags read from the external JSON file.
    """
    # Convert root- arguments to a dict
    root_cache_vars = {}
    for (opt, val) in vars(args).items():
        if opt.startswith("root_OOMPH") and val is not None:
            flag_name = opt.replace("root_OOMPH", "OOMPH")
            root_cache_vars[flag_name] = val

    # Merge in ext_flags first, then user `root_OOMPH_*` flags override if duplicated
    final_cache = dict(ext_flags)
    for (k, v) in root_cache_vars.items():
        final_cache[k] = v

    # Generate root presets
    presets_dict = generate_root_preset(root_build, final_cache, args.generator)
    preset_path = root_dir / "CMakeUserPresets.json"
    write_presets_file(presets_dict, preset_path)

    # Configure with the 'main' preset
    print_progress(">>> Configuring main project", pad_to=60, end="\n" if verbose else "")
    config_cmd = ["cmake", "--preset", "main"]
    start_time = time.perf_counter()
    run_command(config_cmd, root_dir, verbose)
    time_elapsed = time.perf_counter() - start_time
    print_time(time_elapsed, verbose=verbose)

    # Build and install in one go
    print_progress(">>> Building and installing main project", pad_to=60, end="\n" if verbose else "")
    build_cmd = ["cmake", "--build", "--preset", "main", "--target", "install"]
    start_time = time.perf_counter()
    run_command(build_cmd, root_dir, verbose)
    time_elapsed = time.perf_counter() - start_time
    print_time(time_elapsed, verbose=verbose)


def parse_args():
    """Parse and return command-line arguments via argparse."""
    # fmt: off
    parser = ArgumentParser(description="Build automator.")

    parser.add_argument("-v", "--verbose", action="store_true", help="Don't suppress detailed output, show only high-level progress messages.")

    # Only want to build TPLs or root project?
    parser.add_argument("--skip-tpl-build", action="store_true", help="Skip building and installing the external_distributions/ libraries.")
    parser.add_argument("--just-build-tpl", action="store_true", help="Skip configuring/building/installing the root project.")

    # Flags recognised by CMake
    general_group = parser.add_argument_group("general cmake flags")
    general_group.add_argument("-b", "--build-type", default="Release", choices=["Debug", "Release", "RelWithDebInfo", "MinSizeRel"], help="Specify CMAKE_BUILD_TYPE (default: Release). Applies to both external and root.")
    general_group.add_argument("-g", "--generator", default="Ninja", help="CMake generator to use (default: Ninja). For example, 'Unix Makefiles' or 'Ninja'.")

    # Flags common to both external_distributions and the root project
    common_group = parser.add_argument_group("common build flags")
    common_group.add_argument("--OOMPH_ENABLE_MPI", metavar="ON/OFF", help="Enable MPI in both external distributions and root project.")

    # External distributions flags
    ext_group = parser.add_argument_group("external_distributions flags")
    ext_group.add_argument("--ext-OOMPH_BUILD_OPENBLAS", metavar="ON/OFF", help="Build OpenBLAS in third-party libraries")
    ext_group.add_argument("--ext-OOMPH_BUILD_SUPERLU", metavar="ON/OFF", help="Build SuperLU in third-party libraries")
    ext_group.add_argument("--ext-OOMPH_BUILD_SUPERLU_DIST", metavar="ON/OFF", help="Build SuperLU_DIST in third-party libraries")
    ext_group.add_argument("--ext-OOMPH_BUILD_CGAL", metavar="ON/OFF", help="Build CGAL in third-party libraries")
    ext_group.add_argument("--ext-OOMPH_BUILD_MUMPS", metavar="ON/OFF", help="Build MUMPS in third-party libraries")
    ext_group.add_argument("--ext-OOMPH_BUILD_HYPRE", metavar="ON/OFF", help="Build Hypre in third-party libraries")
    ext_group.add_argument("--ext-OOMPH_BUILD_TRILINOS", metavar="ON/OFF", help="Build Trilinos in third-party libraries")
    ext_group.add_argument("--ext-OOMPH_DISABLE_THIRD_PARTY_LIBRARY_TESTS", metavar="ON/OFF", help="Disable tests for third-party libraries")
    ext_group.add_argument("--ext-OOMPH_THIRD_PARTY_INSTALL_DIR", metavar="PATH", help="Custom install directory for third-party libraries")
    ext_group.add_argument("--ext-OOMPH_USE_OPENBLAS_FROM", metavar="PATH", help="Use a preinstalled OpenBLAS from the given path")
    ext_group.add_argument("--ext-OOMPH_USE_GKLIB_FROM", metavar="PATH", help="Use a preinstalled GKlib from the given path")
    ext_group.add_argument("--ext-OOMPH_USE_METIS_FROM", metavar="PATH", help="Use a preinstalled METIS from the given path")
    ext_group.add_argument("--ext-OOMPH_USE_PARMETIS_FROM", metavar="PATH", help="Use a preinstalled ParMETIS from the given path")
    ext_group.add_argument("--ext-OOMPH_USE_BOOST_FROM", metavar="PATH", help="Use a preinstalled Boost from the given path")

    # Any additional flags for the external distributions project
    ext_group.add_argument("--ext-extra-flags", nargs="+", dest="ext_extra_flags", metavar="FLAG", help="Additional raw CMake flags for external_distributions (e.g. -DXYZ=VALUE)")

    # Root project oomph-lib flags
    root_group = parser.add_argument_group("root project oomph-lib flags")
    root_group.add_argument("--root-OOMPH_DONT_SILENCE_USELESS_WARNINGS", metavar="ON/OFF", help="Don't silence certain warnings in oomph-lib")
    root_group.add_argument("--root-OOMPH_ENABLE_MPI_OVERSUBSCRIPTION", metavar="ON/OFF", help="Allow MPI oversubscription in oomph-lib")
    root_group.add_argument("--root-OOMPH_ENABLE_PARANOID", metavar="ON/OFF", help="Enable paranoid checks in oomph-lib")
    root_group.add_argument("--root-OOMPH_ENABLE_RANGE_CHECKING", metavar="ON/OFF", help="Enable range checking in oomph-lib")
    root_group.add_argument("--root-OOMPH_ENABLE_SANITISER_SUPPORT", metavar="ON/OFF", help="Enable sanitizer support in oomph-lib")
    root_group.add_argument("--root-OOMPH_ENABLE_MUMPS_AS_DEFAULT_LINEAR_SOLVER", metavar="ON/OFF", help="Use MUMPS as the default solver in oomph-lib")
    root_group.add_argument("--root-OOMPH_SUPPRESS_TRIANGLE_LIB", metavar="ON/OFF", help="Suppress usage of Triangle library")
    root_group.add_argument("--root-OOMPH_SUPPRESS_TETGEN_LIB", metavar="ON/OFF", help="Suppress usage of TetGen library")
    # fmt: on

    args = parser.parse_args()

    # Set the build type for both external_distributions and the root project
    setattr(args, "ext_CMAKE_BUILD_TYPE", args.build_type)
    setattr(args, "root_CMAKE_BUILD_TYPE", args.build_type)

    # If the user specifies --OOMPH_ENABLE_MPI, override both ext- and root- settings to "ON"
    if args.OOMPH_ENABLE_MPI:
        setattr(args, "ext_OOMPH_ENABLE_MPI", args.OOMPH_ENABLE_MPI)
        setattr(args, "root_OOMPH_ENABLE_MPI", args.OOMPH_ENABLE_MPI)
    return args

if __name__ == "__main__":
    args = parse_args()

    project_root = Path(__file__).resolve().parent
    root_build = project_root / "build"
    external_dist_dir = project_root / "external_distributions"
    external_dist_build_dir = external_dist_dir / "build"

    # Create build directories
    external_dist_build_dir.mkdir(parents=True, exist_ok=True)
    root_build.mkdir(parents=True, exist_ok=True)

    # Let's go
    print_progress("\n>>> Starting project build...")
    start_time = time.perf_counter()

    # 1. Configure, build and install external libraries, retrieving the JSON file path
    if not args.skip_tpl_build:
        configure_build_and_install_external_libs(args, external_dist_dir, external_dist_build_dir, args.verbose)
    else:
        print_progress(">>> Skipping third-party libraries build")

    # Attempt to locate 'cmake_flags_for_oomph_lib.json'
    external_json_path = external_dist_build_dir / "cmake_flags_for_oomph_lib.json"
    if not external_json_path.is_file():
        print("ERROR: Could not find 'cmake_flags_for_oomph_lib.json' after building external_distributions.", file=sys.stderr)
        sys.exit(1)

    # 2. Read external JSON file to get flags to pass to root project
    ext_flags = read_external_json_file(external_json_path)

    # 3. Configure, build, and install root (main) project
    if not args.just_build_tpl:
        configure_build_and_install_root(args, project_root, root_build, ext_flags, args.verbose)
    else:
        print_progress(">>> Skipping root project build")

    print_progress(">>> Build and installation complete!", pad_to=60, end="")
    total_time_elapsed = time.perf_counter() - start_time
    print_time(total_time_elapsed, verbose=args.verbose)

