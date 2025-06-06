#!/usr/bin/env python3

import json
import shutil
import subprocess
import sys
import time
from argparse import BooleanOptionalAction, Namespace, ArgumentParser
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


# ANSI escape sequences for bold green text
@dataclass
class AnsiEscapeCodes:
    BOLD_RED = "\033[1;31m"
    BOLD_YELLOW = "\033[1;33m"
    BOLD_GREEN = "\033[1;32m"
    RESET_COLOR = "\033[0m"


class DirectoryExistsError(Exception):
    pass


def bold_red(text: str):
    return f"{AnsiEscapeCodes.BOLD_RED}{text}{AnsiEscapeCodes.RESET_COLOR}"


def bold_yellow(text: str):
    return f"{AnsiEscapeCodes.BOLD_YELLOW}{text}{AnsiEscapeCodes.RESET_COLOR}"


def bold_green(text: str):
    return f"{AnsiEscapeCodes.BOLD_GREEN}{text}{AnsiEscapeCodes.RESET_COLOR}"


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


def print_progress(message, *args, pad_to: Optional[int] = None, flush=True, **kwargs):
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


def write_presets_file(presets_dict: dict, json_path: Path):
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


def generate_oomph_preset(
    oomph_build_dir: Path,
    oomph_cache_vars: dict,
    generator="Ninja",
) -> dict:
    """
    Create a dictionary representing CMakeUserPresets.json for the oomph-lib project.
    oomph_cache_vars is a dict of {var_name: var_value}.
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
                "binaryDir": str(oomph_build_dir),
                "cacheVariables": oomph_cache_vars,
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
    # Convert ext- arguments to a dict of {OOMPH_*: value, CMAKE_*: value}
    # Example:
    #           ext_OOMPH_ENABLE_MPI -> OOMPH_ENABLE_MPI
    #           ext_CMAKE_BUILD_TYPE -> CMAKE_BUILD_TYPE
    ext_cache = {}
    for (opt, val) in vars(args).items():
        if val is None:
            continue
        if not (opt.startswith("ext_OOMPH") or opt.startswith("ext_CMAKE")):
            continue
        flag_name = opt.replace("ext_OOMPH", "OOMPH").replace("ext_CMAKE", "CMAKE")
        if isinstance(val, Path):
            val = str(val)
        ext_cache[flag_name] = val

    # Specify the build type; 'args.config' will always be set
    ext_cache["CMAKE_BUILD_TYPE"] = args.config

    # Do we need to add 'ccache' support?
    if args.enable_ccache:
        ext_cache["CMAKE_C_COMPILER_LAUNCHER"] = "ccache"
        ext_cache["CMAKE_CXX_COMPILER_LAUNCHER"] = "ccache"

    # Do we need to enable MPI?
    if args.OOMPH_ENABLE_MPI:
        ext_cache["OOMPH_ENABLE_MPI"] = args.OOMPH_ENABLE_MPI

    # Generate CMakeUserPresets.json for external project
    presets_dict = generate_external_dist_preset(external_dist_build_dir, ext_cache, args.generator)
    preset_path = external_dist_dir / "CMakeUserPresets.json"
    write_presets_file(presets_dict, preset_path)

    # The ending of the string
    end = "\n" if verbose else ""

    # Construct configuration command
    print_progress(">>> Configuring third-party libraries", pad_to=60, end=end)
    config_cmd = ["cmake", "--preset", "tpl"]
    if args.ext_extra_flags:
        config_cmd += args.ext_extra_flags
    start_time = time.perf_counter()
    run_command(config_cmd, external_dist_dir, verbose)
    time_elapsed = time.perf_counter() - start_time
    print_time(time_elapsed, verbose=verbose)

    # Build and install external
    print_progress(">>> Building and installing third-party libraries", pad_to=60, end=end)
    build_cmd = ["cmake", "--build", "--preset", "tpl"]
    start_time = time.perf_counter()
    run_command(build_cmd, external_dist_dir, verbose)
    time_elapsed = time.perf_counter() - start_time
    print_time(time_elapsed, verbose=verbose)


def configure_build_and_install_oomph(
    args: Namespace,
    oomph_dir: Path,
    oomph_build_dir: Path,
    ext_flags: dict,
    verbose: bool
):
    """
    Generate oomph-lib project presets, run cmake configure and build+install.
    ext_flags is the dict of flags read from the external JSON file.
    """
    # Convert oomph- arguments to a dict of {OOMPH_*: value, CMAKE_*: value}
    # Example:
    #           oomph_OOMPH_ENABLE_MPI -> OOMPH_ENABLE_MPI
    #           oomph_CMAKE_BUILD_TYPE -> CMAKE_BUILD_TYPE
    oomph_cache_vars = {}
    for (opt, val) in vars(args).items():
        if val is None:
            continue
        if not (opt.startswith("oomph_OOMPH") or opt.startswith("oomph_CMAKE")):
            continue
        flag_name = opt.replace("oomph_OOMPH", "OOMPH").replace("oomph_CMAKE", "CMAKE")
        if isinstance(val, Path):
            val = str(val)
        oomph_cache_vars[flag_name] = val

    # Specify the build type; 'args.config' will always be set
    oomph_cache_vars["CMAKE_BUILD_TYPE"] = args.config

    # Do we need to add 'ccache' support?
    if args.enable_ccache:
        oomph_cache_vars["CMAKE_C_COMPILER_LAUNCHER"] = "ccache"
        oomph_cache_vars["CMAKE_CXX_COMPILER_LAUNCHER"] = "ccache"

    # Do we need to enable MPI?
    if args.OOMPH_ENABLE_MPI:
        oomph_cache_vars["OOMPH_ENABLE_MPI"] = args.OOMPH_ENABLE_MPI

    # Merge in ext_flags first, then user `oomph_OOMPH_*` flags override if duplicated
    final_cache = dict(ext_flags)
    for (k, v) in oomph_cache_vars.items():
        final_cache[k] = v

    # Generate oomph-lib presets
    presets_dict = generate_oomph_preset(oomph_build_dir, final_cache, args.generator)
    preset_path = oomph_dir / "CMakeUserPresets.json"
    write_presets_file(presets_dict, preset_path)

    # Configure with the 'main' preset
    print_progress(">>> Configuring main project", pad_to=60, end="\n" if verbose else "")
    config_cmd = ["cmake", "--preset", "main"]
    start_time = time.perf_counter()
    run_command(config_cmd, oomph_dir, verbose)
    time_elapsed = time.perf_counter() - start_time
    print_time(time_elapsed, verbose=verbose)

    # Build and install in one go
    print_progress(">>> Building and installing main project", pad_to=60, end="\n" if verbose else "")
    build_cmd = ["cmake", "--build", "--preset", "main", "--target", "install"]
    start_time = time.perf_counter()
    run_command(build_cmd, oomph_dir, verbose)
    time_elapsed = time.perf_counter() - start_time
    print_time(time_elapsed, verbose=verbose)


def configure_build_doc(
    args: Namespace,
    doc_dir: Path,
    verbose: bool,
):
    """
    Configure and build the docs.
    """
    # Configure
    print_progress(">>> Configuring docs", pad_to=60, end="\n" if verbose else "")
    config_cmd = ["cmake", "-G", args.generator, "-B", "build"]
    start_time = time.perf_counter()
    run_command(config_cmd, doc_dir, verbose)
    time_elapsed = time.perf_counter() - start_time
    print_time(time_elapsed, verbose=verbose)

    # Build
    print_progress(">>> Building", pad_to=60, end="\n" if verbose else "")
    build_cmd = ["cmake", "--build", "build"]
    start_time = time.perf_counter()
    run_command(build_cmd, doc_dir, verbose)
    time_elapsed = time.perf_counter() - start_time
    print_time(time_elapsed, verbose=verbose)


def is_empty_and_untracked_tracked_dir(path: Path) -> bool:
    """
    Return True if `path` (file or folder) is known/tracked in the current Git index.
    If you pass a directory, Git treats it as "does any tracked item start with path/?"
    """
    is_empty_dir = not any(path.iterdir())

    # Make sure we pass a relative or absolute path in a form Git understands.
    # If you're not at repo root, you may need to prefix with the repo’s root dir
    # or `-C <repo_root>`; here we assume you run this from inside the repo.
    is_untracked = True
    try:
        subprocess.run(
            ["git", "ls-files", "--error-unmatch", path],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True
        )
        is_untracked = False
    except subprocess.CalledProcessError:
        pass
    return is_empty_dir and is_untracked


def wipe_dir_if_found(p: Path):
    if p.exists() and p.is_dir():
        shutil.rmtree(p)


def parse_args():
    """
    Parse and return command-line arguments via argparse.

    NOTE: argparse automatically converts flags with hyphens inbetween to variables with underscores
    in, e.g. '--enable-ccache' --> 'args.enable_ccache'.
    """
    def expanded_path(p):
        if p is None:
            return
        return Path(p).resolve()

    # fmt: off
    parser = ArgumentParser(description="Build automator.")

    parser.add_argument("-v", "--verbose", action="store_true", help="Don't suppress detailed output, show only high-level progress messages.")
    parser.add_argument("-s", "--silence-warnings-about-existing-build-and-install-directories", action="store_true", help="Suppress warnings about existing build/install directories.")

    # Only want to build TPLs or oomph-lib project?
    parser.add_argument("--build-tpl", action=BooleanOptionalAction, default=True, help="Build and install the third-party libraries (default: '--build-tpl').")
    parser.add_argument("--build-oomph", action=BooleanOptionalAction, default=True, help="Configure, build, and install the oomph-lib project (default: '--build-oomph').")
    parser.add_argument("--build-doc", action=BooleanOptionalAction, default=False, help="Configure and build the docs (default: '--no-build-docs').")

    # If the third-party libraries build/install dirs already exist, do we reuse or wipe?
    tpl_group = parser.add_mutually_exclusive_group()
    tpl_group.add_argument("--wipe-tpl", action="store_true", help="Wipe the third-party libraries build/install directories if they already exist.")
    tpl_group.add_argument("--reuse-tpl", action="store_true", help="Don't stop running if the third-party libraries build/install directories already exist.")

    # If the oomph-lib project build/install dirs already exist, do we reuse or wipe?
    oomph_group = parser.add_mutually_exclusive_group()
    oomph_group.add_argument("--wipe-oomph", action="store_true", help="Wipe the oomph-lib project build/install directories if they already exist.")
    oomph_group.add_argument("--reuse-oomph", action="store_true", help="Don't stop running if the oomph-lib project build/install directories already exist.")

    # If the oomph-lib project build/install dirs already exist, do we reuse or wipe?
    doc_group = parser.add_mutually_exclusive_group()
    doc_group.add_argument("--wipe-doc", action="store_true", help="Wipe the doc build directory if it already exists.")
    doc_group.add_argument("--reuse-doc", action="store_true", help="Don't stop running if the doc build directory already exists.")

    # Flags recognised by CMake
    general_group = parser.add_argument_group("general cmake flags")
    general_group.add_argument("-c", "--config", default="Release", choices=["Debug", "Release", "RelWithDebInfo", "MinSizeRel"], help="Specify CMAKE_BUILD_TYPE (default: Release) for both external and oomph-lib.")
    general_group.add_argument("-g", "--generator", default="Ninja", help="CMake generator to use. For example, 'Unix Makefiles' or 'Ninja' (default: Ninja).")
    general_group.add_argument(      "--enable-ccache", action="store_true", help="Enable use of 'ccache' for compilation.")

    # Flags common to both external_distributions and the oomph-lib project
    common_group = parser.add_argument_group("common build flags")
    common_group.add_argument("--OOMPH_ENABLE_MPI", metavar="ON/OFF", choices=["ON", "OFF"], help="Enable MPI in both external distributions and oomph-lib project.")

    # External distributions flags
    ext_group = parser.add_argument_group("external_distributions flags")
    ext_group.add_argument("--ext-CMAKE_INSTALL_PREFIX", type=expanded_path, metavar="PATH", help="Custom install directory for third-party libraries.")
    ext_group.add_argument("--ext-OOMPH_BUILD_OPENBLAS", metavar="ON/OFF", choices=["ON", "OFF"], help="Build OpenBLAS in third-party libraries.")
    ext_group.add_argument("--ext-OOMPH_BUILD_SUPERLU", metavar="ON/OFF", choices=["ON", "OFF"], help="Build SuperLU in third-party libraries.")
    ext_group.add_argument("--ext-OOMPH_BUILD_SUPERLU_DIST", metavar="ON/OFF", choices=["ON", "OFF"], help="Build SuperLU_DIST in third-party libraries.")
    ext_group.add_argument("--ext-OOMPH_BUILD_CGAL", metavar="ON/OFF", choices=["ON", "OFF"], help="Build CGAL in third-party libraries.")
    ext_group.add_argument("--ext-OOMPH_BUILD_MUMPS", metavar="ON/OFF", choices=["ON", "OFF"], help="Build MUMPS in third-party libraries.")
    ext_group.add_argument("--ext-OOMPH_BUILD_HYPRE", metavar="ON/OFF", choices=["ON", "OFF"], help="Build Hypre in third-party libraries.")
    ext_group.add_argument("--ext-OOMPH_BUILD_TRILINOS", metavar="ON/OFF", choices=["ON", "OFF"], help="Build Trilinos in third-party libraries.")
    ext_group.add_argument("--ext-OOMPH_ENABLE_THIRD_PARTY_LIBRARY_TESTS", metavar="ON/OFF", choices=["ON", "OFF"], help="Enable tests for third-party libraries.")
    ext_group.add_argument("--ext-OOMPH_USE_OPENBLAS_FROM", type=expanded_path, metavar="PATH", help="Use a preinstalled OpenBLAS from the given path.")
    ext_group.add_argument("--ext-OOMPH_USE_GKLIB_FROM", type=expanded_path, metavar="PATH", help="Use a preinstalled GKlib from the given path.")
    ext_group.add_argument("--ext-OOMPH_USE_METIS_FROM", type=expanded_path, metavar="PATH", help="Use a preinstalled METIS from the given path.")
    ext_group.add_argument("--ext-OOMPH_USE_PARMETIS_FROM", type=expanded_path, metavar="PATH", help="Use a preinstalled ParMETIS from the given path.")
    ext_group.add_argument("--ext-OOMPH_USE_BOOST_FROM", type=expanded_path, metavar="PATH", help="Use a preinstalled Boost from the given path.")

    # Any additional flags for the external distributions project
    ext_group.add_argument("--ext-extra-flags", nargs="+", dest="ext_extra_flags", metavar="FLAG", help="Additional raw CMake flags for external_distributions (e.g. -DXYZ=VALUE).")

    # Root project oomph-lib flags
    oomph_group = parser.add_argument_group("oomph-lib project flags")
    oomph_group.add_argument("--oomph-CMAKE_INSTALL_PREFIX", type=expanded_path, metavar="PATH", help="Custom installation directory for the main project.")
    oomph_group.add_argument("--oomph-OOMPH_ALLOW_INSTALL_AS_SUPERUSER", metavar="ON/OFF", choices=["ON", "OFF"], help="Allow the user to install to the default system install path (if CMAKE_INSTALL_PREFIX is not set).")
    oomph_group.add_argument("--oomph-OOMPH_INSTALL_HEADERS_AS_SYMLINKS", metavar="ON/OFF", choices=["ON", "OFF"], help="Install symlinks to the oomph-lib headers instead of copying them (default: OFF).")
    oomph_group.add_argument("--oomph-OOMPH_DONT_SILENCE_USELESS_WARNINGS", metavar="ON/OFF", choices=["ON", "OFF"], help="Don't silence certain warnings in oomph-lib.")
    oomph_group.add_argument("--oomph-OOMPH_ENABLE_MPI_OVERSUBSCRIPTION", metavar="ON/OFF", choices=["ON", "OFF"], help="Allow MPI oversubscription in oomph-lib.")
    oomph_group.add_argument("--oomph-OOMPH_ENABLE_PARANOID", metavar="ON/OFF", choices=["ON", "OFF"], help="Enable paranoid checks in oomph-lib.")
    oomph_group.add_argument("--oomph-OOMPH_ENABLE_RANGE_CHECKING", metavar="ON/OFF", choices=["ON", "OFF"], help="Enable range checking in oomph-lib.")
    oomph_group.add_argument("--oomph-OOMPH_ENABLE_SANITISER_SUPPORT", metavar="ON/OFF", choices=["ON", "OFF"], help="Enable sanitizer support in oomph-lib.")
    oomph_group.add_argument("--oomph-OOMPH_ENABLE_MUMPS_AS_DEFAULT_LINEAR_SOLVER", metavar="ON/OFF", choices=["ON", "OFF"], help="Use MUMPS as the default solver in oomph-lib.")
    oomph_group.add_argument("--oomph-OOMPH_SUPPRESS_TRIANGLE_LIB", metavar="ON/OFF", choices=["ON", "OFF"], help="Suppress usage of Triangle library.")
    oomph_group.add_argument("--oomph-OOMPH_SUPPRESS_TETGEN_LIB", metavar="ON/OFF", choices=["ON", "OFF"], help="Suppress usage of TetGen library.")
    # fmt: on
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()

    project_root = Path(__file__).resolve().parent
    external_dist_dir = project_root / "external_distributions"
    doc_dir = project_root / "doc"

    # Third-party libraries working directories
    external_dist_build_dir = external_dist_dir / "build"
    external_dist_install_dir = external_dist_dir / "install"

    # Root project working directories
    oomph_build_dir = project_root / "build"
    oomph_install_dir = project_root / "install"

    # Root project working directories
    doc_build_dir = doc_dir / "build"

    # If the user has provided custom install directories
    if args.ext_CMAKE_INSTALL_PREFIX:
        external_dist_install_dir = args.ext_CMAKE_INSTALL_PREFIX
    if args.oomph_CMAKE_INSTALL_PREFIX:
        oomph_install_dir = args.oomph_CMAKE_INSTALL_PREFIX

    if args.wipe_tpl:
        wipe_dir_if_found(external_dist_build_dir)
        wipe_dir_if_found(external_dist_install_dir)

    if args.wipe_oomph:
        wipe_dir_if_found(oomph_build_dir)
        wipe_dir_if_found(oomph_install_dir)

    if args.wipe_doc:
        # There is a 'clean' target that can be invoked to clean up *most* (but not all) of
        # the generated files.
        if doc_build_dir.exists():
            clean_up_cmd = ["cmake", "--build", "build", "--target", "clean"]
            run_command(clean_up_cmd, doc_dir, args.verbose)

        # Now we need to go to the demo drivers directory and clean up the empty validata/
        # directories (which used to have an index.html file in)
        demo_drivers_dir = project_root / "demo_drivers"
        if not demo_drivers_dir.exists():
            raise FileNotFoundError(f"Unable to locate 'demo_drivers' directory at: {demo_drivers_dir}")

        # Clear out the generated validata/ dirs. We'll make sure they're empty and not tracked
        # by Git before we actually delete them
        for d in demo_drivers_dir.rglob("validata/"):
            if is_empty_and_untracked_tracked_dir(d):
                if args.verbose:
                    print(f"Removing '{d}'")
                d.rmdir()

        # Now kill the build directory
        wipe_dir_if_found(doc_build_dir)

    # Where to inherit the flags output by external_distributions after we've built
    # the third-party libraries that we want
    external_dist_json_path = external_dist_install_dir / "cmake_flags_for_oomph_lib.json"

    # Warn the user that the third-party libraries already exist. We can skip this warning if
    # the user is skipping their build and just building the top level project, or if they've
    # silenced the warnings
    if not args.silence_warnings_about_existing_build_and_install_directories:
        found_dirty_directory = False

        barrier = "=" * 100
        if args.build_tpl:
            for d in (external_dist_build_dir, external_dist_install_dir):
                if d.exists():
                    found_dirty_directory = True
                    if not args.reuse_tpl:
                        raise DirectoryExistsError(
                            f"The directory:\n\n\t{d}\n\nalready exists. If you have changed any settings you should wipe this directory before\nreinstalling oomph-lib. You can bypass this error with the flag '--reuse-tpl' or wipe them with '--wipe-tpl'"
                        )

        # Warn the user that the oomph-lib project build/install dirs already exist. We can skip this
        # warning if the user is only building the third-party libraries, or if they've silenced
        # the warnings
        if args.build_oomph:
            for d in (oomph_build_dir, oomph_install_dir):
                if d.exists():
                    found_dirty_directory = True
                    if not args.reuse_oomph:
                        raise DirectoryExistsError(
                            f"The directory:\n\n\t{d}\n\nalready exists. If you have changed any settings you should wipe this directory before\nreinstalling oomph-lib. You can bypass this error with the flag '--reuse-oomph' or wipe them with '--wipe-oomph'"
                        )

        # Warn the user that the oomph-lib project build/install dirs already exist. We can skip this
        # warning if the user is only building the third-party libraries, or if they've silenced
        # the warnings
        if args.build_doc:
            if doc_build_dir.exists():
                found_dirty_directory = True
                if not args.reuse_doc:
                    raise DirectoryExistsError(
                        f"The directory:\n\n\t{doc_build_dir}\n\nalready exists. If you have changed any settings you should wipe this directory before\nreinstalling oomph-lib. You can bypass this error with the flag '--reuse-doc' or wipe them with '--wipe-doc'"
                    )

        # Pause long enough for the user to see this warning
        if found_dirty_directory:
            time.sleep(2)

    # Let's go
    print_progress("\n>>> Starting project build...")
    start_time = time.perf_counter()

    # Configure, build and install external libraries, retrieving the JSON file path
    if args.build_tpl:
        configure_build_and_install_external_libs(args, external_dist_dir, external_dist_build_dir, args.verbose)
    else:
        print_progress(">>> Skipping build of third-party libraries")

    # If we've *silently* built the third-party libs but the user isn't building the oomph-lib
    # project, we should dump the usage instructions incase they want to do it themselves
    if args.build_tpl and not args.build_oomph and not args.verbose:
        contents = (external_dist_build_dir / "usage.txt").read_text()
        print(f"{contents}\n")

    # Configure, build, and install oomph-lib (main) project
    if args.build_oomph:
        # Attempt to locate 'cmake_flags_for_oomph_lib.json'
        if not external_dist_json_path.is_file():
            print("ERROR: Could not find 'cmake_flags_for_oomph_lib.json' after building external_distributions.", file=sys.stderr)
            sys.exit(1)

        # Read external JSON file to get flags to pass to oomph-lib project
        ext_flags = read_external_json_file(external_dist_json_path)

        configure_build_and_install_oomph(args, project_root, oomph_build_dir, ext_flags, args.verbose)
    else:
        print_progress(">>> Skipping build of oomph-lib project")

    # Configure and build the docs
    if args.build_doc:
        configure_build_doc(args, doc_dir, args.verbose)
    else:
        print_progress(">>> Skipping build of docs")

    print_progress(">>> Build and installation complete!", pad_to=60, end="")
    total_time_elapsed = time.perf_counter() - start_time
    print_time(total_time_elapsed, verbose=args.verbose)
