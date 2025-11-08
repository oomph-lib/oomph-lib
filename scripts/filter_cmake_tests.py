#!/usr/bin/env python3

import math
import os
import re
import sys
from argparse import ArgumentParser, Namespace, REMAINDER
from pathlib import Path
from subprocess import run
from typing import Set, List


def find_matching_directories(root_dir: Path, keywords: List[str], verbose: bool = False) -> Set[Path]:
    """
    Recursively search for C++ files containing any of the given keywords
    and return their parent directories.

    This version uses a single regex to match all keywords (OR logic).
    It reads files line by line to avoid loading entire files into memory.

    Args:
        root_dir (Path): The directory to search within.
        keywords (List[str]): One or more keywords to look for in source files.

    Returns:
        Set[Path]: A set of directories containing at least one C++ file that matches.
    """
    matching_dirs: Set[Path] = set()

    # Create a regex that matches *any* of the given keywords
    # e.g. keywords=["foo","bar"] => pattern: r"(foo|bar)"
    pattern_str = "|".join(re.escape(kw) for kw in keywords)
    keyword_pattern = re.compile(pattern_str)

    # We only look at these file patterns
    cpp_patterns = ["*.h", "*.cc"]

    # Collect all C++ source/header files first to know how many we are processing
    all_files = []
    for patt in cpp_patterns:
        all_files.extend(root_dir.rglob(patt))

    total_files = len(all_files)
    if total_files == 0:
        return matching_dirs

    for i, cpp_file in enumerate(all_files, start=1):
        progress_fraction = i / total_files
        progress_percent = math.floor(progress_fraction * 100)
        if verbose:
            print(f"\rScanning files: {progress_percent}% ({i}/{total_files})", end="", flush=True)

        try:
            # Read line by line to avoid loading entire file into memory
            with open(cpp_file, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if keyword_pattern.search(line):
                        matching_dirs.add(cpp_file.parent.resolve())
                        break  # No need to check more lines in this file
        except Exception as e:
            print(f"\nError reading {cpp_file}: {e}", file=sys.stderr)

    print()  # Move to a new line after the progress bar finishes
    return matching_dirs


def extract_test_names(cmake_file: Path) -> List[str]:
    """
    Extract test names from a CMakeLists.txt file that uses oomph_add_test or oomph_add_pure_cpp_test.

    Args:
        cmake_file (Path): Path to the CMakeLists.txt file.

    Returns:
        List[str]: A list of test names found.
    """
    test_names: List[str] = []
    pattern = re.compile(
        r"oomph_add_(?:pure_cpp_)?test\s*\((?:[^()]*\n)*?[^()]*?TEST_NAME\s+([^\s)]+)",
        re.DOTALL,
    )
    try:
        content = cmake_file.read_text(encoding="utf-8")
        matches = pattern.findall(content)
        test_names.extend(matches)
    except Exception as e:
        print(f"Error reading {cmake_file}: {e}", file=sys.stderr)
    return test_names


def parse_args() -> Namespace:
    """
    Parse command-line arguments.

    Returns:
        Namespace: Parsed arguments.
    """
    # fmt: off
    parser = ArgumentParser(description="Filters tests in demo_drivers/ based on one or more keywords in C++ files.")
    parser.add_argument("--root", required=True, help="Path to the 'demo_drivers' (sub)directory to filter tests from.")
    parser.add_argument("--keywords", nargs="+", help="One or more keywords to search for in C++ source files (OR match).")
    parser.add_argument("--run-from", default=None, type=Path, help="Path to the build directory to run the tests from.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print helpful information to the screen.")
    parser.add_argument('ctest_args', nargs="*", help="Extra arguments to pass to ctest (use '--' to separate them from other arguments).")
    args = parser.parse_args()
    # fmt: on

    # Sanity checks
    if args.run_from is not None:
        if not args.run_from.is_dir():
            parser.error("Expected argument to '--run-from' to be a directory.")
        if not (args.run_from / "CTestTestfile.cmake").exists():
            parser.error("Directory passed to '--run-from' does not contain a 'CTestTestfile.cmake', therefore we cannot run 'ctest' from it.")
    return args


def main() -> None:
    args = parse_args()

    demo_drivers_root = Path(args.root).resolve()
    if not demo_drivers_root.exists():
        raise FileNotFoundError(f"[ERROR]: Directory '{demo_drivers_root}' does not exist.")

    # Show which keywords we're searching for
    keyword_list_str = ", ".join(args.keywords)
    if args.verbose:
        print(f"Searching for keywords: {keyword_list_str}")
        print(f"Within directory: {os.path.relpath(demo_drivers_root, Path.cwd())}\n")

    # 1) Identify all directories that contain a file with the given keywords
    matching_dirs = find_matching_directories(demo_drivers_root, args.keywords, args.verbose)
    if not matching_dirs:
        print("No directories found containing any of the specified keywords.")
        sys.exit(0)

    # 2) From those directories, extract any test names from the local CMakeLists.txt
    all_test_names: List[str] = []
    for d in sorted(matching_dirs):
        cmake_file = d / "CMakeLists.txt"
        if cmake_file.exists():
            tests_found = extract_test_names(cmake_file)
            if tests_found and args.verbose:
                rel_path = os.path.relpath(cmake_file, demo_drivers_root)
                print(f"Found tests in {rel_path}: {tests_found}")
            all_test_names += tests_found

    if all_test_names:
        if args.verbose:
            print(f"Found {len(all_test_names)} matching tests.")
    else:
        print("No tests found in matched directories.")
        sys.exit(0)

    # 3) Generate a single 'ctest -R' command
    # Use a combined OR-regex: ^(testA|testB|...)$
    # Also escape special regex characters in test names
    escaped_names = [re.escape(name) for name in all_test_names]
    or_pattern = "|".join(escaped_names)
    full_pattern = f"^({or_pattern})$"
    run_command = f"ctest -R '{full_pattern}'"

    # Append any extra arguments to the ctest command
    if args.ctest_args:
        extra = " ".join(args.ctest_args)
        run_command += " " + extra

    # [OPTIONAL] 4) Run tests from specified build directory
    if args.run_from:
        try:
            if args.verbose:
                print(f"\nRunning command:\n\t{run_command}\n")
            run(run_command, cwd=args.run_from, shell=True)
        except Exception as _:
            print(f"[ERROR] Failed to run command: {run_command}", file=sys.stderr)
            sys.exit(1)
    else:
        print("\nTo run these tests, run the following command from your demo_drivers/build/ directory:")
        print(f"\n\t{run_command}\n")


if __name__ == "__main__":
    main()
