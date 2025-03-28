#!/usr/bin/env python3

from pathlib import Path
from argparse import ArgumentParser, Namespace
import re
import sys
import math
from typing import Set, List


def find_matching_directories(root_dir: Path, keywords: List[str]) -> Set[Path]:
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

    for (i, cpp_file) in enumerate(all_files, start=1):
        progress_fraction = i / total_files
        progress_percent = math.floor(progress_fraction * 100)
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
        re.DOTALL
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
    parser = ArgumentParser(description="Find OOMPH tests in demo_drivers subdirectories by searching for one or more keywords in C++ files.")
    parser.add_argument("--root", required=True, help="Path to the 'demo_drivers' directory")
    parser.add_argument("keywords", nargs="+", help="One or more keywords to search for in C++ source files (OR match).")
    # fmt: on
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    demo_drivers_root = Path(args.root).resolve()
    if not demo_drivers_root.exists():
        raise FileNotFoundError(f"Error: The specified directory '{demo_drivers_root}' does not exist.")

    # Show which keywords we're searching for
    keyword_list_str = ", ".join(args.keywords)
    print(f"Searching for keywords: {keyword_list_str}")
    print(f"Within directory: {demo_drivers_root}\n")

    # 1) Identify all directories that contain a file with the given keywords
    matching_dirs = find_matching_directories(demo_drivers_root, args.keywords)
    if not matching_dirs:
        print("No directories found containing any of the specified keywords.")
        sys.exit(0)

    # 2) From those directories, extract any test names from the local CMakeLists.txt
    all_test_names: List[str] = []
    for d in sorted(matching_dirs):
        cmake_file = d / "CMakeLists.txt"
        if cmake_file.exists():
            tests_found = extract_test_names(cmake_file)
            if tests_found:
                rel_path = cmake_file.relative_to(demo_drivers_root)
                print(f"Found tests in {rel_path}: {tests_found}")
            all_test_names += tests_found

    if not all_test_names:
        print("No tests found in matched directories.")
        sys.exit(0)

    # 3) Generate a single 'ctest -R' command
    # Use a combined OR-regex: ^(testA|testB|...)$
    # Also escape special regex characters in test names
    escaped_names = [re.escape(name) for name in all_test_names]
    or_pattern = "|".join(escaped_names)
    full_pattern = f"^({or_pattern})$"

    print("\nTo run these tests, use the following command in your build/ directory:\n")
    print(f"\tctest -R '{full_pattern}'\n")


if __name__ == "__main__":
    main()