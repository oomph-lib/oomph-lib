#!/usr/bin/env python3

import argparse
import glob
import os
import re
import sys
from pathlib import Path
from typing import List


def load_file(fname: Path) -> str:
    """Reads in the data stored in the input file.

    :param fname: the path to a file.
    :type fname: Path
    :return: the contents of 'fname'.
    :rtype: str
    """
    contents = None
    with open(fname, "r") as f:
        contents = f.read()
    return contents


def process_doxygen_formula_block(text: str) -> str:
    """Removes newlines and comments from the passed Doxygen formula block.

    :param text: a Doxygen formula block, possibly broken over several lines.
    :type text: str
    :return: the processed formula text, as a single line of text (no newlines).
    :rtype: str
    """
    # Remove triple forward slashes first
    processed_text = text.replace("///", "")

    # Remove any double forward slashes
    processed_text = processed_text.replace("//", "")

    # Replace newlines with spaces
    processed_text = processed_text.replace("\n", " ")

    # Squeeze whitespaces
    processed_text = re.sub(" +", " ", processed_text)

    return processed_text


def find_number_of_leading_whitespaces(text: str, index: int) -> int:
    """Given the index of a position in a text, it finds the corresponding line
    in the text and the number of leading whitespaces on that line.

    :param text: a block of text.
    :type text: str
    :param index: the index in 'text' to determine the number of leading
                  whitespaces from.
    :type index: int
    :return: the number of leading whitespaces.
    :rtype: int
    """
    assert index < len(text)
    line_number = text[:index].count("\n")
    line = text.splitlines()[line_number]
    n_leading_whitespace = len(line) - len(line.lstrip(" "))
    return n_leading_whitespace


def process_file(fname: Path, pattern: str, add_newline: bool, in_place: bool = False) -> None:
    """Performs regex magic to fix broken Doxygen formulae.

    :param fname: the path to the file to fix.
    :type fname: Path
    :param add_newline: flag indicating whether new lines should be added around
                        the processed text, e.g. inline code might not but a
                        formula block might.
    :type add_newline: bool
    :param in_place: [description], defaults to False
    :type in_place: bool, optional
    """

    # Expand the "~/" in the filename and read the file in
    text = load_file(fname)

    # Search for all matching patterns
    regex_iterator = re.finditer(pattern, text)

    # Initialise update text
    updated_text = ""
    current_position = 0

    # If we found matches
    if regex_iterator:
        for match in regex_iterator:
            # Get the (start, end) indices of the match in the original text
            (posn_start, posn_end) = match.span()

            # Add on any text that lies between the end of the previous match
            # and the start of the current match
            updated_text += text[current_position:posn_start]

            # Find out the number of leading whitespaces on the line containing
            # the formula. We can use this to decide how many whitespaces to add
            # if we add any newlines before or after this block has been
            # formatted
            n_leading_whitespace = find_number_of_leading_whitespaces(
                text, index=posn_start)

            # The text required to introduce a new comment line with the correct
            # amount of leading whitespaces
            new_comment_line = "\n" + (" " * n_leading_whitespace) + "///"

            # Extract the match as a string and process it into the desired form
            processed_match = process_doxygen_formula_block(match.group())

            # Make sure the formula starts on its own line if it doesn't already
            if add_newline and (updated_text.rstrip(" ")[-2:] != "//"):
                updated_text += new_comment_line + " "

            # Place the processed formula on a new line
            updated_text += processed_match

            # Make sure there is no text after the formula on the same line
            if add_newline and (not text[posn_end:].lstrip(" ").startswith("\n")):
                updated_text += new_comment_line

            # Continue on from the end of this string
            current_position = posn_end

        # Add on the remaining text from the file
        updated_text += text[current_position:]
    else:
        updated_text = text

    # Replace the original file, if requested. Otherwise, print out the changes
    if in_place:
        with open(os.path.expanduser(fname), "w") as f:
            f.write(updated_text)
    else:
        print("\n===========================================================\n")
        print("File: {}\n".format(fname))
        print(updated_text)
        print("\n===========================================================")


def has_broken_formula(fname: Path, pattern: str) -> bool:
    """Identifies whether the file contains a broken Doxygen formula, i.e. the
    formula is split over more than one line.


    :param fname: the absolute path to a file.
    :type fname: Path
    :return: True if the file is "broken" and False otherwise.
    :rtype: bool
    """
    is_broken = False
    regex_iterator = re.finditer(pattern, load_file(fname))
    if regex_iterator:
        for match in regex_iterator:
            if "\n" in match.group():
                is_broken = True
                break
    return is_broken


def find_all_cpp_files(base_dir: Path) -> List[Path]:
    """Finds all the headers and source files in this directory and below.

    :param base_dir: the path to the directory to search in and below.
    :type base_dir: str
    :return: an array of paths to header/source files.
    :rtype: List[str]
    """
    headers = glob.glob(os.path.join(base_dir, "**/*.h"), recursive=True)
    sources = glob.glob(os.path.join(base_dir, "**/*.cc"), recursive=True)
    headers = [Path(h) for h in headers]
    sources = [Path(s) for s in sources]
    return headers + sources


def run(base_dir: str, in_place: bool = False, just_list_files: bool = False) -> None:
    """Fixes all of the broken Doxygen formulae in the C++ source files of a
    specific directory.

    :param base_dir: the directory in and below which to consider.
    :type base_dir: str
    :param in_place: edits files in-place, defaults to False
    :type in_place: bool, optional
    :param just_list_files: [description], defaults to False
    :type just_list_files: bool, optional
    :return: [description]
    :rtype: [type]
    """
    # A tuple of tuples, with each tuple containing (i) the formula to fix,
    # (ii) the regex pattern to find the formula, and (iii) a flag indicating
    # whether to add surrounding newlines around the formatted formula.
    #                           !!! DO NOT EDIT !!!
    doxygen_formulae_regex_patterns = (
        (r"\f[ ... \f]", r"(\\f\[)([\S\s]*?)(\\f\])", True),
        (r"\f$ ... \f$", r"(\\f\$)([\S\s]*?)(\\f\$)", False)
    )

    # Find all header/source files to process
    cpp_files = find_all_cpp_files(base_dir)

    # Iterate over the different types of Doxygen formulae patterns
    for (formula_type, pattern, add_newline) in doxygen_formulae_regex_patterns:
        def file_has_broken_formula(fname: Path):
            return has_broken_formula(fname, pattern)

        # Find the files in which the pattern is broken over one or more lines
        broken_files = list(filter(file_has_broken_formula, cpp_files))

        # Tell the user which pattern we're fixing
        print("\nFixing formula: {}".format(formula_type))

        if len(broken_files) == 0:
            print("\n   ...but there are no broken files to process!")
            continue

        # Spit out which files are broken
        print("\nBroken files:")
        for file in broken_files:
            print("    * {}".format(file.relative_to(Path.cwd())))
        print()

        # Now process each broken file
        if not just_list_files:
            print("\nProcessing file:")
            for file in broken_files:
                print("    * {}".format(file.relative_to(Path.cwd())))
                process_file(file, pattern=pattern, in_place=in_place,
                             add_newline=add_newline)


if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Fixes Doxygen formulae.')
    parser.add_argument("-d", "--dir", type=lambda p: Path(p).absolute(),
                        default=Path.cwd(), help="The directory to search in and below.")
    parser.add_argument("-i", "--in-place", action='store_true',
                        default=False, help="Edit the files in-place.")
    parser.add_argument("-s", "--silent", action='store_true',
                        default=False, help="Silence the output.")
    parser.add_argument("-l", "--list-files", action='store_true',
                        default=False, help="Just list the broken files.")
    args = parser.parse_args()

    # Extract the commandline arguments
    (base_dir, in_place, silent, list_files) = (
        args.dir, args.in_place, args.silent, args.list_files)

    # Make sure the base directory is a valid directory
    if not base_dir.exists():
        raise ValueError("Requested base directory does not exist!")

    # Write the standard output to /dev/null if silence is requested
    if silent:
        sys.stdout = open(os.devnull, "w")

    # Call the function that does the heavy lifting
    run(base_dir=base_dir, in_place=in_place, just_list_files=list_files)

    # Reset the standard output
    if silent:
        sys.stdout = sys.__stdout__
