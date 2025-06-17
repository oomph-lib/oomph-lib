#!/usr/bin/env python3

import sys
import gzip
import math
import re
from typing import List


def get_type(a: str) -> int:
    """ Distinguish between a number and a string:

        Returns integer 1 if the argument is a number,
                        2 if the argument is a string.
    """
    pattern = re.compile(r"^[+-]?(?:\d*\.\d+|\d+\.\d*|\d+)(?:[EeDd][+-]?\d+)?$")
    return 1 if pattern.match(a) else 2


def append_symbols(string: str, symbol: str, count: int) -> str:
    """Append a given number of symbols to the string followed by a space."""
    return string + (symbol * count) + " "


def read_file(filename: str) -> List[str]:
    """Read file into a list of strings, supporting gzip files."""
    if filename.endswith(".gz"):
        with gzip.open(filename, "rt") as f:
            return f.readlines()
    with open(filename, "r", encoding="utf-8") as f:
        return f.readlines()


def fpdiff_helper(filename1, filename2, relative_error, small, outstream, details_stream):
    """ Calculate the floating-point difference between two data files.
        The idea is to use a looser tolerance than the UNIX diff command,
        so that if two entries have a relative error less than the argument
        relative_error, they are counted as the same.

        Note that the relative error is percentage!

        Information on pass/failure is written to outstream. Details on which
        lines failed are written to details_stream. Warning: if run on
        large files the details_stream may be overwhelmingly long.

        First return value: 0 if the two files are the same, 1 if they are
        different or 5 if the files cannot be opened.

        Second return value: the maximum relative error.

        Third return value: the largest entry that caused an error
        (i.e. what "small" would need to be set as for there to be no
        differences).
    """
    # Load the files, the rearrange so the first file has more lines than the second
    (file1, file2) = (read_file(filename1), read_file(filename2))
    (n1, n2) = (len(file1), len(file2))
    (file1, file2, n) = (file1, file2, n2) if (n1 >= n2) else (file2, file1, n1)

    # Storage for worst case error sizes, number of errors and number of lines with errors
    (max_rel_diff, max_wrong_entry, nerr, nline_error) = (0, 0, 0, 0)

    # Loop over the lines in file1 (the file with the most lines!)
    for (count, line1) in enumerate(file1):
        # If we've run over the end of the file2, issue a warning and end the loop
        if count >= n:
            details_stream.write("\nWarning: files have different numbers of lines")
            details_stream.write(f"\nResults are for first {count} lines of both files\n")
            nerr += 1
            break

        # Read the next line from file2. If the lines are the same, we're done
        line2 = file2[count]
        if line1 == line2:
            continue

        # Split each line into its separate fields then find the number of fields in each
        # line. If the number of fields is not the same, report it as an error
        (fields1, fields2) = (line1.split(), line2.split())
        (nfields1, nfields2) = (len(fields1), len(fields2))
        if nfields1 != nfields2:
            details_stream.write(f"\n =====> line {count+1}: different number of fields\n")
            details_stream.write(f"{nfields1} fields: {line1}")
            details_stream.write(f"{nfields2} fields: {line2}")
            nerr += 1
            continue

        # Otherwise, we now compare field by field

        # Flag to indicate whether there has been a problem in the field
        problem = 0
        # Strings that will hold the output data
        (outputline1, outputline2, outputline3) = ("", "", "")

        # Loop over the fields
        for (f1, f2) in zip(fields1, fields2):
            max_len = max(len(f1), len(f2))
            outputline1 += f1.ljust(max_len) + " "
            outputline3 += f2.ljust(max_len) + " "

            # If the fields are identical, we are fine
            if f1 == f2:
                # Put spaces into the error line
                outputline2 = append_symbols(outputline2, " ", max_len)
                continue

            # Find the type (numeric or string) of each field
            (type1, type2) = (get_type(f1), get_type(f2))

            # If the data-types aren't the same issue an error
            if type1 != type2:
                problem = 1
                nerr += 1
                outputline2 = append_symbols(outputline2, "*", max_len)
                continue

            # If the types are both strings then report the error
            if type1 == 2:
                problem = 1
                nerr += 1
                outputline2 = append_symbols(outputline2, "%", max_len)
                continue

            # Convert strings to floating point number
            x1 = float(f1.lower().replace("d", "e"))
            x2 = float(f2.lower().replace("d", "e"))

            # If both numbers are very small, that's fine
            if (math.fabs(x1) <= small) and (math.fabs(x2) <= small):
                # Put spaces into the error line
                outputline2 = append_symbols(outputline2, " ", max_len)
                continue

            # Find the relative difference based on the largest number
            # Note that this "minimises" the relative error (in some sense)
            # but means that I don't have to separately trap the cases
            # when x1, x2 are zero
            diff = 100.0 * (math.fabs(x1 - x2) / max(math.fabs(x1), math.fabs(x2)))

            # If the relative error is smaller than the tolerance, that's fine
            if diff <= relative_error:
                outputline2 = append_symbols(outputline2, " ", max_len)
                continue

            # Otherwise issue an error
            problem = 1
            nerr += 1
            outputline2 = append_symbols(outputline2, "-", max_len)

            # Track worst case values
            max_rel_diff = max(diff, max_rel_diff)
            max_wrong_entry = max(max_wrong_entry, max(math.fabs(x1), math.fabs(x2)))

        # If there has been any sort of error, print it
        if problem == 1:
            nline_error += 1
            details_stream.write(f"\n =====> line {count+1}\n")
            details_stream.write(f"{outputline1}\n")
            details_stream.write(f"{outputline2}\n")
            details_stream.write(f"{outputline3}\n")

    if nerr > 0:
        outstream.write(f"\n In files {filename1} {filename2}")
        outstream.write(f"\n number of lines processed: {count}")
        outstream.write(f"\n number of lines containing errors: {nline_error}")
        outstream.write(f"\n number of errors: {nerr} ")
        outstream.write(f"\n largest relative error: {max_rel_diff} ")
        outstream.write(f"\n largest abs value of an entry which caused an error: {max_wrong_entry}")
        outstream.write("\n========================================================")
        outstream.write("\n    Parameters used:")
        outstream.write(f"\n        threshold for numerical zero : {small}")
        outstream.write(f"\n        maximum rel. difference [percent] : {relative_error}")
        outstream.write("\n    Legend: ")
        outstream.write("\n        *******  means differences in data type (string vs number)")
        outstream.write("\n        -------  means real data exceeded the relative difference maximum")
        outstream.write("\n        %%%%%%%  means that two strings are different")
        outstream.write("\n========================================================")
        outstream.write("\n\n   [FAILED]\n")
        return (2, max_rel_diff, max_wrong_entry)
    else:
        outstream.write(f"\n\n In files {filename1} {filename2}")
        outstream.write(f"\n   [OK] for fpdiff.py parameters: - max. rel. error = {relative_error} ")
        outstream.write(f"\n                                  - numerical zero  = {small}\n")
        return (0, max_rel_diff, max_wrong_entry)


def fpdiff(filename1, filename2, relative_error=0.1, small=1e-14, outstream=sys.stdout, details_stream=sys.stdout):
    """Wrapper for using fpdiff inside python. Has default args and returns a
    bool."""
    results = fpdiff_helper(filename1, filename2, relative_error, small, outstream, details_stream)
    (ok, max_rel_diff, max_wrong_entry) = results
    return ((ok == 0), max_rel_diff, max_wrong_entry)


def run_as_script(argv):
    """Run fpdiff as a script (handles argument parsing, output as error codes
       and some helpful messages).
    """
    # Note that we shouldn't just put this code this under 'if __name__ ==
    # "__main__":' because variables created there are global. This resulted
    # in some bugs before.

    # Set the defaults
    maxreld = 1.0e-1  # max relative difference in percent
    small = 1.0e-14  # small number -- essentially round-off error

    # Remove the program name from the front of the argument list
    argv.pop(0)

    # Let's find the number of command line arguments
    narg = len(argv)

    # If we're out of range, issue a usage message
    if narg < 2 or narg > 4:
        sys.stdout.write("\n      *********   ERROR   **********\n")
        sys.stdout.write("\nMust specify 2, 3 or 4 keywords on the command line. ")
        sys.stdout.write(f"\nYou have specified {narg}")
        sys.stdout.write("\n   Proper usage:  ")
        sys.stdout.write("\n         fpdiff file1 file2 [max_rel_diff_percent] [small]\n")
        sys.stdout.write("\n      *********  PROGRAM TERMINATING   ***********")
        sys.stdout.write("\n   [FAILED] \n")
        sys.exit(4)

    # Read any optional arguments
    if narg >= 3:
        maxreld = float(argv[2])
        if narg == 4:
            small = float(argv[3])

    # Run the diff and if there is an IO error then fail with a useful message
    try:
        (error_code, _, _) = fpdiff_helper(argv[0], argv[1], maxreld, small, sys.stdout, sys.stdout)
    except IOError as err:
        sys.stdout.write(f"\n   [FAILED] I/O error({err.errno}): {err.strerror} \"{err.filename}\"\n")
        return 5
    return error_code


# What to do if this is run as a script, rather than loaded as a module
if __name__ == "__main__":
    # Run and return whether it succeeded or not
    sys.exit(run_as_script(sys.argv))
