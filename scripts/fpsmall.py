#!/usr/bin/env python


#!/usr/bin/env python

import sys

usage = """Test whether all floating point values in a file are small.

fpsmall [filename] [maximum allowed value]
"""

def pr(*args):
    """Compatible print-like function."""
    sys.stdout.write(" ".join(args))
    sys.stdout.write("\n")

def main(argv):

    if len(argv) < 2 or len(argv) > 3:
        pr(usage)
        return 1

    filename = argv[1]

    if len(argv) > 2:
        small = float(argv[2])
    else:
        small = 1e-3

    f = open(filename)
    lines = f.readlines()
    f.close()
    
    entries = sum([l.strip().split() for l in lines], []) 
    floats = [float(e) for e in entries]

    exit_code = 0
    for fl in floats:
        if fl > small:
            pr("Value too large: ", str(fl))
            exit_code = 2
    
    return exit_code


if __name__ == "__main__":
    sys.exit(main(sys.argv))
