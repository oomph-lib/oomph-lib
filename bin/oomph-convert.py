#!/usr/bin/env python
#
# Copyright (C) 2007- Angelo Simone (TU Delft)
# Licensed under the GNU LGPL Version 2.1
#
# Script for converting from oomph-lib tecplot format to VTK XML
#
#
# Created on 20070905
# Last revision by AS on 20080514
# 
# Report bugs to a.simone@tudelft.nl
#
# The script can deal with 
# - empty lines
# - lines that start with #
# - 2d meshes (plot only 2d elements and skip 1d boundary)
#
# Scan the script for "Assumption" and "Note"
#
#############################################
#

import getopt
import sys
from commands import getoutput

def main(argv):
    "Main function"
    print "* oomph-convert.py, ver. 20080514"

    # Check python version. Bark at user if < 2.3. 
    if sys.version_info<(2,3): 
        print >>sys.stderr, "You need at least Python 2.3 " 
        sys.exit(3) 

    # Get command-line arguments
    try:
        opts, args = getopt.getopt(argv, "hz")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
        
    # Get options
    flag = 0
    for opt, arg in opts:
        if opt in ("-h"):
            usage()
            sys.exit()
        elif opt in ("-z"):
            flag = 1
        
    if len(args) == 1:
        # Get filename and suffix
        ifilename = args[0]
        isuffix = ifilename.split(".")[-1]
        # Set filename and suffix
        extension_len = len(isuffix) + 1
        lenBaseName = len(ifilename) - extension_len
        ofilename = ifilename[:lenBaseName]+".vtu"
        osuffix = ofilename.split(".")[-1]      
    elif len(args) == 2:
        # Get filenames and suffixes
        ifilename = args[0]
        ofilename = args[1]
        isuffix = ifilename.split(".")[-1]
        osuffix = ofilename.split(".")[-1]
    else:
        usage()
        sys.exit(2)

    if isuffix == "dat" and osuffix == "vtu":
        # Convert from oomph-lib Tecplot format to VTK XML format
        tecplot_to_vtkxml(ifilename, ofilename, flag)
    else:   
        error("Sorry, cannot convert between .%s and .%s file formats." % (isuffix, osuffix))

#
#############################################
#

def usage():
    "Display usage"
    print """\

NAME
        oomph-convert.py - script for converting from oomph-lib tecplot format to VTK XML


SYNOPSIS
        oomph-convert.py [OPTION] [input_file.dat] [output_file.vtu]


OPTIONS
        -h         display this help text and exit
        -z         add trailing zeros to the output filename 


TYPICAL USAGE EXAMPLES
        oomph-convert.py -h                             -> display help text  
        oomph-convert.py soln12.dat soln12.vtu          -> generate soln12.vtu
        oomph-convert.py -z soln12.dat soln12.vtu       -> generate soln00012.vtu
        oomph-convert.py soln12.dat                     -> generate soln12.vtu
        oomph-convert.py -z soln12.dat                  -> generate soln00012.vtu
        oomph-convert.py soln12.dat nsol2.vtu           -> generate nsol2.vtu
        oomph-convert.py -z soln12.dat nsol2.vtu        -> generate nsol00002.vtu
        oomph-convert.py soln.dat                       -> generate soln.vtu
        oomph-convert.py -z soln.dat                    -> generate soln00000.vtu    
"""

#
#############################################
#

def error(message):
    "Write an error message"
    for line in message.split("\n"):
        print "*** %s" % line
    sys.exit(2)

#
#############################################
#

def tecplot_to_vtkxml(ifilename, ofilename, flag):
    "Convert from oomph-lib Tecplot format to VTK XML format."

    print "Convert from oomph-lib Tecplot format to VTK XML format."
    
    # -- Open files

    try:
        ifile = open(ifilename, "r")
    except:
        error("input file does not exist")

    ofile = open(ofilename, "w")

    # -- Scan file for known offending keywords
    while 1:
        line = ifile.readline()
        if not line: break
        line = line.strip()
        if "TEXT" in line:
	    print "\n"
	    print line
	    print "\n"
	    error("Cannot process the Tecplot file\n" 
	    "The offending line is above\n"
	    "To continue remove the data related to TEXT")

    # -- Define leading dimension

    ifile.seek(0)
    dim3 = 0
    dim2 = 0
    dim1 = 0
    while 1:
        line = ifile.readline()
        if not line: break
        line = line.strip()
        if "ZONE" in line:
            if "K=" in line:
                dim3 = 1
            if "J=" in line:
                dim2 = 1
            if "I=" in line:
                dim1 = 1
        
        if dim3 == 1:
            dim = 3
        elif dim3 == 0 and dim2 == 1:
            dim = 2
        elif dim3 == 0 and dim2 == 0 and dim1 == 1:
            error("You are better off using Gnuplot for a 1D data set")
            
    print "Dimension of the problem:",dim

    # -- Define the number of ZONE

    ifile.seek(0)
    num_zone = 0
    while 1:
        line = ifile.readline()
        if not line: break
        line = line.strip()
        if "ZONE" in line:
           num_zone += 1

    # Check that we got the number of ZONE
    if num_zone == 0:
        error("Unable to find the number of ZONE.")
        
    #print "Number of ZONE:",num_zone   

    # -- Define number of points and cells

    # Step to beginning of file
    ifile.seek(0)

    number_of_points = 0
    number_of_plot_cells = 0

    for Q in range(1, num_zone+1):
    
        (line) = grab_line(ifile)
        (dimI,dimJ,dimK) = define_zone_dimensions(line)

        if dim == 2:

            Zone_number_of_points=dimI*dimJ
            zone_number_of_plot_cells=(dimI-1)*(dimJ-1)

            for i in range(dimI*dimJ):
                ifile.readline()

            number_of_points = number_of_points + dimI*dimJ
            number_of_plot_cells = (dimI-1)*(dimJ-1) + number_of_plot_cells

        elif dim == 3:
            error("3D not coded yet") 

    print "Plot cells defined"

    # -- Define the number of field values

    #    Assumption: 
    #    - the number of field values is the same for all zones 
    #   (there is a check on the first zone only)
    #    and
    #    - the first zone is _not_ a boundary but an element

    ifile.seek(0)

    line = grab_line(ifile)
    (dimI,dimJ,dimK) = define_zone_dimensions(line)

    line = grab_line(ifile)

    entries = line.split()
    field_variables = 0 

    for token in entries:
        field_variables += 1

    field_variables = field_variables - dim

    print "Field variables = ",field_variables

    # -- Define coordinates
    
    write_header(ofile)
    write_header_coordinates(ofile,number_of_points,number_of_plot_cells)

    ifile.seek(0)

    for Q in range(1, num_zone+1):

        line = grab_line(ifile)
        (dimI,dimJ,dimK) = define_zone_dimensions(line)

        for i in range(dimI*dimJ):
        
            line = grab_line(ifile)

            if dim == 2:
                (x,y) = line.split()[:2]        
                ofile.write(""" %f %f %f \n""" %(float(x), float(y), 0.0))
            elif dim == 3:
                (x,y,z) = line.split()[:3]    
                ofile.write(""" %f %f %f \n""" %(float(x), float(y), float(z)))

    write_footer_coordinates(ofile)

    # -- Define connectivity list

    write_header_connectivities(ofile)

    ifile.seek(0)
    
    first_node_zone = 0
    
    for Q in range(1, num_zone+1):

        line = grab_line(ifile)
        (dimI,dimJ,dimK) = define_zone_dimensions(line)

        # Read coordinates
        for i in range(dimI*dimJ):
            ifile.readline()

        # Write connectivity table
        k = 0
        a = 0
        b = 0
        c = 0 
        d = 0

        for i in range(1, dimJ): # range from 1 to dimJ-1
            for j in range(1, dimI): # range from 1 to dimI-1
                k = k + 1
                a = k
                b = a + 1
                c = b + dimI
                d = c - 1
                ofile.write(""" %d %d %d %d \n""" %(first_node_zone+a-1,first_node_zone+b-1,first_node_zone+c-1,first_node_zone+d-1))

            k = a + 1
        
        first_node_zone = first_node_zone + (c-1) + 1 # This is the first node of a Tecplot ZONE

    write_footer_connectivities(ofile)

    # -- Define element offset

    write_header_offset(ofile)

    ifile.seek(0)

    first_cell_zone = 0
    for Q in range(1, num_zone+1):

        line = grab_line(ifile)
        (dimI,dimJ,dimK) = define_zone_dimensions(line)

        zone_number_of_plot_cells=(dimI-1)*(dimJ-1)
        
        # Read coordinates
        for i in range(dimI*dimJ):
            ifile.readline()

        if dim == 2:
            for i in range(1, zone_number_of_plot_cells+1): # range from 1 to zone_number_of_plot_cells
                print >>ofile, first_cell_zone+i*4
        elif dim == 3:
          error("3D Offset not coded yet.")

        first_cell_zone = first_cell_zone + i*4 

    write_footer_offset(ofile)

    # -- Dump element type

    write_header_types(ofile)

    ifile.seek(0)

    for Q in range(1, num_zone+1):

        line = grab_line(ifile)
        (dimI,dimJ,dimK) = define_zone_dimensions(line)

        zone_number_of_plot_cells=(dimI-1)*(dimJ-1)
        
        # Read coordinates
        for i in range(dimI*dimJ):
            ifile.readline()

        if dim == 2:
            for i in range(1, zone_number_of_plot_cells+1): # range from 1 to zone_number_of_plot_cells
                element_type = 9 # Q4
                print >>ofile, element_type
        elif dim == 3:
            error("3D element type not coded yet.")

    write_footer_types(ofile)

    # -- Dump field value

    write_header_point(ofile)

    for P in range(1,field_variables+1):

        write_header_field(ofile,P)

        ifile.seek(0)

        for Q in range(1, num_zone+1):

            line = grab_line(ifile)
            (dimI,dimJ,dimK) = define_zone_dimensions(line)
        
            for i in range(dimI*dimJ):

                line = grab_line(ifile)
                (V) = line.split()[dim+P-1]
                ofile.write(""" %f \n""" %(float(V)))

        write_footer_field(ofile)

    write_footer_point(ofile)

    write_footer(ofile)

    # -- Close files
    
    ifile.close();
    ofile.close();
    if flag == 0:
        print '* Output file name: %(fn)s ' %{'fn': ofilename}

    if flag == 1:
        # -- Rename output file (adding trailing zeros)

        # 1) Define the number of digits in file name (up to limit_digits)
        
	# Assumption: 
	# - the last four characters in ofilename  
        # define the extension of ofilename(".vtu")
        
	limit_digits = 11
        extension_len = 4
        digits = 0
        for i in range(1, limit_digits+1):
            try:
                dummy = int(ofilename[-extension_len-i:-extension_len])
                digits = i
            except:
                pass

        # 2) Rename the file
        # - digits in ofilename: ofilename[-extension_len-digits:-extension_len]
        # - extension in ofilename: ofilename[-extension_len:]
        # - base file name in ofilename: ofilename[:lenBaseName]
	
        # Note: only filenames with 4 or less digits will have a trailing zero if 
        # the second digits in "%05i.vtu" is 5 in os.rename
        
	import os
        lenBaseName = len(ofilename) - digits - extension_len
        if digits == 0:
            cifer = 0
        else:
            cifer = int(ofilename[-extension_len-digits:-extension_len])
        newofilename = ofilename[:lenBaseName]+"%05i.vtu" % cifer
        os.rename(ofilename, newofilename)
        print '* Output file name: %(fn)s ' %{'fn': newofilename}

#
#############################################
#

# Write field footer
def write_footer_field(ofile):
    ofile.write("""\
</DataArray> 
""")

#
#############################################
#

# Write field header
def write_header_field(ofile,P):
    ofile.write("""\
<DataArray  type="Float32"  Name="V%d"  format="ascii">  
""" %(P))

#
#############################################
#

# Write file header
def write_header(ofile):
    ofile.write("""\
<?xml version=\"1.0\"?>
<VTKFile type="UnstructuredGrid" version="0.1"  byte_order="LittleEndian" compressor="vtkZLibDataCompressor">
<UnstructuredGrid>
""" )
    print "Conversion started"

#
#############################################
#

# Write connectivity header
def write_header_connectivities(ofile):
    ofile.write("""\
<Cells>
<DataArray type="Int32" Name="connectivity" format="ascii">
""")

#
#############################################
#

# Write connectivity footer
def write_footer_connectivities(ofile):
    ofile.write("""\
    </DataArray>
""")
    print "Connectivities defined"

#
#############################################
#

# Write offset header
def write_header_offset(ofile):
    ofile.write("""\
    <DataArray type="Int32" Name="offsets" format="ascii">
    """)

#
#############################################
#

# Write offset footer
def write_footer_offset(ofile):
    ofile.write("""\
</DataArray>
""")
    print "Offset defined"

#
#############################################
#

# Write type header
def write_header_types(ofile):
    ofile.write("""\
<DataArray type="UInt8" Name="types" format="ascii">
    """)

#
#############################################
#

# Write point data header
def write_header_point(ofile):
    ofile.write("""\
<PointData>
    """)

#
#############################################
#

# Write point data footer
def write_footer_point(ofile):
    ofile.write("""\
</PointData>
    """)
    print "Field data defined"

#
#############################################
#

# Write type footer
def write_footer_types(ofile):
    ofile.write("""\
</DataArray>
</Cells>
    """)
    print "Element types defined"

#
#############################################
#

# Write file footer
def write_footer(ofile):
    ofile.write("""\
</Piece>
</UnstructuredGrid>
</VTKFile>
    """)
    print "Conversion done"

#
#############################################
#

# Write coordinate header
def write_header_coordinates(ofile,Points,Cells):
    ofile.write("""\
<Piece NumberOfPoints=" %d "  NumberOfCells=" %d ">
<Points>
<DataArray type="Float32" NumberOfComponents="3" format="ascii">
    """ % (Points,Cells))

#
#############################################
#

# Write coordinate footer
def write_footer_coordinates(ofile):
    ofile.write("""\
</DataArray>
</Points>
    """)
    print "Coordinate defined"

#
#############################################
#

# Define zone dimensions 
def define_zone_dimensions(line):
    dimI = 1
    dimJ = 1
    dimK = 1
    if "K=" in line:
        (zone, dimI, dimJ, dimK) = line.split()[:4]
        dimI = int(dimI.replace("I=", ""))
        dimJ = int(dimJ.replace("J=", ""))
        dimK = int(dimK.replace("K=", ""))
    elif "J=" in line:
        (zone, dimI, dimJ) = line.split()[:3]
        dimI = int(dimI.replace("I=", ""))
        dimJ = int(dimJ.replace("J=", ""))
    elif "I=" in line:
        (zone, dimI) = line.split()[:2]
        dimI = int(dimI.replace("I=", ""))

    return dimI, dimJ, dimK

#
#############################################
#

# Grab line skipping empty lines and comments 
def grab_line(ifile):
    line = ifile.readline()
    line = line.strip()
    line = line.replace(",", "")

    # skip empty line (= "\n") and comment line (beginning with #)
    if not line: 
        line = ifile.readline()
        line = line.strip()
        line = line.replace(",", "")
    elif line[0] == '#':    
        line = ifile.readline()
        line = line.strip()
        line = line.replace(",", "")

    return line

#
#############################################
#

if __name__ == "__main__":
    main(sys.argv[1:])

