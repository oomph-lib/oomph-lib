#!/usr/bin/env python3
#
# Copyright (C) 2007- Angelo Simone (TU Delft)
# Licensed under the GNU LGPL Version 2.1
#
# Script for converting from oomph-lib tecplot format to VTK XML
# 
# Report bugs to a.simone@tudelft.nl
#
# The script can convert
# - 2d and 3d quad and tetra meshes
# - 2d and 3d points from all meshes
#
# The script can NOT convert
# - 1d meshes
#
# Scan the script for "Assumption" and "Note"
#
###############################################################################
# CHANGELOG
#
# 20070905: Angelo Simone
#           - Creation
#
# 20110521: Alexandre Raczynski (alexandre.raczynski@technogerma.fr)
#           - Parsing and conversion algorithm rewritten
#           - 3d quad meshes support added
#           - Trailing zeros option improved (do not rename the file afterward)
#
# 20110610: Angelo Simone
#           - Added support for triangle and tetrahedral meshes
#
# 20110615: Jeremy Chu Van
#           - Added point-only extraction (if there are zones that cannot
#             be handled or an invalid file)
#
# 20130419: David Shepherd
#           - Send errors to stderr not stdout.
#
# 20120714: Matthias Heil
#           - skip files if output files already exist (unless -o is specified
#             on command line). Avoids costly (!) re-generation of existing
#             files.
#           - Fixed vtp file generation with -p flag: Ignore any line
#             the contains non-floating point numbers (rather than
#             ignoring any lines that contains an ascii character which
#             ignores floating point numbers in scientific notation!
#
# 20220610  Andrew Hazel
#           - Converted to python3
#
###############################################################################
import os
import getopt
import string
import sys
import time
import re
from subprocess import getoutput

#
###############################################################################
#

def usage():
    """ Display usage
    """
    print("""\

NAME
        oomph-convert.py - script for converting from oomph-lib tecplot format to VTK XML. 


SYNOPSIS
        oomph-convert.py [OPTION] [input_file.dat] [output_file]


OPTIONS
        -p2 or -p3 outputs only points in 2D or 3D (.vtp)
        -h         display this help text and exit
        -z         add trailing zeros to the output filename
        -o         force overwrite existing files

        By default output files will overwrite old ones only if the input file is
        newer than the output file.


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
        oomph-convert.py -p3 soln12.dat soln12.vtp          -> generate soln12.vtp
        oomph-convert.py -p2 -z soln12.dat soln12.vtp       -> generate soln00012.vtp
        oomph-convert.py -p3 soln12.dat                     -> generate soln12.vtp
        oomph-convert.py -p2 -z soln12.dat                  -> generate soln00012.vtp
        oomph-convert.py -p3 soln12.dat nsol2.vtp           -> generate nsol2.vtp
        oomph-convert.py -p2 -z soln12.dat nsol2.vtp        -> generate nsol00002.vtp
        oomph-convert.py -p3 soln.dat                       -> generate soln.vtp
        oomph-convert.py -p2 -z soln.dat                    -> generate soln00000.vtp""")

#
###############################################################################
#

def main(argv):
    """ Main function
    """
    print("* oomph-convert.py, ver. 20110615")

    # Check python version. Bark at user if major is < 3
    major = sys.version_info[0]
    if major < 3:
        print("Python 2 is not supported by this script. Use oomph-convert.py2 instead", file=sys.stderr)
        sys.exit(3)

    # Get command-line arguments
    try:
        opts, args = getopt.getopt(argv, "hozp:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
        
    # Get options
    zero_pad_name_flag = False
    write_points_flag = False
    overwrite_flag = False
    for opt, arg in opts:
        if opt in ("-h"):
            usage()
            sys.exit()
        if opt in ("-z"):
            zero_pad_name_flag = True
        if opt in ("-p"):
            write_points_flag = True
            argdim = arg
            if argdim not in ("2","3"):    
                usage()
                sys.exit() 
        if opt in ("-o"):
            overwrite_flag = True

          
    if len(args) == 1:
        # Get filename and suffix
        ifilename = args[0]
        isuffix = ifilename.split(".")[-1]
        # Set filename and suffix
        extension_len = len(isuffix) + 1
        lenBaseName = len(ifilename) - extension_len
        
        if write_points_flag:
            ofilename = ifilename[:lenBaseName]+".vtp"
        else:
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


    # Zero pad output names if requested
    if zero_pad_name_flag:
        ofilename = addTrailingZeros(ofilename, osuffix)


    # Check that we are allowed to write to the ofilename (allowed to
    # overwrite or nothing to overwrite).
    if ok_to_write_ofile(overwrite_flag, ifilename, ofilename):
        start = time.time()

        # Convert from oomph-lib Tecplot format to VTK XML format
        if isuffix == "dat" and osuffix == "vtu":
            tecplot_to_vtkxml(ifilename, ofilename)

        # Convert from oomph-lib Tecplot format to VTP XML format
        elif isuffix == "dat" and osuffix == "vtp":
            tecplot_to_vtpxml(ifilename, ofilename,string.atoi(argdim))

        else:   
            error("Sorry, cannot convert between .%s and .%s file formats." % (isuffix, osuffix))

        end = time.time()
        print("* Conversion done in %d seconds" % (end - start))
        print('* Output file name: %(fn)s ' %{'fn': ofilename})


def ok_to_write_ofile(flag, in_file_name, out_file_name):
    """Check if the output file exists, if so decide if we should overwrite it
    and tell the user what we are doing."""

    if not os.path.exists(out_file_name):
        return True

    print("File", out_file_name, "already exists!")

    # Get modification times (in seconds since unix epoch I think)
    in_mtime = os.path.getmtime(in_file_name)
    out_mtime = os.path.getmtime(out_file_name)
    
    if flag:
        print("Overwriting regardless of modification times because of flag.")
        return True
    
    elif in_mtime > out_mtime:
        print("Overwriting because input file is newer than output file")
        return True

    else:
        print("Not overwriting.")
        return False


#
################################################################################
#

def addTrailingZeros(filename,osuffix):

    # -- Rename file by adding trailing zeros

    # 1) Define the number of digits in file name (up to limit_digits)

    # Assumption:
    # - the last four characters in ofilename define the extension of filename(".vtu" or ".vtp")

    limit_digits = 11
    extension_len = 4
    digits = 0
    for i in range(1, limit_digits+1):
        try:
            dummy = int(filename[-extension_len-i:-extension_len])
            digits = i
        except:
            pass

    # 2) Rename the file
    # - digits in filename: filename[-extension_len-digits:-extension_len]
    # - extension in filename: filename[-extension_len:]
    # - base file name in filename: filename[:lenBaseName]

    # Note: only filenames with 4 or less digits will have a trailing zero if
    # the second digits in "%05i.vtu" is 5 in os.rename

    lenBaseName = len(filename) - digits - extension_len
    if digits == 0:
        cifer = 0
    else:
        cifer = int(filename[-extension_len-digits:-extension_len])

    return filename[:lenBaseName]+"%05i." % cifer+"%s" % osuffix

#
################################################################################
#

def tecplot_to_vtkxml(inputFilename, outputFilename):
    """ Converts Tecplot file generates by Oomph into Vtk XML file.
    """

    #---------------------------------------------------------------------------
    # Retrieve the zones from the input Tecplot file
    #---------------------------------------------------------------------------
    try:
        input = open(inputFilename, "r")
    except:
        error("Failed to open input file for reading !")

    sys.stdout.write("Parse input file for Tecplot zones........")
    sys.stdout.flush()
    line = 0
    ignoredlines = 0
    zones = list()
    while 1:
        zone = None
        linetmp = line
        try:
            (zone, line) = TecplotZone.parse(input, line)
        except TecplotParsingError as e:
            input.close()
            error(str(e))
        
        if zone:
            zones.append(zone)
            ignoredlines+=line-linetmp-1-zone.nodesCount()
        else:
            print("done")
            input.close() # Close input file
            break

    nbzones = len(zones)
    sys.stdout.write("* %d lines ignored\n" % ignoredlines)
    if nbzones == 0:

        #---------------------------------------------------------------------------
        # Dummy output
        #---------------------------------------------------------------------------
        try:
            output = open(outputFilename, "w")
        except:
            error("Failed to open output file for writing !")
        
        output.close()
        error("The input file " + inputFilename + "\n does not contain any Tecplot zone! Created an empty file... \n You may want to try converting this file to point \n data with -p2 option if dim == 2 or -p3 option if dim == 3")

    #---------------------------------------------------------------------------
    # Compute global informations
    #---------------------------------------------------------------------------

    # Get the solution dimension (compute the maximum value)
    dimension=0
    for zone in zones:
        if (zone.dimension)[0] > dimension:
            dimension=(zone.dimension)[0]

    if dimension == 1:
        error("1D is not supported. Use GnuPlot to display the data file.")

    # Loop over the zones to get the number of nodes and cells
    nodesCount = 0
    cellsCount = 0
    for zone in zones:
        nodesCount += zone.nodesCount()
        cellsCount += zone.cellsCount[0]

    # Get the field count
    # Assumption: the number of fields is assumed to be constant over all zones
    fieldsCount = len(zones[0].nodes[0].fields)

    #---------------------------------------------------------------------------
    # Write into Vtk XML output file
    #---------------------------------------------------------------------------
    try:
        output = open(outputFilename, "w")
    except:
        error("Failed to open output file for writing !")

    output.write(VtkXml.header)
    output.write(VtkXml.unstructuredGridHeader)
    output.write(VtkXml.pieceHeader % (nodesCount, cellsCount))

    #---------------------------------------------------------------------------
    # Nodes
    #---------------------------------------------------------------------------
    sys.stdout.write("Write nodal coordinates...................")
    sys.stdout.flush()
    output.write(VtkXml.pointsHeader)
    for zone in zones:
        for node in zone.nodes:
            output.write("%e %e %e\n" %(node.coordinates[0], node.coordinates[1], node.coordinates[2]))
    output.write(VtkXml.pointsFooter)
    print("done")

    #---------------------------------------------------------------------------
    # Cells
    #---------------------------------------------------------------------------
    output.write(VtkXml.cellsHeader)

    # Cell connectivity
    #---------------------------------------------------------------------------
    sys.stdout.write("Write cell connectivity...................")
    sys.stdout.flush()
    output.write(VtkXml.connectivityHeader)
    pos = 0 # Current cell origin node index

    zoneCount=0
    for zone in zones:

        zoneCount+=1
        if zone.cellFormat[0] == 1:

            if zoneCount == 1: maxNode=0

            # dump connectivities
            for element in range (0, zone.cellsCount[0]):
                conn=zone.connectivities[element]
                for i in range(0, len(conn)):
                    dum=int(conn[i])
                    # renumber connectivity table if it belongs to zone > 1
                    if zoneCount > 1: dum+=maxNode+1
                    conn[i]=dum-1
                    output.write("%i " %(conn[i]))
                output.write("\n")

            # compute maximum node number for renumbering of next connectivity table
            maxNode=0
            for element in range (0, zone.cellsCount[0]):
                conn=zone.connectivities[element]
                for i in range(0, len(conn)):
                    conn[i]=int(conn[i])
                dum=max(conn)
                if dum > maxNode:
                    maxNode = dum

        if zone.cellFormat[0] == 2:

            dimI = zone.edges[0]

            if zone.dimension[0] == 1: # Line (dim I = 2)
                indexes = [pos, pos+1]
                output.write(" ".join(map(str, indexes)) + "\n")
                pos += dimI

            if zone.dimension[0] == 2: # Quad
                dimJ = zone.edges[1]
                indexes = 4 * [0]
                for j in range(dimJ - 1):
                    for i in range(dimI - 1):
                        # Unique face of the cell
                        indexes[0] = pos            # bottom-left node
                        indexes[1] = pos + 1        # bottom-right node
                        indexes[2] = pos + 1 + dimI # top-right node
                        indexes[3] = pos + dimI     # to-left node
                        output.write(" ".join(map(str, indexes)) + "\n")
                        # Next cell
                        pos += 1
                    # Next row of cells
                    pos += 1
                # Next zone
                pos += dimI
        
            if zone.dimension[0] == 3: # Hexahedron
                indexes = 8 * [0]
                dimJ = zone.edges[1]
                dimK = zone.edges[2]
                for k in range(dimK - 1):
                    for j in range(dimJ - 1):
                        for i in range(dimI - 1):
                            # Front face of the cell
                            indexes[0] = pos            # bottom-left node
                            indexes[1] = pos + 1        # bottom-right node
                            indexes[2] = pos + 1 + dimI # top-right node
                            indexes[3] = pos + dimI     # to-left node
                            # Back face of the cell
                            backPos = pos + dimI * dimJ
                            indexes[4] = backPos
                            indexes[5] = backPos + 1
                            indexes[6] = backPos + 1 + dimI
                            indexes[7] = backPos + dimI
                            output.write(" ".join(map(str, indexes)) + "\n")
                            # Next cell
                            pos += 1
                        # Next row of cells
                        pos += 1
                    # Next k
                    pos += dimI
                # Next zone
                pos += dimI * dimJ

    output.write(VtkXml.connectivityFooter)
    print("done")

    # Cell offset
    #---------------------------------------------------------------------------
    sys.stdout.write("Write cell offsets........................")
    sys.stdout.flush()
    output.write(VtkXml.offsetsHeader)
    offset = 0
    for zone in zones:
        for i in range(1,zone.cellsCount[0]+1):
            if zone.cellType[0] == 3: # VTK_LINE
                offset +=2
            if zone.cellType[0] == 5: # VTK_TRIANGLE
                offset +=3
            if zone.cellType[0] == 9: # VTK_QUAD
                offset +=4
            if zone.cellType[0] == 10: # VTK_TETRAHEDRON
                offset +=4
            if zone.cellType[0] == 12: # VTK_HEXAHEDRON
                offset +=8
            output.write(str(offset) + "\n")
    output.write(VtkXml.offsetsFooter)
    print("done")

    # Cell types
    #---------------------------------------------------------------------------
    sys.stdout.write("Write cell types..........................")
    sys.stdout.flush()

    output.write(VtkXml.typesHeader)    
    
    cellType0 = zones[0].cellType[0]
    warn = 0
    for zone in zones:
        if zone.cellType[0] == 3: # VTK_LINE
            cellType = "3"
        if zone.cellType[0] == 5: # VTK_TRIANGLE
            cellType = "5"
        if zone.cellType[0] == 9: # VTK_QUAD
            cellType = "9"
        if zone.cellType[0] == 10: # VTK_TETRAHEDRON
            cellType = "10"
        if zone.cellType[0] == 12: # VTK_HEXAHEDRON
            cellType = "12"

        if warn == 0 and not zone.cellType[0] == cellType0:#Check if types are differents
            warn = 1
        output.write(zone.cellsCount[0] * (cellType + "\n"))

    output.write(VtkXml.typesFooter)
    print("done")
    if warn == 1:
        sys.stdout.write("Warning: Different types of elements \n")
    output.write(VtkXml.cellsFooter)

    #---------------------------------------------------------------------------
    # Fields
    #---------------------------------------------------------------------------
    output.write(VtkXml.pointDataHeader)
    for fieldIndex in range(fieldsCount):
        sys.stdout.write("Write field %02d/%02d........................." % (fieldIndex + 1, fieldsCount))
        sys.stdout.flush()
        output.write(VtkXml.fieldHeader % (fieldIndex + 1))
        for zone in zones:
            for node in zone.nodes:
                output.write("%e\n" % node.fields[fieldIndex])

        print("done")
        output.write(VtkXml.fieldFooter)

    output.write(VtkXml.pointDataFooter)

    output.write(VtkXml.pieceFooter)
    output.write(VtkXml.unstructuredGridFooter)
    output.write(VtkXml.footer)

    #---------------------------------------------------------------------------
    # Close output file
    #---------------------------------------------------------------------------
    output.close()

#
################################################################################
#

def error(message):
    """ Write an error message
    """
    sys.stderr.write("\nCONVERSION FAILED\n")
    for line in message.split("\n"):
        sys.stderr.write("*** %s\n" % line)
    sys.exit(2)

#
################################################################################
#

class TecplotParsingError(Exception):
    """ Raised when an error occurs while parsing a Tecplot file.
    """
    def __init__(self, title, message, line, reminder=False):
        self.title = title
        self.message = message
        self.line = line
        self.reminder = reminder

    def __str__(self):
        msg = "Parsing failed at line %d: [%s] %s" % (self.line, self.title, self.message)

        return msg

#
################################################################################
#

class TecplotNode:
    """ This class handles Tecplot node informations.
    """
    def __init__(self):
        self.coordinates = [.0, .0, .0]
        """ Node coordinates
        """
        self.fields = list()
        """ Fields values
        """

#
################################################################################
#

class TecplotZone:
    """ This class handles Tecplot zone informations.
    """
    def __init__(self):
        self.edges = list()
        self.nodes = list()
        self.dimension = list()
        self.cellType = list()
        self.cellFormat = list()
        self.cellsCount = list()
        self.connectivities = list()
        """ List of :class:`TecplotZone`
        """

    def nodesCount(self):
        """ :return: The number of nodes in this zone
        """
        return len(self.nodes)

    @staticmethod
    def parse(file, line):
        """ Parse the zone of the given file starting from current position.

            :param file: A file open for reading
            :param line: The current position in the file as line number (used for error messages)
            :return: :class:`TecplotZone` if a zone has been found before the end of the file otherwise :const:`None`
            :return: The updated line position in input file
            :raise:  :class:`TecplotParsingError` on parsing error

            .. note:: This method will seek until it reach a Tecplot zone header.
        """
        #-----------------------------------------------------------------------
        # Seek to the next Tecplot zone
        #-----------------------------------------------------------------------
        
        while 1:
            header = file.readline()
            line += 1

            if not header:
                # We reach the end of the file
                return (None, line)

            if len(header) > 3 and header[0:4] == "ZONE":
                # We got a zone !
                break

        #-----------------------------------------------------------------------
        # Create the zone
        #-----------------------------------------------------------------------
        zone = TecplotZone()

        # Read header informations
        edges = header[len("ZONE"):].strip().split(",")

        # Assumption: Only two zone formats are possible.
        # The first (format 1) is defined with a zone header of the type
        # ZONE N=15, E=16, F=FEPOINT, ET=TRIANGLE
        # In this case, edges = ['N=15', ' E=16', ' F=FEPOINT', ' ET=TRIANGLE']
        # This zone has been coded for two ET keywords: TRIANGLE and TETRAHEDRON
        # The second (format 2) is defined in terms of the IJK indexes and a 
        # typical zone reads
        # ZONE I=5, J=5, K=5
        # One index, I, indicates a line element, 2 indexes, I and J, indicate
        # a 4 node quadrilateral element, and three indexes, I J and K, indicate
        # a tetrahedral element.

        if "TRIANGLE" in header:

            # Define data
            cell = 5 # VTK_TRIANGLE
            dim = 2
            format = 1
            nodesCount = int(edges[0].strip().strip("N="))
            cellsCount = int(edges[1].strip().strip("E="))

            # Append data to zone
            zone.dimension.append(dim)
            zone.cellType.append(cell)
            zone.cellFormat.append(format)
            zone.cellsCount.append(cellsCount)

        elif "TETRAHEDRON" in header:

            # Define data
            cell = 10 # VTK_TETRAHEDRON
            dim = 3
            format = 1
            nodesCount = int(edges[0].strip().strip("N="))
            cellsCount = int(edges[1].strip().strip("E="))

            # Append data to zone
            zone.dimension.append(dim)
            zone.cellType.append(cell)
            zone.cellFormat.append(format)
            zone.cellsCount.append(cellsCount)

        else:

            if not len(edges) in range(1, 4):
                raise TecplotParsingError("Invalid zone header",\
                                          "wrong edges count! Try to convert with -p option",\
                                          line,\
                                          reminder=True)

            format = 2

            # define dimension of zone cell by counting the IJK indexes
            dim = len(edges)

            if dim == 1:
                cell = 3 # VTK_LINE
            if dim == 2:
                cell = 9 # VTK_QUAD
            if dim == 3:
                cell = 12 # VTK_HEXAHEDRON

            # append dimension and cell type/format to the zone
            zone.dimension.append(dim)
            zone.cellType.append(cell)
            zone.cellFormat.append(format)

            # extract information from string
            labels = ["I", "J", "K"]
            for (index, edge) in enumerate(edges):
                edge = edge.strip()
                if len(edge) < 2 or edge[0] != labels[index]:
                    raise TecplotParsingError("Invalid zone header",\
                                              "wrong edges format (%s). Try to convert with -p option" % edges,\
                                              line,\
                                              reminder=True)
                try:
                    zone.edges.append(int(edge[2:]))
                except:
                    raise TecplotParsingError("Invalid zone header",\
                                              "wrong edge format (%s). Try to convert with -p option" % edge,\
                                              line)

            #-------------------------------------------------------------------
            # Node count in cell
            #-------------------------------------------------------------------
            nodesCount = 1
            for edge in zone.edges:
                nodesCount *= edge
                
            #-------------------------------------------------------------------
            # Cell count in zone
            #-------------------------------------------------------------------
            cellsCount = 1
            for edge in zone.edges:
                cellsCount *= edge - 1
            zone.cellsCount.append(cellsCount)

        #-----------------------------------------------------------------------
        # Parse nodes
        #-----------------------------------------------------------------------
        for nodeIndex in range(nodesCount):
            data = file.readline()
            line += 1

            if not data or data.strip() == "":
                data = file.readline() # ...in case of empty line

            data = data.strip().split(" ")
            if len(data) < dim:
                raise TecplotParsingError("Invalid zone", "not enough values for this node ! Try to convert with -p option", line)

            node = TecplotNode()
            for i, value in enumerate(data):
                try:
                    if i < dim:
                        node.coordinates[i] = float(value)
                    else:
                        node.fields.append(float(value))
                except ValueError:
                    raise TecplotParsingError("Invalid zone", "wrong node values ! Try to convert with -p option", line)

            # Append this node to the zone
            zone.nodes.append(node)

        #-----------------------------------------------------------------------
        # Parse connectivities (format == 1)
        #-----------------------------------------------------------------------
        if format == 1:
            for cellIndex in range(cellsCount):
                data = file.readline()
                line += 1

                if not data or data.strip() == "":
                    data = file.readline() # ...in case of empty line

                data = data.strip().split(" ")
                if cell == 3:  cellNodes = 2 # VTK_LINE
                if cell == 5:  cellNodes = 3 # VTK_TRIANGLE
                if cell == 9:  cellNodes = 4 # VTK_QUAD
                if cell == 10: cellNodes = 4 # VTK_TETRAHEDRON
                if cell == 12: cellNodes = 8 # VTK_HEXAHEDRON

                if len(data) != cellNodes:
                    raise TecplotParsingError("Invalid zone", "wrong connectivity list for this element!", line)

                # Append these connectivities to the zone
                zone.connectivities.append(data)

        return (zone, line)

#
################################################################################
#


class VtkXml:
    """ This class acts as a namespace for the definition Vtk XML headers and footers.
    """
    header = '<?xml version=\"1.0\"?>\n' \
             '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n'
    footer = '</VTKFile>\n'

    unstructuredGridHeader = '<UnstructuredGrid>\n'
    unstructuredGridFooter = '</UnstructuredGrid>\n'

    pieceHeader = '<Piece NumberOfPoints="%d" NumberOfCells="%d">\n'
    pieceFooter = '</Piece>'

    pointsHeader = '<Points>\n<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n'
    pointsFooter = '</DataArray>\n</Points>\n'

    cellsHeader = '<Cells>\n'
    cellsFooter = '</Cells>\n'

    connectivityHeader = '<DataArray type="Int32" Name="connectivity" format="ascii">\n'
    connectivityFooter = '</DataArray>\n'

    offsetsHeader = '<DataArray type="Int32" Name="offsets" format="ascii">\n'
    offsetsFooter = '</DataArray>\n'

    typesHeader = '<DataArray type="UInt8" Name="types" format="ascii">\n'
    typesFooter = '</DataArray>\n'

    pointDataHeader = '<PointData>\n'
    pointDataFooter = '</PointData>\n'

    fieldHeader = '<DataArray type="Float32" Name="V%d" format="ascii">\n'
    fieldFooter = '</DataArray>\n'


#
################################################################################
#

def tecplot_to_vtpxml(inputFilename, outputFilename, dim):
    """ Converts Tecplot file generates by Oomph into Vtp XML file.
    """

    #---------------------------------------------------------------------------
    # Retrieve the points from the input Tecplot file
    #---------------------------------------------------------------------------
    try:
        input = open(inputFilename, "r")
    except:
        error("Failed to open input file for reading !")

    sys.stdout.write("Parse input file for points........")
    sys.stdout.flush()
    line = 0
    prev_line=0
    offset_list=[]
    count=0
    nzone=1
    points = list()
    while 1:
        point = None
        try:
            (point, line) = InputPoints.parse(input, line, dim)
            #print "line %d %d is..." % (line, prev_line)
            if line != prev_line+1 :
                #print "...normal point"
                #else:
                #print "...start point"
                if prev_line !=0: 
                    offset_list.append(count)
                    nzone+=1
            count+=1
            prev_line=line

        except TecplotParsingError as e:
            input.close()
            error(str(e))

        if point:
            points.append(point)
        else:
            print("done")
            input.close() # Close input file
            break

    offset_list.append(count-1)
    nbpoints = len(points)

    if nbpoints == 0:
        error("The input file does not contain any point !")


    #---------------------------------------------------------------------------
    # Compute global informations
    #---------------------------------------------------------------------------

    # Get the field count
    # Assumption: the number of fields is assumed to be constant over all points
    fieldsCount = len(points[0].fields)

    #---------------------------------------------------------------------------
    # Write into Vtp XML output file
    #---------------------------------------------------------------------------
    try:
        output = open(outputFilename, "w")
    except:
        error("Failed to open output file for writing !")

    output.write(VtpXml.header)
    output.write(VtpXml.polyDataHeader)
    output.write(VtpXml.pieceHeader % (nbpoints,nzone))

    #---------------------------------------------------------------------------
    # Nodes
    #---------------------------------------------------------------------------
    sys.stdout.write("Write points coordinates...................")
    sys.stdout.flush()
    output.write(VtpXml.pointsHeader)
    for point in points:
        output.write("%e %e %e\n" %(point.coordinates[0], point.coordinates[1], point.coordinates[2]))
    output.write(VtpXml.pointsFooter)
    print("done")
 
    #---------------------------------------------------------------------------
    # Fields
    #---------------------------------------------------------------------------
    output.write(VtpXml.pointDataHeader)
    for fieldIndex in range(fieldsCount):
        sys.stdout.write("Write field %02d/%02d........................." % (fieldIndex + 1, fieldsCount))
        sys.stdout.flush()
        output.write(VtpXml.fieldHeader % (fieldIndex + 1))
        for point in points:
            output.write("%e\n" % point.fields[fieldIndex])

        print("done")
        output.write(VtpXml.fieldFooter)

    output.write(VtpXml.pointDataFooter)

    #---------------------------------------------------------------------------
    # Headers and footers
    #---------------------------------------------------------------------------

    output.write(VtpXml.vertsHeader)
    output.write(VtpXml.vertsFooter)
    output.write(VtpXml.linesHeader)
    #Prepare line information:
    output.write(VtkXml.connectivityHeader)
    for i in range(0,nbpoints):
        output.write("%i "%i)
    output.write(VtkXml.connectivityFooter)
    output.write(VtkXml.offsetsHeader)
    for offset in offset_list:
         output.write("%i "%offset)
    output.write(VtkXml.offsetsFooter)
    #end line information
    output.write(VtpXml.linesFooter)
    output.write(VtpXml.stripsHeader)
    output.write(VtpXml.stripsFooter)
    output.write(VtpXml.polysHeader)
    output.write(VtpXml.polysFooter)
    output.write(VtpXml.pieceFooter)
    output.write(VtpXml.polyDataFooter)
    output.write(VtpXml.footer)

    #---------------------------------------------------------------------------
    # Close output file
    #---------------------------------------------------------------------------
    output.close()

#
################################################################################
#

class InputPoints:
    """ This class handles Tecplot node informations.
    """
    def __init__(self):
        self.coordinates = [.0, .0, .0]
        """ Node coordinates
        """
        self.fields = list()
        """ Fields values
        """

    @staticmethod
    def parse(file, line, dim): 


        #-----------------------------------------------------------------------
        # Create the point
        #-----------------------------------------------------------------------
        point = InputPoints()

        #-----------------------------------------------------------------------
        # Seek to the next point
        #-----------------------------------------------------------------------
        while 1:
            pointline = file.readline()
            line += 1

            if not pointline:
                # We reach the end of the file
                return (None, line)

            # Boolean indicating that the line contains only numerical data
            is_fp_line=1 
                
            #Read the line and split into individual entries
            values = list()
            values = pointline.strip().split(" ")
            
            # Check that all entries are floating points numbers
            for i, value in enumerate(values):
                
                # Define regular expression for floating point number
                float_point_re=re.compile('[-+]?[0-9]*\.?[0-9]*([e,E][-+]?[0-9]+)?') 
                
                match_index=re.match(float_point_re,value)
                if match_index is not None:
                    if len(match_index.group()) != len(value):
                        is_fp_line=0
                        break
                else:
                    is_fp_line=0
                    break
                
            # If it's a line containing only floating point numbers: 
            # Extract coordinates and values
            if is_fp_line > 0:
                for i, value in enumerate(values):
                    try:
                        # Coordinate
                        if i < dim:
                            point.coordinates[i] = float(value)
                        # Value
                        else:
                            point.fields.append(float(value))
                    except ValueError:
                        raise TecplotParsingError("Invalid zone", "wrong node values !", line)
                    
                return (point, line)
            
#
################################################################################
#

class VtpXml:
    """ This class acts as a namespace for the definition Vtk XML headers and footers.
    """
    header = '<?xml version=\"1.0\"?>\n' \
             '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n'
    footer = '</VTKFile>\n'

    polyDataHeader = '<PolyData>\n'
    polyDataFooter = '</PolyData>\n'

    pieceHeader = '<Piece NumberOfPoints="%d" NumberOfVerts="0" NumberOfLines="%d" NumberOfStrips="0" NumberOfPolys="0">\n'
    pieceFooter = '</Piece>'

    pointsHeader = '<Points>\n<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n'
    pointsFooter = '</DataArray>\n</Points>\n'

    cellDataHeader = '<CellData>\n'
    cellDataFooter = '</CellData>\n'

    pointDataHeader = '<PointData>\n'
    pointDataFooter = '</PointData>\n'

    vertsHeader = '<Verts>\n'
    vertsFooter = '</Verts>\n'

    linesHeader = '<Lines>\n'
    linesFooter = '</Lines>\n'

    stripsHeader = '<Strips>\n'
    stripsFooter = '</Strips>\n'

    fieldHeader = '<DataArray type="Float32" Name="V%d" format="ascii">\n'
    fieldFooter = '</DataArray>\n'

    polysHeader = '<Polys>\n'
    polysFooter = '</Polys>\n'

#
################################################################################
#


if __name__ == "__main__":
    main(sys.argv[1:])
