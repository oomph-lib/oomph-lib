#!MC 1200
# Created by Tecplot 360 build 12.0.0.3454

$!VARSET |stem_of_input_file| = "soln"
$!VARSET |nstep| = 100
$!VARSET |first_soln_number| = 0
$!VARSET |soln_skip| = 1 # (1==do every solution)

###############################################################
######################## START OF LOOP ########################
###############################################################

## Start of loop with |nstep| steps
$!LOOP |nstep|

## Set the loop counter
$!VARSET |soln_counter| = ( |first_soln_number| + ((|LOOP| - 1)*|soln_skip|) )

## Compute the counter used in the output (png) file (pad with zeros)
$!VARSET |output_counter| = "|soln_counter|"
$!IF |soln_counter| < 10000
$!VARSET |output_counter| = "0|soln_counter|"
$!ENDIF
$!IF |soln_counter| < 1000
$!VARSET |output_counter| = "00|soln_counter|"
$!ENDIF
$!IF |soln_counter| < 100
$!VARSET |output_counter| = "000|soln_counter|"
$!ENDIF
$!IF |soln_counter| < 10
$!VARSET |output_counter| = "0000|soln_counter|"
$!ENDIF

## Set up Export file type and file name.
$!EXPORTSETUP EXPORTFORMAT = PNG
$!EXPORTSETUP EXPORTFNAME = "spin_up_stream|output_counter|.png"
$!EXPORTSETUP IMAGEWIDTH = 600
$!EXPORTSETUP USESUPERSAMPLEANTIALIASING = YES
$!EXPORTSETUP SUPERSAMPLEFACTOR = 4

#-------------------------------------------------------------#
#------------------------ LEFT FRAME -------------------------#
#-------------------------------------------------------------#

## Change size/position of frame that gets created by default
$!FRAMELAYOUT XYPOS{X = 0.0}
$!FRAMELAYOUT XYPOS{Y = 0.0}
$!FRAMELAYOUT WIDTH = 5.0
$!FRAMELAYOUT HEIGHT = 7.0
$!FRAMELAYOUT SHOWBORDER = NO

###############################################################
########## START OF MACRO FILE GENERATED GRAPHICALLY ##########
######################## (WITH TWEAKS) ########################
###############################################################

## Read in data file
$!READDATASET  '"|stem_of_input_file||soln_counter|.dat" '
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN2D
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5" "V6"'

## Read in data file
$!READDATASET  '"max_min.dat" '
  READDATAOPTION = APPEND
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN2D
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5" "V6"'

## Rename variables
$!RENAMEDATASETVAR
  VAR = 1
  NAME = 'r' 
$!RENAMEDATASETVAR 
  VAR = 2
  NAME = 'z' 

## Make all zones (including dummy zones) active for plotting
$!ACTIVEFIELDZONES  =  [1-|NUMZONES|]

## Set axes to be independent and reset to variable min/max
$!TWODAXIS AXISMODE = INDEPENDENT
$!VIEW AXISFIT
  AXIS = 'X'
  AXISNUM = 1
$!VIEW AXISFIT
  AXIS = 'Y'
  AXISNUM = 1

## Change edge layer properties for all zones
$!FIELDMAP [1-|NUMZONES|]  EDGELAYER{LINETHICKNESS = 0.0875}

## Change font size on the axes
$!TWODAXIS XDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 2.675}}}
$!TWODAXIS YDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 2.675}}}
$!TWODAXIS YDETAIL{TITLE{TEXTSHAPE{HEIGHT = 3.2125}}}
$!TWODAXIS XDETAIL{TITLE{TEXTSHAPE{HEIGHT = 3.2125}}}
$!TWODAXIS XDETAIL{TICKLABEL{OFFSET = 0.875}}
$!TWODAXIS XDETAIL{TITLE{OFFSET = 4.25}}

## Change axis details
$!TWODAXIS XDETAIL{AXISLINE{LINETHICKNESS = 0.3575}}
$!TWODAXIS XDETAIL{TICKS{LENGTH = 1.7875}}
$!TWODAXIS XDETAIL{TICKS{LINETHICKNESS = 0.3575}}
$!TWODAXIS XDETAIL{TICKS{MINORLENGTH = 0.857}}
$!TWODAXIS XDETAIL{TICKS{MINORLINETHICKNESS = 0.0875}}
$!TWODAXIS YDETAIL{AXISLINE{LINETHICKNESS = 0.3575}}
$!TWODAXIS YDETAIL{TICKS{LENGTH = 1.7875}}
$!TWODAXIS YDETAIL{TICKS{LINETHICKNESS = 0.3575}}
$!TWODAXIS YDETAIL{TICKS{MINORLENGTH = 0.857}}
$!TWODAXIS YDETAIL{TICKS{MINORLINETHICKNESS = 0.0875}}

## Turn on contours
$!FIELDLAYERS SHOWCONTOUR = YES
$!GLOBALCONTOUR 1  VAR = 5
$!GLOBALCONTOUR 1  COLORMAPFILTER{COLORMAPDISTRIBUTION = CONTINUOUS}
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 1
  APPROXNUMVALUES = 15

## Show contour legend
#$!GLOBALCONTOUR 1  LEGEND{SHOW = YES}

## Create a field of vectors using variables 3 and 4
$!GLOBALTWODVECTOR UVAR = 3
$!GLOBALTWODVECTOR VVAR = 4
$!RESETVECTORLENGTH 
#$!FIELDLAYERS SHOWVECTOR = YES

$!ATTACHTEXT 
  ANCHORPOS
    {
    X = 0.0
    Y = 94.0
    }
  COLOR = PURPLE
  TEXTSHAPE
    {
    FONT = HELV
    HEIGHT = 16
    }
  TEXT = 'azimuthal velocity distribution' 

## Define termination line for streamlines
$!STREAMTRACE SETTERMINATIONLINE
  RAWDATA
6
0.52 1.0
0.52 1.5
-0.05 1.5
-0.05 -0.1
0.52 -0.1
0.52 0.4
$!GLOBALSTREAM TERMLINE{ISACTIVE = YES}
$!GLOBALSTREAM TERMLINE{SHOW = NO}

## Add streamlines
$!STREAMTRACE ADD
  STREAMTYPE = TWODLINE
  DIRECTION = BOTH
  STARTPOS
    {
    X = 0.1
    Y = 0.6
    }
$!STREAMTRACE ADD
  STREAMTYPE = TWODLINE
  DIRECTION = BOTH
  STARTPOS
    {
    X = 0.1
    Y = 0.8
    }
$!STREAMTRACE ADD
  STREAMTYPE = TWODLINE
  DIRECTION = BOTH
  STARTPOS
    {
    X = 0.25
    Y = 0.9
    }
$!STREAMTRACE ADD
  STREAMTYPE = TWODLINE
  DIRECTION = BOTH
  STARTPOS
    {
    X = 0.25
    Y = 0.5
    }
$!STREAMTRACE ADD
  STREAMTYPE = TWODLINE
  DIRECTION = BOTH
  STARTPOS
    {
    X = 0.4
    Y = 0.3
    }
$!STREAMTRACE ADD
  STREAMTYPE = TWODLINE
  DIRECTION = BOTH
  STARTPOS
    {
    X = 0.4
    Y = 1.1
    }

## Streamline options
$!GLOBALSTREAM 
  STREAMTIMING
    {
    DELTATIME = 0.1
    }
  CELLFRACTION = 0.1
  MINCELLFRACTION = 1E-06
  MAXSTEPS = 5000

## Fit to everything
$!VIEW FIT

###############################################################
########### END OF MACRO FILE GENERATED GRAPHICALLY ###########
###############################################################

#-------------------------------------------------------------#
#-------------------- END OF LEFT FRAME ----------------------#
#-------------------------------------------------------------#

#-------------------------------------------------------------#
#----------------------- RIGHT FRAME -------------------------#
#-------------------------------------------------------------#

## Create a new frame
$!CREATENEWFRAME
  XYPOS
   {
    X = 5.0
    Y = 0.0
   }
   WIDTH = 5.0
   HEIGHT = 7.0
$!FRAMELAYOUT SHOWBORDER = NO

###############################################################
########## START OF MACRO FILE GENERATED GRAPHICALLY ##########
######################## (WITH TWEAKS) ########################
###############################################################

## Read in data file
$!READDATASET  '"|stem_of_input_file||soln_counter|.dat" '
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN2D
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5" "V6"'

## Read in data file
$!READDATASET  '"max_min.dat" '
  READDATAOPTION = APPEND
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN2D
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5" "V6"'

## Rename variables
$!RENAMEDATASETVAR
  VAR = 1
  NAME = 'r' 
$!RENAMEDATASETVAR 
  VAR = 2
  NAME = 'z' 

## Make all zones (including dummy zones) active for plotting
$!ACTIVEFIELDZONES  =  [1-|NUMZONES|]

## Set axes to be independent and reset to variable min/max
$!TWODAXIS AXISMODE = INDEPENDENT
$!VIEW AXISFIT
  AXIS = 'X'
  AXISNUM = 1
$!VIEW AXISFIT
  AXIS = 'Y'
  AXISNUM = 1

## Change edge layer properties for all zones
$!FIELDMAP [1-|NUMZONES|]  EDGELAYER{LINETHICKNESS = 0.0875}

## Change font size on the axes
$!TWODAXIS XDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 2.675}}}
$!TWODAXIS YDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 2.675}}}
$!TWODAXIS YDETAIL{TITLE{TEXTSHAPE{HEIGHT = 3.2125}}}
$!TWODAXIS XDETAIL{TITLE{TEXTSHAPE{HEIGHT = 3.2125}}}
$!TWODAXIS XDETAIL{TICKLABEL{OFFSET = 0.875}}
$!TWODAXIS XDETAIL{TITLE{OFFSET = 4.25}}

## Change axis details
$!TWODAXIS XDETAIL{AXISLINE{LINETHICKNESS = 0.3575}}
$!TWODAXIS XDETAIL{TICKS{LENGTH = 1.7875}}
$!TWODAXIS XDETAIL{TICKS{LINETHICKNESS = 0.3575}}
$!TWODAXIS XDETAIL{TICKS{MINORLENGTH = 0.857}}
$!TWODAXIS XDETAIL{TICKS{MINORLINETHICKNESS = 0.0875}}
$!TWODAXIS YDETAIL{AXISLINE{LINETHICKNESS = 0.3575}}
$!TWODAXIS YDETAIL{TICKS{LENGTH = 1.7875}}
$!TWODAXIS YDETAIL{TICKS{LINETHICKNESS = 0.3575}}
$!TWODAXIS YDETAIL{TICKS{MINORLENGTH = 0.857}}
$!TWODAXIS YDETAIL{TICKS{MINORLINETHICKNESS = 0.0875}}

## Turn on contours
$!FIELDLAYERS SHOWCONTOUR = YES
$!GLOBALCONTOUR 1  VAR = 6
$!GLOBALCONTOUR 1  COLORMAPFILTER{COLORMAPDISTRIBUTION = CONTINUOUS}
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 1
  APPROXNUMVALUES = 15

## Show contour legend
#$!GLOBALCONTOUR 1  LEGEND{SHOW = YES}

## Create a field of vectors using variables 3 and 4
$!GLOBALTWODVECTOR UVAR = 3
$!GLOBALTWODVECTOR VVAR = 4
$!RESETVECTORLENGTH 
#$!FIELDLAYERS SHOWVECTOR = YES

$!ATTACHTEXT 
  ANCHORPOS
    {
    X = 0.0
    Y = 94.0
    }
  COLOR = PURPLE
  TEXTSHAPE
    {
    FONT = HELV
    HEIGHT = 16
    }
  TEXT = 'pressure distribution' 

## Define termination line for streamlines
$!STREAMTRACE SETTERMINATIONLINE
  RAWDATA
6
0.52 1.0
0.52 1.5
-0.05 1.5
-0.05 -0.1
0.52 -0.1
0.52 0.4
$!GLOBALSTREAM TERMLINE{ISACTIVE = YES}
$!GLOBALSTREAM TERMLINE{SHOW = NO}

## Add streamlines
$!STREAMTRACE ADD
  STREAMTYPE = TWODLINE
  DIRECTION = BOTH
  STARTPOS
    {
    X = 0.1
    Y = 0.6
    }
$!STREAMTRACE ADD
  STREAMTYPE = TWODLINE
  DIRECTION = BOTH
  STARTPOS
    {
    X = 0.1
    Y = 0.8
    }
$!STREAMTRACE ADD
  STREAMTYPE = TWODLINE
  DIRECTION = BOTH
  STARTPOS
    {
    X = 0.25
    Y = 0.9
    }
$!STREAMTRACE ADD
  STREAMTYPE = TWODLINE
  DIRECTION = BOTH
  STARTPOS
    {
    X = 0.25
    Y = 0.5
    }
$!STREAMTRACE ADD
  STREAMTYPE = TWODLINE
  DIRECTION = BOTH
  STARTPOS
    {
    X = 0.4
    Y = 0.3
    }
$!STREAMTRACE ADD
  STREAMTYPE = TWODLINE
  DIRECTION = BOTH
  STARTPOS
    {
    X = 0.4
    Y = 1.1
    }

## Streamline options
$!GLOBALSTREAM 
  STREAMTIMING
    {
    DELTATIME = 0.1
    }
  CELLFRACTION = 0.1
  MINCELLFRACTION = 1E-06
  MAXSTEPS = 5000

## Fit to everything
$!VIEW FIT

###############################################################
########### END OF MACRO FILE GENERATED GRAPHICALLY ###########
###############################################################

#-------------------------------------------------------------#
#------------------- END OF RIGHT FRAME ----------------------#
#-------------------------------------------------------------#

## Redraw all frames (must do this here!)
$!REDRAWALL 

## Export all frames to the current .png file
$!EXPORT 
  EXPORTREGION = ALLFRAMES

## End of loop
$!ENDLOOP

###############################################################
######################### END OF LOOP #########################
###############################################################
