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
$!EXPORTSETUP EXPORTFNAME = "spin_up_carpet|output_counter|.png"
$!EXPORTSETUP IMAGEWIDTH = 600
$!EXPORTSETUP USESUPERSAMPLEANTIALIASING = YES
$!EXPORTSETUP SUPERSAMPLEFACTOR = 4

#-------------------------------------------------------------#
#-------------------------- FRAME 1 --------------------------#
#-------------------------------------------------------------#

## Change size/position of frame that gets created by default
$!FRAMELAYOUT XYPOS{X = 0.0}
$!FRAMELAYOUT XYPOS{Y = 0.0}
$!FRAMELAYOUT WIDTH = 5.0
$!FRAMELAYOUT HEIGHT = 4.6
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
  INITIALPLOTTYPE = CARTESIAN3D
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
  INITIALPLOTTYPE = CARTESIAN3D
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

## Choose which field to display on the z-axis
$!THREEDAXIS ZDETAIL{VARNUM = 3}

## Set axes to be independent and reset to variable min/max
$!THREEDAXIS AXISMODE = INDEPENDENT
$!VIEW AXISFIT
  AXIS = 'X'
  AXISNUM = 1
$!VIEW AXISFIT
  AXIS = 'Y'
  AXISNUM = 1
$!VIEW AXISFIT
  AXIS = 'Z'
  AXISNUM = 1

## Set 3D axis preferences
$!THREEDAXIS XDETAIL{SHOWAXIS = YES}
$!THREEDAXIS YDETAIL{SHOWAXIS = YES}
$!THREEDAXIS ZDETAIL{SHOWAXIS = YES}
$!THREEDAXIS XDETAIL{GRIDLINES{SHOW = NO}}
$!THREEDAXIS YDETAIL{GRIDLINES{SHOW = NO}}
$!THREEDAXIS ZDETAIL{GRIDLINES{SHOW = NO}}
$!THREEDAXIS GRIDAREA{USELIGHTSOURCETOFILL = NO}
$!THREEDAXIS XDETAIL{TITLE{OFFSET = 12}}
$!THREEDAXIS XDETAIL{TITLE{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS XDETAIL{TITLE{TEXTSHAPE{HEIGHT = 14}}}
$!THREEDAXIS YDETAIL{TITLE{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS YDETAIL{TITLE{TEXTSHAPE{HEIGHT = 14}}}
$!THREEDAXIS ZDETAIL{TITLE{SHOWONAXISLINE = NO}}
$!THREEDAXIS XDETAIL{TICKLABEL{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS XDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 10}}}
$!THREEDAXIS YDETAIL{TICKLABEL{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS YDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 10}}}
$!THREEDAXIS ZDETAIL{TICKLABEL{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS ZDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 10}}}
$!THREEDAXIS XDETAIL{AUTOGRID = NO}
$!THREEDAXIS XDETAIL{GRSPACING = 0.5}
$!THREEDAXIS YDETAIL{AUTOGRID = NO}
$!THREEDAXIS YDETAIL{GRSPACING = 0.2}
$!THREEDAXIS ZDETAIL{AUTOGRID = NO}
$!THREEDAXIS ZDETAIL{GRSPACING = 0.005}
$!THREEDAXIS XDETAIL{AXISLINE{LINETHICKNESS = 0.1}}
$!THREEDAXIS YDETAIL{AXISLINE{LINETHICKNESS = 0.1}}
$!THREEDAXIS ZDETAIL{AXISLINE{LINETHICKNESS = 0.1}}

## Set 3D view
$!ROTATE3DVIEW PSI
  ANGLE = 10
  ROTATEORIGINLOCATION = DEFINEDORIGIN
$!ROTATE3DVIEW THETA
  ANGLE = 10
  ROTATEORIGINLOCATION = DEFINEDORIGIN

## Change the axis scale factor
## (to make the four plots roughly the same dimensions)
$!GLOBALTHREED AXISSCALEFACT{Z = 70.0}

## Change edge layer properties for all zones
$!FIELDMAP [1-|NUMZONES|]  EDGELAYER{LINETHICKNESS = 0.0875}

## Remove orientation axis (in top-right corner)
$!THREEDAXIS FRAMEAXIS{SHOW = NO}

## Add title text
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
  TEXT = 'radial velocity' 

## Fit to everything
$!VIEW FIT

## Reduce size slightly to fit axis labels on
$!VIEW SETMAGNIFICATION
  MAGNIFICATION = 0.92

###############################################################
########### END OF MACRO FILE GENERATED GRAPHICALLY ###########
###############################################################

#-------------------------------------------------------------#
#---------------------- END OF FRAME 1 -----------------------#
#-------------------------------------------------------------#

#-------------------------------------------------------------#
#-------------------------- FRAME 2 --------------------------#
#-------------------------------------------------------------#

## Create a new frame
$!CREATENEWFRAME
  XYPOS
   {
    X = 5.0
    Y = 0.0
   }
   WIDTH = 5.0
   HEIGHT = 4.6
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
  INITIALPLOTTYPE = CARTESIAN3D
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
  INITIALPLOTTYPE = CARTESIAN3D
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

## Choose which field to display on the z-axis
$!THREEDAXIS ZDETAIL{VARNUM = 4}

## Set axes to be independent and reset to variable min/max
$!THREEDAXIS AXISMODE = INDEPENDENT
$!VIEW AXISFIT
  AXIS = 'X'
  AXISNUM = 1
$!VIEW AXISFIT
  AXIS = 'Y'
  AXISNUM = 1
$!VIEW AXISFIT
  AXIS = 'Z'
  AXISNUM = 1

## Set 3D axis preferences
$!THREEDAXIS XDETAIL{SHOWAXIS = YES}
$!THREEDAXIS YDETAIL{SHOWAXIS = YES}
$!THREEDAXIS ZDETAIL{SHOWAXIS = YES}
$!THREEDAXIS XDETAIL{GRIDLINES{SHOW = NO}}
$!THREEDAXIS YDETAIL{GRIDLINES{SHOW = NO}}
$!THREEDAXIS ZDETAIL{GRIDLINES{SHOW = NO}}
$!THREEDAXIS GRIDAREA{USELIGHTSOURCETOFILL = NO}
$!THREEDAXIS XDETAIL{TITLE{OFFSET = 12}}
$!THREEDAXIS XDETAIL{TITLE{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS XDETAIL{TITLE{TEXTSHAPE{HEIGHT = 14}}}
$!THREEDAXIS YDETAIL{TITLE{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS YDETAIL{TITLE{TEXTSHAPE{HEIGHT = 14}}}
$!THREEDAXIS ZDETAIL{TITLE{SHOWONAXISLINE = NO}}
$!THREEDAXIS XDETAIL{TICKLABEL{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS XDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 10}}}
$!THREEDAXIS YDETAIL{TICKLABEL{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS YDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 10}}}
$!THREEDAXIS ZDETAIL{TICKLABEL{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS ZDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 10}}}
$!THREEDAXIS XDETAIL{AUTOGRID = NO}
$!THREEDAXIS XDETAIL{GRSPACING = 0.5}
$!THREEDAXIS YDETAIL{AUTOGRID = NO}
$!THREEDAXIS YDETAIL{GRSPACING = 0.2}
$!THREEDAXIS ZDETAIL{AUTOGRID = NO}
$!THREEDAXIS ZDETAIL{GRSPACING = 0.005}
$!THREEDAXIS XDETAIL{AXISLINE{LINETHICKNESS = 0.1}}
$!THREEDAXIS YDETAIL{AXISLINE{LINETHICKNESS = 0.1}}
$!THREEDAXIS ZDETAIL{AXISLINE{LINETHICKNESS = 0.1}}

## Set 3D view
$!ROTATE3DVIEW PSI
  ANGLE = 10
  ROTATEORIGINLOCATION = DEFINEDORIGIN
$!ROTATE3DVIEW THETA
  ANGLE = 10
  ROTATEORIGINLOCATION = DEFINEDORIGIN

## Change the axis scale factor
## (to make the four plots roughly the same dimensions)
$!GLOBALTHREED AXISSCALEFACT{Z = 54.0}

## Change edge layer properties for all zones
$!FIELDMAP [1-|NUMZONES|]  EDGELAYER{LINETHICKNESS = 0.0875}

## Remove orientation axis (in top-right corner)
$!THREEDAXIS FRAMEAXIS{SHOW = NO}

## Add title text
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
  TEXT = 'axial velocity' 

## Fit to everything
$!VIEW FIT

## Reduce size slightly to fit axis labels on
$!VIEW SETMAGNIFICATION
  MAGNIFICATION = 0.92

###############################################################
########### END OF MACRO FILE GENERATED GRAPHICALLY ###########
###############################################################

#-------------------------------------------------------------#
#---------------------- END OF FRAME 2 -----------------------#
#-------------------------------------------------------------#

#-------------------------------------------------------------#
#-------------------------- FRAME 3 --------------------------#
#-------------------------------------------------------------#

## Create a new frame
$!CREATENEWFRAME
  XYPOS
   {
    X = 0.0
    Y = 4.6
   }
   WIDTH = 5.0
   HEIGHT = 4.6
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
  INITIALPLOTTYPE = CARTESIAN3D
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
  INITIALPLOTTYPE = CARTESIAN3D
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

## Choose which field to display on the z-axis
$!THREEDAXIS ZDETAIL{VARNUM = 5}

## Set axes to be independent and reset to variable min/max
$!THREEDAXIS AXISMODE = INDEPENDENT
$!VIEW AXISFIT
  AXIS = 'X'
  AXISNUM = 1
$!VIEW AXISFIT
  AXIS = 'Y'
  AXISNUM = 1
$!VIEW AXISFIT
  AXIS = 'Z'
  AXISNUM = 1

## Set 3D axis preferences
$!THREEDAXIS XDETAIL{SHOWAXIS = YES}
$!THREEDAXIS YDETAIL{SHOWAXIS = YES}
$!THREEDAXIS ZDETAIL{SHOWAXIS = YES}
$!THREEDAXIS XDETAIL{GRIDLINES{SHOW = NO}}
$!THREEDAXIS YDETAIL{GRIDLINES{SHOW = NO}}
$!THREEDAXIS ZDETAIL{GRIDLINES{SHOW = NO}}
$!THREEDAXIS GRIDAREA{USELIGHTSOURCETOFILL = NO}
$!THREEDAXIS XDETAIL{TITLE{OFFSET = 12}}
$!THREEDAXIS XDETAIL{TITLE{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS XDETAIL{TITLE{TEXTSHAPE{HEIGHT = 14}}}
$!THREEDAXIS YDETAIL{TITLE{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS YDETAIL{TITLE{TEXTSHAPE{HEIGHT = 14}}}
$!THREEDAXIS ZDETAIL{TITLE{SHOWONAXISLINE = NO}}
$!THREEDAXIS XDETAIL{TICKLABEL{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS XDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 10}}}
$!THREEDAXIS YDETAIL{TICKLABEL{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS YDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 10}}}
$!THREEDAXIS ZDETAIL{TICKLABEL{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS ZDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 10}}}
$!THREEDAXIS XDETAIL{AUTOGRID = NO}
$!THREEDAXIS XDETAIL{GRSPACING = 0.5}
$!THREEDAXIS YDETAIL{AUTOGRID = NO}
$!THREEDAXIS YDETAIL{GRSPACING = 0.2}
$!THREEDAXIS ZDETAIL{AUTOGRID = NO}
$!THREEDAXIS ZDETAIL{GRSPACING = 0.2}
$!THREEDAXIS XDETAIL{AXISLINE{LINETHICKNESS = 0.1}}
$!THREEDAXIS YDETAIL{AXISLINE{LINETHICKNESS = 0.1}}
$!THREEDAXIS ZDETAIL{AXISLINE{LINETHICKNESS = 0.1}}

## Set 3D view
$!ROTATE3DVIEW PSI
  ANGLE = 10
  ROTATEORIGINLOCATION = DEFINEDORIGIN
$!ROTATE3DVIEW THETA
  ANGLE = 10
  ROTATEORIGINLOCATION = DEFINEDORIGIN

## Change the axis scale factor
## (to make the four plots roughly the same dimensions)
$!GLOBALTHREED AXISSCALEFACT{Z = 0.9}

## Change edge layer properties for all zones
$!FIELDMAP [1-|NUMZONES|]  EDGELAYER{LINETHICKNESS = 0.0875}

## Remove orientation axis (in top-right corner)
$!THREEDAXIS FRAMEAXIS{SHOW = NO}

## Add title text
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
  TEXT = 'azimuthal velocity' 

## Fit to everything
$!VIEW FIT

## Reduce size slightly to fit axis labels on
$!VIEW SETMAGNIFICATION
  MAGNIFICATION = 0.92

###############################################################
########### END OF MACRO FILE GENERATED GRAPHICALLY ###########
###############################################################

#-------------------------------------------------------------#
#---------------------- END OF FRAME 3 -----------------------#
#-------------------------------------------------------------#

#-------------------------------------------------------------#
#-------------------------- FRAME 4 --------------------------#
#-------------------------------------------------------------#

## Create a new frame
$!CREATENEWFRAME
  XYPOS
   {
    X = 5.0
    Y = 4.6
   }
   WIDTH = 5.0
   HEIGHT = 4.6
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
  INITIALPLOTTYPE = CARTESIAN3D
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
  INITIALPLOTTYPE = CARTESIAN3D
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

## Choose which field to display on the z-axis
$!THREEDAXIS ZDETAIL{VARNUM = 6}

## Set axes to be independent and reset to variable min/max
$!THREEDAXIS AXISMODE = INDEPENDENT
$!VIEW AXISFIT
  AXIS = 'X'
  AXISNUM = 1
$!VIEW AXISFIT
  AXIS = 'Y'
  AXISNUM = 1
$!VIEW AXISFIT
  AXIS = 'Z'
  AXISNUM = 1

## Set 3D axis preferences
$!THREEDAXIS XDETAIL{SHOWAXIS = YES}
$!THREEDAXIS YDETAIL{SHOWAXIS = YES}
$!THREEDAXIS ZDETAIL{SHOWAXIS = YES}
$!THREEDAXIS XDETAIL{GRIDLINES{SHOW = NO}}
$!THREEDAXIS YDETAIL{GRIDLINES{SHOW = NO}}
$!THREEDAXIS ZDETAIL{GRIDLINES{SHOW = NO}}
$!THREEDAXIS GRIDAREA{USELIGHTSOURCETOFILL = NO}
$!THREEDAXIS XDETAIL{TITLE{OFFSET = 12}}
$!THREEDAXIS XDETAIL{TITLE{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS XDETAIL{TITLE{TEXTSHAPE{HEIGHT = 14}}}
$!THREEDAXIS YDETAIL{TITLE{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS YDETAIL{TITLE{TEXTSHAPE{HEIGHT = 14}}}
$!THREEDAXIS ZDETAIL{TITLE{SHOWONAXISLINE = NO}}
$!THREEDAXIS XDETAIL{TICKLABEL{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS XDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 10}}}
$!THREEDAXIS YDETAIL{TICKLABEL{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS YDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 10}}}
$!THREEDAXIS ZDETAIL{TICKLABEL{TEXTSHAPE{SIZEUNITS = POINT}}}
$!THREEDAXIS ZDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 10}}}
$!THREEDAXIS XDETAIL{AUTOGRID = NO}
$!THREEDAXIS XDETAIL{GRSPACING = 0.5}
$!THREEDAXIS YDETAIL{AUTOGRID = NO}
$!THREEDAXIS YDETAIL{GRSPACING = 0.2}
$!THREEDAXIS ZDETAIL{AUTOGRID = NO}
$!THREEDAXIS ZDETAIL{GRSPACING = 0.5}
$!THREEDAXIS XDETAIL{AXISLINE{LINETHICKNESS = 0.1}}
$!THREEDAXIS YDETAIL{AXISLINE{LINETHICKNESS = 0.1}}
$!THREEDAXIS ZDETAIL{AXISLINE{LINETHICKNESS = 0.1}}

## Set 3D view
$!ROTATE3DVIEW PSI
  ANGLE = 10
  ROTATEORIGINLOCATION = DEFINEDORIGIN
$!ROTATE3DVIEW THETA
  ANGLE = 10
  ROTATEORIGINLOCATION = DEFINEDORIGIN

## Change the axis scale factor
## (to make the four plots roughly the same dimensions)
$!GLOBALTHREED AXISSCALEFACT{Z = 0.37}

## Change edge layer properties for all zones
$!FIELDMAP [1-|NUMZONES|]  EDGELAYER{LINETHICKNESS = 0.0875}

## Remove orientation axis (in top-right corner)
$!THREEDAXIS FRAMEAXIS{SHOW = NO}

## Add title text
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
  TEXT = 'pressure' 

## Fit to everything
$!VIEW FIT

## Reduce size slightly to fit axis labels on
$!VIEW SETMAGNIFICATION
  MAGNIFICATION = 0.92

###############################################################
########### END OF MACRO FILE GENERATED GRAPHICALLY ###########
###############################################################

#-------------------------------------------------------------#
#---------------------- END OF FRAME 4 -----------------------#
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
