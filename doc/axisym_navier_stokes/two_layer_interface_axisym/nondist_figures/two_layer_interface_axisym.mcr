#!MC 1100

# Use png output (otherwise avi)
$!VARSET |PNG| = 0

# Ask user for number of steps to loop
#$!GETUSERINPUT |nstep| INSTRUCTIONS = "Loop: Number of Steps?"
$!VARSET |nstep| = 160 # PATRICKFLAG

# Set number of interface elements (same as number of elements
# in the r direction)
#$!GETUSERINPUT |NumInterfaceElements| INSTRUCTIONS = "Number of interface elements?"
$!VARSET |NumInterfaceElements| = 16 # PATRICKFLAG

# Set first step
$!VARSET |lostep| = 0

# Set step increment
$!VARSET |dlstep| = 1

###############################################################
######################## START OF LOOP ########################
###############################################################

# Start of loop with |nstep| steps
$!LOOP |nstep|

# Set value of the current step (using |LOOP| variable which returns the
# loop counter e.g. |LOOP| = 1 on first time through loop etc.)
$!VARSET |istep| = ( ( |dlstep| * ( |LOOP|-1 ) ) + |lostep| )

# Clear the current layout and start again
$!NEWLAYOUT

# Turn off all graphics drawing
$!DRAWGRAPHICS FALSE

# Set paper to be A4 and portrait
$!PAPER SHOWPAPER = YES
$!PAPER ORIENTPORTRAIT = YES
$!PAPER PAPERSIZE = A4

#-------------------------------------------------------------#
#-------------------------- FRAME 1 --------------------------#
#------------------------ (TOP LEFT) -------------------------#
#-------------------------------------------------------------#

# Change size/position of frame that gets created by default
$!FRAMELAYOUT
  XYPOS
   {
    X = 0.0
    Y = 0.0
   }
   WIDTH = 4.0
   HEIGHT = 5.0

# Turn off current frame's border
$!FRAMELAYOUT SHOWBORDER = NO

###############################################################
########## START OF MACRO FILE GENERATED GRAPHICALLY ##########
######################## (WITH TWEAKS) ########################
###############################################################

# Read in data file
$!READDATASET  '"soln|istep|.dat"' 
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN2D
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5" "V6"'

# Add max/min data to the dataset
$!READDATASET  '"max_min.dat" '
  READDATAOPTION = APPEND
  RESETSTYLE = NO
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN2D
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5" "V6"'

# Rename variables
$!RENAMEDATASETVAR
  VAR = 1
  NAME = 'r' 
$!RENAMEDATASETVAR 
  VAR = 2
  NAME = 'z' 
$!RENAMEDATASETVAR 
  VAR = 3
  NAME = 'u_r(r,z,t)' 
$!RENAMEDATASETVAR 
  VAR = 4
  NAME = 'u_z(r,z,t)' 
$!RENAMEDATASETVAR 
  VAR = 5
  NAME = 'u_azi(r,z,t)' 
$!RENAMEDATASETVAR 
  VAR = 6
  NAME = 'p(r,z,t)' 

# Make all zones active for plotting
$!ACTIVEFIELDZONES  =  [1-|NUMZONES|]

# Set the number of zones to show (not including max/min data)
$!VARSET |NumElementZones| = (|NUMZONES|-1)

# Show edge for *all* zones (including interface elements) except max/min
# zone (this is the last zone)
# Make mesh lines dark green and have linewidth 0.1
$!FIELDMAP [1-|NumElementZones|]
  EDGELAYER
  {
   SHOW = YES
   COLOR = CUSTOM25
   LINETHICKNESS = 0.05
  }

# Make free surface mesh lines thicker PATRICKFLAG WHY DOESN'T THIS WORK?
$!VARSET |FirstFreeSurfaceZone| = (|NUMZONES|-|NumInterfaceElements|)
$!VARSET |LastFreeSurfaceZone| = |NumElementZones|
$!FIELDMAP [|FirstFreeSurfaceZone|-|LastFreeSurfaceZone|]
  EDGELAYER
  {
   SHOW = YES # If this is not included, interface elements' edges not shown
   COLOR = BLACK
   LINETHICKNESS = 0.5
  }

# Turn off mesh
$!FIELDLAYERS SHOWMESH = NO

# Turn on contour
$!FIELDLAYERS SHOWCONTOUR = YES

# Set global attributes associated with contour plot
$!GLOBALCONTOUR 1
  # Change variable associated with contour to be variable number 6 (pressure)
  VAR = 6
  # Set colour map to continuous rather than banded contours
  COLORMAPFILTER { COLORMAPDISTRIBUTION = CONTINUOUS }
  # Assign legend attributes
  LEGEND
   {
    SHOW = YES
    ISVERTICAL = YES
    AUTORESIZE = NO
    NUMBERTEXTSHAPE { HEIGHT = 3 }
    HEADERTEXTSHAPE { HEIGHT = 3 }
    ANCHORALIGNMENT = TOPRIGHT
    BOX { BOXTYPE = NONE }
    XYPOS
     {
      X = 95
      Y = 100
     }
   }

# Reset the contour levels to a set of evenly distributed, nice values
# spanning the entire range of the currently selected contouring variable,
# with a specified number of entries
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 1
  APPROXNUMVALUES = 11

# Use variables 3 and 4 (u_r and u_z) for vector arrows
#$!GLOBALTWODVECTOR 
#  UVAR = 3
#  VVAR = 4
#  RELATIVELENGTH = 0.075

# Assign attributes of the 2 axes
$!TWODAXIS
  XDETAIL
   {
    # Change variable associated with x-axis to be variable number 1
    VARNUM = 1
    SHOWAXIS = YES
    RANGEMIN = 0
    RANGEMAX = 1.0
    TICKLABEL { TEXTSHAPE { HEIGHT = 3 } }
    TITLE
     {
      TEXTSHAPE { HEIGHT = 3.5 }
     }
   }
  YDETAIL
   {
    # Change variable associated with y-axis to be variable number 2
    VARNUM = 2
    SHOWAXIS = YES
    RANGEMIN = 0
    RANGEMAX = 2.0
    TICKLABEL { TEXTSHAPE { HEIGHT = 3 } }
    TITLE
     {
      TEXTSHAPE { HEIGHT = 3.5 }
      OFFSET = 8
     }
   }

# Resize x axis to max/min values (preserve axis scale so that axis
# does not stretch out)
$!TWODAXIS PRESERVEAXISSCALE = YES
$!VIEW AXISFIT
  AXIS = 'X' 
  AXISNUM = 1

# Fit entire plot to grid area
#$!VIEW FIT

# Add text label
$!ATTACHTEXT
  ANCHORPOS
   {
    X = 5
    Y = 95
   }
  COLOR = PURPLE
  TEXTSHAPE
   {
    FONT = HELV
    HEIGHT = 12
   }
  BOX
   {
    MARGIN = 10
    LINETHICKNESS = 0.4
   }
  SCOPE = GLOBAL
  TEXT = 'pressure distribution'

###############################################################
########### END OF MACRO FILE GENERATED GRAPHICALLY ###########
###############################################################

# Fit entire paper to workspace
$!WORKSPACEVIEW FITPAPER

# Redraw all frames (must do this here!)
$!REDRAWALL

#####################################

$!IF |PNG|==1


        $!EXPORTSETUP EXPORTFORMAT = PNG
        $!EXPORTSETUP SUPERSAMPLEFACTOR = 3
        $!EXPORTSETUP USESUPERSAMPLEANTIALIASING = YES
        $!EXPORTSETUP IMAGEWIDTH = 600
        $!EXPORTSETUP EXPORTFNAME = 'two_layer_interface_axisym|istep|.png'
        $!EXPORT
          EXPORTREGION = ALLFRAMES

#        $!EXPORTSETUP EXPORTFORMAT = EPS
#        $!EXPORTSETUP USESUPERSAMPLEANTIALIASING = NO
#        $!EXPORTSETUP IMAGEWIDTH = 600
#        $!EXPORTSETUP EXPORTFNAME = 'two_layer_interface_axisym|istep|.eps'
#        $!EXPORTSETUP PRINTRENDERTYPE = IMAGE
#        $!EXPORTSETUP EXPORTFNAME = 'two_layer_interface_axisym|istep|.img.eps'
#        $!EXPORT
#          EXPORTREGION = ALLFRAMES

$!ELSE
 
        $!IF |LOOP|>1
                $!EXPORTNEXTFRAME
        $!ELSE
                $!EXPORTSETUP
                 EXPORTFORMAT = AVI
                 EXPORTFNAME = "two_layer_interface_axisym.avi"
                 # Tell it to export ALL the frames!
                 EXPORTREGION = ALLFRAMES
                $!EXPORTSETUP IMAGEWIDTH = 460
                $!EXPORTSETUP USESUPERSAMPLEANTIALIASING = NO
                $!EXPORTSTART
        $!ENDIF

$!ENDIF


$!VARSET |LAST_STEP|=|istep|

#####################################

# End of loop
$!ENDLOOP

###############################################################
######################### END OF LOOP #########################
###############################################################

$!IF |PNG|==0
        $!EXPORTFINISH
$!ENDIF
