#!MC 1200
# Created by Tecplot 360 build 12.0.0.3454

$!VARSET |stem_of_input_file| = "soln"
$!VARSET |nstep| = 240
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

#-------------------------------------------------------------#
#------------------------ MAIN FRAME -------------------------#
#-------------------------------------------------------------#

## Change size/position of frame that gets created by default
$!FRAMELAYOUT XYPOS{X = 0.0}
$!FRAMELAYOUT XYPOS{Y = 0.0}
$!FRAMELAYOUT WIDTH = 7.0
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
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5"'

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
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5"'

## Rename variables
$!RENAMEDATASETVAR
  VAR = 1
  NAME = 'x' 
$!RENAMEDATASETVAR 
  VAR = 2
  NAME = 'y' 

## Make all zones (including dummy zones) active for plotting
$!ACTIVEFIELDZONES  =  [1-|NUMZONES|]

## Set up Export file type and file name.
$!EXPORTSETUP EXPORTFORMAT = PNG
$!EXPORTSETUP EXPORTFNAME = "soln|output_counter|.png"
$!EXPORTSETUP IMAGEWIDTH = 500
$!EXPORTSETUP USESUPERSAMPLEANTIALIASING = YES
$!EXPORTSETUP SUPERSAMPLEFACTOR = 4

## Set axes to be independent and reset to variable min/max
$!TWODAXIS AXISMODE = INDEPENDENT
$!VIEW AXISFIT
  AXIS = 'X'
  AXISNUM = 1
$!VIEW AXISFIT
  AXIS = 'Y'
  AXISNUM = 1

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
$!FIELDLAYERS SHOWVECTOR = YES

## Fit to everything
$!VIEW FIT

###############################################################
########### END OF MACRO FILE GENERATED GRAPHICALLY ###########
###############################################################

#-------------------------------------------------------------#
#-------------------- END OF MAIN FRAME ----------------------#
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
