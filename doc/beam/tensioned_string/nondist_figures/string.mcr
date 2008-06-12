#!MC 1000

#$!GETUSERINPUT |lostep| INSTRUCTIONS = "Loop. First Step??"
$!VARSET  |lostep|=1
#$!GETUSERINPUT |dlstep| INSTRUCTIONS = "Loop. Step Increment?"
$!VARSET  |dlstep|=1
$!GETUSERINPUT |nstep| INSTRUCTIONS = "Loop. Number of Steps??"

$!LOOP |nstep|
$!VarSet |nnstep| = |LOOP|
$!VarSet |nnstep| -= 1
$!VarSet |iistep| = |dlstep|
$!VarSet |iistep| *= |nnstep|
$!VarSet |iistep| += |lostep|
$!NEWLAYOUT
$!DRAWGRAPHICS FALSE
#    $!IF |iistep| < 10 
#      $!VARSET |istep|='00|iistep|'
#    $!ENDIF
#    $!IF |iistep| > 9 
#      $!VARSET |istep|='0|iistep|'
#    $!ENDIF
#    $!IF |iistep| > 99 
#      $!VARSET |istep|=|iistep|
#    $!ENDIF
$!VARSET |istep|=|iistep|
$!READDATASET '"beam|istep|.dat "'
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INITIALPLOTTYPE = CARTESIAN2D

$!TWODAXIS GRIDAREA{DRAWBORDER = YES}
$!TWODAXIS XDETAIL{TITLE{TITLEMODE = USETEXT}}
$!TWODAXIS XDETAIL{TITLE{TEXT = 'x'}}
$!TWODAXIS YDETAIL{TITLE{TITLEMODE = USETEXT}}
$!TWODAXIS YDETAIL{TITLE{TEXT = 'y'}}
$!REDRAW 
$!FIELD [1-|NUMZONES|]  BOUNDARY{COLOR = BLUE}
$!FIELD [1-|NUMZONES|]  BOUNDARY{SHOW = YES}
$!FIELDLAYERS SHOWMESH = NO
$!REDRAWALL 
$!TWODAXIS YDETAIL{RANGEMIN = -12}
$!TWODAXIS XDETAIL{RANGEMIN = -3}
$!TWODAXIS XDETAIL{RANGEMAX = 13}



$!IF |LOOP|<2 
  $!EXPORTSETUP EXPORTFORMAT = AVI
  $!EXPORTSETUP IMAGEWIDTH = 600
  $!EXPORTSETUP EXPORTFNAME = 'string.avi'
  $!EXPORTSTART 
$!ENDIF
$!IF |LOOP|>1 
  $!EXPORTNEXTFRAME
$!ENDIF




$!EndLoop

$!EXPORTFINISH 
