#!MC 1000

#$!GETUSERINPUT |lostep| INSTRUCTIONS = "Loop. First Step??"
$!VARSET  |lostep|=0
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
#$!VARSET |istep|*=3
$!READDATASET "soln|istep|.dat " 
  READDATAOPTION = NEW
  RESETSTYLE = YES


$!FIELDLAYERS SHOWMESH = NO
$!FIELD [1-|NUMZONES|]  BOUNDARY{COLOR = BLACK}
$!FIELD [1-|NUMZONES|]  BOUNDARY{LINETHICKNESS = 0.1}
$!TWODAXIS GRIDAREA{DRAWBORDER = YES}
$!REDRAWALL 

$!TWODAXIS YDETAIL{RANGEMIN = -1.4}
$!TWODAXIS YDETAIL{RANGEMAX = 1.4}
$!VIEW TRANSLATE
  X = 10
  Y = 0
$!VIEW TRANSLATE
  X = 10
  Y = 0
$!VIEW TRANSLATE
  X = -5
  Y = 0

$!TWODAXIS XDETAIL{TITLE{TITLEMODE = USETEXT}}
$!TWODAXIS XDETAIL{TITLE{TEXT = 'x'}}
$!TWODAXIS YDETAIL{TITLE{TITLEMODE = USETEXT}}
$!TWODAXIS YDETAIL{TITLE{TEXT = 'y'}}


$!DRAWGRAPHICS TRUE
$!REDRAW 




$!EXPORTSETUP EXPORTFORMAT = PNG
$!EXPORTSETUP IMAGEWIDTH = 600
$!EXPORTSETUP EXPORTFNAME = 'static_fish|LOOP|.png'
$!EXPORT
     EXPORTREGION = ALLFRAMES


$!EndLoop


