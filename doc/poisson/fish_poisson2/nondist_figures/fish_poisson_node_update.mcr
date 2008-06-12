#!MC 900

$!VARSET |PNG|=1

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


$!CREATERECTANGULARZONE
  IMAX = 2
  JMAX = 2
  KMAX = 2
  X1 = -0.4
  Y1 = -1.0
  Z1 = 0
  X2 = 1.7
  Y2 = 1.0
  Z2 = 0.2
  XVAR = 1
  YVAR = 2
  ZVAR = 3


$!FRAMEMODE = THREED
$!FIELDLAYERS SHOWMESH = NO
$!FIELDLAYERS SHOWSHADE= YES
$!FIELD [1-|NUMZONES|]  BOUNDARY{COLOR = BLACK}
$!FIELD [1-|NUMZONES|]  BOUNDARY{LINETHICKNESS = 0.02}

$!THREEDVIEW 
  PSIANGLE = 41.6842
  THETAANGLE = 341.835
  ALPHAANGLE = 0
  VIEWERPOSITION
    {
    X = 4.9132659192
    Y = -12.189088911
    Z = 15.0121092049
    }
$!VIEW PUSH



$!VIEW FIT
$!REDRAW 


$!THREEDAXIS YDETAIL{TITLE{TITLEMODE = USETEXT}}
$!THREEDAXIS YDETAIL{TITLE{TEXT = 'x_2'}}
$!THREEDAXIS XDETAIL{TITLE{TITLEMODE = USETEXT}}
$!THREEDAXIS XDETAIL{TITLE{TEXT = 'x_1'}}
$!THREEDAXIS ZDETAIL{TITLE{SHOWONAXISLINE = NO}}
$!THREEDAXIS AXISMODE = XYDEPENDENT
$!GLOBALTHREED AXISSCALEFACT{Z = 2.03513577767}
$!THREEDVIEW
  PSIANGLE = 47.0266
  THETAANGLE = 25.1503
  ALPHAANGLE = 0
  VIEWERPOSITION
    {
    X = -4.58552359496
    Y = -11.1970530788
    Z = 5.77108303458
    }


$!VIEW FIT

$!FIELD [|NUMZONES|]  SHADE{SHOW = NO}
$!FIELD [|NUMZONES|]  BOUNDARY{SHOW = NO}

$!DRAWGRAPHICS TRUE
$!REDRAW 

$!IF |PNG|==1


        $!EXPORTSETUP EXPORTFORMAT = PNG
        $!EXPORTSETUP IMAGEWIDTH = 600
        $!EXPORTSETUP EXPORTFNAME = 'soln|istep|.png'
        $!EXPORT
          EXPORTREGION = CURRENTFRAME

$!ELSE

        $!IF |LOOP|>1
                $!EXPORTNEXTFRAME
        $!ELSE

                $!EXPORTSETUP
                 EXPORTFORMAT = AVI
                 EXPORTFNAME = "soln.avi"
                $!EXPORTSETUP IMAGEWIDTH = 829
                $!EXPORTSTART
        $!ENDIF

$!ENDIF


$!VARSET |LAST_STEP|=|istep|

$!EndLoop


$!IF |PNG|==0
        $!EXPORTFINISH
$!ENDIF



