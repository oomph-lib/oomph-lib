#!MC 1000

#!MC 1000

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
#$!VARSET |istep|+=1
$!VARSET |istep|*=5

$!READDATASET  '"RESLT_full/Wall|istep|.dat" ' 
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  INITIALPLOTTYPE = CARTESIAN2D
  VARNAMELIST = '"V1" "V2"' 
$!TWODAXIS GRIDAREA{DRAWBORDER = YES}
$!TWODAXIS XDETAIL{RANGEMAX = 1.4}
$!TWODAXIS YDETAIL{RANGEMAX = 1.3}
$!TWODAXIS XDETAIL{TITLE{TITLEMODE = USETEXT}}
$!TWODAXIS XDETAIL{TITLE{TEXT = 'x_1'}}

$!TWODAXIS YDETAIL{TITLE{TITLEMODE = USETEXT}}
$!TWODAXIS YDETAIL{TITLE{TEXT = 'x_2'}}
$!REDRAWALL 
$!FIELD [1]  MESH{LINETHICKNESS = 0.4}


$!ATTACHTEXT 
  ANCHORPOS
    {
    X = 22.4347127
    Y = 13.9807747
    }
  COLOR = BLUE
  TEXTSHAPE
    {
    HEIGHT = 20
    }
  TEXT = 'Neumann boundary' 
$!PICK ADD
  X = 1.98163636364
  Y = 5.65344467935
$!ATTACHTEXT 
  ANCHORPOS
    {
    X = 17.892929
    Y = 25.8225588
    }
  COLOR = RED
  TEXTSHAPE
    {
    HEIGHT = 19
    }
  ANGLE = 90
  TEXT = 'Dirichlet boundary' 
$!ATTACHTEXT 
  ANCHORPOS
    {
    X = 39.29373
    Y = 63.40889
    }
  COLOR = RED
  TEXTSHAPE
    {
    HEIGHT = 19
    }
  ANGLE = -45
  TEXT = 'Dirichlet boundary' 
############################


$!REDRAWALL


$!IF |PNG|==1


        $!EXPORTSETUP EXPORTFORMAT = PNG
        $!EXPORTSETUP USESUPERSAMPLEANTIALIASING = YES
        $!EXPORTSETUP IMAGEWIDTH = 600
        $!EXPORTSETUP EXPORTFNAME = 'domain|istep|.png'
        $!EXPORT
          EXPORTREGION = CURRENTFRAME

        $!EXPORTSETUP EXPORTFORMAT = EPS
        $!EXPORTSETUP IMAGEWIDTH = 1423
        $!EXPORTSETUP EXPORTFNAME = 'domain|istep|.eps'

        $!EXPORTSETUP PRINTRENDERTYPE = IMAGE
        $!EXPORTSETUP EXPORTFNAME = 'domain|istep|.img.eps'
        $!EXPORT
          EXPORTREGION = CURRENTFRAME

$!ELSE

        $!IF |LOOP|>1
                $!EXPORTNEXTFRAME
        $!ELSE

                $!EXPORTSETUP
                 EXPORTFORMAT = AVI
                 USESUPERSAMPLEANTIALIASING = YES
                 EXPORTFNAME = "domain.avi"
                $!EXPORTSETUP IMAGEWIDTH = 829
                $!EXPORTSTART
        $!ENDIF

$!ENDIF


$!VARSET |LAST_STEP|=|istep|

$!EndLoop


$!IF |PNG|==0
        $!EXPORTFINISH
$!ENDIF

