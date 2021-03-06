#!MC 1100


$!VARSET |PNG|=1

$!LOOP 5


$!DRAWGRAPHICS FALSE

$!VARSET |STEP|=(|LOOP|-1)




$!READDATASET  '"RESLT/soln|STEP|.dat" ' 
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN2D

$!FIELDLAYERS SHOWMESH = NO
$!FIELDMAP [1-|NUMZONES|]  EDGELAYER{COLOR = RED}
$!PAPER SHOWPAPER = YES
$!TWODAXIS GRIDAREA{DRAWBORDER = YES}
$!TWODAXIS YDETAIL{RANGEMAX = 1.6}



$!GLOBALTWODVECTOR UVAR = 3
$!GLOBALTWODVECTOR VVAR = 4
$!FIELDLAYERS SHOWVECTOR = YES
$!FIELDMAP [1-|NUMZONES|]  VECTOR{SHOW = NO}

$!VARSET |NLAGRSTART|=(|NUMZONES|+1)


$!READDATASET  '"RESLT/lagr|STEP|.dat" ' 
  READDATAOPTION = APPEND
  RESETSTYLE = NO
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN2D

$!TWODAXIS VIEWPORTPOSITION{Y1 = 11.01}
$!TWODAXIS VIEWPORTPOSITION{Y2 = 88.01}
$!TWODAXIS VIEWPORTPOSITION{X2 = 79.7127215}
$!TWODAXIS VIEWPORTPOSITION{X1 = 13.256142544}

$!FRAMELAYOUT XYPOS{X = 0.393696850394}
$!FRAMELAYOUT XYPOS{Y = 0.393696850394}
$!FRAMELAYOUT WIDTH = 5.81843307087
$!FRAMELAYOUT HEIGHT = 7.4802992126


$!ACTIVEFIELDMAPS += [|NLAGRSTART|-|NUMZONES|]
$!FIELDMAP [|NLAGRSTART|-|NUMZONES|]  VECTOR{SHOW = YES}
$!FIELDMAP [|NLAGRSTART|-|NUMZONES|]  VECTOR{COLOR = BLUE}
$!FIELDMAP [|NLAGRSTART|-|NUMZONES|]  POINTS{IJKSKIP{I = 5}}
$!RESETVECTORLENGTH 
$!GLOBALTWODVECTOR RELATIVELENGTH = 0.5

$!TWODAXIS XDETAIL{TITLE{TITLEMODE = USETEXT}}
$!TWODAXIS XDETAIL{TITLE{TEXT = 'x'}}
$!TWODAXIS YDETAIL{TITLE{TITLEMODE = USETEXT}}
$!TWODAXIS YDETAIL{TITLE{TEXT = 'y'}}

$!VIEW FIT
$!DRAWGRAPHICS TRUE
$!REDRAWALL


$!IF |PNG|==1

        $!EXPORTSETUP EXPORTFORMAT = PNG
        $!EXPORTSETUP IMAGEWIDTH = 500
        $!EXPORTSETUP EXPORTFNAME = 'lagr_distort|STEP|.png'
        $!EXPORT 
          EXPORTREGION = ALLFRAMES

$!ELSE


        $!IF |LOOP|>1
                $!EXPORTNEXTFRAME
        $!ELSE

                $!EXPORTSETUP
                  EXPORTREGION = ALLFRAMES
                  EXPORTFORMAT = AVI
                  EXPORTFNAME = "lagr_distort.avi"
                $!EXPORTSETUP IMAGEWIDTH = 750
                $!EXPORTSTART
        $!ENDIF
$!ENDIF


$!ENDLOOP

$!IF |PNG|==0
        $!EXPORTFINISH
$!ENDIF





