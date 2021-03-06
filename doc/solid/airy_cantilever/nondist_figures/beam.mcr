#!MC 1000

# Use png output (otherwise avi)
$!VARSET |PNG|=1



$!VARSET |lostep|=0
$!VARSET |dlstep|=1
$!VARSET |nstep|=5

$!LOOP |nstep|


$!DRAWGRAPHICS FALSE

$!NEWLAYOUT
$!VARSET |istep|=(|lostep|+(|LOOP|-1)*|dlstep|)


$!READDATASET  '"RESLT/soln|istep|.dat" ' 
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN2D
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5"' 
$!FIELDLAYERS SHOWMESH = NO
$!GLOBALCONTOUR 1  VAR = 3
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 1
  APPROXNUMVALUES = 15
$!FIELDLAYERS SHOWCONTOUR = YES
$!FIELDMAP [1-|NUMZONES|]  EDGELAYER{COLOR = BLACK}
$!FIELDMAP [1-|NUMZONES|]  EDGELAYER{LINETHICKNESS = 0.02}
$!TWODAXIS GRIDAREA{DRAWBORDER = YES}
$!TWODAXIS YDETAIL{RANGEMAX = 1}
$!TWODAXIS YDETAIL{RANGEMAX = 1}
$!TWODAXIS YDETAIL{RANGEMIN = -7}
$!TWODAXIS XDETAIL{RANGEMAX = 12}
$!TWODAXIS XDETAIL{TITLE{TITLEMODE = USETEXT}}
$!TWODAXIS XDETAIL{TITLE{TEXT = 'x<SUB>1</SUB>'}}
$!TWODAXIS XDETAIL{TITLE{TEXTSHAPE{HEIGHT = 4.6}}}
$!TWODAXIS XDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 4}}}
$!TWODAXIS YDETAIL{TICKLABEL{TEXTSHAPE{HEIGHT = 4}}}
$!TWODAXIS YDETAIL{TITLE{TEXTSHAPE{HEIGHT = 4.6}}}
$!TWODAXIS YDETAIL{TITLE{TITLEMODE = USETEXT}}
$!TWODAXIS YDETAIL{TITLE{TEXT = 'X<SUB>2</SUB>'}}
$!GLOBALCONTOUR 1  LEGEND{SHOW = YES}
$!GLOBALCONTOUR 1  LEGEND{XYPOS{X = 85}}
$!RENAMEDATASETVAR 
  VAR = 3
  NAME = '<GREEK>s</GREEK><SUB>11</SUB>' 


$!PAPER SHOWPAPER = YES

$!DRAWGRAPHICS TRUE
$!REDRAWALL



$!IF |PNG|==1


        $!EXPORTSETUP EXPORTFORMAT = PNG
        $!EXPORTSETUP IMAGEWIDTH = 600
        $!EXPORTSETUP EXPORTFNAME = 'beam|LOOP|.png'
        $!EXPORT
          EXPORTREGION = ALLFRAMES

$!ELSE

        $!IF |LOOP|>1
                $!EXPORTNEXTFRAME
        $!ELSE

                $!EXPORTSETUP
                  EXPORTREGION = ALLFRAMES
                  EXPORTFORMAT = AVI
                  EXPORTFNAME = "beam.avi"
                $!EXPORTSETUP IMAGEWIDTH = 750
                $!EXPORTSTART
        $!ENDIF

$!ENDIF


$!ENDLOOP


$!IF |PNG|==0
        $!EXPORTFINISH
$!ENDIF
