#!MC 1000

$!VARSET |PNG|=0


#$!GETUSERINPUT |lostep| INSTRUCTIONS = "Loop. First Step??"
$!VARSET  |lostep|=0
#$!GETUSERINPUT |dlstep| INSTRUCTIONS = "Loop. Step Increment?"
$!VARSET  |dlstep|=1
$!GETUSERINPUT |nstep| INSTRUCTIONS = "Loop. Number of Steps??"
#$!VARSET |nstep| = 300

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
#$!VARSET |istep|*=10

$!VarSet |LFDSFN1| = '"RESLT/soln|istep|.dat"'
$!VarSet |LFDSVL1| = '"V1" "V2" "V3"'
$!VarSet |LFDSFN2| = '"RESLT/Wall|istep|.dat"'
$!VarSet |LFDSVL2| = '"V1" "V2" "V3"'
$!SETSTYLEBASE FACTORY
$!PAPER 
  BACKGROUNDCOLOR = WHITE
  ISTRANSPARENT = YES
  ORIENTPORTRAIT = NO
  SHOWGRID = YES
  SHOWRULER = YES
  SHOWPAPER = YES
  PAPERSIZE = A4
  PAPERSIZEINFO
    {
    A3
      {
      WIDTH = 11.693
      HEIGHT = 16.535
      }
    A4
      {
      WIDTH = 8.2677
      HEIGHT = 11.693
      LEFTHARDCLIPOFFSET = 0.125
      RIGHTHARDCLIPOFFSET = 0.125
      TOPHARDCLIPOFFSET = 0.125
      BOTTOMHARDCLIPOFFSET = 0.125
      }
    }
  RULERSPACING = ONECENTIMETER
  PAPERGRIDSPACING = ONETENTHCENTIMETER
  REGIONINWORKAREA
    {
    X1 = -0.05
    Y1 = -0.05
    X2 = 11.74
    Y2 = 8.318
    }
$!COLORMAP 
  CONTOURCOLORMAP = SMRAINBOW
$!COLORMAPCONTROL RESETTOFACTORY
### Frame Number 1 ###
$!READDATASET  '|LFDSFN1|' 
  INITIALPLOTTYPE = CARTESIAN2D
  INCLUDETEXT = YES
  INCLUDEGEOM = YES
  VARLOADMODE = BYNAME
  VARNAMELIST = '|LFDSVL1|' 
$!REMOVEVAR |LFDSVL1|
$!REMOVEVAR |LFDSFN1|

$!VARSET |FIRST_ZONES| = |NUMZONES|


$!READDATASET  '|LFDSFN2|' 
  INITIALPLOTTYPE = CARTESIAN2D
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  READDATAOPTION = APPEND
  RESETSTYLE = NO
  VARLOADMODE = BYNAME
  VARNAMELIST = '|LFDSVL2|' 
$!REMOVEVAR |LFDSVL2|
$!REMOVEVAR |LFDSFN2|
$!FRAMELAYOUT 
  SHOWBORDER = NO
  HEADERCOLOR = RED
  XYPOS
    {
    X = 0.3937
    Y = 0.3937
    }
  WIDTH = 10.906
  HEIGHT = 7.4803
$!PLOTTYPE  = CARTESIAN2D
$!FRAMENAME  = 'Frame 001' 
$!ACTIVEFIELDZONES  =  [1-|NUMZONES|]
$!GLOBALRGB 
  RANGEMIN = 0
  RANGEMAX = 1
$!GLOBALCONTOUR  1
  VAR = 3
  LEGEND
    {
    SHOW = YES
    SHOWHEADER = NO
    XYPOS
      {
      X = 95
      }
    }
  COLORCUTOFF
    {
    RANGEMIN = -0.504052549601
    RANGEMAX = 0.50600245595
    }
  COLORMAPFILTER
    {
    CONTINUOUSCOLOR
      {
      CMIN = -1.00908005238
      CMAX = 1.01102995872
      }
    }
$!CONTOURLEVELS NEW
  CONTOURGROUP = 1
  RAWDATA
19
-0.9
-0.8
-0.7
-0.6
-0.5
-0.4
-0.3
-0.2
-0.1
0
0.1
0.2
0.3
0.4
0.5
0.6
0.7
0.8
0.9
$!GLOBALCONTOUR  2
  LEGEND
    {
    XYPOS
      {
      X = 95
      }
    }
  COLORMAPFILTER
    {
    CONTINUOUSCOLOR
      {
      CMIN = 0
      CMAX = 1
      }
    }
$!GLOBALCONTOUR  3
  LEGEND
    {
    XYPOS
      {
      X = 95
      }
    }
  COLORMAPFILTER
    {
    CONTINUOUSCOLOR
      {
      CMIN = 0
      CMAX = 1
      }
    }
$!GLOBALCONTOUR  4
  LEGEND
    {
    XYPOS
      {
      X = 95
      }
    }
  COLORMAPFILTER
    {
    CONTINUOUSCOLOR
      {
      CMIN = 0
      CMAX = 1
      }
    }
$!GLOBALSCATTER 
  LEGEND
    {
    XYPOS
      {
      X = 95
      }
    }
  REFSCATSYMBOL
    {
    COLOR = RED
    FILLCOLOR = RED
    }
$!FIELD  [1-|FIRST_ZONES|]
  MESH
    {
    COLOR = RED
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = RED
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = RED
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = WHITE
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = WHITE
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    }
  SURFACES
    {
    SURFACESTOPLOT = KPLANES
    }
  VOLUMEMODE
    {
    VOLUMEOBJECTSTOPLOT
      {
      SHOWISOSURFACES = NO
      SHOWSLICES = NO
      SHOWSTREAMTRACES = NO
      }
    }



$!FIELD  [|NUMZONES|]
  MESH
    {
    COLOR = RED
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = RED
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = RED
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = WHITE
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.8
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    }
  SURFACES
    {
    SURFACESTOPLOT = JPLANES
    }
  VOLUMEMODE
    {
    VOLUMEOBJECTSTOPLOT
      {
      SHOWISOSURFACES = NO
      SHOWSLICES = NO
      SHOWSTREAMTRACES = NO
      }
    }



$!TWODAXIS 
  XDETAIL
    {
    VARNUM = 1
    }
  YDETAIL
    {
    VARNUM = 2
    }
$!VIEW FIT
$!TWODAXIS 
  DEPXTOYRATIO = 1
$!TWODAXIS 
  XDETAIL
    {
    RANGEMIN = 0
    RANGEMAX = 1.6331078116701687186
    GRSPACING = 0.2
    TITLE
      {
      TITLEMODE = USETEXT
      TEXT = 'x' 
      }
    }
$!TWODAXIS 
  YDETAIL
    {
    RANGEMIN = 0
    RANGEMAX = 1.1499999999999999112
    GRSPACING = 0.2
    TITLE
      {
      TITLEMODE = USETEXT
      TEXT = 'y' 
      }
    }
$!GLOBALISOSURFACE 
  ISOVALUE1 = -0.504052549601
  ISOVALUE2 = 0.000974953174591
  ISOVALUE3 = 0.50600245595
  MARCHINGCUBEALGORITHM = CLASSICPLUS
$!GLOBALSLICE 
  BOUNDARY
    {
    SHOW = NO
    }
$!FIELDLAYERS 
  SHOWMESH = NO
  SHOWCONTOUR = YES
$!ATTACHTEXT 
  ANCHORPOS
    {
    X = 32.7473533265
    Y = 92.5454133066
    }
  TEXTSHAPE
    {
    FONT = HELV
    HEIGHT = 26
    }
  BOX
    {
    MARGIN = 10
    LINETHICKNESS = 0.4
    }
  SCOPE = GLOBAL
  MACROFUNCTIONCOMMAND = '' 
  TEXT = 'Contours of the solution' 

$!SETSTYLEBASE CONFIG


#------------------------
# don't show dummy zones
#------------------------
$!VARSET |LAST_BUT_ONE|=(|NUMZONES|-1)
$!VARSET |LAST_BUT_TWO|=(|NUMZONES|-2)

$!ACTIVEFIELDZONES -= [|LAST_BUT_TWO|-|LAST_BUT_ONE|]


############################



$!REDRAWALL


$!IF |PNG|==1


        $!EXPORTSETUP EXPORTFORMAT = PNG
        $!EXPORTSETUP IMAGEWIDTH = 600
        $!EXPORTSETUP EXPORTFNAME = 'unsteady_heat_contour|istep|.png'
        $!EXPORT
          EXPORTREGION = CURRENTFRAME

        $!EXPORTSETUP EXPORTFORMAT = EPS
        $!EXPORTSETUP IMAGEWIDTH = 1423
        $!EXPORTSETUP EXPORTFNAME = 'unsteady_heat_contour|istep|.eps'

        $!EXPORTSETUP PRINTRENDERTYPE = IMAGE
        $!EXPORTSETUP EXPORTFNAME = 'unsteady_heat_contour|istep|.img.eps'
        $!EXPORT
          EXPORTREGION = CURRENTFRAME

$!ELSE

        $!IF |LOOP|>1
                $!EXPORTNEXTFRAME
        $!ELSE

                $!EXPORTSETUP
                 EXPORTFORMAT = AVI
                 EXPORTFNAME = "unsteady_heat_contour.avi"
                $!EXPORTSETUP IMAGEWIDTH = 829
                $!EXPORTSTART
        $!ENDIF

$!ENDIF


$!EndLoop




$!IF |PNG|==0
        $!EXPORTFINISH
$!ENDIF


#$!QUIT

