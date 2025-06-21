#!MC 1000

$!VARSET |PNG|=0


#$!GETUSERINPUT |lostep| INSTRUCTIONS = "Loop. First Step??"
$!VARSET  |lostep|=0
#$!GETUSERINPUT |dlstep| INSTRUCTIONS = "Loop. Step Increment?"
$!VARSET  |dlstep|=1
$!GETUSERINPUT |nstep| INSTRUCTIONS = "Loop. Number of Steps??"
#$!VARSET |nstep| = 75

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

$!DRAWGRAPHICS FALSE


$!VarSet |LFDSFN1| = '"RESLT_CR/soln|istep|.dat"'
$!VarSet |LFDSVL1| = '"V1" "V2" "V3" "V4" "V5"'
$!VarSet |LFDSFN2| = '"RESLT_CR/exact_soln|istep|.dat"'
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
$!RENAMEDATASETVAR 
  VAR = 1
  NAME = 'x' 
$!RENAMEDATASETVAR 
  VAR = 2
  NAME = 'y' 
$!RENAMEDATASETVAR 
  VAR = 3
  NAME = 'u(y,t)' 
$!RENAMEDATASETVAR 
  VAR = 4
  NAME = 'v(y,t)' 
$!RENAMEDATASETVAR 
  VAR = 5
  NAME = 'p(y,t)' 
$!READDATASET  '|LFDSFN2|' 
  INITIALPLOTTYPE = CARTESIAN2D
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  READDATAOPTION = APPEND
  RESETSTYLE = NO
$!REMOVEVAR |LFDSFN1|
$!REMOVEVAR |LFDSFN2|
$!REMOVEVAR |LFDSVL1|
$!FRAMELAYOUT 
  SHOWBORDER = NO
  SHOWHEADER = NO
  HEADERCOLOR = RED
  XYPOS
    {
    X = 0.125
    Y = 0.20954
    }
  WIDTH = 11.443
  HEIGHT = 7.8486
$!PLOTTYPE  = CARTESIAN2D
$!FRAMENAME  = 'Frame 001' 
$!ACTIVEFIELDZONES  =  [1-100]
$!GLOBALRGB 
  RANGEMIN = 0
  RANGEMAX = 1
$!GLOBALCONTOUR  1
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
$!GLOBALTHREEDVECTOR 
  UVAR = 1
  VVAR = 2
  WVAR = 3
$!GLOBALTWODVECTOR 
  UVAR = 3
  VVAR = 4
  RELATIVELENGTH = 0.4
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
$!FIELD  [1]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [2]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [3]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [4]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [5]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [6]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [7]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [8]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [9]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [10]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [11]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [12]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [13]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [14]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [15]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [16]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [17]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [18]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [19]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [20]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [21]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [22]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [23]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [24]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [25]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [26]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [27]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [28]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [29]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [30]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [31]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [32]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [33]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [34]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [35]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [36]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [37]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [38]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [39]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [40]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [41]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [42]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [43]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [44]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [45]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [46]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [47]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [48]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [49]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [50]
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
    ARROWHEADSTYLE = FILLED
    COLOR = RED
    LINETHICKNESS = 0.6
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
$!FIELD  [51]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [52]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [53]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [54]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [55]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [56]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [57]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [58]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [59]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [60]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [61]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [62]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [63]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [64]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [65]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [66]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [67]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [68]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [69]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [70]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [71]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [72]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [73]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [74]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [75]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [76]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [77]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [78]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [79]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [80]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [81]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [82]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [83]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [84]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [85]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [86]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [87]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [88]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [89]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [90]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [91]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [92]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [93]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [94]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [95]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [96]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [97]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [98]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [99]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
$!FIELD  [100]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLACK
    LINETHICKNESS = 0.02
    }
  POINTS
    {
    POINTSTOPLOT = SURFACENODES
    IJKSKIP
      {
      I = 2
      J = 2
      }
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
  GROUP = 2
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
    RANGEMIN = -0.5
    RANGEMAX = 1.5
    GRSPACING = 0.5
    }
$!TWODAXIS 
  YDETAIL
    {
    RANGEMIN = -0.10000000000000000555
    RANGEMAX = 1.3083537533863498492
    GRSPACING = 0.5
    }
$!GLOBALISOSURFACE 
  MARCHINGCUBEALGORITHM = CLASSICPLUS
$!GLOBALSLICE 
  BOUNDARY
    {
    SHOW = NO
    }
$!FIELDLAYERS 
  SHOWMESH = NO
  SHOWVECTOR = YES

### Frame Number 2 ###
$!CREATENEWFRAME 
$!FRAMELAYOUT 
  SHOWBORDER = NO
  SHOWHEADER = NO
  ISTRANSPARENT = YES
  HEADERCOLOR = RED
  XYPOS
    {
    X = 0.125
    Y = 0.20954
    }
  WIDTH = 11.443
  HEIGHT = 7.8486
$!PLOTTYPE  = CARTESIAN3D
$!FRAMENAME  = 'Frame 001' 
$!ACTIVEFIELDZONES  =  [1-100]
$!GLOBALRGB 
  RANGEMIN = 0
  RANGEMAX = 1
$!GLOBALCONTOUR  1
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
$!GLOBALTHREEDVECTOR 
  UVAR = 1
  VVAR = 2
  WVAR = 3
$!GLOBALTWODVECTOR 
  UVAR = 3
  VVAR = 4
  RELATIVELENGTH = 0.4
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
$!FIELD  [1]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [2]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [3]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [4]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [5]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [6]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [7]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [8]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [9]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [10]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [11]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [12]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [13]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [14]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [15]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [16]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [17]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [18]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [19]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [20]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [21]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [22]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [23]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [24]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [25]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [26]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [27]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [28]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [29]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [30]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [31]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [32]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [33]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [34]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [35]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [36]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [37]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [38]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [39]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [40]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [41]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [42]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [43]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [44]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [45]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [46]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [47]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [48]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [49]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [50]
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
    SHOW = NO
    COLOR = RED
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = RED
    }
  SHADE
    {
    COLOR = CUSTOM10
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = RED
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
$!FIELD  [51]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [52]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [53]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [54]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [55]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [56]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [57]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [58]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [59]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [60]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [61]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [62]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [63]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [64]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [65]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [66]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [67]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [68]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [69]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [70]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [71]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [72]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [73]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [74]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [75]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [76]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [77]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [78]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [79]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [80]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [81]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [82]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [83]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [84]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [85]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [86]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [87]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [88]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [89]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [90]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [91]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [92]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [93]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [94]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [95]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [96]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [97]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [98]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [99]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!FIELD  [100]
  MESH
    {
    COLOR = GREEN
    }
  CONTOUR
    {
    CONTOURTYPE = FLOOD
    COLOR = GREEN
    USELIGHTINGEFFECT = YES
    }
  VECTOR
    {
    SHOW = NO
    COLOR = BLUE
    LINETHICKNESS = 0.4
    }
  SCATTER
    {
    COLOR = GREEN
    }
  SHADE
    {
    COLOR = CUSTOM34
    }
  BOUNDARY
    {
    SHOW = YES
    COLOR = BLUE
    LINETHICKNESS = 0.8
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
  GROUP = 2
$!THREEDAXIS 
  XDETAIL
    {
    VARNUM = 1
    }
  YDETAIL
    {
    VARNUM = 2
    }
  ZDETAIL
    {
    VARNUM = 3
    }
$!VIEW FIT
$!THREEDAXIS 
  AXISMODE = XYDEPENDENT
  XYDEPXTOYRATIO = 1
  DEPXTOYRATIO = 1
  DEPXTOZRATIO = 1
  FRAMEAXIS
    {
    SHOW = NO
    }
$!THREEDAXIS 
  XDETAIL
    {
    SHOWAXIS = NO
    RANGEMIN = -0.10000000000000000555
    RANGEMAX = 1.3000000000000000444
    GRSPACING = 0.5
    AXISLINE
      {
      EDGE = 3
      }
    }
$!THREEDAXIS 
  YDETAIL
    {
    SHOWAXIS = NO
    RANGEMIN = -0.10000000000000000555
    RANGEMAX = 1.3000000000000000444
    GRSPACING = 0.2
    AXISLINE
      {
      EDGE = 3
      }
    }
$!THREEDAXIS 
  ZDETAIL
    {
    SHOWAXIS = NO
    RANGEMIN = -1
    RANGEMAX = 1
    GRSPACING = 0.5
    GRIDLINES
      {
      SHOW = NO
      }
    AXISLINE
      {
      EDGE = 2
      }
    }
$!GLOBALISOSURFACE 
  MARCHINGCUBEALGORITHM = CLASSICPLUS
$!GLOBALSLICE 
  BOUNDARY
    {
    SHOW = NO
    }
$!GLOBALTHREED 
  AXISSCALEFACT
    {
    X = 1.09149543728
    Y = 1.09149543728
    Z = 0.5423
    }
  ROTATEORIGIN
    {
    X = 1
    Y = 0.6
    Z = 0
    }
$!THREEDVIEW 
  PSIANGLE = 90
  THETAANGLE = 90
  ALPHAANGLE = -90
  VIEWERPOSITION
    {
    X = -8.24
    Y = 0.604038410778
    Z = -0.0298791329236
    }
  VIEWWIDTH = 2.87992
$!FIELDLAYERS 
  SHOWMESH = NO
$!ATTACHTEXT 
  ANCHORPOS
    {
    X = 77.5737914277
    Y = 58.5725919455
    }
  COLOR = RED
  TEXTSHAPE
    {
    HEIGHT = 16
    }
  BOX
    {
    MARGIN = 10
    LINETHICKNESS = 0.4
    }
  SCOPE = GLOBAL
  MACROFUNCTIONCOMMAND = '' 
  TEXT = 'Computed Solution' 
$!ATTACHTEXT 
  ANCHORPOS
    {
    X = 77.5737914277
    Y = 62.1453546645
    }
  COLOR = BLUE
  TEXTSHAPE
    {
    HEIGHT = 16
    }
  BOX
    {
    MARGIN = 10
    LINETHICKNESS = 0.4
    }
  SCOPE = GLOBAL
  MACROFUNCTIONCOMMAND = '' 
  TEXT = 'Exact Solution' 
$!SETSTYLEBASE CONFIG



############################







$!IF |PNG|==1


        $!EXPORTSETUP EXPORTFORMAT = PNG
        $!EXPORTSETUP SUPERSAMPLEFACTOR = 3
        $!EXPORTSETUP USESUPERSAMPLEANTIALIASING = YES
        $!EXPORTSETUP IMAGEWIDTH = 600
        $!EXPORTSETUP EXPORTFNAME = 'velocity_vectors_CR|istep|.png'
        $!EXPORT
          EXPORTREGION = CURRENTFRAME

        $!EXPORTSETUP EXPORTFORMAT = EPS
        $!EXPORTSETUP USESUPERSAMPLEANTIALIASING = NO
        $!EXPORTSETUP IMAGEWIDTH = 1423
        $!EXPORTSETUP EXPORTFNAME = 'velocity_vectors_CR|istep|.eps'

        $!EXPORTSETUP PRINTRENDERTYPE = IMAGE
        $!EXPORTSETUP EXPORTFNAME = 'velocity_vectors_CR|istep|.img.eps'
        $!EXPORT
          EXPORTREGION = CURRENTFRAME

$!ELSE

        $!IF |LOOP|>1
                $!EXPORTNEXTFRAME
        $!ELSE

                $!EXPORTSETUP
                 EXPORTFORMAT = AVI
                 EXPORTFNAME = "velocity_vectors_CR.avi"
                $!EXPORTSETUP IMAGEWIDTH = 829
                $!EXPORTSETUP USESUPERSAMPLEANTIALIASING = NO
                $!EXPORTSTART
        $!ENDIF

$!ENDIF


$!VARSET |LAST_STEP|=|istep|

$!EndLoop

$!FIELDLAYERS SHOWBOUNDARY = NO
$!REDRAWALL





$!IF |PNG|==0
        $!EXPORTFINISH
$!ENDIF


$!QUIT

