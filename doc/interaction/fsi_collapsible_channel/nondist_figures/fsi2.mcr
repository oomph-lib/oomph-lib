#!MC 1100

#=========================================================
# Produce animation that shows the wall degrees of 
# freedom that affect a given fluid node
#=========================================================


# Use png output (otherwise avi)
$!VARSET |PNG|=0


#################################################
# Overall step
#################################################
$!VARSET |STEP|=3


#$!GETUSERINPUT |lostep| INSTRUCTIONS = "First Step/wall element?"
$!VARSET |lostep|=100
#$!GETUSERINPUT |dlstep| INSTRUCTIONS = "Step/wall element number increment?"
$!VARSET |dlstep|=1
#$!GETUSERINPUT |nstep| INSTRUCTIONS = "Number of steps?"
$!VARSET |nstep|=4 #20



$!LOOP |nstep|

$!VarSet |nnstep| = |LOOP|
$!VarSet |nnstep| -= 1
$!VarSet |iistep| = |dlstep|
$!VarSet |iistep| *= |nnstep|
$!VarSet |iistep| += |lostep|
$!NEWLAYOUT

$!VARSET |istep|=|iistep|



$!DRAWGRAPHICS FALSE
$!READDATASET  '"RESLT/fsi_doc_fluid_mesh|STEP|.dat" "RESLT/fsi_doc_fluid_element|STEP|-|istep|.dat" ' 
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYPOSITION
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN2D
$!FIELDLAYERS SHOWMESH = NO
$!TWODAXIS XDETAIL{AXISLINE{SHOW = NO}}
$!TWODAXIS XDETAIL{AXISLINE{SHOW = YES}}
$!TWODAXIS AXISMODE = INDEPENDENT
$!VIEW FIT
$!FIELDLAYERS SHOWSCATTER = YES
$!FIELDMAP [1-|NUMZONES|]  SCATTER{SHOW = NO}
$!FIELDMAP [|NUMZONES|]  SCATTER{SHOW = YES}
$!REDRAWALL 
$!FIELDMAP [|NUMZONES|]  SCATTER{COLOR = RED}
$!FIELDMAP [|NUMZONES|]  SCATTER{FILLMODE = USELINECOLOR}
$!FIELDMAP [1-|NUMZONES|]  EDGELAYER{COLOR = BLUE}
$!FIELDMAP [1-|NUMZONES|]  EDGELAYER{LINETHICKNESS = 0.02}
$!REDRAWALL 

$!TWODAXIS XDETAIL{TITLE{TITLEMODE = USETEXT}}
$!TWODAXIS XDETAIL{TITLE{TEXT = 'x'}}
$!TWODAXIS YDETAIL{TITLE{TITLEMODE = USETEXT}}
$!TWODAXIS YDETAIL{TITLE{TEXT = 'y'}}
$!TWODAXIS GRIDAREA{DRAWBORDER = YES}
$!VIEW FIT
$!TWODAXIS XDETAIL{RANGEMIN = -1}
$!TWODAXIS XDETAIL{RANGEMAX = 26}
$!TWODAXIS YDETAIL{RANGEMIN = -0.1}
$!TWODAXIS YDETAIL{RANGEMAX = 1.1}

$!DRAWGRAPHICS TRUE
$!REDRAWALL



$!IF |PNG|==1


        $!EXPORTSETUP EXPORTFORMAT = PNG
        $!EXPORTSETUP IMAGEWIDTH = 750
#        $!EXPORTSETUP USESUPERSAMPLEANTIALIASING = YES
        $!EXPORTSETUP EXPORTFNAME = 'fsi2-|LOOP|.png'
        $!EXPORT
          EXPORTREGION = ALLFRAMES

#        $!EXPORTSETUP EXPORTFORMAT = EPS
#        $!EXPORTSETUP IMAGEWIDTH = 1423
#        $!EXPORTSETUP EXPORTFNAME = 'fsi|LOOP|.eps'
#        $!EXPORT
#          EXPORTREGION = ALLFRAMES

$!ELSE

        $!IF |LOOP|>1
                $!EXPORTNEXTFRAME
        $!ELSE

                $!EXPORTSETUP
                  EXPORTREGION = ALLFRAMES
                  EXPORTFORMAT = AVI
                  EXPORTFNAME = "fsi2.avi"
                $!EXPORTSETUP IMAGEWIDTH = 750
#                $!EXPORTSETUP USESUPERSAMPLEANTIALIASING = YES
                $!EXPORTSTART
        $!ENDIF

$!ENDIF


$!ENDLOOP


$!IF |PNG|==0
        $!EXPORTFINISH
$!ENDIF

# $!QUIT
