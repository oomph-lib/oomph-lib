#!MC 1100



#=========================================================
# Produce animation that shows the fluid (and solid!) degrees of 
# freedom that affect a given wall element. Also show
# the coincidence between the wall Gauss points and their
# equivalents in the fluid domain.
#=========================================================


# Use png output (otherwise avi)
$!VARSET |PNG|=0


#################################################
# Overall step
#################################################
$!VARSET |STEP|=3




#$!GETUSERINPUT |lostep| INSTRUCTIONS = "First Step/wall element?"
$!VARSET |lostep|=0
#$!GETUSERINPUT |dlstep| INSTRUCTIONS = "Step/wall element number increment?"
$!VARSET |dlstep|=1
#$!GETUSERINPUT |nstep| INSTRUCTIONS = "Number of steps?"
$!VARSET |nstep|=7


$!LOOP |nstep|

$!VarSet |nnstep| = |LOOP|
$!VarSet |nnstep| -= 1
$!VarSet |iistep| = |dlstep|
$!VarSet |iistep| *= |nnstep|
$!VarSet |iistep| += |lostep|
$!NEWLAYOUT

$!VARSET |istep|=|iistep|



$!DRAWGRAPHICS FALSE
$!READDATASET  '"RESLT/fsi_doc_fluid_mesh|STEP|.dat" "RESLT/fsi_doc_wall_element|STEP|-|istep|.dat" ' 

$!VARSET |LAST_BUT_ONE|=(|NUMZONES|-1)

  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYPOSITION
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN2D
$!FIELDLAYERS SHOWMESH = NO
$!FIELDMAP [1-|NUMZONES|]  EDGELAYER{COLOR = BLUE}
$!FIELDMAP [1-|NUMZONES|]  EDGELAYER{LINETHICKNESS = 0.02}
$!TWODAXIS AXISMODE = INDEPENDENT
$!FIELDLAYERS SHOWSCATTER = YES
$!FIELDMAP [1-|NUMZONES|]  SCATTER{SHOW = NO}
$!FIELDMAP [|LAST_BUT_ONE|-|NUMZONES|]  SCATTER{SHOW = YES}
$!FIELDMAP [|LAST_BUT_ONE|-|NUMZONES|]  SCATTER{SYMBOLSHAPE{GEOMSHAPE = DEL}}
$!FIELDMAP [|NUMZONES|]  SCATTER{SYMBOLSHAPE{GEOMSHAPE = GRAD}}
$!FIELDMAP [|LAST_BUT_ONE|-|NUMZONES|]  SCATTER{COLOR = RED}
$!FIELDMAP [|NUMZONES|]  SCATTER{COLOR = GREEN}
$!VIEW FIT

$!TWODAXIS YDETAIL{TITLE{TITLEMODE = USETEXT}}
$!TWODAXIS YDETAIL{TITLE{TEXT = 'y'}}
$!TWODAXIS XDETAIL{TITLE{TITLEMODE = USETEXT}}
$!TWODAXIS XDETAIL{TITLE{TEXT = 'x'}}
$!TWODAXIS GRIDAREA{DRAWBORDER = YES}

$!VIEW ZOOM
  X1 = 4.4054948597
  Y1 = 0.963508958095
  X2 = 15.5439226049
  Y2 = 1.00247269242

$!VIEW ZOOM
  X1 = 3.98032541608
  Y1 = 0.681881201582
  X2 = 15.7165917827
  Y2 = 1.03063714746


$!DRAWGRAPHICS TRUE
$!REDRAWALL



$!IF |PNG|==1


        $!EXPORTSETUP EXPORTFORMAT = PNG
        $!EXPORTSETUP IMAGEWIDTH = 600
        $!EXPORTSETUP USESUPERSAMPLEANTIALIASING = YES
        $!EXPORTSETUP EXPORTFNAME = 'fsi|LOOP|.png'
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
                  EXPORTFNAME = "fsi.avi"
                $!EXPORTSETUP IMAGEWIDTH = 750
#                $!EXPORTSETUP USESUPERSAMPLEANTIALIASING = YES
                $!EXPORTSTART
        $!ENDIF

$!ENDIF


$!ENDLOOP


$!IF |PNG|==0
        $!EXPORTFINISH
$!ENDIF

