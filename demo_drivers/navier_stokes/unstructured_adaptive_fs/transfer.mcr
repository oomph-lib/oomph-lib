#!MC 1120


$!VarSet |dir| = "RESLT"
$!VarSet |orig| = 79
$!VarSet |next|=(|orig|+1)



$!READDATASET  '"|dir|/soln|orig|.dat" '
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"x" "y" "u" "v" "p" "du/dt" "dv/dt" "u_m" "v_m" "x_h1" "y_h1" "x_h2" "y_h2" "u_h1" "v_h1" "u_h2" "v_h2" "du/dx" "du/dy" "dv/dx" "dv/dy" "error" "size"'

$!VARSET |FIRST|=1
$!VARSET |LAST|=|NUMZONES|

$!FIELDLAYERS SHOWMESH = NO
$!FIELDMAP [|FIRST|-|LAST|]  EDGELAYER{COLOR = RED}
$!FIELDMAP [|FIRST|-|LAST|]  EDGELAYER{LINETHICKNESS = 0.0200000000000000004}
$!REDRAWALL 




$!VARSET |FIRST|=(|NUMZONES|+1)
$!READDATASET  '"|dir|/soln|orig|_adapt.dat" '
  READDATAOPTION = APPEND
  RESETSTYLE = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"x" "y" "u" "v" "p" "du/dt" "dv/dt" "u_m" "v_m" "x_h1" "y_h1" "x_h2" "y_h2" "u_h1" "v_h1" "u_h2" "v_h2" "du/dx" "du/dy" "dv/dx" "dv/dy" "error" "size"'
$!VARSET |LAST|=|NUMZONES|

$!ACTIVEFIELDMAPS += [|FIRST|-|LAST|]
$!FIELDMAP [|FIRST|-|LAST|]  GROUP = 2
$!FIELDMAP [|FIRST|-|LAST|]  EDGELAYER{COLOR = BLUE}
$!FIELDMAP [|FIRST|-|LAST|]  EDGELAYER{LINETHICKNESS = 0.0200000000000000004}
$!REDRAWALL 

$!VARSET |FIRST|=(|NUMZONES|+1)
$!READDATASET  '"|dir|/soln|next|.dat" '
  READDATAOPTION = APPEND
  RESETSTYLE = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN3D
  VARNAMELIST = '"x" "y" "u" "v" "p" "du/dt" "dv/dt" "u_m" "v_m" "x_h1" "y_h1" "x_h2" "y_h2" "u_h1" "v_h1" "u_h2" "v_h2" "du/dx" "du/dy" "dv/dx" "dv/dy" "error" "size"'
$!VARSET |LAST|=|NUMZONES|

$!ACTIVEFIELDMAPS += [|FIRST|-|LAST|]
$!FIELDMAP [|FIRST|-|LAST|]  GROUP = 3
$!FIELDMAP [|FIRST|-|LAST|]  EDGELAYER{COLOR = GREEN}
$!REDRAWALL 
$!FIELDMAP [|FIRST|-|LAST|]  EDGELAYER{LINETHICKNESS = 0.100000000000000006}
$!FIELDMAP [|FIRST|-|LAST|]  EDGELAYER{LINETHICKNESS = 0.0200000000000000004}

$!THREEDAXIS ZDETAIL{VARNUM = 6}
$!REDRAWALL 
$!VIEW FIT
$!VIEW ZOOM
  X1 = -0.66815024235
  Y1 = -0.259272486435
  X2 = 0.120903377187
  Y2 = 0.112351410789

