#!MC 1120

$!VARSET |file1|="Validation/RESLT/soln30_on_proc0.dat"
$!VARSET |file2|="Validation/RESLT/soln30_on_proc1.dat"

$!VARSET |file1|="Validation/RESLT_doubly_adaptive_load_balanced_for_restart/soln7_on_proc0.dat"
$!VARSET |file2|="Validation/RESLT_doubly_adaptive_load_balanced_restarted_from_load_balanced/soln7_on_proc0.dat"

$!READDATASET '|file1|'
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDECUSTOMLABELS = NO

$!FIELDLAYERS SHOWMESH = NO
$!FIELDMAP [1-|NUMZONES|]  EDGELAYER{COLOR = RED}

$!VARSET |next_first|=(|NUMZONES|+1)
$!READDATASET  '|file2|'
  READDATAOPTION = APPEND
  RESETSTYLE = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  ASSIGNSTRANDIDS = YES
  INITIALPLOTTYPE = CARTESIAN2D
$!ACTIVEFIELDMAPS += [|next_first|-|NUMZONES|]
$!FIELDMAP [|next_first|-|NUMZONES|]  EDGELAYER{COLOR = GREEN}
$!REDRAWALL 

