#!MC 1000
$!DRAWGRAPHICS FALSE
$!READDATASET  '"RESLT1/neighbours0.dat" ' 
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  INITIALPLOTTYPE = CARTESIAN3D

$!FIELDLAYERS SHOWMESH = NO
$!FIELDLAYERS SHOWSCATTER = YES

$!FIELDLAYERS SHOWSHADE = YES
$!FIELDLAYERS USETRANSLUCENCY = YES
$!FIELD [1-|NUMZONES|]  SURFACEEFFECTS{SURFACETRANSLUCENCY = 95}
$!FIELD [1]  SHADE{SHOW = NO}

$!GLOBALSCATTER VAR = 4
$!FIELD [1-|NUMZONES|]  SCATTER{SIZEBYVARIABLE = YES}
$!FIELD [1-|NUMZONES|]  SCATTER{FILLMODE = USELINECOLOR}

$!THREEDAXIS GRIDAREA{ISFILLED = NO}


# switch off all zones
$!ACTIVEFIELDZONES -= [2-|NUMZONES|]

$!GETUSERINPUT |nelements| INSTRUCTIONS = "Loop. Number of elements??"


$!VARSET |firstdone|=0

# initialise zone number for first element
$!VARSET |centralelement|=-12

$!FIELD [1]  BOUNDARY{SHOW = YES}

#=======================
# LOOP OVER ALL ELEMENTS
#=======================
$!LOOP |nelements|



$!VARSET |centralelement|+=13
$!ACTIVEFIELDZONES += [|centralelement|]
$!FIELD [|centralelement|]  SCATTER{SHOW = NO}


$!THREEDAXIS ZDETAIL{TITLE{TITLEMODE = USETEXT}}
$!THREEDAXIS YDETAIL{TITLE{TITLEMODE = USETEXT}}
$!THREEDAXIS XDETAIL{TITLE{TITLEMODE = USETEXT}}
$!THREEDAXIS YDETAIL{TITLE{TEXT = 'down  <---> up '}}
$!THREEDAXIS XDETAIL{TITLE{TEXT = 'left  <---> right '}}
$!THREEDAXIS ZDETAIL{TITLE{TEXT = 'back <---> front '}}

$!THREEDVIEW
  PSIANGLE = 54.9806
  THETAANGLE = -40.5871
  ALPHAANGLE = 0
  VIEWERPOSITION
    {
    X = 20.9031931957
    Y = -20.8190051758
    Z = 21.2610787102
    }


$!VIEW PUSH


#=======================
# NEIGHBOUR1
#=======================
$!VARSET |neighbourzone|=|centralelement|
$!VARSET |neighbourzone|+=1
$!ACTIVEFIELDZONES += [|neighbourzone|]
$!FIELD [|neighbourzone|]  SCATTER{SHOW = NO}


$!ATTACHTEXT 
  ANCHORPOS
    {
    X = 7.52364676016
    Y = 70.3414869946
    }
  TEXT = 'Element: |LOOP| \nNeighbour: Left\n' 



#=======================
# MATCHING POINT
#=======================
$!VARSET |matchzone|=|centralelement|
$!VARSET |matchzone|+=2
$!ACTIVEFIELDZONES += [|matchzone|]
$!FIELD [|matchzone|]  BOUNDARY{SHOW = NO}

$!DRAWGRAPHICS TRUE
$!VIEW FIT
$!REDRAW

$!IF  |firstdone|==1
  $!EXPORTNEXTFRAME
$!ENDIF
$!IF |firstdone|==0
  $!EXPORTSETUP EXPORTFORMAT = AVI
  $!EXPORTSETUP BITDUMPREGION = ALLFRAMES
  $!EXPORTSETUP IMAGEWIDTH = 829
  $!EXPORTSETUP EXPORTFNAME = 'neighbours_checking.avi'
  $!EXPORTSTART
  $!VARSET |firstdone|=1
$!ENDIF

$!PICK ADDALL
  SELECTTEXT = YES
$!PICK CUT


#=======================
# FILL IN THE REMAINING 5 NEIGHBOURS AND MATCHING POINTS 
#=======================

#=======================
# NEIGHBOUR2
#=======================
$!VARSET |neighbourzone|=|centralelement|
$!VARSET |neighbourzone|+=3
$!ACTIVEFIELDZONES += [|neighbourzone|]
$!FIELD [|neighbourzone|]  SCATTER{SHOW = NO}
$!FIELD [|matchzone|]  SCATTER{SHOW = NO}

$!ATTACHTEXT 
  ANCHORPOS
    {
    X = 7.52364676016
    Y = 70.3414869946
    }
  TEXT = 'Element: |LOOP| \nNeighbour: Right\n'


#=======================
# MATCHING POINT
#=======================
$!VARSET |matchzone|=|centralelement|
$!VARSET |matchzone|+=4
$!ACTIVEFIELDZONES += [|matchzone|]
$!FIELD [|matchzone|]  BOUNDARY{SHOW = NO}

$!DRAWGRAPHICS TRUE
$!VIEW FIT
$!REDRAW

$!IF  |firstdone|==1
  $!EXPORTNEXTFRAME
$!ENDIF
$!IF |firstdone|==0
  $!EXPORTSETUP EXPORTFORMAT = AVI
  $!EXPORTSETUP BITDUMPREGION = ALLFRAMES
  $!EXPORTSETUP IMAGEWIDTH = 829
  $!EXPORTSETUP EXPORTFNAME = 'neighbours_checking.avi'
  $!EXPORTSTART
  $!VARSET |firstdone|=1
$!ENDIF

$!PICK ADDALL
  SELECTTEXT = YES
$!PICK CUT


#=======================
# NEIGHBOUR3
#=======================
$!VARSET |neighbourzone|=|centralelement|
$!VARSET |neighbourzone|+=5
$!ACTIVEFIELDZONES += [|neighbourzone|]
$!FIELD [|neighbourzone|]  SCATTER{SHOW = NO}
$!FIELD [|matchzone|]  SCATTER{SHOW = NO}

$!ATTACHTEXT 
  ANCHORPOS
    {
    X = 7.52364676016
    Y = 70.3414869946
    }
  TEXT = 'Element: |LOOP| \nNeighbour: Down\n'

#=======================
# MATCHING POINT
#=======================
$!VARSET |matchzone|=|centralelement|
$!VARSET |matchzone|+=6
$!ACTIVEFIELDZONES += [|matchzone|]
$!FIELD [|matchzone|]  BOUNDARY{SHOW = NO}

$!DRAWGRAPHICS TRUE
$!VIEW FIT
$!REDRAW

$!IF  |firstdone|==1
  $!EXPORTNEXTFRAME
$!ENDIF
$!IF |firstdone|==0
  $!EXPORTSETUP EXPORTFORMAT = AVI
  $!EXPORTSETUP BITDUMPREGION = ALLFRAMES
  $!EXPORTSETUP IMAGEWIDTH = 829
  $!EXPORTSETUP EXPORTFNAME = 'neighbours_checking.avi'
  $!EXPORTSTART
  $!VARSET |firstdone|=1
$!ENDIF

$!PICK ADDALL
  SELECTTEXT = YES
$!PICK CUT


#=======================
# NEIGHBOUR4
#=======================
$!VARSET |neighbourzone|=|centralelement|
$!VARSET |neighbourzone|+=7
$!ACTIVEFIELDZONES += [|neighbourzone|]
$!FIELD [|neighbourzone|]  SCATTER{SHOW = NO}
$!FIELD [|matchzone|]  SCATTER{SHOW = NO}

$!ATTACHTEXT 
  ANCHORPOS
    {
    X = 7.52364676016
    Y = 70.3414869946
    }
  TEXT = 'Element: |LOOP| \nNeighbour: Up\n'


#=======================
# MATCHING POINT
#=======================
$!VARSET |matchzone|=|centralelement|
$!VARSET |matchzone|+=8
$!ACTIVEFIELDZONES += [|matchzone|]
$!FIELD [|matchzone|]  BOUNDARY{SHOW = NO}

$!DRAWGRAPHICS TRUE
$!VIEW FIT
$!REDRAW

$!IF  |firstdone|==1
  $!EXPORTNEXTFRAME
$!ENDIF
$!IF |firstdone|==0
  $!EXPORTSETUP EXPORTFORMAT = AVI
  $!EXPORTSETUP BITDUMPREGION = ALLFRAMES
  $!EXPORTSETUP IMAGEWIDTH = 829
  $!EXPORTSETUP EXPORTFNAME = 'neighbours_checking.avi'
  $!EXPORTSTART
  $!VARSET |firstdone|=1
$!ENDIF

$!PICK ADDALL
  SELECTTEXT = YES
$!PICK CUT


#=======================
# NEIGHBOUR5
#=======================
$!VARSET |neighbourzone|=|centralelement|
$!VARSET |neighbourzone|+=9
$!ACTIVEFIELDZONES += [|neighbourzone|]
$!FIELD [|neighbourzone|]  SCATTER{SHOW = NO}
$!FIELD [|matchzone|]  SCATTER{SHOW = NO}

$!ATTACHTEXT 
  ANCHORPOS
    {
    X = 7.52364676016
    Y = 70.3414869946
    }
  TEXT = 'Element: |LOOP| \nNeighbour: Back\n'



#=======================
# MATCHING POINT
#=======================
$!VARSET |matchzone|=|centralelement|
$!VARSET |matchzone|+=10
$!ACTIVEFIELDZONES += [|matchzone|]
$!FIELD [|matchzone|]  BOUNDARY{SHOW = NO}

$!DRAWGRAPHICS TRUE
$!VIEW FIT
$!REDRAW

$!IF  |firstdone|==1
  $!EXPORTNEXTFRAME
$!ENDIF
$!IF |firstdone|==0
  $!EXPORTSETUP EXPORTFORMAT = AVI
  $!EXPORTSETUP BITDUMPREGION = ALLFRAMES
  $!EXPORTSETUP IMAGEWIDTH = 829
  $!EXPORTSETUP EXPORTFNAME = 'neighbours_checking.avi'
  $!EXPORTSTART
  $!VARSET |firstdone|=1
$!ENDIF

$!PICK ADDALL
  SELECTTEXT = YES
$!PICK CUT



#=======================
# NEIGHBOUR6
#=======================
$!VARSET |neighbourzone|=|centralelement|
$!VARSET |neighbourzone|+=11
$!ACTIVEFIELDZONES += [|neighbourzone|]
$!FIELD [|neighbourzone|]  SCATTER{SHOW = NO}
$!FIELD [|matchzone|]  SCATTER{SHOW = NO}

$!ATTACHTEXT 
  ANCHORPOS
    {
    X = 7.52364676016
    Y = 70.3414869946
    }
  TEXT = 'Element: |LOOP| \nNeighbour: Front\n'


#=======================
# MATCHING POINT
#=======================
$!VARSET |matchzone|=|centralelement|
$!VARSET |matchzone|+=12
$!ACTIVEFIELDZONES += [|matchzone|]
$!FIELD [|matchzone|]  BOUNDARY{SHOW = NO}

$!DRAWGRAPHICS TRUE
$!VIEW FIT
$!REDRAW

$!IF  |firstdone|==1
  $!EXPORTNEXTFRAME
$!ENDIF
$!IF |firstdone|==0
  $!EXPORTSETUP EXPORTFORMAT = AVI
  $!EXPORTSETUP BITDUMPREGION = ALLFRAMES
  $!EXPORTSETUP IMAGEWIDTH = 829
  $!EXPORTSETUP EXPORTFNAME = 'neighbours_checking.avi'
  $!EXPORTSTART
  $!VARSET |firstdone|=1
$!ENDIF

$!PICK ADDALL
  SELECTTEXT = YES
$!PICK CUT


#$!ACTIVEFIELDZONES -= [|neighbourzone|-|matchzone|]


#====================================================
# FROM NOW ON DON'T SHOW FIRST ELEMENT ANY MORE 
# BUT KEEP IT ALIVE SO THERE'S ALWAYS AT LEAST ONE
# ELEMENT THERE...
#====================================================
$!FIELD [1]  BOUNDARY{SHOW = NO}


$!IF |LOOP|<> |nelements|
        $!ACTIVEFIELDZONES -= [2-|NUMZONES|]
$!ENDIF



#=======================
# END OF LOOP OVER ALL ELEMENTS
#=======================
$!ENDLOOP



$!EXPORTFINISH 

