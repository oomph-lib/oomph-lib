try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

soln_pvd = PVDReader( FileName='RESLT/soln.pvd' )

AnimationScene1 = GetAnimationScene()
#AnimationScene1.EndTime = 50.0
AnimationScene1.PlayMode = 'Snap To TimeSteps'

RenderView1 = GetRenderView()
DataRepresentation1 = Show()
DataRepresentation1.ScaleFactor = 0.1
DataRepresentation1.ScalarOpacityUnitDistance = 0.19193831036664846
DataRepresentation1.SelectionPointFieldDataArrayName = 'V1'
DataRepresentation1.EdgeColor = [0.0, 0.0, 0.5000076295109483]

WarpByScalar1 = WarpByScalar()

RenderView1.CameraViewUp = [0.03948755513634052, 0.20355959016033767, 0.9782659281826747]
RenderView1.CameraPosition = [0.26859028408374336, -2.1632681672113696, 0.5635191461130061]
RenderView1.CameraClippingRange = [1.6505022867005104, 4.100087435032]
RenderView1.CameraFocalPoint = [0.49999999999999994, 0.49999999999999994, 1.0732287550786024e-17]
RenderView1.CameraParallelScale = 0.7071067811865476
RenderView1.CenterOfRotation = [0.5, 0.5, 0.0]

WarpByScalar1.Scalars = ['POINTS', 'V1']

DataRepresentation2 = Show()
DataRepresentation2.ScaleFactor = 0.1
DataRepresentation2.ScalarOpacityUnitDistance = 0.19193831036664846
DataRepresentation2.SelectionPointFieldDataArrayName = 'V1'
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]

a1_V1_PVLookupTable = GetLookupTableForArray( "V1", 1, NanColor=[0.25, 0.0, 0.0], RGBPoints=[-0.10000000149011612, 0.23, 0.299, 0.754, -0.10000000149011612, 0.706, 0.016, 0.15], VectorMode='Magnitude', ColorSpace='Diverging', ScalarRangeInitialized=1.0 )

a1_V1_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )

AnimationScene1.AnimationTime = 1.0

RenderView1.CameraViewUp = [-0.26660789808928614, -0.28415692315242297, 0.9209642076112192]
RenderView1.CacheKey = 1.0
RenderView1.CameraPosition = [1.042955148812281, -2.098791141516002, -0.6446595413814319]
RenderView1.CameraClippingRange = [1.5253628356644837, 4.230934312248583]
RenderView1.ViewTime = 1.0
RenderView1.UseCache = 1
RenderView1.CameraFocalPoint = [0.4999999999999999, 0.49999999999999956, -2.210168370816743e-16]

DataRepresentation1.ScalarOpacityFunction = a1_V1_PiecewiseFunction
DataRepresentation1.ColorArrayName = 'V1'
DataRepresentation1.LookupTable = a1_V1_PVLookupTable
DataRepresentation1.Visibility = 1

DataRepresentation2.Representation = 'Surface With Edges'


RenderView1.CameraViewUp = [0.3047361523061217, 0.20642580945980588, 0.9297979687365014]
RenderView1.CameraPosition = [-1.6813012470711173, -1.3372511744291389, 0.7632237760387917]
RenderView1.CameraClippingRange = [1.467671158890028, 5.083109160489926]
RenderView1.UseCache = 0
RenderView1.CameraFocalPoint = [0.5, 0.5, -0.3595764935016632]
RenderView1.CameraParallelScale = 0.793281321271938
RenderView1.CenterOfRotation = [0.5, 0.5, -0.3595764935016632]


flux_pvd = PVDReader( FileName='RESLT/flux.pvd' )

RenderView1.CameraClippingRange = [1.467671158890028, 5.083109160489926]
RenderView1.UseCache = 0

DataRepresentation3 = Show()
DataRepresentation3.ScaleFactor = 0.1
DataRepresentation3.SelectionPointFieldDataArrayName = 'V1'
DataRepresentation3.EdgeColor = [0.0, 0.0, 0.5000076295109483]

Calculator1 = Calculator()

Calculator1.AttributeMode = 'point_data'

Calculator1.Function = '(V2*V3*iHat+V2*V4*jHat)'
Calculator1.ResultArrayName = 'melt'

DataRepresentation4 = Show()
DataRepresentation4.ScaleFactor = 0.1
DataRepresentation4.SelectionPointFieldDataArrayName = 'melt'
DataRepresentation4.EdgeColor = [0.0, 0.0, 0.5000076295109483]

Glyph1 = Glyph( GlyphType="Arrow", GlyphTransform="Transform2" )

DataRepresentation3.Visibility = 0

Glyph1.Scalars = ['POINTS', 'V1']
Glyph1.SetScaleFactor = 0.1
Glyph1.Vectors = ['POINTS', 'melt']
Glyph1.GlyphTransform = "Transform2"
Glyph1.GlyphType = "Arrow"

DataRepresentation5 = Show()
DataRepresentation5.ScaleFactor = 0.10000000000100001
DataRepresentation5.SelectionPointFieldDataArrayName = 'V1'
DataRepresentation5.EdgeColor = [0.0, 0.0, 0.5000076295109483]

RenderView1.CameraViewUp = [0.4471618071411335, 0.3323053002561054, 0.8304333240278687]
RenderView1.CameraPosition = [-1.4047564251721567, -1.1928495088321383, 1.3434827993321454]
RenderView1.CameraClippingRange = [1.4688565474321054, 5.081617000529558]
RenderView1.CameraFocalPoint = [0.5, 0.5, -0.3595764935016632]
RenderView1.UseCache = 0


DataRepresentation1.Opacity = 0.7000000000000001

my_representation4 = GetDisplayProperties(Glyph1)
a1_V2_PVLookupTable = GetLookupTableForArray( "V2", 1, NanColor=[0.25, 0.0, 0.0], RGBPoints=[-0.18805700540542603, 0.23, 0.299, 0.754, 2.5486900806427, 0.706, 0.016, 0.15], VectorMode='Magnitude', ColorSpace='Diverging', ScalarRangeInitialized=1.0 )

a1_V2_PiecewiseFunction = CreatePiecewiseFunction( Points=[0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0] )
my_representation4.ColorArrayName = ''
my_representation4.DiffuseColor = [1.0, 1.0, 0.0]

a1_V1_PVLookupTable = GetLookupTableForArray( "V1", 1, RGBPoints=[-1.179229974746704, 0.7058823529411765, 0.01568627450980392, 0.14901960784313725] )





Render()
