try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

# Load the state
servermanager.LoadState("fpsi.pvsm")

#Turn background black for all windows
RenderView1 = GetRenderView()
RenderView2 = GetRenderViews()[1]
RenderView3 = GetRenderViews()[2]
RenderView4 = GetRenderViews()[3]
RenderView5 = GetRenderViews()[4]
#RenderView6 = GetRenderViews()[5]
#RenderView7 = GetRenderViews()[6]
RenderView1.Background = [0.0, 0.0, 0.0]
RenderView2.Background = [0.0, 0.0, 0.0]
RenderView3.Background = [0.0, 0.0, 0.0]
RenderView4.Background = [0.0, 0.0, 0.0]
RenderView5.Background = [0.0, 0.0, 0.0]
#RenderView6.Background = [0.0, 0.0, 0.0]
#RenderView7.Background = [0.0, 0.0, 0.0]


# Make sure that the view in the state is the active one
SetActiveView(GetRenderView())

# Now render
Render()

# Write animation
WriteAnimation('animation.png', Magnification=1, Quality=2, FrameRate=1.000000)
