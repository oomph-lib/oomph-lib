/**

\mainpage Example problem: Steady flow in a curved tube 

 The problem of steady flow in a curved tube is 
 considered with a prescribed Poiseuille flow at the inlet 
 and a traction-free outlet condition. It is
 not clear that the latter is appropriate, but the main aim of this example is
 to check that the \c TubeMesh works correctly. 

 A detailed comparison between the flow field and the Dean solution should
 be performed for validation purposes, but the qualitative features seem
 reasonable.
 
\image html curve_geom.gif "Sketch of the problem with pressure contours. " 
\image latex curve_geom.eps "Sketch of the problem with pressure contours. " width=0.75\textwidth

\image html curve_sec.gif "Contours of axial velocity and secondary streamlines. " 
\image latex curve_sec.eps "Contours of axial velocity and secondary streamlines. " width=0.75\textwidth

Detailed documentation to be written. Here's the driver code...


\include curved_pipe.cc



<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

