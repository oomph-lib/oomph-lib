/**

\mainpage The spatially-adaptive solution of the azimuthally Fourier-decomposed equations of 3D time-harmonic linear elasticity on unstructured meshes

In this tutorial we re-visit the solution of the
time-harmonic equations of 3D linear elasticity in cylindrical polar
coordinates, using a Fourier decomposition of the solution in the 
azimuthal direction. The driver code is very similar to the one discussed in 
<a href="../../cylinder/html/index.html">another tutorial
</a> -- the main purpose of the current tutorial is to demonstrate
the use of spatial adaptivity on unstructured meshes. Compared to
the test case considered in the <a href="../../cylinder/html/index.html">
other tutorial</a> we study a slightly less contrived test problem:
the forced time-harmonic oscillations of a finite-length, hollow cylinder,
loaded by a time-periodic pressure load on its inner surface.



<CENTER>
<TABLE BORDER=1, WIDTH=500px>
<TR>
<TD bgcolor="cornsilk">
<CENTER>
<B>Acknowledgement:</B>
This implementation of the equations and the documentation 
were developed jointly with Robert Harter (Thales Underwater Systems
Ltd) with financial support from a KTA Secondment grant from
University of Manchester's EPSRC-funded Knowledge Transfer Account.
</CENTER>
</TD>
</TR>
</TABLE>
</CENTER>

<HR>
<HR>
 
\section test The test problem
The figure below shows the problem considered in this tutorial:
an annular elastic body that occupies the region \f$ r_{\rm min}\leq
r\leq r_{\rm max}, z_{\rm min}\leq
z\leq z_{\rm max} \f$ is loaded by 
a time-harmonic pressure load acting on its inner surface 
(at \f$ r = r_{\rm min} \f$). The upper and lower ends of the
hollow cylinder (at \f$ z = z_{\rm min} \f$ and \f$ z = z_{\rm min} \f$)
are held at a fixed position.

Here is an animation of the resulting displacement field
for \f$ r_{\rm min} = 0.1 \f$ and  \f$ r_{\rm max} = 1.1. \f$


\image html anim.gif "Forced oscillations of a thick-walled, hollow cylinder, subject to a pressure load on its inner surface. The pink shape in the background shows the cylinder's undeformed shape (in a radial plane); the mesh plotted in the region r < 0 illustrates how spatial adaptivity refines the mesh in regions of sharp displacement gradients (near the loaded surface and the supports). " 
\image latex anim.eps "Forced oscillations of a thick-walled, hollow cylinder, subject to a pressure load on its inner surface. The pink shape in the background shows the cylinder's undeformed shape (in a radial plane); the mesh plotted in the region r < 0 illustrates how spatial adaptivity refines the mesh in regions of sharp displacement gradients (near the loaded surface and the supports). " width=0.6\textwidth

<HR>
<HR>
  
\section num The numerical solution
The driver code for this problem is very similar to the one discussed in
<a href=../../cylinder/html/index.html>another tutorial</a>.
Running \c sdiff on the two driver codes
<CENTER>
<A HREF="../../../../demo_drivers/time_harmonic_fourier_decomposed_linear_elasticity/cylinder/cylinder.cc">
demo_drivers/time_harmonic_fourier_decomposed_linear_elasticity/cylinder/cylinder.cc
</A>
</CENTER>
and
<CENTER>
<A HREF="../../../../demo_drivers/time_harmonic_fourier_decomposed_linear_elasticity/cylinder/pressure_loaded_cylinder.cc">
demo_drivers/time_harmonic_fourier_decomposed_linear_elasticity/cylinder/pressure_loaded_cylinder.cc
</A>
</CENTER>
shows you the differences, the most important of which are:
- The change of the forcing to a spatially constant 
  pressure load on the inside boundary.
- The provision of the \c actions_before/after_adapt() functions and a 
  helper function \c complete_problem_setup()
  which rebuilds the elements (by passing the problem parameters
  to the elements) following the unstructured mesh adaptation.
  (The need/rationale for such a function is discussed in 
  <a href="../../../meshes/mesh_from_inline_triangle/html/index.html">
  another tutorial.</a>)
- The mesh generation and the application of boundary conditions at
  the upper and lower boundaries of the hollow cylinder.
.  
All of this is reasonably straightforward and provides a
powerful code that automatically adapts the mesh in regions of 
large displacement gradients. Have a look through the driver code 
and play with it.
 
<HR>
<HR> 

\section code Code listing
Here's a listing of the complete driver code:

\include pressure_loaded_cylinder.cc



<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/time_harmonic_fourier_decomposed_linear_elasticity/cylinder/">
demo_drivers/time_harmonic_fourier_decomposed_linear_elasticity/cylinder/
</A>
</CENTER>
- The driver code is:
<CENTER>
<A HREF="../../../../demo_drivers/time_harmonic_fourier_decomposed_linear_elasticity/cylinder/pressure_loaded_cylinder.cc">
demo_drivers/time_harmonic_fourier_decomposed_linear_elasticity/cylinder/pressure_loaded_cylinder.cc
</A>
</CENTER>
.




<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

