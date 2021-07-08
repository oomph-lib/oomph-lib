/**

\mainpage The adaptive solution of the equations of time-harmonic linear elasticity on unstructured meshes

The aim of this tutorial is to demonstrate the adaptive solution of the
time-harmonic equations of linear elasticity in cartesian
coordinates on unstructured meshes.


  The driver code is very similar to the one presented in
<a href=../../elastic_annulus/html/index.html>another tutorial</a>
and we only discuss the changes necessary to deal with the generation
of the unstructured mesh and the assignment of different material
properties to different parts of the domain.

<HR>
<HR> 

\section test A test problem: Time-harmonic oscillations of an elastic annulus reinforced by a T-rib

The sketch below shows the problem setup: A 2D elastic annulus
which is reinforced with two T-ribs is 
subjected to a time-periodic pressure load of magnitude
\f[
{\bf t} = P ( \exp(\alpha(\varphi-\pi/4)^2)  + \exp(\alpha(\varphi-3\pi/4)^2) )
\f]
(where \f$ \varphi \f$ is the polar angle) along its outer edge.
The parameter \f$ \alpha \f$ controls the "sharpness" of the
pressure load. For  \f$ \alpha=0 \f$ we obtain a uniform, axisymmetric load;
the sketch below shows the pressure distribution (red vectors
indicating the traction) for \f$ \alpha = 1. \f$

\image html setup.gif "Sketch of the problem setup. " 
\image latex setup.eps "Sketch of the problem setup. " width=0.8\textwidth

The structure is symmetric and we only discretise the right half
(\f$ x_1 > 0 \f$) of the domain and apply symmetry conditions
(zero horizontal displacement) on the \f$ x_2-\f$ axis. 


<HR>
<HR>
 

\section results Results
The figure below shows an animation of the structure's time-harmonic
oscillation. The blue shaded region shows the shape of the
oscillating structure while the pink region shows its initial configuration.
The left half of the plot is used to show the (mirror image of
the) adaptive unstructured mesh on which the solution was computed:

\image html all.gif "Animation of the time-harmonic deformation. " 
\image latex all.eps "Animation of the time-harmonic deformation. " width=0.8\textwidth

This looks very pretty and shows that we can compute in non-trivial
geometries but should you believe the results? Here's an
attempt to convince you: If we make the rib much softer than the
annulus, the rib will not offer much structural resistance
and the annular region will deform as if the rib was not present.
If we then set \f$ \alpha = 0 \f$ we apply an axisymmetric forcing
onto the structure and would expect the resulting displacement
field (at least in the annular region) to be axisymmetric. 
For this case it is easy to find an analytical solution to the
problem. The next two figures show a comparison between
the analytical (green spheres) and computed solutions (shaded)
for the real part of the horizontal and vertical displacements,
respectively.


\image html validate_real.gif "Real part of the horizontal displacement, computed (shaded) and exact (spheres). " 
\image latex validate_real.eps "Real part of the horizontal displacement, computed (shaded) and exact (spheres). " width=0.8\textwidth

\image html validate_imag.gif "Real part of the vertical displacement, computed (shaded) and exact (spheres). " 
\image latex validate_imag.eps "Real part of the vertical displacement, computed (shaded) and exact (spheres). " width=0.8\textwidth

Convinced?


<HR> 
<HR>

\section num The numerical solution
The driver code for this problem is very similar to the one discussed in
<a href=../../elastic_annulus/html/index.html>another tutorial</a>.
Running \c sdiff on the two driver codes
<CENTER>
<A HREF="../../../../demo_drivers/time_harmonic_linear_elasticity/elastic_annulus/time_harmonic_elastic_annulus.cc">
demo_drivers/time_harmonic_linear_elasticity/elastic_annulus/time_harmonic_elastic_annulus.cc
</A>
</CENTER>
and
<CENTER>
<A HREF="../../../../demo_drivers/time_harmonic_linear_elasticity/elastic_annulus/unstructured_time_harmonic_elastic_annulus.cc">
demo_drivers/time_harmonic_linear_elasticity/elastic_annulus/unstructured_time_harmonic_elastic_annulus.cc
</A>
</CENTER>
shows you the differences, the most important of which are:
- The provision of multiple elasticity tensors and frequency
  parameters for the two different regions (the rib and the 
  annulus).
  \n\n
- The provision of a helper function \c complete_problem_setup()
  which rebuilds the elements (by passing the problem parameters
  to the elements) following the unstructured mesh adaptation.
  (The need/rationale for such a function is discussed in 
  <a href="../../../meshes/mesh_from_inline_triangle/html/index.html">
  another tutorial.</a>)
  \n\n
- The mesh generation -- the specification of the curvilinear
  boundaries and the geometry of the rib is somewhat tedious.
  We refer to 
  <a href="../../../meshes/mesh_from_inline_triangle_internal_boundaries/html/index.html">
  another tutorial</a> for a discussion of how to define the
  internal mesh boundary that separates the two regions
  (the rib and the annular region) so that we can assign 
  different material properties to them. 
  \n\n
.  
All of this is reasonably straightforward and provides a
powerful code that automatically adapts the mesh while
respecting the curvilinear boundaries of the domain. 
Have a look through the driver code and play with it.

<HR>
<HR> 

\section code Code listing
Here's a listing of the complete driver code:

\include unstructured_time_harmonic_elastic_annulus.cc

<HR>
<HR> 


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/time_harmonic_linear_elasticity/elastic_annulus/">
demo_drivers/time_harmonic_linear_elasticity/elastic_annulus
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/time_harmonic_linear_elasticity/elastic_annulus/unstructured_time_harmonic_elastic_annulus.cc">
demo_drivers/time_harmonic_linear_elasticity/elastic_annulus/unstructured_time_harmonic_elastic_annulus.cc
</A>
</CENTER>
.




<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

