/**

\mainpage Demo problem: A static interface between two viscous fluids

\section overview Overview of the problem

 We consider a closed rectangular container of unit width that
contains two immiscible fluids at rest. The lower fluid is of a
prescribed volume \f$ \mathcal{V} \f$ and the interface between the
two fluids meets the wall of the container at a contact angle \f$
\theta_{c} \f$. In the absence of any body forces or external forces,
a static solution is obtained in which the velocity field is zero, 
the fluid pressure in each layer is constant,
 and free surface is 
 of constant curvature (an arc of a circle in
 two-dimensions), set by the
 contact angle and the geometry of the domain. From simple geometry, the
 mean curvature of the interface in the present problem 
 is \f$ \kappa = 1/r = 2\cos\theta_{c} \f$. 

 We shall be rather brief in our discussion of this problem because it is
 extremely similar to the 
 <a href="../../static_single_layer/html/index.html"> static free surface
 bounding a single layer of viscous fluid. </a> In fact, the main
 difference is that the region of upper fluid is no longer treated as
 a single external pressure, but must be meshed so that the
 Navier--Stokes equations can be solved within it. In fact, the most
 significant difference between the two problems is that we need to 
 <a href= #mesh> construct a two-layer mesh. </a>
 Dynamic two-fluid problems
 are introduced in <a href="../../two_layer_interface/html/index.html"> 
 another tutorial,
 </a> but the static problem discussed here is again complicated by the need to
 enforce a constant volume constraint. 

\section vol_const Enforcing the volume constraint
 
 Unlike the <a href="../../static_single_layer/html/index.html">
 equivalent single-fluid problem </a>, there is no external pressure,
 so the volume constraint must be associated with an internal pressure
 degree of freedom. Thus, we must hijack a pressure variable and we
 choose to do so in the upper fluid.

 \dontinclude static_two_layer.cc
 \skipline Hijack
 \until hijack_internal

 In addition, we must fix another fluid pressure at a fixed reference value
  so that the problem is non-degenerate. It is conceptually appealing
  to fix the reference pressure in the other (lower) fluid because the
  (constant) pressure in the upper fluid is already "being set" by the
  volume constraint.

 \skipline Pin a single
 \until fix_pre

\section contact_angle Enforcing the contact angle constraint

 The method of enforcing the contact angle constraint is exactly the
 same as discussed in the 
 <a href="../../static_single_layer/html/index.html#contact_angle">
 single-layer tutorial. </a>

\section mesh Constructing the Two-Layer Elastic Mesh

 We must create a two-layer elastic mesh that will allow us access to
 the elements in each fluid. We will also need to set the volume
 constraint by adding \c ElasticLineVolumeConstraintBoundingElements
 to the boundaries that surround one of the fluids and to add the \c
 ElasticLineFluidInterfaceElements along the interface. Thus, we need
 to change the boundaries of the existing mesh.
 
 We begin by inheriting from the \c
 ElasticRectangularQuadMesh. 
 \dontinclude static_two_layer.cc 
 \skipline start_of_specific_mesh_class
 \until public ElasticRectangularQuadMesh<ELEMENT>

 The arguments to the constructor specify
 the number of elements in the
 horizontal direction and in each layer and also the width of the
 container and the height of each layer. The remaining arguments
 determine whether the mesh should be made periodic in the x direction
 and are the \c TimeStepper object.

 \skipline ElasticTwoLayerMesh(
 \until Default_TimeStepper

 We provide
 separate storage for elements above and below the interface 

 \skipline Set up the pointers
 \until end of upper and lower

 and the
 elements adjacent to the interface in the upper and lower
 fluid. 
 
 \skipline Set the elements 
 \until } //end of bulk elements

 We will use these elements adjacent to the interface to construct the
 \c ElasticLineFluidInterfaceElements and it is important that we only
 add one layer of interface elements, as discussed in another tutorial.

 We next change the number of boundaries 

 \skipline Reset the number
 \until set_nboundary(7)
  
 and then reassign the existing boundary nodes to the new numbering
 scheme. This is tedious and not terribly instructive, so is not
 shown, but it's all in the code if you want to see how it's done.

 Finally, we add the nodes to the new interface boundary and setup the
 lookup schemes for the bulk elements adjacent to the new boundaries.

\dontinclude static_two_layer.cc
\skipline Add the nodes to the interface
\until setup_boundary_element_info();

\section problem The problem class

 The problem class is extremely similar to that in the 
 <a href="../../static_single_layer/html/index.html"> singer-layer
 problem. </a> The main differences are that the \c
 ElasticTwoLayerMesh is used instead of the \c
 ElasticRectangularQuadMesh and the boundary conditions are modified
 to take the new boundary numbering into account. In addition, the
 free surface elements and volume constraint elements 
 are created using the lookup schemes in the \c
 ElasticTwoLayerMesh to ensure that only a single layer of elements
 are included on the free surface. If we used the lookup schemes
 assigned by the generic function \c
 Mesh::setup_boundary_element_info() bulk elements on
 both sides of internal boundaries will be included, so we would
 construct twice as many interface elements as required. 
  <HR>
<HR>

\section comments Comments and Exercises

\subsection com Comments

- The driver code also contains a formulation in which \c SpineElements
  are used. Happily, the answers produced by both formulations are the same.

- An axisymmetric version of the problem is also included in the
  library. The problem is so similar that it is not described in any
  further detail. For a list of the differences between the
  two-dimensional and axisymmetric formulations of the single-layer
  problem see 
  <a href="../../../axisym_navier_stokes/axi_static_cap/html/index.html">
  another tutorial. </a>
.

\subsection exercises Exercises

-# Confirm that the computed pressure difference across the interface
   agrees with the analytic expression and the single-layer
   computations.
-# What happens if you do not fix a pressure?
-# What happens if you fix the pressure that is traded for the volume
   constraint?
-# Can you fix the reference pressure and choose the traded pressure
   value from the same fluid?
-# Use the generic \c Mesh::boundary_element_pt() function to construct the 
   interface elements. What happens? Why?
-# Modify the problem to include a non-zero gravitational body force?
   What happens to the interface in this case?
.

<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/static_cap/">
demo_drivers/navier_stokes/static_cap/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/static_cap/static_two_layer.cc">
demo_drivers/navier_stokes/static_cap/static_two_layer.cc
</A>
</CENTER>\n\n
- Source files for the axisymmetric version of the problem are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/axisym_navier_stokes/two_fluid_spherical_cap/">
demo_drivers/axisym_navier_stokes/two_fluid_spherical_cap/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/axisym_navier_stokes/two_fluid_spherical_cap/two_fluid_spherical_cap.cc">
demo_drivers/axisym_navier_stokes/two_fluid_spherical_cap/two_fluid_spherical_cap.cc
</A>
</CENTER>
.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

