/**

\mainpage Demo problem: Compression of 2D circular disk

In this example we study the compression of a 2D circular disk, loaded
by an external pressure. We also demonstrate:
- how to "upgrade" a \c Mesh to a \c SolidMesh
- why it is necessary to use "undeformed MacroElements" to ensure that
  the numerical results converge to the correct solution under 
  mesh refinement if the domain has curvilinear boundaries.
- how to switch between different constitutive equations
- how to incorporate isotropic growth into the model
.
We validate the numerical results by comparing them 
against the analytical solution of the equations of linear
elasticity which are valid for small deflections. 

<HR>
<HR>


\section problem The problem
The figure below shows a sketch of the basic problem: A 2D circular disk
of radius \f$ a \f$ is loaded by a uniform pressure \f$ p_0^*\f$.
We wish to compute the disk's deformation for a variety of
constitutive equations.

\image html disk_compression.gif "Sketch of the problem. " 
\image latex disk_compression.eps "Sketch of the problem. " width=0.75\textwidth

The next sketch shows a variant of the problem: We assume that 
the material undergoes isotropic growth (e.g. via a
biological growth process or thermal expansion, say) with a constant growth
factor \f$ \Gamma \f$. We refer to the the document
<A HREF="../../solid_theory/html/index.html#isotropic_growth">
"Solid mechanics: Theory and implementation"</A> for a detailed 
discussion of the theory of isotropic growth. Briefly, the growth 
factor defines the relative increase in the volume of an infinitesimal material
element, relative to its volume in the stress-free reference 
configuration. If the growth factor is spatially uniform, isotropic 
growth leads to a uniform expansion of the material. For a circular disk, 
uniform growth increases the disk's
radius from \f$ a_0 \f$ to \f$ a = a_0 \sqrt{\Gamma} \f$ without
inducing any internal stresses. This uniformly expanded
disk may then be regarded as the stress-free reference configuration
upon which the external pressure acts.


\image html disk_compression_with_growth.gif "Sketch of the problem. " 
\image latex disk_compression_with_growth.eps "Sketch of the problem. " width=0.75\textwidth

      
<HR>
<HR>

\section results Results

The animation shows the disk's deformation when subjected to uniform
growth of \f$ \Gamma = 1.1 \f$ and loaded by a pressure that ranges
from negative to positive values. All lengths were scaled on the 
disks initial radius (i.e. its radius in the absence of growth and 
without any external load).

\image html disk_compression_anim.gif "Deformation of the uniformly grown disk when subjected to an external pressure. " 
\image latex disk_compression_anim.eps "Deformation of the uniformly grown disk when subjected to an external pressure. " width=0.75\textwidth


The figure below illustrates the disk's load-displacement characteristics 
by plotting the disk's non-dimensional radius as function of the 
non-dimensional pressure, \f$ p_0 = p_0^*/{\cal S} \f$ , where 
\f$ {\cal S} \f$  is the characteristic stiffness of the material, for
a variety of constitutive equations.  


\image html trace.gif "Load displacement characteristics for various constitutive equations. " 
\image latex trace.eps "Load displacement characteristics for various constitutive equations. " width=0.75\textwidth

\subsection hooke Generalised Hooke's law

The blue, dash-dotted line corresponds to \c oomph-lib's
generalisation of Hooke's law (with Young's modulus \f$ E \f$ and
Poisson ratio \f$ \nu \f$ ) in which the dimensionless second 
Piola Kirchhoff stress tensor (non-dimensionalised with the
material's Young's modulus \f$ E \f$, so that \f$ {\cal S} = E \f$ ) 
is given by
\f[
\sigma^{ij} = \frac{1}{2(1+\nu)}
\left( 
G^{ik} G^{jl} + G^{il} G^{jk}  +
\frac{2\nu}{1-2\nu} G^{ij} G^{kl} 
\right) \gamma_{kl}.
\f]
Here \f$ \gamma_{ij}= 1/2 (G_{ij}-g_{ij})\f$ is Green's strain
tensor, formed from the difference between the deformed and
undeformed metric tensors, \f$ G_{ij} \f$  and \f$ g_{ij} \f$, 
respectively. The three different markers identify the
results obtained with the two forms of the
principle of virtual displacement, employing the displacement
formulation (squares), and a pressure/displacement formulation
with a continuous (delta) and a discontinuous (nabla) pressure interpolation.

For zero pressure the disk's non-dimensional radius is equal to the
uniformly grown radius \f$ \sqrt{\Gamma} = 1.0488. \f$ For
small pressures the load-displacement curve follows
the linear approximation
\f[
r =  \sqrt{\Gamma} \big( 1-p_0 (1+\nu)(1-2\nu) \big).
\f]
We note that the generalised Hooke's law leads to strain
softening behaviour under compression (the pressure required
to reduce the disk's radius to a given value increases more 
rapidly than predicted by the linear approximation) whereas
under expansion (for negative external pressures) the behaviour is 
strain softening.

\subsection mr Generalised Mooney-Rivlin law

The red, dashed line illustrates the behaviour when Fung \& Tong's
generalisation of the Mooney-Rivlin law (with Young's modulus, \f$ E \f$,
Poisson ratio \f$ \nu \f$ and Mooney-Rivlin parameter \f$ C_1\f$)
is used as the constitutive equation. For this constitutive law, 
the non-dimensional strain energy function \f$ W = W^*/{\cal S}\f$,
where the characteristic stress is given by Young's modulus, i.e.  \f$ {\cal
S} =E \f$, is given by
\f[
W =  \frac{1}{2} (I_1-3) + (G-C_1)(I_2-3) + (C_1-2G)(I_3-1)
+ (I_3-1)^2 \ \frac{G(1-\nu)}{2(1-2\nu)} \ ,
\f]
where 
\f[
G = \frac{E}{2(1+\nu)}
\f]
is the shear modulus, and \f$ I_1, I_2 \f$ and \f$ I_3 \f$ are the
three invariants of Green's strain tensor. See 
<A HREF="../../solid_theory/html/index.html#strain-energy">
"Solid mechanics: Theory and implementation"</A> for a detailed 
discussion of strain energy functions. The figure shows that
for small deflections, the disk's behaviour is again well approximated
by linear elasticity. However, in the large-displacement regime
the Mooney-Rivlin is strain hardening under extension and softening
under compression when compared to the linear elastic behaviour.
      
<HR>
<HR>

\section namespace Global parameters and functions
As usual we define the global problem parameters in a namespace.
We provide pointers to the constitutive equations and strain energy
functions to be explored, and define the associated constitutive
parameters.

\dontinclude disk_compression.cc
\skipline namespace_for
\until C1=1.3;


Next we define the pressure load, using the general interface
defined in the \c SolidTractionElement class.  The arguments
of the function reflect that  the load on 
a solid may be a function of the Lagrangian and Eulerian coordinates, 
and the external unit normal on the solid. Here we apply a spatially
constant external pressure of magnitude \c P which acts in the
direction of the negative outer unit normal on the solid.

\until end of pressure load

Finally, we define the growth function and impose a spatially uniform
expansion that (in the absence of any external load) would 
increase the disk's volume by 10%. 

\until end namespace

<HR>
<HR>

\section main The driver code

The driver code is very short: We store the command line arguments
(as usual, we use a non-zero number of command line arguments 
as an indication that the code is run in self-test mode and reduce
the number of steps performed in the parameter study)
and create a strain-energy-based constitutive equation: Fung \& Tong's
generalisation of the Mooney-Rivlin law. 

\skipline start_of_main
\until Strain_energy_function_pt);

We build a problem object, using the displacement-based 
\c RefineableQPVDElements to discretise the domain, and perform
a parameter study, exploring the disk's deformation for a range
of external pressures. 

\until done case 0

We repeat the exercise with elements from the 
\c RefineableQPVDElementWithContinuousPressure family
which discretise the principle of virtual displacements (PVD) 
in the pressure/displacement formulation, using continuous pressures
(Q2Q1; Taylor Hood).

\until done case 1

The next computation employs
\c RefineableQPVDElementWithPressure elements
in which the pressure is interpolated by piecewise linear
but globally discontinuous basis functions (Q2Q-1; Crouzeiux-Raviart).

\until done case 2

Next, we change the constitutive equation to \c oomph-lib's
generalised Hooke's law,

\until  Global_Physical_Variables::E);

before repeating the parameter studies with the same three 
element types:

\until end of main




 


<HR>
<HR>

\section mesh The mesh
We formulate the problem in cartesian coordinates (ignoring the
problem's axisymmetry) but discretise only one quarter of 
the domain, applying appropriate symmetry conditions along the 
x and y axes. The computational domain may be discretised
with the \c RefineableQuarterCircleSectorMesh that we 
already used in many previous examples. To use the mesh in 
this solid mechanics problem we must first "upgrade" 
it to a \c SolidMesh. This is easily done by multiple inheritance:

\dontinclude disk_compression.cc
\skipline start_mesh
\until {

The constructor calls the constructor of the underlying \c
RefineableQuarterCircleSectorMesh and sets the Lagrangian
coordinates of the nodes to their current Eulerian positions,
making the initial configuration stress-free. 

\until }

We also provide a helper function that creates a mesh of 
\c SolidTractionElements which are attached to the curved domain 
boundary (boundary 1). These elements will be used to apply the 
external pressure load.

\until };


<HR>
<HR>

\section problem_class The Problem class
 
The definition of the Problem class is very straightforward. 
In addition to the
constructor and the (empty) \c actions_before_newton_solve() and 
\c actions_after_newton_solve() functions, we provide the function
\c parameter_study(...) which performs a parameter study, computing the disk's 
deformation for a range of external pressures. The member data
includes pointers to the mesh of "bulk" solid elements, and the
mesh of \c SolidTractionElements that apply the pressure load.
The trace file is used to document the disk's load-displacement
characteristics by plotting the radial displacement of the nodes on the
curvilinear boundary, pointers to which are stored in the vector
\c Trace_node_pt.

\until };

<HR>
<HR>

\section constructor The Constructor

We start by constructing the mesh of "bulk" \c SolidElements, using
the \c Ellipse object to specify the shape of the curvilinear
domain boundary. 

\until xi_hi);

Next we choose the nodes on the curvilinear domain boundary (boundary
1) as the nodes whose displacement we document in the trace file.

\until boundary_node_pt(1,j)

The \c QuarterCircleSectorMesh that forms the basis of the
"bulk" mesh contains only three elements -- not enough to 
expect the solution to be accurate. Therefore we apply 
one round of uniform mesh refinement before attaching the
\c SolidTractionElements to the mesh boundary 1, using the
function \c make_traction_element_mesh() in the 
\c ElasticRefineableQuarterCircleSectorMesh. 

\until make_traction

We add both meshes to the \c Problem and build a combined
global mesh:

\until build_global_mesh();

Symmetry boundary conditions along the horizontal and vertical 
symmetry lines require that the nodes' vertical position is pinned
along boundary 0, while their horizontal position is pinned 
along boundary 2. 

\until pin_position(1)

Since we are using refineable solid elements, we pin any "redundant"
pressure degrees of freedom in the "bulk" solid mesh (see
<A HREF="../../airy_cantilever/html/index.html#ex">
the exercises in another tutorial</A> for a more detailed 
discussion of this issue). 

\until Solid_mesh_pt

Next, we complete the build of the elements in the "bulk" solid mesh
by passing the pointer to the constitutive equation and the pointer
to the isotropic-growth function to the elements:

\until }

We repeat this exercise for the \c SolidTractionElements which must be
given a pointer to the function that applies the pressure load

\until }

Finally, we set up the equation numbering scheme and report the number
of unknowns.

\until }


<HR>
<HR>

\section doc_solution Post-processing 
The post-processing function outputs the shape of the deformed disk.
We use the trace file to record how the
disk's volume (area) and the radii of the control nodes on the
curvilinear domain boundary vary with the applied pressure. To 
facilitate the validation of the results against the analytical
solution, we also add the radius predicted by the linear theory 
to the trace file.
 
\until end of doc_solution

<HR>
<HR>

\section run_it Performing the parameter study

The function \c parameter_study(...) computes the disk's deformation
for a range of external pressures and outputs the results. 
The output directory is labelled by the \c unsigned 
function argument. This ensures that parameter studies performed
with different constitutive equations are written into different
directories.

\until end of parameter study

<HR>
<HR>

\section comm_and_ex Comments and Exercises

\subsection macro The use of MacroElements in solid mechanics problems
<A HREF="../../../poisson/fish_poisson2/html/index.html">Recall</A>
how \c oomph-lib employs \c MacroElements to represent the 
exact domain shapes in adaptive computations involving problems 
with curvilinear boundaries. When an element is refined, the
(Eulerian) position of any newly-created nodes is based on the
element's \c MacroElement counterpart, rather than being determined
by finite-element interpolation from the "father element". This ensures
that <b>(i)</b> newly-created nodes on curvilinear domain boundaries are placed
exactly onto those boundaries and <b>(ii)</b> that newly-created nodes in the
interior are placed at positions that match smoothly onto the boundary.

This strategy is adapted slightly for solid mechanics problems:
-# The Eulerian position of newly-created \c SolidNodes is
   determined by finite element interpolation from the 
   "father element", \b unless the newly-created \c SolidNode is located
   on a domain boundary and its position is pinned by displacement
   boundary conditions. 
   \n\n
-# The same procedure is employed to determine the Lagrangian 
   coordinates of newly-created \c SolidNodes.
.
These modifications ensure that, as before, newly-created nodes 
on curvilinear domain boundaries are placed exactly onto those
boundaries if their positions are pinned by displacement boundary
conditions. (If the nodal positions are not pinned, the node's
Eulerian position will be determined as part of the
solution.) The use of finite-element interpolation from the "father element"
in the interior of the domain for both Lagrangian and Eulerian
coordinates ensures that the creation of new nodes does not induce
any stresses into a previously computed solution.

   
\subsection ex Exercises

-# Our discretisation of the problem in cartesian coordinates 
   did not exploit the problem's axisymmetry. Examine the trace file to 
   assess to which extent the computation retained the axisymmetry. 
   \n\n



<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

