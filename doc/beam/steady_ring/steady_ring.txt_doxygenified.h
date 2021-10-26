/**

\mainpage Demo problem: Large-displacement post-buckling of a pressure-loaded thin-walled elastic ring -- displacement-control techniques
 
In this example we study a more challenging beam problem: The
non-axisymmetric buckling of a thin-walled elastic 
ring, loaded by a spatially-constant external pressure, \f$ p_{ext}\f$ .
For sufficiently small positive (or arbitrarily large negative)
values of  \f$ p_{ext}\f$ the ring deforms axisymmetrically and in
this mode it is very stiff, implying that large changes in pressure
are required to change its radius. However, if  \f$ p_{ext}\f$ exceeds
a critical threshold, the axisymmetric configuration becomes unstable,
causing the ring to buckle non-axisymmetrically. Once buckled,
the ring is much more flexible and small changes in \f$
p_{ext}\f$ are sufficient to induce dramatic changes in its shape. 

The rapid change in stiffness following the buckling makes it 
difficult to compute strongly buckled solutions by continuation 
methods as the solution computed for a certain value of  \f$ p_{ext}\f$ 
may represent a poor approximation of the solution at a 
slightly larger pressure. In extreme cases, this can 
cause the Newton method (whose convergence relies on the provision of
good initial guesses for the unknowns) to diverge. 

We will demonstrate the use of so-called "displacement control" techniques
to overcome these problems. Displacement-control techniques are 
useful in problems in which some \e a-priori knowledge about the
expected displacement field allows us to re-formulate the
problem. Rather than prescribing the load on the elastic solid 
and computing the displacement field, we prescribe the displacement 
of a carefully-selected material "control" point and regard the 
load required to achieve this displacement as an unknown. 

<CENTER>
<TABLE>
<TR>
<TD>
\n
<CENTER>
<B>Non-axisymmetric buckling of a pressure-loaded, thin-walled elastic ring</B>
</CENTER> 

We wish to compute the deformation of a linearly-elastic, circular 
ring of undeformed radius \f$ R_0 \f$ and wall thickness \f$ h^* \f$,
subject to a spatially-constant external pressure \f$ p_{ext}^* \f$. 
We choose the undeformed radius \f$ R_0 \f$ as the lengthscale
for the non-dimensionalisation of the problem and assume that 
\f$ h = h^*/R_0 \ll 1 \f$ , justifying the use of beam theory.
As discussed in <A HREF="../../tensioned_string/html/index.html">
the previous example</A> we scale the pressure on the ring's 
effective Young's modulus, \f$ p_{ext} = p_{ext}^*/E_{eff} \f$ , 
where \f$ E_{eff} = E/(1-\nu^2) \f$ . We use the non-dimensional
arclength  \f$ \xi \in [0,2\pi] \f$ along the ring's undeformed 
centreline as the Lagrangian coordinate and parametrise the position vector 
to the ring's undeformed centreline as
\f[
{\bf r}_w(\xi) = \left( 
\begin{array}{c}
\cos(\xi) \\
\sin(\xi) 
\end{array}
\right).
\f]
We wish to compute the position vector \f$ {\bf R}_w(\xi) \f$ 
to the deformed ring's centreline.


Here is a sketch of the problem:

\image html ring_sketch.gif "Sketch of the buckling ring. " 
\image latex ring_sketch.eps "Sketch of the buckling ring. " width=0.75\textwidth

Note that we have chosen the Eulerian coordinate axes so that they 
coincide with the ring's lines of symmetry. 

</TD> 
</TR>
</TABLE>  
</CENTER>  

<HR>

\section theory The expected deformation and the displacement-control formulation of the problem 
Standard linear stability analysis [see,
e.g., G.J. Simitses "An introduction to the elastic stability of 
structures", Prentice-Hall, (1976)] predicts the ring to become 
unstable to non-axisymmetric perturbations with an azimuthal
wavenumber of \f$ N=2 \f$  at a pressure
of 
\f[ 
p_{ext}^{[buckl]} = 3 \times \frac{1}{12} h^3.
\f]
For \f$ p_{ext} > p_{ext}^{[buckl]} \f$ we therefore expect the ring to 
deform into a shape similar to the one indicated by the dashed line
in the above sketch. As the ring buckles non-axisymmetrically, 
the material point on 
the vertical symmetry line (i.e. the material point with Lagrangian coordinate 
\f$ \xi=\pi/2\f$ ) moves radially inwards. This makes it 
an excellent control point for this problem as we can "sweep" through 
the entire range of the ring's post-buckling deformation by
varying its \f$ x_2 \f$ -- coordinate 
from \f$ R_{ctrl} =1 \f$ (corresponding to the undeformed,
axisymmetric configuration) to \f$ R_{ctrl} = 0 \f$
(corresponding to a configuration in which the ring is collapsed 
to the point of opposite wall contact). We note that Flaherty, 
Keller & Rubinow's  analysis [SIAM J. Appl. Math \b 23, 446--455
(1972)], based on an inextensible beam model, predicts opposite 
wall contact to occur at an external pressure of
\f[
p_{ext}^{[owc]} = 5.22 \times \frac{1}{12} h^3.
\f]

To apply displacement control we change the formulation of
the problem from 
<CENTER>
<TABLE>
<TR>
<TD>

<CENTER>
<B> Original formulation of the problem</B>
</CENTER> 

Determine the position vector to the centreline of the 
deformed ring, \f$ {\bf R}_w(\xi) \f$ , in terms of the Lagragian 
coordinate \f$ \xi \in [0,2\pi] \f$ , for a given value of the 
external pressure \f$ p_{ext} \in [0, 5.22 \times h^3/12] \f$ .
\n
</TD> 
</TR>
</TABLE>  
</CENTER>  


to 

<CENTER>
<TABLE>
<TR>
<TD>

<CENTER>
<B> Displacement-control formulation of the problem</B>
</CENTER> 

Determine the position vector to the centreline of the 
deformed ring, \f$ {\bf R}_w(\xi) \f$ , in terms of the Lagragian 
coordinate \f$ \xi \in [0,2\pi] \f$ , for an external pressure 
\f$ p_{ext} \f$ that is such that the vertical coordinate of the
control point is equal to \f$ R_{ctrl} \f$, i.e. 
\f[
{\bf R}_w(\xi=\pi/2) \cdot {\bf e}_2 = R_{ctrl},
\f]
where \f$ R_{ctrl} \in [0,1] \f$ is given.

</TD> 
</TR>
</TABLE>  
</CENTER>  


<HR>

\section results Results

The figure below shows computed ring shapes for \f$ h^*/R_0 =
1/20 \f$. They were obtained in two phases: During the first phase 
of the computation we subjected the ring to a load of the form
\f[
{\bf f} = \bigg( p_{ext} - p_{cos} \cos(2 \xi) \bigg) \, {\bf N}
\f]
where \f$ {\bf N} \f$ is the outer unit normal on the deformed ring
and \f$ p_{cos} \ll 1 \f$ is a small cosinusoidal perturbation
to the external pressure which forces the ring to buckle
"in the right direction". The undeformed configuration was 
used as the initial guess for an initial computation with 
\f$ p_{ext} = 0\f$ [or, in the case of displacement control, 
\f$ R_{ctrl} = 1\f$ ]. Subsequently we increased \f$ p_{ext} \f$
[or decreased \f$ R_{ctrl}\f$ ] in small steps, using the previously
computed solutions for the displacement field [and \f$ p_{ext} \f$ ]
as initial guesses for the unknowns. This procedure was continued until
the ring was collapsed up to the point of opposite wall contact. 
During the second phase, we set \f$ p_{cos} = 0\f$ and reversed the
procedure to re-trace the deformation to the axisymmetric state. 

\image html ring.gif "Computed ring shapes. " 
\image latex ring.eps "Computed ring shapes. " width=0.75\textwidth


The figure below illustrates the load/displacement characteristics
computed by this procedure. The graph shows the radii of two 
material points on the ring: The green line shows the radius \f$ R_1 =
R_{ctrl} \f$ of the control point; the red line the radius \f$ R_2 \f$ 
of the material point located at \f$ \xi = 0\f$, i.e. the
material point located on the ring's second line of symmetry. 
(Because of the symmetry of the buckling pattern this line may also be
interpreted as the load-displacement curve for the control point when
the ring buckles in the "opposite" direction). The dash-dotted blue line
shows the load/displacement curve when the ring deforms axisymmetrically.
In this mode, the radii of the two material points are obviously 
identical. Finally,  the dashed lines shows the load/displacement 
path when the ring is subjected to 
a non-zero perturbation pressure, \f$ p_{cos} \ne 0\f$ , which deforms
it non-axisymmetrically so that \f$ R_1 \ne R_2 \f$ even for 
\f$ p_{ext} < p_{ext}^{[buckl]}\f$ . 

\image html trace.gif "Bifurcation diagram for a buckling ring. " 
\image latex trace.eps "Bifurcation diagram for a buckling ring. " width=0.75\textwidth

The diagram clearly illustrates the enormous change in stiffness
when the ring changes from the axisymmetric to the non-axisymmetric 
regime. The non-axisymmetric regime emanates from the axisymmetric
state (via a supercritical bifurcation at a pressure of \f$ p_{ext} =
3.0 \times h^3/12 \f$ , as predicted by the linear stability analysis)
and opposite wall contact occurs at  \f$ p_{ext} = 5.22 \times h^3/12
\f$ , in perfect agreement with Flaherty, Keller & Rubinow's 
theoretical analysis.



<HR>

\section displ_control Applying displacement control in oomph-lib

To facilitate the solution of solid mechanics problems with 
displacement-control techniques, \c oomph-lib provides a
\c DisplacementControlElement, which may be used to
add the displacement control constraint to the global system of
equations that is solved by Newton's method. Applying displacement
control in a solid mechanics problem involves the following 
straightforward steps:
-# <B>Formulate the solid mechanics problems in its "standard
   form":</B> \n\n In this formulation the (scalar) load level is 
   prescribed and the displacements are regarded as unknowns. \n\n
-# <B>Change the representation of the (scalar) load level to allow
   it to be an unknown in the Problem:</B> \n\n In the original 
   formulation, the load level
   is likely to have been represented by a double precision number. 
   To allow the load level to be regarded as an unknown, it must be 
   represented by a value in a \c Data object. Currently, 
   \c oomph-lib's \c DisplacementControlElement expects the load level 
   to be the one-and-only value in the load \c Data object. We note 
   that computations with a prescribed load are still possible and 
   simply require pinning of the value that represents the load. \n\n
-# <B>"Tell" the solid mechanics elements that the load level is a (potential) 
   unknown:</B> \n\n
   Since the load on the solid mechanics elements affects their
   residuals, their elemental Jacobian matrices (which contain 
   the derivatives of the entries in the elemental residual vector
   with respect to \e all unknowns that affect it) must take the
   dependence on the (potentially) unknown load level into account. This 
   may be achieved by adding the \c Data object that stores the load
   level to the element's external \c Data (i.e. \c Data whose
   values \e affect the element's residual vector, but whose values are
   not determined \e by the element). External \c Data is
   automatically included in an element's equation numbering
   procedure. Furthermore, since the elemental Jacobian matrices of 
   \c SolidFiniteElements are generated by finite-differencing,
   the derivatives of the element's residual vector with respect 
   to the load level are computed automatically. Consequently, the 
   application of displacement control does not require a
   re-implementation of the solid mechanics elements. \n\n
-# <B>Identify a material point in the solid that can serve as a 
   "control point":</B> \n \n Ideally, this control point should move
   monotonically along one of the (Eulerian) coordinate directions as the
   deformation of the solid body increases.  \n\n
-# <B>Create a \c DisplacementControlElement and add it to the 
   \c Mesh.</B> \n\n
   The \c DisplacementControlElement adds the displacement constraint
   to the global system of equations and thus provides the 
   additional equation required to determine the unknown load
   level. We note that the \c DisplacementControlElement has 
   two constructors: \n\n
   - The first version expects a pointer to the \c Data object
     whose one-and-only value contains the unknown load level. 
     This version of the constructor is appropriate in cases where
     the load \c Data has already been created elsewhere (as described 
     above).  \n\n
   - The second version of the constructor creates the required 
     load \c Data object internally and provides access to it via
     a pointer-based access function. \n\n
   .
   The section \ref global_data provides a more detailed discussion 
   of when to use which version of the constructor. 
   \n\n 
-# <B>Done!</B> \n\n Set the desired value of the control displacement
   and solve the problem. \n\n
.
 
The driver code discussed below illustrates these steps.

<HR>
<HR>



\section namespace Global parameters and functions
The namespace \c Global_Physical_Variables, used to define the problem
parameters and the load function, is very similar to that
used in <A HREF="../../tensioned_string/html/index.html">the previous
example</A> without displacement control. They main difference is that
the adjustable load (the external pressure \f$ p_{ext}\f$ ) is now
defined as a \c Data value, rather than a double precision number.
This allows it to become an unknown in the problem when displacement
control is used.

\dontinclude steady_ring.cc
\skipline start_of_namespace
\until end of namespace


<HR>


\section main The driver code
The main function builds two different versions of the problem, 
demonstrating the use of displacement control in the cases when
the \c Data object that contains the adjustable load already exists,
or has to be created by the \c DisplacementControlElement,
respectively. The two versions only differ in the  \ref constructor. 

\dontinclude steady_ring.cc
\skipline start_of_main
\until end of main

<HR>

\section problem_class The problem class

The problem class is very similar to the one used  
in <A HREF="../../tensioned_string/html/index.html">the previous
example</A> without displacement control. Here we store the
number of solid mechanics elements in the Problem's private member
data since we will add the \c DisplacementControlElement to the mesh. 

\dontinclude steady_ring.cc
\skipline start_of_problem_class
\until end of problem class

<HR>

\section constructor The constructor

The first half of the constructor is similar to the one used
in <A HREF="../../tensioned_string/html/index.html">the previous
example</A> without displacement control. We create a \c GeomObject
(an \c Ellipse with unit half axes) to define the ring's undeformed
geometry, and build the 1D Lagrangian mesh. Note that, because of 
the symmetry of the buckling pattern, we only discretise one 
quarter of the domain.

\dontinclude steady_ring.cc
\skipline start_of_constructor
\until OneDLagrangianMesh

Next we apply the symmetry boundary conditions: Zero vertical [horizontal]
displacements and infinite [zero] slope at \f$ \xi =0
 \f$  [at \f$ \xi =\pi/2 \f$ ]. (See 
<A HREF="../../tensioned_string/html/index.html">
the previous example</A> for a more detailed discussion of the
boundary conditions for beam elements.)

\skipline Boundary condition:
\until pin_position(1,1)

We now distinguish between the cases with and without displacement
control:

- <B>Case 1: No displacement control</B> \n\n
  If we don't use displacement control, we create the \c Data object
  whose one-and-only value stores the adjustable load level (the
  external pressure), and store the pointer to the newly created
  \c Data object in \c Global_Physical_Variables::Pext_data_pt. 
  Since the value of the external pressure is prescribed (i.e. not 
  an unknown in the problem), we pin the value.  \n\n
  \skipline Normal load incrementation
  \until Global_Physical_Variables::Pext_data_pt->pin(0);
  \n\n
- <B>Case 2: Displacement control</B> \n\n
  If displacement control is used, we identify the control point
  (the material point at the "end" of the last element in the mesh)
  and the (Eulerian) coordinate direction (the vertical coordinate) 
  in which its position is to be controlled:\n\n
  \skipline Displacement control
  \until    DisplacementControlElement* displ_control_el_pt;
  \n\n
  - <B>Case 2a: Load \c Data does not yet exist</B> \n\n
    We can now call the constructor for the \c
    DisplacementControlElement: \n\n
    \skipline  Build displacement control element
    \until &Global_Physical_Variables::Xprescr);
    \n\n
    This version of the constructor creates the load \c Data 
    object whose one-and-only value stores the adjustable load. 
    We obtain a pointer to the newly-created \c Data object from 
    the access function
    \c DisplacementControlElement::displacement_control_load_pt() and
    store it in \c Global_Physical_Variables::Pext_data_pt to make it
    accessible to the load function 
    \c Global_Physical_Variables::press_load(...) \n\n
    \skipline The constructor of the
    \until displacement_control_load_pt();
    \n\n
  - <B>Case 2b: The load \c Data already exists</B> \n\n
    In some applications, the \c Data object that specifies the adjustable
    load level might already have been created elsewhere in the
    \c Problem. For such cases, \c oomph-lib provides an alternative 
    constructor whose argument list includes the pointer to the 
    already-existing load \c Data object. To demonstrate its use
    we create a suitable \c Data object and store it in the \c
    Problem's "global" \c Data (see \ref global_data for a more
    detailed discussion of this step) and pass it to the constructor: \n\n
    \skipline // Demonstrate use of displacement control with some
    \until }
  .
  \n\n
  In both  cases, we add the newly created \c
  DisplacementControlElement to the mesh to ensure that it is included
  in the \c Problem. \n\n
  \skipline Add the displacement-control
  \until }
.
Next, we execute the usual loop over the elements to pass the
pointers to the problem's non-dimensional parameters the elements. 
If displacement control is used, we also pass the pointer to the 
load \c Data object to the elements' external \c Data to indicate 
that the element's residual vectors depends on the unknown 
load. Finally, we set up the equation numbering scheme. 

\skipline Loop over the elements
\until end of constructor

<HR>

\section doc Post-processing
The post-processing function documents the ring shapes and 
adds the load/displacement characteristics of the two material
points on the ring's symmetry lines to the trace file.
\skipline start_of_doc
\until end of doc


<HR>

\section parameter The parameter study
We start by opening a trace file to record the load/displacement
characteristics, and output the initial configuration. 

\skipline start_of_run
\until doc_solution

Next we set up the increment of the control parameter, choosing 
the displacement or load increments such that the ring's deformation
is increased from the axisymmetric initial state to total collapse
with opposite wall contact in 11 steps. 

\skipline Number of steps
\until Global_Physical_Variables::Pext_data_pt->set_value

Without displacement-control the Newton method can require a large
number of iterations, therefore we increment the maximum number of 
iterations. 

\skipline Without
\until }

We start the parameter study by increasing the ring's compression
(either by increasing the external pressure or by reducing the
\f$ x_2\f$ - coordinate of the control point) with \f$ p_{cos} > 0 \f$
to induce buckling in the desired direction.
 
\skipline // Downward loop over parameter incrementation with pcos>0 
\until end of downward


Then we reset the perturbation pressure to zero and reduce the ring's
collapse by decreasing the external pressure (or by increasing the
\f$ x_2\f$ - coordinate of the control point). 

\skipline // Reset
\until end of run

Done!

<HR>
<HR>
 
\section further Further comments and Exercises
\subsection global_data "Internal" and "external" Data in elements and "global" Data in Problems

In the section \ref constructor we encountered two 
different constructors for the \c DisplacementControlElement. 
- The first version (in which the \c Data object that contains the 
  unknown load value is created by the constructor) is most natural 
  for the problem considered here as the load \c Data is
  created specifically for the purpose of allowing displacement
  control to be applied. It therefore natural to regard the
  \c DisplacementControlElement as being "in charge of" the load 
  \c Data , and storing it in the element's "internal Data". 
  (<A HREF="../../../the_data_structure/html/index.html">Recall</A>
  that once \c Data is stored in an element's internal \c Data,
  the element becomes responsible for performing tasks such as equation
  numbering, timestepping, etc. that must be performed exactly once
  for each \c Data value in the \c Problem.) \n\n
- The second version of the constructor is appropriate for problems 
  in which the load \c Data has already been created by another 
  element, implying that the other element is already "in charge" 
  of the \c Data object and performs the equation numbering,
  timestepping etc. for its values. In that case, the load \c Data
  is regarded as "external Data" for the \c
  DisplacementControlElement. 
.
In our example code, we simulated the second scenario
by creating the load \c Data object before calling the
second version of the constructor. While this ensures that
the load \e can be regarded as an unknown in the problem,
the \c Problem remains "unaware" of the additional unknown, 
as none of the elements in the \c Problem's mesh is "in charge" 
of it. While this is clearly a 
somewhat artificial scenario, \c oomph-lib provides a 
mechanism for handling such cases: Adding a \c Data object
to the \c Problem's "global Data" by calling

\code
Problem::add_global_data(...)
\endcode

puts the \c Problem itself "in charge" of this \c Data.
 
<HR>

\subsection exercises Exercises
-# Run the code without displacement control (e.g. by setting the
   boolean flag \c displ_control in the main function to \c false)
   and confirm that a non-zero perturbation pressure, \f$ p_{cos} \ne 0\f$
   is required to induce the ring's non-axisymmetric collapse
   This shows that the numerical model is not as sensitive to 
   non-axisymmetric perturbations as the theory suggests -- 
   roundoff error alone is not sufficient to initiate
   non-axisymmetric buckling. Use this version of the code to 
   compute the load-displacement characteristics of the ring in 
   its axisymmetric state, i.e. the dash-dotted blue line in 
   the bifurcation diagram shown at the beginning of this document. \n\n
-# Run the code without displacement control and explain why 
   during the second phase of the parameter study (when 
   \f$ p_{cos} = 0 \f$ ),  the Newton method converges very 
   slowly and provides very inaccurate results 
   when \f$ p_{ext} = 3\f$ .\n\n
-# We claimed that the load/displacement curves for the post-buckling
   regime emanate from the axisymmetric branch at \f$ p_{ext} = 3 \f$ , yet
   closer inspection of the bifurcation diagram shows a small kink
   in the post-buckling curves near the bifurcation. Explain what
   causes this kink and why it is practically impossible to avoid
   its occurrence. [Hint: The load/displacement curves contain 
   individual data points, each one of which corresponds to 
   a solution of the governing equations. The solutions were obtained
   by Newton's method which tends to converge to a solution in which
   the unknowns are "close" to their values at the beginning of the
   iteration. Is it obvious that an initial guess that 
   corresponds to a non-axisymmetric configuration is necessarily "closer" 
   to another non-axisymmetric solution than to a nearby axisymmetric 
   solution?] \n\n
-# The computation shown above was performed with the default
   non-dimensional wall thickness of \f$ h^*/R_0 = 1/20 \f$ .
   Repeat the computation with smaller wall-thicknesses (the 
   <A HREF="../../tensioned_string/html/index.html">
   previous example</A> shows how to change the default value)
   and confirm the theoretical predictions for the
   dependence of \f$ p_{ext}^{[buckl]} \f$ and \f$ p_{ext}^{[owc]}
   \f$ on \f$ h^*/R_0 \f$. \n\n
-# Comment out the command that stores the load \c Data as
   global \c Data and explain why the code fails with a segmentation
   fault. Use the \c Problem::self_test() function to identify the problem. 
   [Hint: What is the default value for a \c Data object's global 
   equation numbers? What happens if this default is not overwritten?
   Why is the default assignment not overwritten?] 





<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="
../../../../
demo_drivers/beam/steady_ring/
">
demo_drivers/beam/steady_ring/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="
../../../../
demo_drivers/beam/steady_ring/steady_ring.cc
">
demo_drivers/beam/steady_ring/steady_ring.cc
</A>
</CENTER>
.




<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

