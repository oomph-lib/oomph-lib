/**

\mainpage Demo problem: Adaptive solution of Poisson's equation in a fish-shaped domain.

In this document, we discuss the solution of a 2D Poisson problem
using \c oomph-lib's powerful mesh adaptation routines:

<CENTER>
<TABLE>
<TR> 
<TD>
<CENTER>
<B>Two-dimensional model Poisson problem in a non-trivial domain</B>
</CENTER> 
Solve
\f[
\sum_{i=1}^2 \frac{\partial^2u}{\partial x_i^2} = -1,
 \ \ \ \ \ \ \ \ \ \ (1)
\f]
in the fish-shaped domain \f$D_{fish} \f$, with homogeneous
Dirichlet boundary conditions
\f[
\left. u\right|_{\partial D_{fish}}=0.
\ \ \ \ \ \ \ \ \ \ (2)
\f]
</TD>
</TR>
</TABLE>  
</CENTER>  

\image html fish_poisson_soln.gif "Plot of the solution " 
\image latex fish_poisson_soln.eps "Plot of the solution " width=0.75\textwidth

The sharp corners in the domain create singularities 
in the solution (its derivatives are unbounded) and so 
accurate results can only be obtained if we use 
a sufficiently fine discretisation. Implementing this by
uniform mesh refinement would create a huge number of elements in the
interior of the domain where the fine discretisation 
is not required.

To avoid this problem, \c oomph-lib provides
mesh adaptation routines that 
automatically adapt the mesh, based on \e a \e posteriori error estimates.
Regions in which an error estimator indicates that the
solution is not resolved to the required accuracy are refined;
automatic unrefinement is performed in regions where the 
discretisation is unnecessarily fine. 

We provide a detailed discussion of the driver code 
<A HREF="../../../../demo_drivers/poisson/fish_poisson/fish_poisson.cc">
fish_poisson.cc</A> which illustrates a variety of mesh refinement
procedures. [The alternative driver code 
<A HREF="../../../../demo_drivers/poisson/fish_poisson/fish_poisson_no_adapt.cc">
fish_poisson_no_adapt.cc</A> solves the same problem without
mesh adaptation. Its structure is very similar to that in 
the <A HREF="../../../poisson/two_d_poisson/html/index.html">
2D Poisson problem considered earlier</A>. It is provided mainly
to illustrate how easy to it is incorporate adaptivity into 
a \c Problem.]

In the current example we demonstrate how to \e use existing 
refineable meshes and elements. Two further examples will demonstrate
how easy it is to \e create refineable meshes in 
<A HREF="../../../poisson/two_d_poisson_adapt/html/index.html">
domains with polygonal boundaries</A> and in
<A HREF="../../../poisson/fish_poisson2/html/index.html">
domains with curvilinear boundaries.</A>


<HR>   
<HR>


\section global Global parameters and functions
The namespace \c ConstSourceForPoisson only contains
the constant source function \f$ f(x)=-1 \f$. 

\dontinclude fish_poisson.cc
\skipline start_of_namespace
\until end of namespace


<HR>
<HR>

\section main The driver code

The main code is very short and calls two functions that
illustrate two different adaptation strategies: 
- A black-box approach in which the adaptation cycle 
  -# solve the problem on the initial, coarse mesh
  -# compute an error estimate
  -# adapt the mesh
  -# solve again
  .
  is performed automatically until the solution satisfies the 
  required error bounds (or until the maximum permitted number of 
  adaptation steps has been reached). 
- In the second approach we start by performing a number of uniform 
  mesh refinement steps, and then use incremental adaptations, allowing
  us to document how the refinement proceeds.
.

\dontinclude fish_poisson.cc
\skipline start_of_main
\until end of main

<HR>

\subsection black_box Black-box adaptation

We start by creating the Problem object, using the refineable
equivalent of the \c QPoissonElement -- the 
\c RefineableQPoissonElement, which is templated by the dimension and
the number
of nodes along the element's edges; the \c
RefineableQPoissonElement<2,3> is a nine-node (bi-quadratic) quad element.

\dontinclude fish_poisson.cc
\skipline start_black_box
\until RefineableFish


 After creating the \c DocInfo object, we document the (default)
adaptivity targets:

\skipline Setup labels
\until doc_ad

These include
- <B>The target for the maximum error:</B>
  Any elements whose error estimate exceed this value will be split
  into four "sons".
- <B>The target for the minimum error: </B> Any elements whose error
  estimate lies below this value are deemed to be unnecessarily
  small and are scheduled for (possible) unrefinement. [Elements can
  only be unrefined (i.e. merged with their "brothers") if their
  "brothers" are also scheduled for unrefinement.]
- <B>The minimum refinement level:</B> In problems with curvilinear
  domain boundaries it is often necessary to retain a reasonably accurate 
  representation of the boundary (e.g. for postprocessing
  purposes), even if the error estimate suggests that the 
  mesh could be unrefined further.
- <B>The maximum refinement level:</B> In problems where the solution
  has singularities, the refinement process would continue indefinitely,
  therefore an upper bound on the refinement level must be imposed.
- Finally, because unrefinement is done purely to speed up the computation,
  it would not make sense to adapt the mesh if this process would only
  remove a few elements, while forcing the re-computation of the solution
  on an only slightly coarsened mesh. Therefore, no mesh adaptation is
  performed if
  - the adaptation would only perform unrefinements
  - <B>and</B> the number of elements scheduled for unrefinement is below
    a certain threshold.
  .
.
These default parameters can be changed by the user; see \ref comments.

The fully-adaptive solution of the problem is very simple. We simply pass the
maximum number of adaptations to the Newton solver and document the
results. Done! 

\dontinclude fish_poisson.cc
\skipline with fully automatic adaptation
\until end black box

<HR>

\subsection incremental Incremental adaptation

To allow the user more control over the mesh adaptation process,
\c oomph-lib provides a number of functions that perform individual
adaptation steps without re-computing the solution immediately. This allows
the user to 
- perform uniform mesh refinement and unrefinement, 
- impose a specific refinement pattern,
- monitor/document the progress of the automatic adaptation. 
.

The second driver function illustrates some of these functions. We
start by setting up the problem, create the \c DocInfo object and document
the adaptivity targets, exactly as before:

\dontinclude fish_poisson.cc
\skipline start_of_incremental
\until doc_ad

Next, we solve the problem on the original, very coarse mesh
and document the result:

\dontinclude fish_poisson.cc
\skipline on the initial, very coarse
\until doc_info.number()++;


We know that the result is unlikely to be very accurate, 
so we apply three levels of uniform refinement, increasing the number
of elements from 4 to 256, and re-compute:

\dontinclude fish_poisson.cc
\skipline Do three rounds
\until doc_info.number()++;

The solution looks much smoother but we suspect that the corner
regions are still under-resolved. Therefore, we call the
\c Problem::adapt() function which computes an error estimate
for all elements and automatically performs a single 
mesh adaptation (refinement/unrefinement) step. If this adaptation changes the 
mesh, we recompute the solution, using the
"normal" Newton solver without automatic adaptation. 
We document the solution and  continue the adaptation 
cycle until \c Problem::adapt() ceases to change the mesh:



\dontinclude fish_poisson.cc
\skipline Now do (up to)
\until end of incremental

The progress of the adaptation is illustrated in the animated gif
at the beginning of this document. The first frame displays
the solution on the original four-element mesh; the next frame
shows the solution on the uniformly refined mesh; the final two
frames show the progress of the subsequent, error-estimate-driven
mesh adaptation.

<HR>
<HR>
 

\section problem The problem class
The problem class is virtually identical to that used in the 
<A HREF="../../../poisson/two_d_poisson/html/index.html">2D Poisson problem
without mesh refinement</A>. In the present problem, we leave the function 
\c Problem::actions_before_newton_solve() empty because the boundary 
conditions do not change. The function 
\c RefineableFishPoissonProblem::mesh_pt() overloads the
(virtual) function \c Problem::mesh_pt() since it returns a pointer
to a generic \c Mesh object, rather than a pointer to the specific
mesh used in this problem. This avoids explicit re-casts
in the rest of the code where member functions of the
specific mesh need to be accessed.



\dontinclude fish_poisson.cc
\skipline start_of_problem_class
\until end of problem class
 
<HR>
<HR>

\section constructor The Problem constructor
We start by creating the mesh, using \c oomph-lib's
\c RefineableFishMesh object:

\skipline start_of_constructor
\until RefineableFishMesh

Next, we create an error estimator for the problem. The \c
Z2ErrorEstimator is based on Zhu and Zienkiewicz's flux recovery
technique and can be used with all elements that 
are derived from the \c ElementWithZ2ErrorEstimator base class
(or with functions that implement the pure virtual functions that are defined
in this class) -- the \c RefineableQPoissonElement is an element of 
this type.


\skipline Create/set error estimator
\until new Z2Error

Next we pin the nodal values on all boundaries, apply the 
homogeneous Dirichlet boundary conditions,
pass the pointer to the source function to the elements,
and set up the equation
numbering scheme. 

\skipline Set the 
\until end of constructor


<HR>
<HR>

\section doc Post-processing
The post-processing routine writes the computed result to 
an output file, labeled with the identifiers specified
in the \c DocInfo object. 

\skipline start_of_doc
\until end of doc


<HR>
<HR>

\section comments Comments and Exercises

The purpose of this example was to provide a high-level overview
of \c oomph-lib's mesh adaptation procedures. 
We demonstrated that the implementation of full adaptivity 
only required us to
- replace the \c FishMesh and the 
  \c QPoissonElement objects by their refineable equivalents, 
  \c RefineableFishMesh  and \c RefineableQPoissonElement, respectively
- specify the error estimator, and
- specify the maximum number of adaptations for the black-box adaptive
  Newton solver.
. 

(Compare the Problem specification for the current problem
to that of its non-refineable equivalent, contained in the
alternative driver code 
<A HREF="../../../../demo_drivers/poisson/fish_poisson/fish_poisson_no_adapt.cc">
fish_poisson_no_adapt.cc.</A>)


Since most of the "hard work" involved in the mesh adaptation is
"hidden" from the user, we highlight some important aspects
of the procedure:

\subsection adapt_bc Automatic transfer of the solution/boundary conditions during the mesh adaptation
The \c Problem::adapt() function automatically determines the correct
boundary conditions for newly created nodes on the Mesh boundary; 
it automatically updates the equation numbering scheme, and 
interpolates the solution from the original mesh onto the adapted 
mesh. This is important in nonlinear problems where the provision 
of a good initial guess for the Newton iteration is vital;
and in time-dependent problems where the solution at one timestep
provides initial conditions for the next one. See the discussion of
<A HREF="../../../unsteady_heat/two_d_unsteady_heat_adapt/html/index.html">
the adaptive solution of the unsteady heat equation</A>
for more details.
Furthermore, the source function pointers are automatically
passed to an element's four "son" elements when the element
is subdivided. This allows the adaptation to proceed completely
automatically, without any intervention by the "user". On return
from \c Problem::adapt() the problem can immediately be re-solved.

In some special cases, certain actions may need to be performed 
before or after the mesh adaptation (e.g. if flux boundary
conditions are applied by \c  FaceElements; this is explained in
<A HREF="../../../poisson/two_d_poisson_flux_bc_adapt/html/index.html"> 
another example</A>). To ensure that these steps are performed when 
the adaptation is controlled by the "black-box" adaptive Newton solver, 
the \c Problem class provides the two empty virtual functions 
\code
Problem::actions_before_adapt()
\endcode
and 
\code
Problem::actions_after_adapt()
\endcode
which are called automatically before and after the adaptation.
The "user" can overload these in his/her specific \c Problem class
to implement such actions.


\subsection how_to Automatic mesh adaptation in domains with curvilinear boundaries
The mesh adaptation not only increases the number of elements 
but also produces a more accurate representation of the
curvilinear domain boundary -- new boundary nodes are
placed exactly onto the analytically-defined, curvilinear boundary, 
rather than on the boundaries of the "father" element, which only provides
an approximate representation of the exact domain boundary. 
This is achieved by employing a \c MacroElement-based representation
of the \c Domain -- we will discuss this in more detail in
<A HREF="../../../poisson/fish_poisson2/html/index.html">
another example.</A>

\subsection problem_vs_mesh Problem adaptation vs. Mesh adaptation
Many adaptation routines in the \c Problem class have 
equivalents in the \c RefineableMesh class. It is important
to appreciate the important differences between them: If adaptation
is performed at the \c Problem level, the adapted \c Problem is
fully functional, i.e. boundary conditions will have been
assigned for newly created nodes on the mesh boundary,
the equation numbering scheme will have been updated, etc.
The adapted \c Problem can therefore be re-solved immediately.
Conversely, if a mesh is refined directly, using the member functions
of the \c RefineableMesh class, many of these additional
tasks need to be performed "by hand" before the adapted
\c Problem can be resolved. 


<HR>

\subsection exercises Exercises
To familiarise yourself with \c oomph-lib's mesh adaptation procedures 
we suggest the following exercises:
-# When the Poisson problem is solved with the default refinement
   targets, no elements are unrefined. Increase the minimum
   permitted error from its default value of \f$ 10^{-5} \f$ to 
   \f$ 10^{-4} \f$ by adding the statement 
   \code 
   problem.mesh_pt()->min_permitted_error()=1.0e-4;
   \endcode
   before
   \code
   problem.mesh_pt()->doc_adaptivity_targets(cout);
   \endcode
   This value forces an unrefinement of several elements in the mesh:
\image html fish_poisson_soln2.gif "Plot of the solutions obtained with the modified adaptivity targets. " 
\image latex fish_poisson_soln2.eps "Plot of the solutions obtained with the modified adaptivity targets. " width=0.75\textwidth
-# Convince yourself that \c Problem::adapt() does indeed interpolate the
   solution from the coarse mesh to the fine mesh -- call
   \c Problem::doc_solution(...) before and after its execution.  
-# The \c Problem::refine_uniformly() function has a counterpart
    \c Problem::unrefine_uniformly(). Why does this function 
    not simply unrefine every single element in the mesh? 
    Explore the action of \c
    Problem::unrefine_uniformly() by plotting the solution
    before and after a few executions of this function.
\image html fish_poisson_soln3.gif "Uniform unrefinement " 
\image latex fish_poisson_soln3.eps "Uniform unrefinement " width=0.75\textwidth
-# Impose a "user-defined" refinement pattern
   by calling the function \c Problem::refine_selected_elements(...). 
.


<HR>
<HR>

\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:
<CENTER>
<A HREF="../../../../demo_drivers/poisson/fish_poisson/">
demo_drivers/poisson/fish_poisson/
</A>
</CENTER>
- The driver code is: 
<CENTER>
<A HREF="../../../../demo_drivers/poisson/fish_poisson/fish_poisson.cc">
demo_drivers/poisson/fish_poisson/fish_poisson.cc
</A>
</CENTER>
.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

