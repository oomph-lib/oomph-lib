/**

\mainpage The deformation of a thin-shell material with a small strain, using the Kirchhoff-Love shell theory 

In this document, we discuss the solution of a simple two-dimensional thin elastic problem, using oomph-lib's Kirchhoff-Love shell elements. 

Specifically, in this document we demonstrate

- the descriptions of the governing equation of a thin-shell deformation with a small strain
.
and
- how to implement a shell problem with the \c BellElement and the \c C1CurvedElement in \c oomph-lib 
.

The reader is referred to <A HREF="../../../c1_element/bell_element/html/index.html">
the Bell triangular finite element</A> and <A HREF="../../../c1_element/curved_element/html/index.html">
the \f$ C^1 \f$-curved triangular finite element tutorials</A>, 
for more detailed descriptions on the \c BellElement and \c C1CurvedElement.

<HR>   
<HR>

\section overview Overview of a thin shell

A shell is defined as a thin three-dimensional elastic body where the thickness, \f$ h \f$, is smaller compared to the other two dimensions. Many analyses of thin shells neglect the effect of the transverse shear and follow the theory of Kirchhoff-Love. This theory states that a normal vector of the undeformed mid-surface remains normal to the deformed mid-surface throughout the deformation and it deforms inextensionally. 

Employing the Kirchhoff-Love assumption in a shell theory intends to reduce the dimension of a shell problem from the three-dimensional to the two-dimensional theory. Therefore, the shell governing equation which is derived from the principle of virtual displacement can be reduced to two-dimensional space. This is a result of allowing the integration in the coordinate perpendicular to the mid-surface to be carried out analytically. Therefore, all quantities can be expressed only on the two Lagrangian coordinates of the mid-surface.

Since the elastic material under consideration is a thin shell, the structure can experience large deflections and rotations, although strains and stresses may remain small. In such thin bodies, an assumption for small deformations (strains) is employed to simplify excessive deflections. This has been done by employing the linear (or Hookean) constitutive equation to represent the stress components as a linear function of all strain components.

<HR>

\subsection theory A thin-shell governing equation with a small strain

In this section, we illustrate the linearisation of the governing equation of a static shell in a general geometry. It will be seen that a thin-elastic shell can be modelled by equations defined on a two-dimensional domain. 

To develop a mathematical description of the deformation of a three-dimensional shell, we consider a shell geometry illustrated in the following figure. From the thinness feature of a shell, its geometry is allowed to specify by the two-dimensional reference surface and its thickness. We can choose the coordinate \f$ \xi_3 = 0 \f$ to be the shell's mid-surface and the reference surface. Therefore, we have that two coordinates \f$ \xi_1 \f$ and \f$ \xi_2 \f$ are located on the reference surface where the third coordinate \f$ \xi_3 \f$ is normal to the reference surface. The upper and lower surfaces of the shell have the coordinates \f$ \xi_3 = h/2 \f$ and \f$ \xi_3 = -h/2 \f$, respectively, where \f$ h \f$ is the thickness of the shell.

\image html shell_geometry.gif "A shell geometry. " 
\image latex shell_geometry.eps "A shell geometry. " width=0.45\textwidth

The governing equation of a static shell with zero pre-stress is parametrised by the two coordinates \f$ \xi_1, \xi_2 \f$ on the mid-surface of the shell and can be obtained as 

\f[
 \int_{\Omega}{\left\{ \tilde{E}^{\alpha\beta\gamma\delta}\left( \gamma_{\alpha\beta}\delta\gamma_{\gamma\delta} + \frac{1}{12}h^2\kappa_{\alpha\beta}\delta\kappa_{\gamma\delta}\right) -
\frac{1}{h}\sqrt{\frac{A}{a}}\mathbf{f} \cdotp \delta\mathbf{R}\right\}\sqrt{a}}d\xi_1d\xi_2=0, \ \ \ \ (1) 
\f]

where \f$ \tilde{E}^{\alpha\beta\gamma\delta} \f$ represents the plane stress stiffness tensor defined on the mid-surface and can be computed by 

\f[
\tilde{E}^{\alpha\beta\gamma\delta} = \frac{E}{2(1+\nu)}\left( t^{\alpha\gamma}t^{\beta\delta} + t^{\alpha\delta}t^{\beta\gamma} + \frac{2\nu}{1-\nu}t^{\alpha\beta}t^{\gamma\delta}\right), \ \ \ \ (2) 
\f]

where \f$ \nu \f$ denotes the Poisson's ratio. This equation is presented with a non-dimensional form and is described in details in <A HREF="../../../beam/tensioned_string/html/index.html"> another tutorial</A>. Non-dimensionalisation of quantities is also explained here.

Note that in the small strain regime, an area element between the undeformed and deformed configurations are indistinguishable. Hence, the external force is preferable to express on the area element of the undeformed midplane which is related with that of the deformed midplane as
\f[
  \hat{\mathbf{f}}\sqrt{a} = \mathbf{f}\sqrt{A},
\f]
where \f$ \bf{f}, \hat{\bf{f}} \f$ denote the external force expressed on the area element of the undeformed and the deformed midplane, respectively. 

In order to obtain a linearised version of the governing equation of a thin-shell deformation, all terms in the equation (1) have to be linearised. A small strain is then assumed. Therefore, we have that \f$ \forall\epsilon << 1 \f$,

\f[ \mathbf{f} = \epsilon\mathbf{\tilde{f}}, \ \ \ \mathbf{u} = \epsilon\mathbf{\tilde{u}}, \f]

where a displacement \f$ \bf{u} \f$ will be considered in the tangential and normal directions (rather than the Cartesian system) in this study. Therefore, a displacement \f$ \mathbf{u} \f$ on the midplane which is parametrised by the two-dimensional coordinates, \f$ \xi_1, \xi_2 \f$ can be decomposed into two tangential and normal components as

\f[ \mathbf{u} = u^{j}\mathbf{t}_{j}, \f]

where the covariant base vectors \f$ \mathbf{t}_{1}, \mathbf{t}_{2} \f$ are tangent in direction of coordinate lines \f$ \xi_{1}, \xi_{2}\f$, respectively, and \f$ \mathbf{\hat{t}}_{3} \f$ is a unit normal vector to the undeformed mid-surface. Coefficients \f$ u^{j}; j=1,2,3\f$ are associated components of a displacement \f$ \mathbf{u} \f$ in two tangential and one normal directions. 

The linear version of the undeformed covariant metric tensor of the mid-surface is expressed as

\f[ a_{\alpha\beta} = \mathbf{t}_{\alpha}\cdotp \mathbf{t}_{\beta}, \f]
and, the linear version of the deformed metric tensor can be considered as 

\f[ A_{\alpha\beta} = a_{\alpha\beta} + u^{\beta}|_{\alpha} + u^{\alpha}|_{\beta}, 
\f]
where the quantity \f$ u^{j}|_{\alpha} = (u^{j}_{,\alpha}+u^{i}\Gamma^{j}_{i\alpha}) \f$ is the \f$ j \f$-component of \f$ \alpha \f$-derivative of a displacement \f$ \mathbf{u} \f$ in the tangential and normal directions. Furthermore, \f$ \Gamma^{j}_{i\alpha} = \left(\frac{\partial\mathbf{t}_{i}}{\partial\xi_{\alpha}}\cdotp\mathbf{t}_{j}\right) \f$ denotes Christoffel symbol of the second kind. Note that the comma preceding the subscript \f$ j \f$ signifies partial differentiation with respect to the coordinate line \f$ \xi_j. \f$ Also, a Latin index represents any of the numbers \f$ 1,2,3 \f$ and a Greek index represents the numbers \f$ 1,2. \f$

The determinants of the covariant metric tensor in undeformed and deformed shell are \f$ a \f$ and \f$ A, \f$ respectively, and can be calculated by the determinant of the associated metric tensors. 

Then, the linearised strain tensor is expressed as 

\f[
	\gamma_{\alpha\beta} = \frac{1}{2}\left(A_{\alpha\beta} - a_{\alpha\beta} \right) = \frac{1}{2}\left( u^{\beta}|_{\alpha} + u^{\alpha}|_{\beta} \right).
\f]

Next, the linear version of the curvature tensor in the undeformed and deformed configurations can be specified as, respectively,

\f[ b_{\alpha\beta} = \mathbf{\hat{t}}_3\cdotp\mathbf{r}_{,\alpha\beta}, \f]

and,

\f[ B_{\alpha\beta} = b_{\alpha\beta} + \left[ \frac{1}{a}L_{j}\Gamma^j_{\alpha\beta} -\frac{L_{3}}{a^2}\Gamma^3_{\alpha\beta}\right] +  \left( u^{i}|_{\alpha}\Gamma^{3}_{i\beta} + u^{3}_{,\alpha\beta} + u^{k}\frac{\partial \Gamma^{3}_{k\alpha}}{\partial\xi_{\beta}} + \frac{\partial u^{k}}{\partial\xi_{\beta}}\Gamma^{3}_{k\alpha}\right). \f]

Then, the linearised bending tensor can be obtained as 
\f[
	\kappa_{\alpha\beta} = b_{\alpha\beta} - B_{\alpha\beta} = - \left[ \frac{1}{a}L_{j}\Gamma^{j}_{\alpha\beta}-\frac{L_{3}}{a^2}\Gamma^{3}_{\alpha\beta}\right] -  \left( u^{i}|_{\alpha}\Gamma^{3}_{i\beta} + u^{3}_{,\alpha\beta} + u^{k}\frac{\partial \Gamma^{3}_{k\alpha}}{\partial\xi_{\beta}} + \frac{\partial u^{k}}{\partial\xi_{\beta}}\Gamma^{3}_{k\alpha}\right),
\f]

where \f$ L_{j} \f$ is defined as

\f[ L = \left( a_{12} u^{3}|_{2} - a_{22} u^{3}|_{1}, -a_{11} u^{3}|_{2} + a_{21} u^{3}|_{1}, a_{11} u^{2}|_{2} - a_{12} u^{1}|_{2} + a_{22} u^{1}|_{1} - a_{21} u^{2}|_{1} \right)^T.
\f]

<HR>

\subsection repr A finite element representation of the displacements

When we substitute all linearised terms of the strain tensor, \f$ \gamma_{\alpha\beta} \f$ and the bending tensor, \f$ \kappa_{\alpha\beta} \f$ and their variations back into the shell governing equation (1), it can be seen that the governing equation contains a second-order derivative of a normal displacement and a first-order derivative of tangential displacements in both directions. Therefore, we have that the normal displacement requires \f$ C^1 \f$-continuity while the tangential displacements in both directions require only \f$ C^0 \f$-continuity.

To approximate the normal displacement, a \f$ C^1 \f$-continuous interpolation function has to be considered in order to ensure the continuity of its derivatives in the finite element method. In \c oomph-lib, the Bell shape functions can be employed to provide the \f$ C^1 \f$-continuity when the straight boundary domain is concerned. Alternatively, the \f$ C^1 \f$-curved triangular shape functions can be used when the curvilinear boundary is concerned. The Bell and the \f$ C^1 \f$-curved triangular shape functions can be overloaded from \c BellElementShape<> and \c C1CurvedElementShape<>, respectively.

Unlike the normal displacement, an interpolation function approximating the solution of the tangential displacements does not require a continuity for its derivatives. Therefore, a Lagrange-type interpolation function will be employed to approximate the tangential displacements, \f$ u_1, u_2. \f$ These interpolation functions will be overloaded from \c TElementShape<>.

\c Oomph-lib's \c BellShellElement and \c oomph-lib's \c C1CurvedShellElement provide a discretisation of the variational principle (1) with two-dimensional, subparametric, triangular finite elements on a straight-sided and curvilinear boundaries, respectively. In these elements, the displacements are regarded as unknowns, and the \f$ C^1 \f$-interpolation is used to interpolate the normal displacement  while the Lagrange polynomials is used to interpolate the tangential displacements. Furthermore, the geometry is approximated by the linear Lagrange interpolations for a straight-sided boundary domain while the cubic polynomial is employed to approximate the curved boundary domain.

<HR>
<HR>

\section reslt Numerical results

In this section, implementations of the linearised governing equations for thin-shell problems will be illustrated. There are three problems considered in this document; square- and circular-plate, and circular tube bending. The implementations will be based on the governing equations (1) that we derived in section 1.1.

In all cases, the problems will be solved with the assumption that the thickness is thin so that the linear theory can be applied. Our choice of thickness is 0.01. Also, applied forces will be applied in normal direction to the undeformed surface with no initial stress.

Furthermore, in order to perform the finite element implementations, the domains of interest will be discretised by triangular elements with an unstructured mesh. In this study, \c Oomph-lib's \c BellShellElement is employed when straight-sided boundaries are concerned while \c Oomph-lib's \c C1CurvedShellElement is employed when curvilinear boundaries are involved.

Note that figures of the displacements that we will illustrate throughout this document will be in the Cartesian coordinate system. Since the displacements obtained as the solution of (1) are in the tangential and normal coordinates, the transformation to the Cartesian coordinate system has to be done by

\f[ u_x = u^1t^1_1 + u^2t^1_2 + u^3\hat{t}^1_3, \ \ \ \ u_y = u^1t^2_1 + u^2t^2_2 + u^3\hat{t}^2_3, \ \ \ \  u_z = u^1t^3_1 + u^2t^3_2 + u^3\hat{t}^3_3, \f]

where \f$ t^j_i; j=1,2,3, \f$ denote components of the tangent base vector \f$ \mathbf{t}_i; i=1,2 \f$, and the unit normal vector, \f$ \mathbf{\hat{t}}_3 \f$, respectively. 

<HR>
\subsection plate The square-plate bending problem

Here, we will consider a deformation of a flat plate which is subjected to a pressure loading on its upper surface as illustrated in the following figure. The length of the plate in both directions is 1. The boundary conditions in this problem are two clamped boundaries, \f$ \xi_1=0, \xi_1=1 \f$, and two free edges, \f$ \xi_2=0, \xi_2=1 \f$. Therefore, the displacement and its rotational degrees of freedom are pinned at the boundaries \f$ \xi_1=0, \xi_1=1 \f$ while all degrees of freedom are set to be free at \f$ \xi_2=0, \xi_2=1 \f$. 
 
\image html plate_geometry.gif "The geometry of the square plate with two clamped edges and two free edges. Forces applied to a body are uniform in the outward normal direction. " 
\image latex plate_geometry.eps "The geometry of the square plate with two clamped edges and two free edges. Forces applied to a body are uniform in the outward normal direction. " width=0.55\textwidth

In order to implement the finite element method of a two-dimensional space in this study, the domains of interest which is the unit square will be discretised by triangles. Note that the unstructured mesh contains 150 elements. 

Here are figures illustrate displacements in all directions in Cartesian coordinates system for the flat plate problem stated above with the applied loads in the normal direction equal to \f$ 1.0\times 10^{-8} \f$. 

\image html Plate_x_direction.gif "The displacement in x-direction. " 
\image latex Plate_x_direction.eps "The displacement in x-direction. " width=0.65\textwidth

\image html Plate_y_direction.gif "The displacement in y-direction. " 
\image latex Plate_y_direction.eps "The displacement in y-direction. " width=0.65\textwidth

\image html Plate_z_direction.gif "The displacement in z-direction. " 
\image latex Plate_z_direction.eps "The displacement in z-direction. " width=0.65\textwidth

Regarding the displacements in both tangential directions, it can be seen from Figures 1.3 and 1.4 that no deformation occurs in the \f$ x \f$- and \f$ y \f$-directions. The underlying reason is that the forces are applied in the normal direction to the surface of the plate which correspond to the \f$ z \f$-direction. Hence, there is no force applied in both tangential directions which correspond to the \f$ x \f$- and \f$ y \f$-directions. Therefore, there is no contribution to make the body deforms in those directions as the linear governing equations are not coupled between displacements in each direction. To see this, the reader have to expand (1) after substituting \f$ \gamma_{\alpha\beta} \f$, \f$ \kappa_{\alpha\beta} \f$, and their variations in. 

<HR>

\subsection tube The circular-tube problem

In this section, we consider a deformation of a circular tube which is subjected to a pressure loading on its surface as illustrated in the following figure. The loads applied on the tube are uniformly distributed in the normal direction. 

\image html circular_tube.gif "The geometry of the square plate with two clamped edges and two free edges. " 
\image latex circular_tube.eps "The geometry of the square plate with two clamped edges and two free edges. " width=0.55\textwidth

To implement the deformation of the circular tube in this
study, the quarter of a circular tube will be implemented and symmetric conditions are assumed along the tube. The boundary conditions determined in this problem are clamped supports at both ends of the tube.

Similar to the plate bending problem, the domains of interest will be discretised by an unstructured triangular mesh with 248 elements. 

Here are figures illustrate displacements in all directions in Cartesian coordinates system for the circular tube problem stated above with the applied loads in the normal direction equal to \f$ 1.0\times 10^{-6} \f$. 

\image html Tube_x_direction.gif "The displacement in x-direction. " 
\image latex Tube_x_direction.eps "The displacement in x-direction. " width=0.65\textwidth

\image html Tube_y_direction.gif "The displacement in y-direction. " 
\image latex Tube_y_direction.eps "The displacement in y-direction. " width=0.65\textwidth

\image html Tube_z_direction.gif "The displacement in z-direction. " 
\image latex Tube_z_direction.eps "The displacement in z-direction. " width=0.65\textwidth

It can be seen that the displacement in \f$ x \f$- and \f$ y \f$-directions are symmetry. This behaviour depicts that the thin-circular tube deforms axisymmetrically within a small-strain regime. 

<HR>

\subsection plate2 The circular-plate bending problem

Here, we will consider a deformation of a flat plate which is
subjected to a pressure loading on its upper surface like in section 1.2.1. However, the domain of interest considered here will be curved. Therefore, the unit circular plate is considered with clamped boundaries. Note that only one quarter of the unit circular plate is analysed and symmetric conditions are applied in this problem.

In order to implement the finite element method of a two-dimensional space in this study, the domains of interest which is the unit circular plate will be discretised by triangles. Note that the unstructured mesh contains 84 elements. 

Here are figures illustrate displacements in all directions in Cartesian coordinates system for the flat plate problem stated above with the applied loads in the normal direction equal to \f$ 1.0\times 10^{-7} \f$. 

\image html Cir_Plate_x_direction.gif "The displacement in x-direction. " 
\image latex Cir_Plate_x_direction.eps "The displacement in x-direction. " width=0.65\textwidth

\image html Cir_Plate_y_direction.gif "The displacement in y-direction. " 
\image latex Cir_Plate_y_direction.eps "The displacement in y-direction. " width=0.65\textwidth

\image html Cir_Plate_z_direction.gif "The displacement in z-direction. " 
\image latex Cir_Plate_z_direction.eps "The displacement in z-direction. " width=0.65\textwidth

Similarly, there is no deformation in both x- and y-directions as explained in section 1.2.1.

<HR>
<HR>

\section impl Implementation in oomph-lib

In the following, we illustrate the driver codes for the square-plate
bending problem while other problems can be determined similarly.

<HR>

\subsection global Global parameters and functions
The namespace \c Physical_Variables is where the source
function and the exact solution are defined. The source function can
be specified via \c Physical_Variables::source_function() while the exact solutions are defined via \c Physical_Variables::get_exact_u(). Note that the six exact solutions correspond to the six degrees of freedom defined on each node. Furthermore, the applied source functions are required to be in the direction of the unit normal vector.

\dontinclude unstructured_clamped_square_plate.cc
\skipline start_of_namespace
\until end of namespace 


<HR>

\subsection main The driver code

The driver code is very simple and short. It is where the problem is
defined. In this study, the problem is constructed using the
unstructured mesh with triangular elements in 2D. A number of nodes in an element has to be specific as a template parameter in
the problem set up. This is crucial in order to take care of element
nodes when the number of nodes on the element defined differently for each interpolation functions. 

Normally, in the problem that \f$ C^1 \f$-shape functions are the only interpolation functions used to approximate the variable space, the number of \c NNODE_1D remains 2 as required by the \f$ C^1 \f$-shape functions (see <A HREF="../../../c1_element/bell_element/html/index.html">
the Bell triangular finite element</A> and <A HREF="../../../c1_element/curved_element/html/index.html">
the \f$ C^1 \f$-curved triangular finite element tutorials</A>). However, in a linear shell problem, different shape functions are used to interpolate displacements in different directions as different order of continuity is required in their governing equations. 

In this study, the displacements in tangential directions are chosen to approximate by the quadratic Lagrange interpolations which provide only \f$ C^0 \f$ continuity and are defined to have 3 nodes per side on a triangle. On the other hand, the approximation of the normal displacement employs the \f$ C^1 \f$-interpolation functions which are defined to have 2 nodes per side on a triangle. 

Therefore, the number of \c NNODE_1D has been modified
in the \c BellShellElement<DIM,NNODE_1D> (and \c C1CurvedShellElement when dealing with the curvilinear boundary domain) to be 3 in this problem. Consequencely, the extra nodes that are not necessary for the \f$ C^1 \f$-shape functions have to be taken care. 

Following the usual self-test, we call the function \c MyLinearisedShellProblem::parameter_study() to compute the solution of the problem within a range of external pressures.

\dontinclude unstructured_clamped_square_plate.cc
\skipline start_of_main
\until end of main

<HR>

\subsection problem The problem class

The problem class has five member functions, illustrated as follows:

- The problem constructor
- \c action_before_newton_solve() : Update the problem specifications before solve. Boundary conditions maybe set here.
- \c action_after_newton_solve() : Update the problem specifications after solve.
- \c doc_solution() : Pass the number of the case considered, so that output files can be distinguished.
- \c parameter_study() : Computes the shell's deformation for a range of external pressures.
.

From the above mentioned functions, only the problem constructor is non-trivial. The reader is referred to <A HREF="../../../beam/tensioned_string/html/index.html"> another tutorial</A>  for a description on \c parameter_study.

In the present problem, the function \c Problem::actions_after_newton_solve() is not required, so it remains empty. Also, the class includes two private data members which store pointers to a source function and to the geometric object that specifies the shell's undeformed shape.

\dontinclude unstructured_clamped_square_plate.cc
\skipline start_of_problem_class
\until end of problem class

<HR>

\subsection constructor The Problem constructor

The problem constructor starts by overloading the function \c
Problem::mesh_pt() and set to the specific mesh used in the
problem. In this tutorial, we implement the problem with 2D triangular unstructured mesh which is externally created by \c Triangle. The generated output will be used to build \c oomph-lib mesh. The reader may refer to <A HREF="../../../meshes/mesh_from_inline_triangle/html/index.html"> another tutorial</A> to create an unstructured triangular mesh internally.

\dontinclude unstructured_clamped_square_plate.cc
\skipline start_of_constructor
\until TriangleMesh

We then create the undeformed centreline of the shell to one of an oomph-lib's standard shell geometric objects.

\skipline set the undeformed shell
\until Plate

Prior to consider the boundary conditions, we will illustrate how to take care extra nodes on the element for the \f$ C^1 \f$-interpolations. Since there is no degree of freedom of the \f$ C^1 \f$-interpolations defined on the mid-sided nodes, they have to be pinned.

In order to do this, firstly, we loop over the element and pin all degrees of freedom on non-vertex nodes that associated with the \f$ C^1 \f$-interpolations. Since nodes in the element situate anticlockwise and the midside nodes fill in progressing along the consecutive edges, the pinning procedure is easily done by starting to pin from the third node and so on. 

\skipline pinning the middle nodes
\until end of the middle node pinning

Next, the boundary conditions of the problem will be taken care. We
pin the nodal values on the boundaries where the boundary conditions applied. Note that, at the clamped boundaries, all second-order derivatives degrees of freedom are also pinned in order to reduce the number of degrees of freedom in the problem. These second-order derivatives are the derived boundary conditions that can be taken care by the natural boundary conditions.

\skipline start_of_boundary_conditions
\until end of boundary conditions

We then loop over the elements and set the pointer to the physical
parameters (if any), the function pointer to the source function, and the pointer to the geometric object that specifies the undeformed surface of the shell.

\skipline Loop over elements
\until end of loop over elements

We finish the constructor by assigning the equation numbering scheme.
\skipline Setup equation numbering scheme
\until end of constructor

<HR>

\subsection actionbefore Action before newton solve

In the \c action_before_newton_solve(), the problem specifications
 will be updated before performing the newton solve. The boundary
 values will be (re)set from the exact solutions. 

\dontinclude unstructured_clamped_square_plate.cc
\skipline start_of_actions_before_newton_solve
\until end of actions before solve

<HR>
<HR>

\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/shell/">
demo_drivers/shell/
</A>
</CENTER>\n

- The driver code for the square plate bending problem is: \n\n
<CENTER>
<A HREF="
../../../../
demo_drivers/shell/plate/unstructured_clamped_square_plate.cc
">
demo_drivers/shell/plate/unstructured_clamped_square_plate.cc
</A>
</CENTER>\n

- The driver code for the circular tube problem is: \n\n
<CENTER>
<A HREF="
../../../../
demo_drivers/shell/clamped_shell/unstructured_clamped_curved_shell.cc
">
demo_drivers/shell/clamped_shell/unstructured_clamped_curved_shell.cc
</A>
</CENTER>\n

- The driver code for the circular plate bending problem is: \n\n
<CENTER>
<A HREF="
../../../../
demo_drivers/shell/circular_plate/unstructured_clamped_circular_plate.cc
">
demo_drivers/shell/circular_plate/unstructured_clamped_circular_plate.cc
</A>
</CENTER>
.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

