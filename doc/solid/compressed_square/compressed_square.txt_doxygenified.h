/**

\mainpage Demo problem: Compressible and incompressible behaviour

The purpose of this tutorial is to discuss the use of compressible
and incompressible constitutive equations for solid mechanics
problems. As discussed in the 
<a href="../../solid_theory/html/index.html#elastic_constitutive">
Solid Mechanics Theory Tutorial</a>, problems in which the solid
is truly incompressible require a solution based on the 
pressure/displacement formulation of the principle of 
virtual displacements (PVD). Mathematically, this is because 
the pressure acts as the Lagrange multiplier that enforces 
incompressibility. The pressure/displacement form of the
PVD may be discretised with finite elements that employ
continuous (e.g. Taylor-Hood-type) or discontinuous
(e.g. Crouzeix-Raviart-type) elements. 
 
  Some constitutive equations allow for compressible and
incompressible behaviour, depending on their parameters.
Problems involving such constitutive equations may be solved 
with large number of formulations. As an example, consider \c
oomph-lib's generalised Hookean constitutive equation
(see the <a href="#disclaimer">disclaimer</a> below).


<b> I. The displacement form </b>

\c oomph-lib's generalised Hookean constitutive law assumes that 
the 2nd Piola Kirchhoff stress \f$ \sigma^{ij} \f$ (non-dimensionalised on 
Young's modulus \f$ E \f$) is related to Green's strain tensor 
\f$ \gamma_{ij} = 1/2 \ (G_{ij} - g_{ij}) \f$ via
\f[
\sigma^{ij} = \frac{1}{2(1+\nu)}
\left( 
G^{ik} G^{jl} + G^{il} G^{jk}  +
\frac{2\nu}{1-2\nu} G^{ij} G^{kl} 
\right) \gamma_{kl}.
\hspace{3cm} (1)
\f]
Here \f$ g_{ij} \f$ and \f$ G_{ij} \f$ are the metric tensors 
associated with the undeformed and deformed configurations. 
This constitutive law reduces to the
the classical version of Hooke's law for small strains when
\f$ G_{ij} \to g_{ij}. \f$ In the above form, the constitutive
law can be used directly in the displacement-based form of the PVD, 
unless \f$ \nu = 1/2 \f$ -- the case that corresponds to 
incompressible behaviour in the small-displacement regime. 
(While the above formulation only breaks down completely
when \f$ \nu = 1/2, \f$ numerical solutions based on the
above form of the constitutive equation become badly-behaved
as \f$ \nu \f$ approaches that value.)


<b> II. The pressure-displacement form for compressible and
incompressible behaviour</b>

To avoid the problems that arise as \f$ \nu \to 1/2, \f$ the constitutive 
equation may be rewritten in the form
\f[
\sigma^{ij} =  \overline{\sigma}^{ij} - p \ G^{ij} \ , 
\f]
where the deviatoric stress is given by
\f[ 
\overline{\sigma}^{ij} = 
\frac{1}{2(1+\nu)}
\bigg( 
G^{ik} G^{jl} + G^{il} G^{jk} \bigg) \gamma_{kl}.
\f]
To remain consistent with (1), the pressure
\f$ p \f$ must then be determined via the equation
\f[
p \ \frac{1}{\kappa} + d = 0, \hspace{4cm} (2)
\f]
where 
\f[
d = G^{ij} \gamma_{ij}
\f]
is the generalised dilatation (which reduces to the actual dilatation
in the small-strain limit). The inverse bulk modulus is defined as
\f[
\frac{1}{\kappa} =  \frac{(1-2\nu)(1+\nu)}{\nu}
\f]
and tends to zero as \f$ \nu \to 0. \f$ The alternative 
form of the constitutive equation can therefore be used for
any value of \f$ \nu \f$. If \f$ \nu = 1/2, \f$ the 
constraint (2) enforces \f$ d=0, \f$ i.e 
incompressible behaviour (at least in the
small displacement regime -- for large deflections the constraint
 (2) simply ensures that the generalised dilatation 
vanishes; this may not be a physically meaningful constraint). 

<b> III. Truly incompressible behaviour</b>

Finally, we may retain the decomposition of the stress into its
deviatoric and non-deviatoric parts but determine the pressure
from the actual incompressibility constraint, i.e. replace
(2) by
\f[
\det G_{ij} - \det g_{ij} = 0. \hspace{5cm} (3)
\f]
In this case, the deviatoric part of the stress is determined
by the material's constitutive parameters, i.e. its Young's modulus
and its Poisson ratio while the pressure (acting as the Lagrange
multiplier for the constraint (3)) 
enforces true incompressibility.

<a name="disclaimer">

<table border=1>
<tr>
<td bgcolor="cornsilk">
<center>
<b>Disclaimer</b>
</center>
We wish to stress that not all the combinations discussed above
are necessarily physically meaningful. For instance, as already
acknowledged,  setting \f$ \nu = 1/2 \f$ in formulation II does not 
ensure true incompressibility (in the sense of (3)) 
if the solid undergoes large deflections. Similarly, setting 
\f$ \nu \f$ to a value that differs from \f$ 1/2 \f$ while 
enforcing true incompressibility via formulation III would be 
extremely odd. Furthermore, it is not clear (to us) if our
generalisation of Hooke's law (obtained by replacing the undeformed
metric tensor \f$ g_{ij} \f$ by its deformed counterpart \f$ G_{ij}
\f$ yields a constitutive equation that is particularly useful 
for any specific material. However, the same is true for any other
constitutive equation -- not all of them are useful for all materials!
The important point is that all constitutive equations that allow 
for compressible and incompressible behaviour contain combinations 
of constitutive parameters that blow up as incompressibility
is approached. Formulation II shows how to express such 
constitutive laws in a way that can be used for any value of 
the constitutive parameters. 
</td>
</tr>
</TABLE>

 
With this disclaimer in mind, we will now demonstrate the use 
of the various combinations outlined above in a simple test problem
for which an exact (linearised) solution is available. 

<HR>
<HR>
 
\section problem The problem 
Here is a sketch of the problem: Three faces of a square elastic body are
surrounded by "slippery" rigid walls that allow the body to slide freely
along them. The body is subject to a vertical, gravitational
body force, acting in the negative \f$ y \f$ -direction.

\image html sketch.gif "Sketch of the problem. " 
\image latex sketch.eps "Sketch of the problem. " width=0.75\textwidth

We choose the height of the square as the reference length for
the non-dimensionalisation by setting 
\f$ {\cal L} = L, \f$ and use the characteristic stiffness associated with the
body's constitutive equation to scale the stresses and the body
forces. For instance, for linear elastic behaviour, we choose
the reference stress to be the solid's Young's modulus, thus setting
\f$ {\cal S} = E. \f$ The non-dimensional body force 
is then given by \f$ {\bf f} = - g \ {\bf e}_y, \f$ where
\f[
g = \frac{\rho g^* L}{E}
\f]
indicates the magnitude of the gravitational load relative to the
body's stiffness;
see the <a href="../../solid_theory/html/index.html#non-dim_solid">
Solid Mechanics Theory Tutorial</a> for full details on the 
non-dimensionalisation of the governing equations. 



 
<HR>
<HR> 

\section exact An exact solution

Assuming weak loading, i.e. \f$ g \ll 1 \f$, the body will undergo
small deflections and its deformation will be governed by the equations
of linear elasticity. Given that the walls bounding the solid 
are slippery, it is easy to show that a displacement field has
the form \f$ {\bf u} = v(y)\ {\bf e}_y \f$. Inserting this
ansatz into the Navier-Lame equations shows that the non-dimensional
displacement in the vertical direction is given by
\f[
v(y) = g \ \frac{(1+\nu)(1-2\nu)}{(1-\nu)} 
\left( \frac{1}{2} y^2 - y \right).
\f]
 The stresses are given by
\f[
\sigma_{xx} = \frac{g \ \nu}{(1-\nu)} (y-1) \ , \ \ \ \
\sigma_{xy} = 0 \ ,  
\mbox{ \ \ \ and \ \ \ }
\sigma_{yy} = g \ (y-1) \ .
\f]
This shows that, as the solid approaches incompressibility,
\f$ \nu \to 1/2, \f$ the displacements are suppressed and
the stress distribution becomes hydrostatic with \f$ \sigma_{xx} \to 
\sigma_{yy}. \f$

<HR>
<HR>

\section results Results

Here is a plot of the three non-zero quantities (the vertical 
displacement \f$ v \f$, and the horizontal and vertical stresses,
\f$ \sigma_{xx} \f$ and \f$ \sigma_{yy} \f$, respectively) for
\f$ g = 10^{-2}  \f$ and \f$ \nu =0.45. \f$ The shaded surface 
shows the exact solution while the mesh shows the finite element
solution obtained with the displacement formulation of the problem.
The solutions obtained with formulations II and III
are graphically indistinguishable. All other quantities
(the transverse displacement and the shear stress \f$ \sigma_{xy} \f$)
are of order \f$ 10^{-12}. \f$ 

\image html validate_vs_lin.gif "Plot of the vertical displacement and the horizontal and vertical stresses, respectively. " 
\image latex validate_vs_lin.eps "Plot of the vertical displacement and the horizontal and vertical stresses, respectively. " width=0.75\textwidth
 
The driver code listed below also computes solutions for various 
other values of \f$ \nu \f$ (including the incompressible case
\f$ \nu=1/2 \f$),  and  demonstrates how to enforce
"true" incompressibility via formulation III. Furthermore, it
computes a solution based on the incompressible Mooney-Rivlin
constitutive law. The agreement between analytical and
computed results is as good as the one shown in the above plot. 
In particular in all cases where incompressibility is enforced, the vertical
displacement is suppressed and the stress state becomes
hydrostatic, as predicted by the analytical solution. 

<HR>
<HR>

\section wrapper Modifying the element's output function
Since we wish to compare the stresses and displacements
against the analytical solution, we modify the \c SolidElement's
output function by means of a wrapper class that overloads the 
\c SolidElement::output(...).
 
\dontinclude compressed_square.cc
\skipline start_wrapper
\until };
<HR>
<HR>


\section namespace Problem Parameters
As usual we define the various problem parameters in a 
global namespace. We prepare a pointer to a constitutive equation
and define Poisson's ratio for use with the generalised Hookean 
constitutive law. 

\skipline start_namespace
\until Nu

We also prepare a pointer to a strain energy function and 
define the coefficients for the Mooney Rivlin law:

\until C2

Finally, we define the gravitational body force.

\until end namespace


<HR>
<HR>

\section main The driver code

The driver code solves the problem with a large number of
different formulations and constitutive equations. We start
with the generalised Hookean constitutive equation and
consider three different values of Poisson's ratio, corresponding
to compressible, near-incompressible and incompressible behaviour:


\dontinclude compressed_square.cc
\skipline start_of_main
\until ::Nu);

First, we solve the problem in the pure displacement formulation,
using (the wrapped version of the) displacement-based 
\c QPVDElements. (As discussed above, the displacement formulation
cannot be used for \f$ \nu=1/2 \f$.)

\until }

Next we consider the pressure/displacement formulation with
continuous pressures (Taylor-Hood), using the  
\c QPVDElementWithContinuousPressure element.

\until }

We suppress the listing of the remaining combinations (see the driver code
<A HREF="../../../../demo_drivers/solid/compressed_square/compressed_square.cc">
compressed_square.cc
</A>for details):
- Pressure/displacement formulation with discontinuous pressures
  (Crouzeix-Raviart), using the \c QPVDElementWithPressure element. 
- Pressure/displacement formulation with continuous pressures
  (Taylor-Hood), using the \c QPVDElementWithContinuousPressure
  element, with true incompressibility enforced via formulation
  III. 
- Pressure/displacement formulation with discontinuous pressures
  (Crouzeix-Raviart), using the \c QPVDElementWithPressure element, 
  with true incompressibility enforced via formulation
  III.
.

Before the end of the loop over the different \f$ \nu \f$ values we 
delete the constitutive equation, allowing it to be re-built 
with a different Poisson ratio.

\skipline Clean up
\until end generalised Hookean

Next we build the strain-energy-based Mooney-Rivlin constitutive law

\until Global_Physical_Variables::Strain_energy_function_pt);

and solve the problem with \c QPVDElementWithContinuousPressure
and \c QPVDElementWithPressure elements, enforcing true incompressibility 
enforced via formulation III by setting the \c incompressible flag
to true.

\until end of main
<HR>
<HR>

\section class The Problem class
The \c Problem class has the usual member functions. The \c i_case
label is used to distinguish different cases while the boolean
\c incompressible indicates if we wish to enforce incompressibility
via formulation III. 

\dontinclude compressed_square.cc
\skipline begin_problem
\until };

<HR>
<HR>

\section constructor The Problem constructor
We start by building the mesh -- the \c SolidMesh version of the
\c ElasticRectangularQuadMesh. 

\skipline start_of_constructor
\until new Elastic

We complete the build of the elements by specifying the constitutive
equation and the gravitational body force. 

\until ::gravity;

If the element is based on the pressure/displacement form of the principle of
virtual displacements we enforce (true) incompressibility if
required. Note that, by default, \c oomph-lib's
pressure/displacement-based solid mechanics elements do not 
assume incompressibility, i.e. the default is to use formulation II.

\until end compressibility 

Finally, we pick a control node to document the solid's load-displacement
characteristics and apply the boundary conditions (no displacements
normal to the "slippery" walls) before setting up the equation
numbering scheme.

\until end of constructor


<HR>
<HR>

\section doc Post-processing

The post-processing routine documents the load-displacement
characteristics in a trace file and outputs the deformed domain shape.

\until endl;

We then output the exact solution of the linearised equations
using the same format as in the element's overloaded output function.
 

\until end doc


<HR>
<HR>

\section comm_ex Comments and Exercises

As usual, \c oomph-lib provides self-tests that assess if the 
enforcement incompressibility (or the lack thereof) is consistent:
- The compiler will not allow the user to enforce incompressibility
  on elements that are based on the displacement form of the
  principle of virtual displacements. This is because
  the displacement-based elements do not have a member function
  \c incompressible().
  \n\n
- Certain constitutive laws, such as the Mooney-Rivlin law used in
  the present example require an incompressible formulation. If \c oomph-lib
  is compiled with the \c PARANOID flag, an error is thrown if
  such a constitutive law is used by an element for which
  incompressibility has not been requested. 
  \n\n
  <center>
  <b>Recall that the default setting is not to enforce
  incompressibility!</b>
  </center>
  \n
  If the library is
  compiled without the \c PARANOID flag no warning will be issued
  but the results will be "wrong" at least in the sense that
  the material does not behave like an incompressible Mooney-Rivlin
  solid. In fact, it is likely that the Newton solver will
  diverge.  Anyway, as we keep saying, without the \c PARANOID flag, 
  you're on your own!
  \n\n
.
You should experiment with different combinations of constitutive laws and
element types to familiarise yourself with these issues.

<HR>
<HR>

\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER> 
<A HREF="../../../../demo_drivers/solid/compressed_square/">
demo_drivers/solid/compressed_square/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/solid/compressed_square/compressed_square.cc">
demo_drivers/solid/compressed_square/compressed_square.cc
</A>
</CENTER>
.











<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

