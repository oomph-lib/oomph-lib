/**

\mainpage Interfaces, Free Surfaces and Surface Transport: Theory and implementation

This document provides the theoretical background underlying \c oomph-lib's 
free surface, interface and surface transport capabilities. 
We begin with a review of the relevant
theory to establish the overall framework and notation, and then
discuss the implementation of the methodology in \c oomph-lib. 

Here is an overview of the structure of this document:
- \ref theory
  - \ref surface_representation
  - \ref surface_gradient
  - \ref boundary_conditions
  - \ref surface_transport
  .
- \ref fs_implementation
  - \ref fluid_interface
  - \ref line
  - \ref axi
  - \ref spine_formulation
  - \ref elastic_formulation
  .

.
If you don't need all the gory details, you may prefer to start by exploring the
free-surface flow tutorials in \c oomph-lib's 
<A HREF="../../../example_code_list/html/index.html#free_surface_nst">list of example driver codes.</A>

<HR>
<HR>

\section theory Theory

A complete theoretical treatment of moving surfaces in space requires the use of
differential geometry, described in detail by Aris
(1962) <em> Vectors,  Tensors  and  the Basic  Equations of Fluid
Mechanics </em>, Wetherburn (1955) <em> Differential Geometry of Three
Dimensions, Volume I </em> and Green \& Zerna (1968) <em> Theoretical
Elasticity </em> among others. Here, we shall present only the essential
information and the reader is referred to these texts to fill in
additional background.

<HR>

\subsection surface_representation Geometry of Surfaces

For our purposes, a surface is a object of dimension \f$n-1\f$ embedded in
an \f$n\!\f$-dimensional  Euclidean space, where \f$n\f$ is
two or three. We can always represent such a surface by a position
vector from our chosen origin \f$ \mathbf{R}^{*}(\zeta^{\alpha},t)
\f$, parametrised by time \f$t\f$ and intrinsic (surface) coordinates
\f$\zeta^{\alpha}\f$, where the Greek index \f$\alpha = 1,..,n-1\f$. 

Throughout this document we will use the summation convention 
that repeated Roman indices are to be summed over the range from 1 to
\f$n\f$ and repeated Greek indices are to be summed over the range
from 1 to \f$n-1\f$. We will retain the summation signs
for all other sums, such as sums over the nodes etc.

 The covariant base vectors of the surface are defined to be the partial
 derivatives of the position vector with respect to the surface
 coordinates. For a two-dimensional surface,
\f[ \mbox{\boldmath$g$}_{1} =
\frac{\partial\mbox{\boldmath$R$}}{\partial \zeta^{1}}, \quad \quad \quad
\mbox{\boldmath$g$}_{2} =
\frac{\partial\mbox{\boldmath$R$}}{\partial \zeta^{2}}. \f]
 The covariant metric tensor of the surface is formed by taking the dot product
of the covariant base vectors
\f[ g_{\alpha\beta} =  \mbox{\boldmath$g$}_{\alpha}\mbox{\boldmath$\cdot$}
 \mbox{\boldmath$g$}_{\beta}, \f] 
and the determinant of the covariant metric tensor is denoted by 
\f[ g = g_{11}g_{22} - g_{12}g_{21}. \f]
The contravariant metric tensor is the inverse of the covariant metric
 tensor and is denoted by \f$ g^{\alpha\beta} \f$. Thus,
 \f[ g^{11} = g_{22}/g, \quad g^{12} = -g_{12}/g, \quad g^{21} = -
 g_{21}/g, \quad g^{22} = g_{11}/g. \f]

 For a one-dimensional surface, there is a single covariant base
 vector that coincides with the tangent vector
\f[ \mbox{\boldmath$g$}_{1} = \frac{\partial
 \mbox{\boldmath$R$}}{\partial \zeta^1}.\f]
 Here, the covariant metric tensor is simply the inner product of the
 tangent vector with itself \f$ g_{11} = \mbox{\boldmath$g$}_{1}
 \mbox{\boldmath$\cdot g$}_{1}
 \f$ with determinant \f$ g = g_{11} \f$ 
 and the contravariant metric tensor is \f$ g^{11} = 1/g_{11}
 \f$.

 We will need to integrate quantities over the surface, which uses 
 the result that an infinitesimal unit of area is
 \f[ \mbox{ d}S = \sqrt{g} \mbox{ d}\zeta^{\alpha}, \f]
 in terms of the intrinsic coordinates.

\image html surface_sketch.gif "Sketch of section of a two-dimensional surface, S , in a three-dimensional space. The covariant base vectors are tangent to the surface and in the direction of the intrinsic surface coordinates. The intrinsic coordinates parametrise the position vector to material points in the surface. The surface is bounded by the curve B and the vector m is perpendicular to the outer unit normal of the surface and tangent to the bounding curve. " 
\image latex surface_sketch.eps "Sketch of section of a two-dimensional surface, S , in a three-dimensional space. The covariant base vectors are tangent to the surface and in the direction of the intrinsic surface coordinates. The intrinsic coordinates parametrise the position vector to material points in the surface. The surface is bounded by the curve B and the vector m is perpendicular to the outer unit normal of the surface and tangent to the bounding curve. " width=0.5\textwidth
 
<HR>

\subsection surface_gradient Differential Operators On A Surface

 The formulation of transport equations within the surface requires the
 rates of change of surface quantities. The appropriate derivative 
 is the surface gradient, often
 written as \f$ \mbox{\boldmath$\nabla$}_{\!\!_S} \f$ and defined
 unhelpfully in many papers and textbooks 
 as the gradient operator restricted to the surface or
 \f[ \mbox{\boldmath$\nabla$}_{\!\!S} \phi = \mbox{\boldmath$\nabla$}
 \phi - (\mbox{\boldmath$\nabla$}\phi \mbox{\boldmath$\cdot n$})
 \mbox{\boldmath$n$}.\f]
 The problem with the above definition 
 is that if a quantity is defined only on the surface it is
 impossible to take its gradient \f$\mbox{\boldmath$\nabla$}\f$.

 A more helpful definition in terms of the intrinsic coordinates is
 that 
\f[ \mbox{\boldmath$\nabla$}_{\!\!_S} \phi =
 g^{\alpha\beta} \mbox{\boldmath$g$}_{\alpha} \frac{\partial \phi}{\partial
 \zeta^{\beta}}. \f]
 For a one dimensional surface parametrised by the arc-length,
 \f$ s\f$, this reduces to \f$ \mbox{\boldmath$\nabla$}_{\!\!_S} \phi =
 \mbox{\boldmath$t$} \partial \phi / \partial s,\f$ where \f$ t\f$ is a
 unit tangent vector to the surface. By definition the surface
 gradient is tangent to the surface and has no normal component.

 The surface divergence of an \f$ n\!\f$-dimensional vector quantity 
defined on the surface is given by
\f[  \mbox{\boldmath$\nabla$}_{\!\!_S} \mbox{\boldmath$\cdot v$} =
 g^{\alpha\beta} \mbox{\boldmath$g$}_{\alpha} \mbox{\boldmath$\cdot$}
\frac{\partial \mbox{\boldmath$v$}}{\partial \zeta^{\beta}}. \f]
The divergence theorem applied over the surface (Aris, 1955) is
\f[ \int\!\!\!\int_{S} \mbox{\boldmath$\nabla$}_{\!\!_S} \mbox{\boldmath$\cdot
v$}_{t}\,\mbox{ d}S = \int_{B} \mbox{\boldmath$v$}_{t}
\mbox{\boldmath$\cdot m$}\,\mbox{d} l, \f]
where \f$ \mbox{\boldmath$v$}_{t}\f$ is a vector tangential to the
surface; and \f$ \mbox{\boldmath$m$}\f$ is the unit vector tangent 
to the surface,
but perpendicular to the tangent of the bounding curve \f$ B\f$, see the Figure
above. 

If a general vector is decomposed into normal and tangential
components, \f$ \mbox{\boldmath$v$} = \mbox{\boldmath$v$}_{t} + v_{n}
\mbox{\boldmath$n$} \f$, we can write
\f[ \int\!\!\!\int_{S} \mbox{\boldmath$\nabla$}_{\!\!_S}
\mbox{\boldmath$\cdot v$}\,\mbox{ d}S = 
 \int\!\!\!\int_{S} \mbox{\boldmath$\nabla$}_{\!\!_S}
\mbox{\boldmath$\cdot$}
\left(\mbox{\boldmath$v$}_{t} + v_{n} \mbox{\boldmath$n$}\right)\,\mbox{ d}S = 
\int\!\!\!\int_{S} \mbox{\boldmath$\nabla$}_{\!\!_S} \mbox{\boldmath$\cdot
v$}_{t}\,\mbox{ d}S + \int\!\!\!\int_{S}
\mbox{\boldmath$\nabla$}_{\!\!_S} \mbox{\boldmath$\cdot$} (v_{n}
\mbox{\boldmath$n$})\,\mbox{ d}S \f]
\f[ = \int\!\!\!\int_{S} \mbox{\boldmath$\nabla$}_{\!\!_S} \mbox{\boldmath$\cdot
v$}_{t}\,\mbox{ d}S\  + \int\!\!\!\int_{S}
 v_{n} \mbox{\boldmath$\nabla$}_{\!\!_S} \mbox{\boldmath$\cdot$}
 \mbox{\boldmath$n$} + \mbox{\boldmath$n\cdot\nabla$}_{\!\!_S}
 v_{n}\,\mbox{ d}S =  
\int\!\!\!\int_{S} \mbox{\boldmath$\nabla$}_{\!\!_S} \mbox{\boldmath$\cdot
v$}_{t}\,\mbox{ d}S\  - \int\!\!\!\int_{S}
 v_{n} \kappa\,\mbox{ d}S, \f]
where \f$\kappa\f$ is twice the mean curvature of
the surface and equal to minus the surface divergence of the normal. The
term \f$ \mbox{\boldmath$n\cdot\nabla$}_{\!\!_S} v_{n} = 0 \f$
because the normal and surface divergence of any scalar quantity 
are orthogonal. 

Hence using the divergence theorem on the first term on the right-hand
side, we obtain
\f[  \int\!\!\!\int_{S} \mbox{\boldmath$\nabla$}_{\!\!_S}
\mbox{\boldmath$\cdot v$}\,\mbox{ d}S
 =  - \int\!\!\!\int_{S}
 v_{n} \kappa\,\mbox{ d}S
+ \int_{B} \mbox{\boldmath$v$}_{t} \mbox{\boldmath$\cdot m$}\,\mbox{d}
l =  - \int\!\!\!\int_{S}
 v_{n} \kappa\,\mbox{ d}S +
\int_{B} \mbox{\boldmath$v\cdot
m$}\,\mbox{d} l, \quad (1) \f]
 because the normal component of \f$ \mbox{\boldmath$v$}\f$ is
 perpendicular to \f$ \mbox{\boldmath$m$}\f$ and therefore contributes
 nothing to the line integral.

If we now choose \f$ \mbox{\boldmath$v$} = \phi\,
\mbox{\boldmath$e$}_{i}\f$, where \f$ \mbox{\boldmath$e$}_{i}\f$ is
the unit base vector in the \f$ i\!\f$-th Cartesian coordinate
direction, then we obtain
\f[\int\!\!\!\int_{S} \mbox{\boldmath$\nabla$}_{\!\!_S} \phi\,
\mbox{\boldmath$\cdot e$}_{i}\,\mbox{ d}S\  = - \int\!\!\!\int_{S}
 \phi (\mbox{\boldmath$e$}_{i} \mbox{\boldmath$\cdot n$}) \kappa\,\mbox{ d}S +
\int_{B} \phi\, \mbox{\boldmath$e$}_{i} \mbox{\boldmath$\cdot m$}\,\mbox{d} l, \f]
 which is equivalent to the surface divergence theorem for a scalar
 field described by 
 Wetherburn (1955; p 240, Eqn 26) 
\f[
\int\!\!\!\int_{S} \kappa \phi \mbox{\boldmath$n$}\mbox{ d}S =
 \int_{B} \phi \mbox{\boldmath$m$} \mbox{ d} l -
\int\!\!\!\int_{S}\mbox{\boldmath$\nabla$}_{\!\!_{S}}\, \phi\mbox{ d}S.
\quad \quad (2) \f]


<HR> 

\subsection boundary_conditions Free Surface and Interface Boundary Conditions 

\subsubsection dyn_con Dynamic condition

\image html free_surface_sketch.gif "Sketch of the interface between two fluids. " 
\image latex free_surface_sketch.eps "Sketch of the interface between two fluids. " width=0.6\textwidth

 The presence of an interface with non-constant surface tension 
\f$ \sigma^{*} \f$
contributes to the
overall force balance, as also described in <a href="../../single_layer_free_surface/html/index.html"> another
tutorial </a>. The surface tension acts a line force bounding the
interface and acting in the direction \f$ \mbox{\boldmath$m$}\f$. Thus
using the surface divergence theorem in the form (2),
\f[
 \int_{B} \sigma^{*} \mbox{\boldmath$m$} \mbox{ d} l = 
\int\!\!\!\int_{S} \sigma^{*} \kappa + 
 \mbox{\boldmath$\nabla$}_{\!\!_{S}}\, \sigma^{*} \mbox{
 d}S, \f]
which gives the terms to be included in the force balance.

If we 
define the lower fluid in the sketch above to be fluid 1 and the upper fluid to
be fluid 2. The traction exerted by fluid 1 onto fluid 2, \f$
\mathbf{t}^{[1]*} \f$ and that exerted by fluid
2 onto fluid 1, \f$ \mathbf{t}^{[2]*} \f$. Then, balance of forces requires that
\f[ \mathbf{t}^{[1]*} - \mathbf{t}^{[2]*} \equiv 
\left[\left[ \mbox{\boldmath$\tau$}^{*}
\mbox{\boldmath$\cdot n$}^{[1]} \right]\right] = 
\sigma^{*} \, \kappa^{*} \, \mbox{\boldmath$n$}^{[1]} +
 \mbox{\boldmath$\nabla$}_{\!\!_{S}}\, \sigma^{*} \f]
where we have been explicit about the fact that the curvature is
dimensional and 
\f$ \kappa^{*} > 0 \f$ if the centre of curvature lies inside fluid 1. 

After using the non-dimensionalisation described in  <a href="../../single_layer_free_surface/html/index.html"> another
tutorial, </a> the boundary condition becomes
\f[\left[\left[ \mbox{\boldmath$\tau$}
\mbox{\boldmath$\cdot n$}^{[1]} \right]\right]
    = \frac{1}{Ca} \left[\sigma \, \kappa \, \mbox{\boldmath$n$} +
 \mbox{\boldmath$\nabla$}_{\!\!_{S}}\, \sigma  \right] \f]
where \f$ \sigma = \sigma^{*} / \sigma_{ref} \f$ and 
\f$ Ca = \mu_{ref}\, \mathcal{U} / \sigma_{ref},\f$
is the capillary number based on a reference viscosity and surface tension. 
This condition can be incorporated directly into the weak form of the
momentum equations because the surface integral term in these
equations, as described in <a href="../../rayleigh_traction_channel/html/index.html#traction_theory">
another tutorial </a>, 
is \f[ \int\!\!\!\int_{S} \mbox{\boldmath$n\cdot \tau\cdot  \psi$}^{(F)}\,\mbox{ d}S, \f]
where \f$ \mbox{\boldmath$\psi$}^{(F)}\f$ 
are the vector test functions associated with the
momentum equations. Note that it is more compact to work with vector
test functions rather than the Cartesian components of a vector test
function in this case. 
The terms in the integral are
exactly those on the left-hand side of the boundary condition
multiplied by the test function. Hence, they become
\f[ \int\!\!\!\int_{S} 
\frac{1}{Ca} \left[\sigma \, \kappa \, \mbox{\boldmath$n$} +
 \mbox{\boldmath$\nabla$}_{\!\!_{S}}\, \sigma
 \right]\mbox{\boldmath$\cdot\psi$}^{(F)}\, \mbox{ d}S. \f]

Early contributors computed these terms directly, but an undesirable
feature is that computation of the curvature requires taking second
derivatives of the position vector, which requires a higher degree of
smoothness than previously demanded. We can use the surface divergence
theorem again to weaken this requirement. Firstly we must use the
product rule to bring the test function into the surface divergence
\f[ \int\!\!\!\int_{S} 
\frac{1}{Ca} \left[\sigma \mbox{\boldmath$\psi$}^{(F)} \mbox{\boldmath$\cdot$}\, \kappa \, \mbox{\boldmath$n$}  +
 \mbox{\boldmath$\nabla$}_{\!\!_{S}}\, \mbox{\boldmath$\cdot$}\left(\sigma 
 \mbox{\boldmath$\psi$}^{(F)} \right)  - \sigma\,
 \mbox{\boldmath$\nabla$}_{\!\!_{S}}\,\mbox{\boldmath$\cdot\psi$}^{(F)} \right]\, \mbox{ d}S; \f]
and then we can use the surface divergence theorem (1)
on the first two terms to obtain
 \f[ \int_{B} \frac{1}{Ca} \sigma \mbox{\boldmath$\psi$}^{(F)}
 \mbox{\boldmath$\cdot m$}\,\mbox{ d}l 
 - \int\!\!\!\int_{S} \frac{1}{Ca}
\sigma
 \mbox{\boldmath$\nabla$}_{\!\!_{S}}\, \mbox{\boldmath$\cdot\psi$}^{(F)} \, \mbox{ d}S. \f]

In index notation these terms become
 \f[ \int_{B} \frac{1}{Ca} \sigma \psi_{i}^{(F)} m_{i}\,\mbox{ d}l 
 - \int\!\!\!\int_{S} \frac{1}{Ca}
\sigma g^{\alpha\beta} \Big[\mbox{\boldmath$g$}_{\alpha}\Big]_{i}
 \left[\mbox{\boldmath$\psi$}^{(F)}_{,\beta}\right]_{i} \, \mbox{ d}S, \f]
 where \f$ [\mbox{\boldmath$g$}_{\alpha}]_{i} \f$ denotes the \f$
 i\!\f$-th component of the covariant base vector \f$
 \mbox{\boldmath$g$}_{\alpha} \f$ and \f$ \sigma_{,\beta}\f$
 represents the derivative \f$ \partial \sigma / \partial
 \zeta^{\beta} \quad (3) \f$.
These are the terms that are implemented in \c oomph-lib's free
surface elements. Note that
variations in surface tension are taken into account in this
formulation without the explicit need to take its surface
 derivative. An important observation is that the intrinsic surface
 coordinates can be chosen to be the local coordinates of each element
 so that we do not need to introduce another set of coordinates.


\subsubsection kin_con Kinematic condition

 The kinematic condition is that "particles on the surface must
 remain on the surface". In other words the normal velocity of the
 surface must equal the normal rate of change of its position with
 time. The condition is compactly expressed in non-dimensional form as
 \f[ \left(\mbox{\boldmath$U$} - St \frac{\partial \mathbf{R}}{\partial
 t}\right) \mbox{\boldmath$\cdot n$} = 0, \f]
 where \f$ St\f$ is the Strouhal number, see <a href="../../single_layer_free_surface/html/index.html"> another
tutorial </a> for details. Note that, in general \f$ \mbox{\boldmath$U$}\f$ the
 velocity of the surface need <em> not </em> coincide with the
 velocity of the fluid.

 Although this equation must be satisfied the details of how exactly
 it is implemented depend crucially on the mesh-update strategy chosen
 and we make use of the C++ features of inheritance and templating to
 avoid code duplication, see below for details.

<HR> 

\subsection surface_transport Surface Transport Equations

 Consider a chemical species with surface concentration \f$ \Gamma\f$ 
 that is only present on the surface. 
 It can be transported within the moving surface by the
 standard mechanisms of advection and diffusion, but changes in the
 surface area can also induce changes in its concentration. The
 formulation of surface transport equations has been discussed many
 times in the literature and the main confusion surrounds how material
 derivatives are taken. Unlike the
 conventional bulk equations the surface does not occupy every point
 in the domain, so one cannot simply use the standard form of the
 material derivative.

 Rather than covering the literature, here we shall simply state,
 <em> ab initio </em>, the governing equations formulated by
 Huang, Lai \& Tseng as well as  Cermelli et al (2005), who claim it was
 established by Slattery (1972). We shall demonstrate that it is
 equivalent to the form stated by Wong \&
 Rumshitski and used by Campana et al (2004), but that it leads to a
 simpler formulation which avoids explicit calculation of the curvature
 (and hence second derivatives). 

 The dimensionless 
 governing equations in weak form governing the transport of a
 scalar quantity \f$ \Gamma \f$ are
 \f[ \int\!\!\!\int_{S} \left[ St \left(\frac{\partial \Gamma}{\partial t} -
   \dot{\mbox{\boldmath$R$}} \mbox{\boldmath$\cdot$} 
   \mbox{\boldmath$\nabla$}_{\!\!_{S}} \Gamma \right) +
   \mbox{\boldmath$\nabla$}_{\!\!_{S}} \mbox{\boldmath$\cdot$} 
   \left(\Gamma \mbox{\boldmath$U$}\right) -
   \frac{1}{Pe_{s}}
 \mbox{\boldmath$\nabla$}_{\!\!_{S}}\mbox{\boldmath$\cdot$} 
 \mbox{\boldmath$\nabla$}_{\!\!_{S}}
   \Gamma \right]\, \phi\, \mbox{d}S = 0. \f]
 In the above equation the time derivative is taken at fixed ``nodes''
 in the finite element formulation and the ALE-like term compensates
 for tangential movement of these nodes along the surface. The normal
 movement is enforced to be exactly the same as the surface velocity by
 the kinematic condition. 
 Note that the velocity in the third term of the governing equation 
 is the full (bulk) fluid velocity. The dimensionless quantity \f$
 Pe_s = \mathcal{U}\mathcal{L} / D_{s} \f$ is the surface Peclet
 number. 

 The formulation of Campana et al is found by decomposing the
 velocity in this term into normal and tangential components.
 \f[ \mbox{\boldmath$\nabla$}_{S} 
 \mbox{\boldmath$\cdot$} \left(\Gamma \mbox{\boldmath$U$}\right) = 
 \mbox{\boldmath$\nabla$}_{\!\!_{S}} \mbox{\boldmath$\cdot$}
 \left(\Gamma \mbox{\boldmath$U$}_{t} + \Gamma U_{n}
 \mbox{\boldmath$n$} \right). \f]
 The surface gradient yields a vector that is tangential to the
 surface so that its inner product with the unit normal,
 \f$ \mbox{\boldmath$n$}\f$ is
 zero. Thus,
 \f[ \mbox{\boldmath$\nabla$}_{\!\!_{S}} \mbox{\boldmath$\cdot$} 
 \left(\Gamma \mbox{\boldmath$U$}\right) = 
  \mbox{\boldmath$\nabla$}_{\!\!_{S}} \mbox{\boldmath$\cdot$}
  \left(\Gamma \mbox{\boldmath$U$}_{t}\right)  
  + \Gamma U_{n} \mbox{\boldmath$\nabla$}_{\!\!_{S}}
 \mbox{\boldmath$\cdot n$}, \f]
  which is the starting point for Campana <em> et al</em>'s
  derivation because the surface divergence of the normal may be
  replaced by the curvature, as described above.

 Returning to our formulation we use the surface divergence theorem
 (1) to
 integrate the diffusion term and the product rule to handle the third
 term:
 \f[ \int\!\!\!\int_{S} \left[ St \left(\frac{\partial \Gamma}{\partial t} -
   \dot{\mbox{\boldmath$R$}} \mbox{\boldmath$\cdot$} \mbox{\boldmath$\nabla$}_{\!\!_{S}} \Gamma \right) +
   \Gamma \mbox{\boldmath$\nabla$}_{\!\!_{S}} \mbox{\boldmath$\cdot$} \mbox{\boldmath$U$} + \mbox{\boldmath$U$} \mbox{\boldmath$\cdot$}
   \mbox{\boldmath$\nabla$}_{\!\!_{S}} \Gamma \right] \phi +
   \frac{1}{Pe_{s}} \mbox{\boldmath$\nabla$}_{\!\!_{S}} \Gamma \mbox{\boldmath$\cdot$} \mbox{\boldmath$\nabla$}_{\!\!_{S}}
   \phi\, \mbox{d}S  -
 \int_{B} \frac{1}{Pe_{s}} \mbox{\boldmath$\nabla$}_{\!\!_{S}} \Gamma
 \mbox{\boldmath$\cdot m$}\,  \phi \,\mbox{ d} l
= 0,\f]
 \f[ \Rightarrow \quad \int\!\!\!\int_{S} \left[ St\, \frac{\partial
       \Gamma}{\partial t} 
 + (\mbox{\boldmath$U$} -
   St\, \dot{\mbox{\boldmath$R$}}) \mbox{\boldmath$\cdot$} \mbox{\boldmath$\nabla$}_{\!\!_{S}} \Gamma +
   \Gamma \mbox{\boldmath$\nabla$}_{\!\!_{S}} \mbox{\boldmath$\cdot$} \mbox{\boldmath$U$} \right] \phi +
   \frac{1}{Pe_{s}} \mbox{\boldmath$\nabla$}_{\!\!_{S}} \Gamma \mbox{\boldmath$\cdot$} \mbox{\boldmath$\nabla$}_{\!\!_{S}}
   \phi\, \mbox{d}S 
 -
 \int_{B} \frac{1}{Pe_{s}} \mbox{\boldmath$\nabla$}_{\!\!_{S}} \Gamma
 \mbox{\boldmath$\cdot m$}\,  \phi \,\mbox{ d} l
= 0,\f]
The line term represents no-diffusive flux out of the system.
In index notation the equations are
 \f[ \int\!\!\!\int_{S} \left[ St\, \frac{\partial
       \Gamma}{\partial t} 
+ g^{\alpha\beta}
 [\mbox{\boldmath$g$}_{\alpha}]_{i}\left\{
 (U_{i} -  St\, \dot{R}_{i})  \Gamma_{,\beta}
 + \Gamma \left[\mbox{\boldmath$U$}_{,\beta}\right]_{i}\right\}
  \right] \phi +
   \frac{1}{Pe_{s}} g^{\alpha\beta} \Gamma_{,\alpha} \phi_{,\beta}\, \mbox{d}S 
 -
 \int_{B} \frac{1}{Pe_{s}} g^{\alpha\beta}
 [\mbox{\boldmath$g$}_{\alpha}]_{i} m_{i} \Gamma_{,\beta}\,  \phi \,\mbox{ d} l
= 0,\f]

These are the equations implemented in \c oomph-lib using the
definitions of surface derivatives given in \ref
surface_gradient. 

<HR>
<HR>

\section fs_implementation Implementation

We will now discuss how the discrete versions of the 
equations derived above are actually implemented in \c
oomph-lib. The basic idea is that the equations should be implemented
independently of the specific element type and mesh-update strategy and a 
base class oomph::FluidInterfaceElement defines the 
generic functionality for all fluid interface elements. 
The only difference between the different surface geometries are in
the definitions of the surface derivative operators and these are
defined in specific classes oomph::LineDerivatives (1D surface),
oomph::AxisymmetricDerivatives and
oomph::SurfaceDerivatives (2D surface).
The final specific element is created by using a special templated
class that determines the node-update strategy
and takes the base class, derivative class and bulk element
as template arguments.

<HR>

\subsection fluid_interface The FluidInterfaceElement class

The template-free oomph::FluidInterfaceElement class provides storage and member
functions that are common to all free-surface and interface
elements. The most important functions to be aware of are:
- Storage and access functions for (pointers to) the capillary and
Strouhal numbers.
- Storage for a (pointer to) an external pressure degree of freedom if
the boundary is a free surface, rather than an interface.
- The virtual function
\code
 double FluidInterfaceElement::compute_surface_derivatives(...)
 \endcode
  specifies how the surface gradient operators are computed.
- The function 
\code
 double FluidInterfaceElement::sigma(const Vector<double> &s_local) 
\endcode
 that returns the surface tension at the given local coordinate;
 default implementation returns 1.0.
- The function
\code
 FluidInterfaceElement::fill_in_generic_residual_contribution_interface(..)
\endcode
that is responsible for assembling the residual and jacobian
contributions corresponding to the dynamic and kinematic boundary conditions.
- The function
\code
 FluidInterfaceElement::add_additional_residual_contributions_interface(...)
\endcode
 which is called from <em> within </em> the integration loop and is
 used to assemble any additional surface transport equations or
 equations arising from different node update strategies. This
 function is virtual so that it can be overloaded in derived
 classes.
- The function 
\code
 virtual int FluidInterfaceElement::kinematic_local_eqn(...)
\endcode
 that is used to specify the local equation number used for the
 kinematic condition, which depends on the mesh-update strategy chosen.
.


\subsection line The LineDerivatives class

The class oomph::LineDerivatives
implements the specific surface derivatives
associated with a one-dimensional surface in a two-dimensional
domain.The global coordinate system is Cartesian so its base
vectors do not vary with the surface coordinates and 
\f$ [\mbox{\boldmath$\psi$}_{,\beta}^{(F)}]_{i} =
\psi^{(F)}_{i,\beta}\f$;
and the contribution to each component of the momentum equation is
found by setting the appropriate component \f$ \psi^{(F)}_{i} = 0\f$ for
\f$ i = 1,2 \f$.

\subsection axi The AxisymmetricDerivatives class

The class oomph::AxisymmetricDerivatives implements the specific residuals associated
with a two-dimensional surface in a three-dimensional domain, under
the assumption of axisymmetry. Thus, the coordinate system is
cylindrical polar \f$(r,z,\theta)\f$, but it is assumed that there are
no variations in the \f$ \theta \f$ direction.

It is worthwhile including the required mathematics here because the
terms are not the same as in the \c LineDerivatives
class. Specifically, if the surface coordinates are \f$ (\zeta^{1},
\zeta^{2}) = (s,\theta)\f$, the 
position vector to the surface is given by
 \f[ \mbox{\boldmath$R$} = \left(\begin{array}{c} r(s)\cos\theta \\
 r(s)\sin\theta \\ z(s) \end{array}\right) \quad\Rightarrow\quad
  \mbox{\boldmath$g$}_{1} = \left(\begin{array}{c} r'(s)\cos\theta \\
 r'(s)\sin\theta \\ z'(s) \end{array}\right), \quad
  \mbox{\boldmath$g$}_{2} = \left(\begin{array}{c} -r(s)\sin\theta \\
 r(s)\cos\theta \\ 0 \end{array}\right), \f]
 where \f$ r'(s) = \partial r / \partial s\f$ and \f$ z'(s) = \partial
 z / \partial s\f$. Hence,
 \f[ g_{11} = (r')^{2} + (z')^{2}, \quad g_{12} = g_{21} = 0, \quad
 g_{22} = r^{2}, \quad\mbox{and}\quad g = r^{2} \left[(r')^{2} +
 (z')^{2}\right].\f]
 In our standard formulation, the vector test function is given by
 \f[ \mbox{\boldmath$\psi$}^{(F)} = \left(\begin{array}{c}
 \psi_{r}^{(F)}(s) \cos\theta \\ \psi_{r}^{(F)}(s) \sin\theta \\
 \psi_{z}^{(F)}(s) \end{array} \right) \quad\Rightarrow\quad
\mbox{\boldmath$\psi$}^{(F)}_{,1} = \left(\begin{array}{c}
 (\psi')_{r}^{(F)}(s) \cos\theta \\ (\psi')_{r}^{(F)}(s) \sin\theta \\
 (\psi')_{z}^{(F)}(s) \end{array} \right) \quad\mbox{and}\quad 
 \mbox{\boldmath$\psi$}^{(F)}_{,2} = \left(\begin{array}{c}
 -\psi_{r}^{(F)}(s) \sin\theta \\ \psi_{r}^{(F)}(s) \cos\theta \\
 0 \end{array} \right).\f]
 
 Thus the contribution to the momentum equation (3) are
\f[ \int_{B} \frac{1}{Ca} \sigma \psi_{i}^{(F)} m_{i}\,\mbox{ d}l 
 - \int\!\!\!\int_{S} \frac{1}{Ca}
\sigma  \left[\frac{r' (\psi')_{r}^{(F)} + z' (\psi')_{z}^{(F)}}{(r')^{2} +
(z')^2} + \frac{1}{r} \psi_{r}^{(F)}\right]\, \mbox{ d}S, \f]
Two separate contributions are then derived from the cases \f$
\psi_{r}^{(F)} = 0 \f$ and \f$ \psi_{z}^{(F)} = 0\f$. The difference
from the \c LineDerivatives is the final term that accounts for the
azimuthal curvature. In addition, we must also multiply all terms by
the additional factor of \f$ r\f$ in the square-root of the determinant of
the surface metric tensor.

\subsection surface The SurfaceDerivatives class

The class oomph::SurfaceDerivatives implements the specific residuals associated
with a general two-dimensional surface in a three-dimensional domain. 
Once again, the global coordinate is Cartesian, so the
contribution to each momentum equation is found by setting the
two other components of the test function to be zero.

<HR>

\section spine_formulation The SpineLine/Axi/SurfaceFluidInterfaceElement classes

We shall discuss the "line" version of the elements, but
the others are essentially the same. 

The class oomph::SpineLineFluidInterfaceElement is templated by the
bulk element type, \c ELEMENT,  and inherits from \c
FluidInterfaceElement, LineDerivatives and \c
Hijacked<SpineElement<FaceGeometry<ELEMENT> > >. The hijacking is only
required for imposition of contact angle boundary conditions, see 
<a href="../../static_single_layer/html/index.html"> another tutorial
</a> for more details. 

Note that the use of templates to make the code generic makes it
hard to read. A simplified constructor is given below and
simply builds the element based on
the \c FaceGeometry of the bulk element and sets the indices
associated with the bulk fluid velocity components from the bulk
element.

\dontinclude specific_node_update_interface_elements.h
\skipline SpineLineFluidInterfaceElement(
\until {}


 If a spine method is used to update the nodal positions then the
spine height is the unknown associated with the kinematic
condition. Thus, the function \c kinematic_local_eqn(...) is
overloaded accordingly
\dontinclude specific_node_update_interface_elements.h
\skipline kinematic_local_eqn(
\until spine_local_eqn(n)

Finally, the element calculates the geometric contributions to the
jacobian using the generic functionality in \c ElementWithMovingNodes
\skipline fill_in_contribution_to_jacobian
\until End of jacobian contribution

There are no additional contributions to the residuals or jacobian.

<HR>

\section elastic_formulation The ElasticLine/Axi/SurfaceFluidInterfaceElement classes

We shall discuss the "line" version of the elements, but
the others are essentially the same. 

The class oomph::ElasticLineFluidInterfaceElement is templated by the
bulk element type, \c ELEMENT,  and inherits from \c
LineDerivatives, FluidInterfaceElement and \c
Hijacked<FaceGeometry<ELEMENT>  >. The hijacking is again only
required for imposition of contact angle boundary conditions.

The "elastic" versions of the elements are more complicated than the
"spine" versions because the kinematic condition is imposed using
Lagrange multipliers, following the method described by Cairncross et
al A finite element method for free surface flows of incompressible
fluids in three dimensions. Part I. Boundary fitted mesh motion'
(2000). These Lagrange multipliers must be added to the \c Nodes on
the free surface and their introduction adds additional terms to the
equations governing the bulk mesh motion.

The constructor is given below and in addition to building the element based on
the \c FaceGeometry of the bulk element and setting the indices
associated with the bulk fluid velocity components from the bulk
element, it also adds the additional storage required for the Lagrange
multipliers 

\dontinclude specific_node_update_interface_elements.h
\skipline ElasticLineFluidInterfaceElement(
\until //End of constructor

The kinematic boundary condition is associated with the Lagrange multiplier
and \c kinematic_local_eqn(...) is
overloaded accordingly
\dontinclude specific_node_update_interface_elements.h
\skipline Lagrange multiplier) 
\until }

The element calculates the geometric contributions to the
jacobian using the generic functionality in \c SolidElements
\skipline fill_in_contribution_to_jacobian
\until }

The additional contributions to the residuals and jacobian arise from
the Lagrange multiplier contributions to the equations that determine
the position of the nodes. The essential loop is the contribution
below which adds the normal traction to the governing equations of
solid mechanics:
\skipline Loop over the shape functions to assemble contributions
\until End of loop over shape functions 

\section surface Surface Transport Implementation

The recommended strategy for implementing surface transport equations
is to inherit from the appropriate \c FluidInterfaceElement, 
and include the equations independently of the mesh-update
strategy by overloading the function
\code
 FluidInterfaceElement::add_additional_residual_contributions_interface(...)
\endcode
The mesh-update specialisations can be added by further specialisation
as required.
Alternatively, one could simply inherit
directly from the specialised element. This is the approach taken in
this 
<a href="../../../multi_physics/rayleigh_instability_surfactant/html/index.html">
example code. </a>

An important point is that the axisymmetric formulation of the
surface divergence of a vector is not entirely trivial,
as discussed in \ref axi. In the context of surface transport
the term that is always present is \f$
\mbox{\boldmath$\nabla$}_{\!\!_{S}}\mbox{\boldmath$\cdot$}
\mbox{\boldmath$U$} \f$. Using the same surface coordinates as above,
\f$ (\zeta^{1},
\zeta^{2}) = (s,\theta)\f$, the velocity vector is
 \f[ \mbox{\boldmath$U$} = \left(\begin{array}{c} U_r(s)\cos\theta \\
 U_r(s)\sin\theta \\ U_z(s) \end{array}\right)
\quad\Rightarrow\quad
  \mbox{\boldmath$U$}_{,1} = \left(\begin{array}{c} \frac{\partial
  U_{r}}{\partial s}\cos\theta \\
 \frac{\partial U_{r}}{\partial s} \sin\theta \\ \frac{\partial
 U_z}{\partial s} \end{array}\right), \quad
  \mbox{\boldmath$U$}_{,2} = \left(\begin{array}{c} -U_r \sin\theta \\
 U_r\cos\theta \\ 0 \end{array}\right). \f]
Thus, the surface divergence term is
\f[ \mbox{\boldmath$\nabla$}_{\!\!_{S}} \mbox{\boldmath$\cdot$}
\mbox{\boldmath$U$}
 =  \frac{\left(\frac{\partial U_{r}}{\partial s} \frac{\partial
    r}{\partial s} + \frac{\partial U_{z}}{\partial s} \frac{\partial
   z}{\partial s}\right)}{(r')^{2} + (z')^{2}} + \frac{U_{r}}{r},\f]
 which is exactly the same form as the surface divergence of the
 vector test function derived in the section \ref axi, as expected.



<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

