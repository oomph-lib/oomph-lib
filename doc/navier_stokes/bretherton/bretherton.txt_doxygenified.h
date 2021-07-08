/**

\mainpage The Bretherton problem: An air finger propagates into a 2D channel

This tutorial discusses the solution of the classical
Bretherton problem in which an air-finger is driven steadily into a 2D
fluid-filled, rigid-walled channel. 


The driver code listed below was originally developed as a 
test-bed for the paper 

<CENTER>
Hazel, A.L. & Heil, M. (2008) The influence of gravity on the steady 
propagation of a semi-infinite bubble into a flexible channel. 
Physics of Fluids 20, 092109
<a href="http://www.maths.manchester.ac.uk/~mheil/MATTHIAS/PDF/HazelHeilPoF2008.pdf">(pdf reprint)</a>
</CENTER>

in which we studied the elastic-walled equivalent, paying particular
attention to the effect of transverse gravity which plays quite an important
role in this problem.

Compared to other tutorials, there is no detailed discussion of the
driver code itself because the implementation is somewhat messy.
However, given that the driver code is (of course!) very well documented
you should be able to figure out what's going on once you've read
through the discussion of the problem formulation and our brief
comments on the \ref implement.
<a href="../../../contact/html/index.html">Get in touch</a>
if you need more information.


\section problem The problem


The sketch below shows the problem setup: An (inviscid) air
finger propagates at a steady speed \f$ U \f$ into a 2D rigid-walled,
fluid-filled channel of half-width \f$ H_0^* \f$. 
In the process it displaces some of the viscous fluid (of
density \f$ \rho \f$, viscosity
\f$ \mu \f$ and surface tension \f$ \sigma^* \f$)
and deposits a film of thickness \f$ h_0^* \f$ on the channel walls.
[Note that, for the sake of simplicity, we ignore the effect of 
transverse gravity in this discussion; the driver code listed 
below allows the specification of a gravitational body force which will break 
the up-down symmetry of the solution.]


\image html bretherton_sketch.gif "Sketch of the problem setup " 
\image latex bretherton_sketch.eps "Sketch of the problem setup " width=0.75\textwidth


We non-dimensionalise the velocities on \f$ U \f$, the lengths
on \f$ H_0^* \f$ and the pressure on the viscous scale, \f$ \mu U /
H_0^* \f$ and solve the problem in a frame moving with the 
velocity of the finger. The flow is then governed by the
non-dimensional Navier-Stokes equations
\f[
Re \ u_j\frac{\partial u_i}{\partial x_j} = -\frac{\partial p}{\partial x_i} + \
\frac{\partial }{\partial x_j }\left( \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}  \right)
\f]
and the continuity equation
\f[ 
\frac{\partial u_i}{\partial x_i}=0, 
\f] 
where the Reynolds number is defined as \f$ Re=U H_0 \rho/\mu \f$. 
On the free fluid surface, whose outer unit normal we 
denote by \f$ {\bf n} \f$, the fluid normal velocity vanishes,
\f[
{\bf u} \cdot {\bf n} = 0 \mbox{\ \ \ on the air-liquid interface}.
\f]
We measure the fluid pressure \f$ p \f$ relative to the bubble
pressure 
\f$ p_b \f$. Then
the dynamic boundary condition implies that
\f[
-p n_i + \left( \frac{\partial u_i}{\partial x_j} + 
\frac{\partial u_j}{\partial x_i} \right) n_j + \frac{1}{Ca} \
\kappa_f n_i  = 0
\mbox{\ \ \ on the air-liquid interface,}
\f]
where \f$ \kappa_f = \kappa_f^* H_0 \f$ is the non-dimensional interface
curvature and the capillary number \f$ Ca = U\mu/\sigma^* \f$ is 
a non-dimensional measure
of the bubble velocity.
The no-slip condition on the channel walls requires that \f$ {\bf
u}=(-1,0) \f$ at \f$ x_2=\pm 1 \f$. Far behind the bubble tip, the fluid
remains stationary on the walls. In the moving frame of reference, this 
corresponds to the plug flow profile \f$ {\bf u}=(-1,0) \f$ as \f$
x_1\to -\infty \f$.
The residual film thickness \f$ h_0 \f$ has to be determined as part of the
solution and it determines the volume flux through the system.
The Poiseuille velocity profile far ahead of the bubble tip, 
\f[
{\bf u}=(-1+\frac{3}{2}(1-h_0)(1-x_2^2),0) \mbox{\ \ \ as $ x_1\to \infty$,}
\ \ \ \ \ \ \ \ (1)
\f] 
is then determined by overall continuity.  



<HR>
<HR>

\section results Some results

  Here are some pretty pictures of the computed flow field,


\image html bretherton_flow_with_streamlines.gif "Velocity and pressure fields " 
\image latex bretherton_flow_with_streamlines.eps "Velocity and pressure fields " width=0.75\textwidth

and here is a comparison of the computational predictions for the
bubble pressure and film thickness against Bretherton's famous
asymptotic predictions (valid only for small capillary numbers!):


\image html bretherton_trace.gif "Bubble pressure and film thickness far behind the finger tip vs. capillary number; computed solution and Bretherton's predictions for small Ca. " 
\image latex bretherton_trace.eps "Bubble pressure and film thickness far behind the finger tip vs. capillary number; computed solution and Bretherton's predictions for small Ca. " width=0.75\textwidth


<HR>
<HR>

\section implement Implementation

The discretisation of the problem is reasonably straightforward,
the most (human-)time-consuming part being the creation of the spine mesh,
discussed in <a href="../../../navier_stokes/spine_channel/html/index.html">
another tutorial</a>. The real complication arises from the fact
that the application of the "inflow condition" (1) at
the right end of the domain -- superficially a Dirichlet condition for the
velocity -- depends on the non-dimensional film thickness 
\f$ h_0 \f$ at the left end of the domain. This introduces a
non-local interaction between these degrees of freedom.
We handle this by providing a templated wrapper class 
\c BrethertonElement (templated by the type of the "wrapped" fluid
element) which allows the specification of the film thickness as
external \c Data (i.e. \c Data whose values affect the element's residuals but
are not determined by the element; see 
the discussion of \c oomph-lib's various types of elemental \c Data 
in <a href="../../../the_data_structure/html/index.html">the "bottom up"
discussion of the library's data structure</a>). Read the driver
code listed below to see how it works!


<HR>
<HR>

\section driver The driver code

\include bretherton.cc

\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/bretherton">
demo_drivers/navier_stokes/bretherton
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/navier_stokes/bretherton/bretherton.cc">
demo_drivers/navier_stokes/bretherton/bretherton.cc</A>
</CENTER>


<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

