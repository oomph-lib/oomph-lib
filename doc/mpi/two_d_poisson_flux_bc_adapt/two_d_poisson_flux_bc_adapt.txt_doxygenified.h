/**

\mainpage Parallel solution of a 2D Poisson problem with flux boundary conditions
 
This document provides an overview of how to distribute the 
<a href="../../../poisson/two_d_poisson_flux_bc_adapt/html/index.html">2D 
Poisson problem with flux boundary conditions</a>. 
It is part of a
<a href="../../../example_code_list/html/index.html#distributed">
series of tutorials</a>
that discuss how to modify existing serial driver codes so that the
\c Problem object can be distributed across multiple processors.

A feature
of this problem is that the flux boundary conditions are applied by
attaching "flux elements" (derived from the \c FaceElement 
base class) to the "bulk elements" adjacent to the appropriate mesh 
boundary.  As discussed in the 
<a href="../../../poisson/two_d_poisson_flux_bc_adapt/html/index.html">
tutorial for the serial driver code</a>, the \c FaceElements are 
not involved in any adaptation within the bulk mesh. 
Instead, they are detached 
before the bulk mesh is adapted and re-attached afterwards, which
 ensures that the \c FaceElements are only attached
to bulk elements present in the adapted mesh.

The same issue arises when the \c Problem is distributed:
all \c FaceElements must be attached 
before the problem is distributed to allow \c METIS to 
analyse the interaction between face and bulk elements correctly.
However, after the \c Problem has been distributed, some of the bulk 
elements on each processor will have been deleted, leaving the
corresponding \c FaceElement dangling. To deal with such problems,
\c oomph-lib provides the empty virtual functions

\code
Problem::actions_before_distribute() 
\endcode

and

\code
Problem::actions_after_distribute() 
\endcode

which are called automatically by \c Problem::distribute(...).
Specifically, \c Problem::actions_before_distribute() is called 
\b after the problem distribution has been determined by 
\c METIS but \b before the actual distribution (during which
elements are deleted) takes place. \c
Problem::actions_after_distribute() is called after the problem
distribution is complete.


In the present problem we 
overload the functions \c Problem::actions_before_distribute() and \c 
Problem::actions_after_distribute() to perform the same functions as \c 
actions_before_adapt() (i.e. delete the flux elements) and \c
actions_after_adapt() (i.e. re-attach the flux elements).  
We note that any \c FaceElement that is attached to a halo element in the bulk
mesh becomes a halo element itself; see the
<a href="../../general_mpi/html/index.html#face_elements">general 
MPI tutorial</a> for further details. 


Most of driver code is identical to its serial counterpart and
we only discuss the changes required to distribute the problem. Please
 refer to 
<a href="../../../poisson/two_d_poisson_flux_bc_adapt/html/index.html">
another tutorial</a> for a more detailed discussion of the
problem and its (serial) implementation. 


<HR>

\section main_body The main function

The only changes required to the main function are the usual calls to 
initialise and finalise \c oomph-lib's MPI routines and a single call
to \c Problem::distribute() after the problem has been
constructed. The source code is actually slightly more complicated because the
distribution is read in from a file so that the driver can be used as a
self-test. Note that the file must specify 
the partition for \b all elements, 
including the \c FaceElements. (We refer to 
<a href="../../adaptive_driven_cavity/html/index.html#no_disk">
another tutorial</a>
for details on how to create the distribution file.) 

<HR> 

\section problem_class The problem class

The only additions to the problem class are the functions \c
actions_before_distribute() and \c actions_after_distribute().  As
explained above, these perform exactly the same functions as \c
actions_before_adapt() and \c actions_after_adapt(), respectively.

\dontinclude two_d_poisson_flux_bc_adapt.cc
\skipline Actions before distribute:
\until }

\skipline Actions after distribute:
\until }

<HR> 

\section doc_solution The doc_solution() function

As with other driver codes, the output files are modified to allow
each processor to output its elements into files that include the
processor number.

\skipline start_of_doc
\until end of doc

<HR>

The remainder of this driver code is unchanged from the 
<a href="../../../poisson/two_d_poisson_flux_bc_adapt/html/index.html">serial
version</a>.



<HR>
<HR>

\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/mpi/distribution/two_d_poisson_flux_bc_adapt/">
demo_drivers/mpi/distribution/two_d_poisson_flux_bc_adapt/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/mpi/distribution/two_d_poisson_flux_bc_adapt/two_d_poisson_flux_bc_adapt.cc">
demo_drivers/mpi/distribution/two_d_poisson_flux_bc_adapt/two_d_poisson_flux_bc_adapt.cc
</A>
</CENTER>
.





<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

