/**

\mainpage Parallel solution of the FSI channel with leaflet problem: Distributing problems with algebraic node updates

This document provides an overview of how to distribute 
fluid-structure interaction problems in which algebraic 
node update methods are used to deform the fluid mesh in response
to changes in the shape of the domain boundary. It is part of a
<a href="../../../example_code_list/html/index.html#distributed">
series of tutorials</a>
that discuss how to modify existing serial driver codes so that the
\c Problem object can be distributed across multiple processors.

As discussed in the 
<a href="../../general_mpi/html/index.html#alg_node_update">general 
MPI tutorial</a>, the parallel implementation of algebraic node update methods
for specific meshes is greatly facilitated if the \c GeomObject
that describes the motion of the moving domain boundary is available
on all processors. In FSI problems, the moving boundary 
is typically represented by a \c MeshAsGeomObject -- a compound
\c GeomObject formed from the lower-dimensional mesh of 
\c SolidElements that define the moving boundary of the fluid domain.
To ensure that the \c MeshAsGeomObject remains available on all
processors when the underlying mesh is distributed, we use the
function \c Mesh::keep_all_elements_as_halos() to indicate that
all elements should be retained on all processors. 


<CENTER>
<TABLE BORDER=1, WIDTH=500px>
<TR>
<TD bgcolor="cornsilk">
<CENTER>
\b Note
</CENTER>
 The procedure described here is appropriate \b only 
if the solid mesh used to
create the \c MeshAsGeomObject is \b not \b adapted during the solution of
the problem. We refer to 
<a href="../../turek_flag/html/index.html">another tutorial</a> 
for instructions on how to deal with the case when the solid mesh is
itself adapted.
</TD>
</TR>
</TABLE>
</CENTER>

<HR>

\section fsi_channel_with_leaflet Revisiting the FSI channel with leaflet problem

We demonstrate the methodology for the problem of 
<a href="../../../interaction/fsi_channel_with_leaflet/html/index.html">
flow past an elastic leaflet in a 2D channel</a>.
Most of driver code is identical to its serial counterpart
and we only discuss the changes required to distribute the problem.
Please refer to 
<a href="../../../interaction/fsi_channel_with_leaflet/html/index.html">
another tutorial</a> for a more detailed discussion of the
problem and its (serial) implementation.  

<HR>

\subsection main_body The main function

As with other parallel driver codes, the changes from the
<a href="../../../interaction/fsi_channel_with_leaflet/html/index.html">
serial driver code </a>
are extremely modest, essentially including the initialisation and
shutdown of MPI.
Once the problem is set up, we call \c Problem::distribute(...),
using the boolean flag \c report_stats to indicate that
the statistics of the distribution should be reported on screen.

\dontinclude fsi_channel_with_leaflet.cc
\skipline Distribute problem using METIS
\until problem.distribute(

<HR>

\subsection problem_class The problem class

The only additional function is \c actions_after_distribute(),
described <a href="#actions_routines">later in this tutorial</a>.

<HR>

\subsection constructor The problem constructor

 The only difference from the serial counterpart is that we  
insist that all beam elements are kept as halos when the mesh 
is distributed.

\dontinclude fsi_channel_with_leaflet.cc
\skipline Discretise leaflet
\until keep_all_elements_as_halos()

<HR>

\subsection action_after_adapt Actions after adaptation

The \c actions_after_adapt() function requires only a minor change
from the serial version --- a simple re-ordering of the sequence
in which the various steps are performed. In serial it does not matter
in which order the (re-)assignment of the auxiliary node update
functions and (re-)setup of the fluid-structure interaction are
performed. In the parallel version, however, 
the assignment of the auxiliary node update functions must take place
\b after the (re-)setup of the fluid-structure interaction. This is
because the function 
\c FSI_functions::setup_fluid_load_info_for_solid_elements(...)
creates halo copies of the fluid elements (and their nodes!)
if the required fluid element
is not present on the current processor. 
Any newly-created halo fluid nodes on the FSI boundary 
are accessible via the 
usual boundary lookup schemes and 
must be told about the auxiliary node update function which applies
the no-slip condition. 

\skipline start_of_actions_after_adapt
\until end_of_actions_after_adapt

<HR>

\subsection action_after_distribution Actions after distribution

 The actions required after distribution are extremely similar to those
 required after adaptation because distribution deletes certain elements
(or replaces them by halo copies). The only difference is that the
 redundant pressure degrees of freedom do not have to be re-pinned.

\skipline start_of_actions_after_distribute
\until end_of_actions_after_distribute

<HR>

\subsection doc_solution The doc_solution() function

As usual, we add the processor number to the end of the filename for 
each output file to make sure that the different processors
don't over-write each other's output.

The trace file documents the time trace of the imposed influx
and the displacement of the node at the tip of the leaflet. It could be
written by any processor since all solid elements are retained
everywhere. We only write to the trace file from processor 0.


\skipline Write trace file
\until }


<HR>
<HR>

\section sources Source files for this tutorial

The source files for this tutorial can be found in

<CENTER>
<A HREF="../../../../demo_drivers/mpi/multi_domain/fsi_channel_with_leaflet">
demo_drivers/mpi/multi_domain/fsi_channel_with_leaflet
</A>
</CENTER>

Similar examples of modified driver codes for FSI problems for a
channel with a collapsible wall and an oscillating ring can be found
in
<CENTER>
<A HREF="../../../../demo_drivers/mpi/multi_domain/">
demo_drivers/mpi/multi_domain
</A>
</CENTER>


<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

