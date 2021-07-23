/**

\mainpage Creation of internal boundaries for Triangle based meshes

In this document we show how to add internal boundaries and create 
connections among them in 2D unstructured meshes with \c oomph-lib. 
The mesh generation procedure is performed by using 
<A HREF="http://www.cs.berkeley.edu/~jrs/">Jonathan Shewchuk's</A>
<A HREF="http://www.cs.cmu.edu/~quake/triangle.html"> Triangle</A>. 
This tutorial complements the discussion in 
<A HREF="../../../../doc/meshes/mesh_from_inline_triangle/html/index.html">
another tutorial</a>
where we illustrated how to build \c oomph-lib meshes using the 
\c TriangleMeshPolygon and \c TriangleMeshClosedCurve objects.

<HR>

 Additionally, we show via examples 
- how to connect \c TriangleMeshCurveSection objects
- how to add addition holes and extra regions, i.e. regions not defined 
  via \c TriangleMeshPolygon or \c TriangleMeshClosedCurve objects

<HR>
<HR>

\section overview Adding internal boundaries

\subsection open_curve The TriangleMeshOpenCurve object

An internal boundary is defined in \c oomph-lib by the object 
\c TriangleMeshOpenCurve. We can create and add as many internal 
boundaries as we need for defining our domain, see figure below.

\image html figure_extravagant.gif "The \c TriangleMeshOpenCurve object allow us to define domains with internal boundaries. The internal boundaries are shown in red. Notice the connections of the internal boundaries with other boundaries. The connections are shown with filled small circles. " 
\image latex figure_extravagant.eps "The \c TriangleMeshOpenCurve object allow us to define domains with internal boundaries. The internal boundaries are shown in red. Notice the connections of the internal boundaries with other boundaries. The connections are shown with filled small circles. " width=0.75\textwidth

The \c TriangleMeshOpenCurve object allows the creation of complex internal 
boundaries as a combination of \c TriangleMeshPolyLine and/or 
\c TriangleMeshCurviLine objects. In order to create a 
\c TriangleMeshOpenCurve it is necessary to create a 
\c Vector of \c TriangleMeshCurveSection objects that define
the \c TriangleMeshOpenCurve, see figure below.

\image html different_open_curves.gif "Different types of \c TriangleMeshOpenCurve, a) A \c TriangleMeshOpenCurve defined as a combination of two \c TriangleMeshPolyLines (in blue) and a \c TriangleMeshCurviLine (in red). b) A \c TriangleMeshOpenCurve can also be defined as just one \c TriangleMeshCurveSection, in this case a \c TriangleMeshPolyLine. c) A \c TriangleMeshOpenCurve defined as a combination of a \c TriangleMeshPolyLine and a \c TrianglesMeshCurviLine. " 
\image latex different_open_curves.eps "Different types of \c TriangleMeshOpenCurve, a) A \c TriangleMeshOpenCurve defined as a combination of two \c TriangleMeshPolyLines (in blue) and a \c TriangleMeshCurviLine (in red). b) A \c TriangleMeshOpenCurve can also be defined as just one \c TriangleMeshCurveSection, in this case a \c TriangleMeshPolyLine. c) A \c TriangleMeshOpenCurve defined as a combination of a \c TriangleMeshPolyLine and a \c TrianglesMeshCurviLine. " width=0.75\textwidth

<HR>

\section connections Creating connections among current boundaries in the domain and with internal boundaries

Any internal boundary can be connected (via its ends) to any other 
boundary currently in the domain.

It is important to note that there are four different types of
connections:
- Connection between a \c TriangleMeshPolyLine and a \c TriangleMeshPolyLine
- Connection between a \c TriangleMeshCurviLine and a \c TriangleMeshPolyLine
- Connection between a \c TriangleMeshPolyLine and a \c TriangleMeshCurviLine
- Connection between a \c TriangleMeshCurviLine and a \c TriangleMeshCurviLine

The next figure shows an sketch of the different types of connections.

\image html connections.gif "Different types of connections. a) Connection between a \c TriangleMeshPolyLine and a \c TriangleMeshPolyLine. b) Connection between a \c TriangleMeshCurviLine and a \c TriangleMeshPolyLine. c) Connection between a \c TriangleMeshPolyLine and a \c TriangleMeshCurviLine. d) Connection between a \c TriangleMeshCurviLine and a \c TriangleMeshCurviLine. " 
\image latex connections.eps "Different types of connections. a) Connection between a \c TriangleMeshPolyLine and a \c TriangleMeshPolyLine. b) Connection between a \c TriangleMeshCurviLine and a \c TriangleMeshPolyLine. c) Connection between a \c TriangleMeshPolyLine and a \c TriangleMeshCurviLine. d) Connection between a \c TriangleMeshCurviLine and a \c TriangleMeshCurviLine. " width=0.75\textwidth

Connections to \c TriangleMeshPolyLine can only be made by using existing vertices
of the destination \c TriangleMeshPolyLine. Conversely, connections to \c TriangleMeshCurviLine
can be made by specifying the position along the \c TriangleMeshCurviLine. The  
connection point is identified by the \c TriangleMeshCurviLine's intrinsic coordinates. The 
figure below shows the case when a connection with a \c TriangleMeshPolyLine is possible.

\image html figure_poly.gif "In order to connect a boundary to a \c TriangleMeshPolyLine it is necessary that the destination boundary has a vertex on the connection point. a) There is no vertex for connecting to the destination \c TriangleMeshPolyLine. b) There is a vertex in the destination \c TriangleMeshPolyLine to receive the initial end of the source boundary. " 
\image latex figure_poly.eps "In order to connect a boundary to a \c TriangleMeshPolyLine it is necessary that the destination boundary has a vertex on the connection point. a) There is no vertex for connecting to the destination \c TriangleMeshPolyLine. b) There is a vertex in the destination \c TriangleMeshPolyLine to receive the initial end of the source boundary. " width=0.75\textwidth

<HR>
<HR>

\section example_poisson Example: A domain with multiple internal boundaries and connections

We will show how to add internal boundaries and how to create connections among them by
using the same example as in 
<A HREF="../../../../doc/meshes/mesh_from_triangle/html/index.html">
another tutorial</A>.

<HR>

A sketch of the domain is shown on the next figure. The domain has seven boundaries, 
two of which define the outer boundary and the other five define internal boundaries.

\image html sketch_domain.gif "The sketch of the domain; note the different internal boundaries and their connection points. " 
\image latex sketch_domain.eps "The sketch of the domain; note the different internal boundaries and their connection points. " width=0.75\textwidth

<HR>

\subsection def_domain Definition of the domain/mesh

The following code is extracted from the Problem constructor and shows how to create the domain/mesh.

\subsubsection def_out Defining the outer boundaries

The first task is to define the outer boundary which we represent by a \c TriangleMeshClosedCurve
comprising two \c TriangleMeshPolyLines.

\dontinclude mesh_from_inline_triangle_internal_boundaries.cc
\skipline Begin - Outer boundary
\until End - Outer boundary

Note that we defined additional vertices along the two boundaries to allow the creation
of connections to the outer boundary, see figure below.

\image html outer_boundary_vertices.gif "Vertex 1 on boundary 0 and vertices 1 and 3 on boundary 1 are created to allow the subsequent connections with internal boundaries. " 
\image latex outer_boundary_vertices.eps "Vertex 1 on boundary 0 and vertices 1 and 3 on boundary 1 are created to allow the subsequent connections with internal boundaries. " width=0.75\textwidth

\subsubsection def_inopen Defining internal boundaries

We create a \c Vector of five \c TriangleMeshOpenCurves which will represent the 
internal boundaries.

\dontinclude mesh_from_inline_triangle_internal_boundaries.cc
\skipline // Internal open boundaries
\until inner_open_boundaries_pt(n_open_curves);

<HR>

<b> Simple internal boundaries without connections </b>

Our domain has two unconnected internal boundaries: a straight line, 
represented by a \c TriangleMeshPolyLine and a curve line represented
by a \c TriangleMeshCurviLine.

\dontinclude mesh_from_inline_triangle_internal_boundaries.cc
\skipline We start by creating the internal boundaries
\until new TriangleMeshOpenCurve(internal_curve_section2_pt);

<HR>

<b> Adding an internal boundary connected to an outer boundary </b>

We need to add three more internal boundaries which have connections with the outer 
boundary and with some internal boundaries; see the sketch of the domain.

In order to connect an internal boundary with any other boundary in the domain
it is necessary to specify the following:

- Is the initial or the final end of the boundary the one that we want to 
connect?

- Do we want to connect the boundary with a \c TriangleMeshPolyLine or with a 
\c TriangleMeshCurviLine?

- If we want to connect with a \c TriangleMeshPolyLine then we need to 
identify the corresponding vertex number to which we want to connect. If we
want to connect with a \c TriangleMeshCurviLine we need to identify the
corresponding intrinsic coordinate \f$ s \f$ on the curve for connection.

<HR>

Consider the creation of the internal boundary 4 shown on the sketch of the domain. 
This is an example of connecting a \c TriangleMeshCurviLine with a \c TriangleMeshPolyLine 
(the outer boundary). In this case both ends of the \c TriangleMeshCurviLine are 
connected so we need to use the interface for connecting the initial and the final 
ends of the boundary to a \c TriangleMeshPolyLine. There is only one thing left to 
identify, the vertex number of the destination boundary to which we want to connect. 
By looking at figure 1.5 we see that these values are vertex #3 for the initial end 
and vertex #1 for the final end. Once identified, we perform the connection as follow:

- Create the internal boundary as usual 

\dontinclude mesh_from_inline_triangle_internal_boundaries.cc
\skipline Open curve 3
\until boundary_id);

- Specify connections using the destination boundary and the vertex number

\dontinclude mesh_from_inline_triangle_internal_boundaries.cc
\skipline Do the connection with
\until vertex_to_connect_final);

- Finally, create the \c TriangleMeshOpenCurve object as usual

\skipline The open curve 
\until new TriangleMeshOpenCurve(internal_curve_section3_pt);

<HR>

<b> Adding an internal boundary connected with another internal boundary </b>

Boundary 5, a straight line connected to the outer boundary and to an internal boundary 
(the one that we created on the previous step).

The initial end of the boundary is connected to the outer boundary therefore we need to 
identify the vertex number to which it is connected; by looking at the sketch of the domain
we see that the destination vertex number is 1. The final end is connected with a 
\c TriangleMeshCurviLine therefore we need to specify the intrinsic coordinate in the
object where we want to connect (the intrinsic coordinate is the arc-length along the circle).

\dontinclude mesh_from_inline_triangle_internal_boundaries.cc
\skipline This boundary is connected to the outer
\until new TriangleMeshOpenCurve(internal_curve_section4_pt);

<b> Adding an internal boundary connected with another internal boundary 
(the last one from our example) </b>

The fifth internal boundary is connected to two internal boundaries. The initial end 
of the new internal boundary is connected to a \c TriangleMeshPolyLine and the final
end is connected to a \c TriangleMeshCurviLine. The vertex number for 
connection with the \c TriangleMeshPolyLine was established in the definition
of the previous internal boundary, vertex #1. The circles defining boundaries 4 and 6 
meet each other at \f$(x, y) = (\frac{3}{4}, \frac{\sqrt(7)}{4})\f$ which corresponds to
the intrinsic coordinate \f$ s = arctan(\frac{\sqrt{7}}{3}) \f$ on the 
\c TriangleMeshCurviLine object representing boundary 4.

\dontinclude mesh_from_inline_triangle_internal_boundaries.cc
\skipline Open curve 5
\until new TriangleMeshOpenCurve(internal_curve_section5_pt);

\subsubsection create_mesh Build the mesh

We build the mesh specific parameters using a \c TriangleMeshParameters object.

\dontinclude mesh_from_inline_triangle_internal_boundaries.cc
\skipline Create mesh
\until new TriangleMesh<ELEMENT>(triangle_mesh_parameters);

<HR>
<HR>

\subsection def_extra Defining additional holes and regions

It is possible to define holes in the domain using the methods discussed in 
<A HREF="../../../../doc/meshes/mesh_from_inline_triangle/html/index.html">
another tutorial</a>; where we explained how to create holes by using closed
curves. Holes could also be created by specifying additional parameters on the
\c TriangleMeshParameters object. The specification of regions or layers
follows the same principle.

\subsubsection def_extra_holes Defining additional holes

For defining holes it is only necessary to specify a point that lies inside a 
closed boundary (created by the connection of internal open boundaries). For 
example, the region bounded by boundaries 4, 5 and 6 can be turned into a hole
by specifying the coordinate \f$ ( 1.5, 0.75) \f$, see figure below.

\image html sketch_domain_hole.gif "The sketch of the domain with a hole created by the connection of internal lines. " 
\image latex sketch_domain_hole.eps "The sketch of the domain with a hole created by the connection of internal lines. " width=0.75\textwidth

\dontinclude mesh_from_inline_triangle_internal_boundaries_extra.cc
\skipline // Adding a hole on the domain
\until triangle_mesh_parameters.extra_holes_coordinates() = additional_hole;

\skipline // Pass the TriangleMeshParameters object to the TriangleMesh one
\until Problem::mesh_pt() =

\subsubsection def_extra_regions Defining regions

The definition of regions follows the same principle that defining holes, it means,
one only need to specify the coordinates of a point that lies inside a region.
If we would like to define the region showed on light grey on the figure below 
we use the method \c add_region_coordinates of the \c TriangleMeshParameters
object. You can specify any region id for the defined region.

\image html sketch_domain_hole_region.gif "The sketch of the domain with a hole and a defined region. " 
\image latex sketch_domain_hole_region.eps "The sketch of the domain with a hole and a defined region. " width=0.50\textwidth

\dontinclude mesh_from_inline_triangle_internal_boundaries_extra.cc
\skipline // Adding a region on the domain
\until Problem::mesh_pt()

By default the complete domain belongs to region zero. Therefore, one
could leave parts of the domain without a defined region or explicitly
define all the regions on the domain.
NOTE: Be sure to use the region ids when documenting the solution by
its specified regions.

<HR>
<HR>


\section ex Comments and exercises
-# Create two new internal boundaries without connections, one using a 
\c TriangleMeshPolyLine object and other using a \c TriangleMeshCurviLine object.
   \n\n
-# Define a new vertex on the new \c TriangleMeshPolyLine and create a connection
   to that vertex. 
   \n\n
-# Once one performs a connection to a \c TriangleMeshPolyLine or to a 
   \c TriangleMeshCurviLine there is an internal checking on the connection 
   vertices coordinates for the \c TriangleMeshPolyLine and on the connection 
   intrinsic coordinate for the \c TriangleMesheCurviline case, they must be 
   close enough in order to allow the connection. The tolerance values are
   \c tolerance_for_vertex_connection = \f$ 1.0e-14 \f$ for the 
   \c TriangleMeshPolyline case and \c tolerance_for_zeta_connection = \f$ 1.0e-14 \f$
   for the \c TriangleMeshCurviline case. One can explicitly define the allowed 
   tolerance by adding an extra argument when using the methods 
   \c connect_initial_vertex_to_polyline, \c connect_final_vertex_to_polyline 
   \c connect_initial_vertex_to_curviline or \c connect_final_vertex_to_curviline.
   As an exercise, use different values for the allowed connection tolerance.



<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="../../../../demo_drivers/meshing/mesh_from_triangle/">
demo_drivers/meshing/mesh_from_inline_triangle_internal_boundaries/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="../../../../demo_drivers/meshing/mesh_from_inline_triangle_internal_boundaries/mesh_from_inline_triangle_internal_boundaries.cc">
demo_drivers/meshing/mesh_from_inline_triangle/mesh_from_inline_triangle_internal_boundaries.cc
</A> 
</CENTER>
\n\n
- The additional driver code \n\n
<CENTER>
<A HREF="../../../../demo_drivers/meshing/mesh_from_inline_triangle_internal_boundaries/mesh_from_inline_triangle_internal_boundaries_extra.cc">
demo_drivers/meshing/mesh_from_inline_triangle_internal_boundaries/mesh_from_inline_triangle_internal_boundaries_extra.cc
</A>
</CENTER>
\n
shows how to define a hole by connecting internal boundaries.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

