/**

\mainpage Existing (unstructured) meshes: Wrappers to third-party mesh generators

\c oomph-lib does not provide its own unstructured mesh generator
but has several mesh classes that generate unstructured meshes 
from the output of third-party unstructured mesh generators.

<B>Notes:</B>
-# The unstructured tet and triangle meshes listed below 
   can \b not be used with \c oomph-lib's mesh adaptation or 
   node-update procedures. A suitably fine mesh has to
   be generated offline by the third-party mesh generator.
   If required, node-updates (in response to changes in the domain
   boundaries) have to be performed manually. \n\n
-# For some element types, the mesh generation process is not
   particularly efficient (yet!). A suitable warning message
   is issued in such cases. \n\n
-# Since the third-party mesh generators tend to triangulate the domain
   with simplex elements, curvilinear boundaries are not
   resolved more accurately by using higher-order elements unless
   some post-processing is performed. \n\n
-# The meshes have not been tested as extensively as \c oomph-lib's
   structured meshes, described <A HREF="../../mesh_list/html/index.html">
   elsewhere.</A> \n\n
.  

<HR>
<HR>

\section mesh_list Mesh list

<CENTER>
<TABLE>
<TR>
<TD>
<B>Mesh</B>
</TD>
<TD>
<B>Representative Mesh plot</B>
</TD>
</TR>
<TR>
<TD>
<A HREF="../../mesh_from_triangle/html/index.html"><B>
TriangleMesh<ELEMENT>
</B></A> \n\n
- This class creates \c oomph-lib meshes based on the output from 
  <A HREF="http://www.cs.cmu.edu/~jrs">J.R.Shewchuk's</A>
  Delaunay mesh generator 
  <A HREF="http://www.cs.cmu.edu/~quake/triangle.html">
  <TT>Triangle</TT></A> 
- The mesh can be used with all \c FiniteElements that are derived from
  the geometric finite element \c TElement<2,NNODE_1D>.
. 
<B>Example driver codes:</B> \n
- The use of  <A HREF="http://www.cs.cmu.edu/~quake/triangle.html">
  <TT>Triangle</TT></A> and the \c TriangleMesh class
  are explained in a 
  <A HREF="../../mesh_from_triangle/html/index.html">
  separate tutorial.</A> \n\n
- In <A HREF="../../mesh_from_xfig/html/index.html">
  another tutorial</A> we demonstrate how the code 
  <A  HREF="../../../../demo_drivers/meshing/mesh_from_xfig_triangle/fig2poly.cc">fig2poly.cc</A> may be used to generate input files 
  for <A HREF="http://www.cs.cmu.edu/~quake/triangle.html">
  <TT>Triangle</TT></A> based on the output from the open-source 
  drawing program <A HREF="http://en.wikipedia.org/wiki/Xfig">xfig.</A> 
.
</TD>
<TD>
\image html oomph_mesh.gif " " 
\image latex oomph_mesh.eps " " width=0.75\textwidth
</TD>
</TR>
<TR>
<TD>
<A HREF="../../mesh_from_tetgen/html/index.html"><B>
TetgenMesh<ELEMENT>
</B></A> \n\n
- This class creates \c oomph-lib meshes based on the output from 
  <A HREF="http://www.wias-berlin.de/~si ">Hang Si's</A>
  open-source mesh generator 
  <A HREF="http://wias-berlin.de/software/tetgen//"> \c Tetgen </A>.
- The mesh can be used with all \c FiniteElements that are derived from
  the geometric finite element \c TElement<3,NNODE_1D>.
. 
<B>Example driver codes:</B> \n
- The use of   <A HREF="http://wias-berlin.de/software/tetgen//"> \c Tetgen </A>
  and the \c TetgenMesh class
  are explained in a <A HREF="../../mesh_from_tetgen/html/index.html">
  separate tutorial.</A> 
.
</TD>
<TD>
\image html tetgen_mesh.gif " " 
\image latex tetgen_mesh.eps " " width=0.75\textwidth
</TD>
</TR>
<TR>
<TD>
<A HREF="../../mesh_from_vmtk/html/index.html"><B>
Generating meshes from medical scans with VMTK
</B></A> \n\n
- We provide the option to generate tetgen-based meshes
  for physiological fluid-structure interaction problems, 
  using the  <a href="http://www.vmtk.org">Vascular 
  Modeling Toolkit (VMTK).</a>
. 
<B>Example driver codes and tutorials:</B> \n
- We provide a <a href="../../mesh_from_vmtk/html/index.html">separate 
  tutorial</a> that shows how to generate \c oomph-lib meshes
  from medical images.
- The methodology is used in the following driver codes:
  \n\n
  - <a href="../../../solid/vmtk_solid/html/index.html">The inflation of a 
    blood vessel.</a>
    \n\n
  - <a href="../../../navier_stokes/vmtk_fluid/html/index.html">
    Finite Reynolds number flow through a (rigid) iliac bifurcation.</a>
    \n\n
  - <a href="../../../interaction/vmtk_fsi/html/index.html">
    Finite Reynolds number flow through an elastic iliac bifurcation.</a>
    \n\n
  .
.
</TD>
<TD>
\image html iliac.gif " " 
\image latex iliac.eps " " width=0.75\textwidth
</TD>
</TR>
<TR>
<TD>
<A HREF="../../mesh_from_geompack/html/index.html"><B>
GeompackQuadMesh<ELEMENT>
</B></A> \n\n
- This class creates \c oomph-lib meshes based on the output from 
  Barry Joe's mesh generator 
  <A HREF="http://members.shaw.ca/bjoe/">\c Geompack++, </A>
  available as freeware at 
  <A HREF="http://members.shaw.ca/bjoe/">
  http://members.shaw.ca/bjoe/.</A>
- The mesh can be used with all \c FiniteElements that are derived from
  the geometric finite element \c QElement<2,2>.
. 
<B>Example driver codes:</B> \n
- The use of 
  <A HREF="http://members.shaw.ca/bjoe/">\c Geompack++ </A>
  and the \c GeompackQuadMesh class
  are explained in a <A HREF="../../mesh_from_geompack/html/index.html">
  separate tutorial.</A> 
.
</TD>
<TD>
\image html geompack_mesh.gif " " 
\image latex geompack_mesh.eps " " width=0.75\textwidth
</TD>
</TR>
</TABLE> 
</CENTER>    

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

