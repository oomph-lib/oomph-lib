/**

\mainpage Visualising oomph-lib's output files with Paraview

All of \c oomph-lib's existing elements implement the
\c GeneralisedElement::output(...) functions, allowing
the computed solution to be documented via a simple
call to the \c Mesh::output(...) function, e.g.

\code
// Open output file
ofstream output_file("output.dat")

// Call the mesh's output function which loops over the 
// element and calls theirs...
Problem::mesh_pt()->output(output_file);
\endcode

By default, the output is written in a format that is suitable
for displaying the data with <A HREF="http://www.tecplot.com">
tecplot,</A> a powerful and easy-to-use commercial plotting package
-- possibly a somewhat odd choice for a an open-source library. 

We also provide the capability to output data in 
a format that is suitable for display with 
<a href="http://www.paraview.org">paraview</a>,
an open-source 3D plotting package.
For elements for which the relevant output functions are
implemented (they are defined as broken virtual functions in the
\c FiniteElement base class) output files for all the elements
in a certain mesh (here the one pointed t by \c Bulk_mesh_pt) can 
be written as

\dontinclude rayleigh_instability.cc
\skipline Output solution to file in paraview format
\skipline some_file.open(
\until some_file.close();

where the unsigned integer \c npts  controls the number of plot points per
element (just as in the tecplot-based output functions). If \c npts
is set to 2, the solution is output at the elements' vertices. For 
larger values of \c npts the solution is sampled at greater number of (equally
spaced) plot points within the element -- this makes sense for 
higher-order elements, i.e. elements in which the finite-element solution is not
interpolated linearly between the vertex nodes. It is important to note that
when displaying such a solution in paraview's "Surface with Edges" mode, the 
"mesh" that is displayed does not represent the actual finite element mesh but is a
finer auxiliary mesh that is created merely to establish the connectivity
between the plot points.

Paraview makes it possible to animate sequences of plots from
time-dependent simulations. To correctly animate results from
temporally adaptive simulations (where the timestep varies)
paraview can operate on pvd files which provide a list of 
filenames and the associated time. These can be written automatically
from within \c oomph-lib, using the functions in the \c ParaviewHelper
namespace:

\dontinclude mesh.h
\skipline paraview_helper
\until }

Once the pvd file is opened, call \c
ParaviewHelper::write_pvd_header(...) to write the header information
required by paraview; then add the name of each output file 
and the value of the associated value of the continuous time, using \c
ParaviewHelper::write_pvd_information(...). When the simulation is
complete write the footer
information using \c
ParaviewHelper::write_pvd_footer(...), then close to the pvd file.

Currently, the paraview output functions are only implemented for a
relatively small number of elements but it is straightforward
to implement them for others.

The <a href="../../FAQ/html/index.html">FAQ</a> contain
an entry that discusses how to display \c oomph-lib's output
with <A HREF="http://www.gnuplot.info">gnuplot</A> and 
<a href="../../FAQ/html/index.html#tecplot">
how to adjust \c oomph-lib's output functions to different 
formats.</a>
  


Angelo 
Simone has written a python script that converts \c oomph-lib's
output to the vtu format that can be read by 
<a href="http://www.paraview.org">paraview</a>. 
This has since been improved
and extended significantly with input from Alexandre Raczynski and 
Jeremy van Chu. The conversion script
can currently deal with output from meshes that are composed
of 2D triangles and quad and 3D brick and tet elements. 

The \c oomph-lib distribution contains three scripts:

- <a href="../../../bin/oomph-convert.py"><code>bin/oomph-convert.py
  </code>:</a> The python conversion script itself.
  \n\n
- <a href="../../../bin/oomph-convert"><code>bin/oomph-convert
  </code>:</a> A shell script wrapper that allows the processing of
  multiple files.
  \n\n
- <a href="../../../bin/makePvd"><code>bin/makePvd</code>:</a>
  A shell script the creates the \c * \c .pvd files required by
  paraview to produce animations.
. 
 
<hr>
<hr>



\section py The oomph-convert.py script for single files

\subsection py_sample An example session


-# Add \c oomph-lib's bin directory to your path (in the example
   shown here, \c oomph-lib
   is installed in the directory \c /home/mheil/version185/oomph): 
   \n\n
   \code
   biowulf: 10:31:50$ PATH=$PATH:/home/mheil/version185/oomph/bin
   \endcode
   \n
-# Here is what's in the current directory at the moment: \c curved_pipe.dat is
   the \c oomph-lib output produced from 
   <a href="../../navier_stokes/curved_pipe/html/index.html">a
   simulation of steady flow through a curved pipe.</a>
   \n\n
   \code
   biowulf: 11:05:10$ ll
   total 824
   -rw-r--r--    1 mheil    users        2292 May 21 09:19 curved_pipe.dat
   \endcode
   \n
-# Run \c oomph-convert.py 
   \n\n
   \code
    biowulf: 11:16:18$  oomph-convert curved_pipe.dat

    -- Processing curved_pipe.dat
    * oomph-convert.py, ver. 20110531
    Convert from oomph-lib Tecplot format to VTK XML format.
    Dimension of the problem: 3
    Plot cells defined
    Field variables =  4
    Conversion started
    Coordinate defined
    Connectivities defined
    Offset defined
    Element types defined
    Field data defined
    Conversion done
    * Output file name: curved_pipe.vtu
    \endcode
   \n
-# We now have the corresponding \c * \c .vtu file
   \n\n
   \code
   biowulf: 11:32:08$ ll
   total 1024
   -rw-r--r--    1 mheil    users      329874 Jun 21 09:19 curved_pipe.dat
   -rw-r--r--    1 mheil    users      705294 Jun 21 11:16 curved_pipe.vtu
   \endcode
   \n
-# ...which we can visualise with paraview:
   \n\n
   \code
   biowulf: 11:34:08$ paraview --data=curved_pipe.vtu
   \endcode


<hr>
<hr>

If your output file is invalid or contains elements
that cannot currently be converted, you can use the \c -p option (followed by
2 or 3 to indicate the spatial dimension of the problem) to 
extract points only:

   \code
   biowulf: 11:16:13$ oomph-convert.py -p2 soln0.dat
   \endcode
   \n

The output is now a .vtp data file (Visualization Toolkit Polygonal)
which is also supported by Paraview. To display your .vtp data, use the
\image html glyph_button.gif " " 
\image latex glyph_button.eps " " width=0.03\textwidth
\c Glyph filter (displaying the points as crosses, say). Here 
is a representative plot in which 
<a href="../../poisson/fish_poisson/html/index.html"> the adaptive 
solution of a 2D Poisson equation in a fish-shaped domain</a>
is displayed with points. 


\image html paraview00.gif " " 
\image latex paraview00.eps " " width=\textwidth


\subsection display Display your data

Here are a few screenshots from a paraview session to get you
started. When paraview starts up, you have to select
the arrays of values you want to load and click on 
\image html apply_button.gif " " 
\image latex apply_button.eps " " width=0.07\textwidth
\c Apply:
 
\image html paraview01.gif " " 
\image latex paraview01.eps " " width=\textwidth


Select the array of values you want to display in the active window
(V1, V2, V3...). You can only display one at a time. It is applied on
the data set selected in the pipeline:

\image html paraview02.gif " " 
\image latex paraview02.eps " " width=\textwidth


Now choose the plot style of your data. \c Outline display a box containing
the data but not the data itself (it's not entirely clear to us
why you would want to do this, but...). \c Points and \c Wireframe 
best suited for 3D computations because they allow you to "see 
through" the data set. 
\c Surface and
\c Surface \c With \c Edges is best suited for 2D computations
because only the surface is displayed. Here is a view of the data 
in \c Wireframe mode:

\image html paraview03.gif " " 
\image latex paraview03.eps " " width=\textwidth


You can move the figure with buttons
\image html move_button.gif " " 
\image latex move_button.eps " " width=0.30\textwidth
in the toolbar or simply with the
mouse: "Left click + Move" to rotate, "Middle click + Move" to move, zoom in
with scroll up or "Right click + Move up" and zoom out with scroll down
or "Right click + Move down".

\image html paraview04.gif " " 
\image latex paraview04.eps " " width=\textwidth


You can also display the colour legend by clicking on 
\image html legend_button.gif " " 
\image latex legend_button.eps " " width=0.03\textwidth
, change the display colours (HSV,
RGB, Diverging...) by clicking on
\image html color_button.gif " " 
\image latex color_button.eps " " width=0.03\textwidth
, and rescale the colour scale to data by clicking on
\image html scale_button.gif " " 
\image latex scale_button.eps " " width=0.03\textwidth
. Modifying the colour scale is easy: add
points by clicking on the colour bar to change the distribution of colours
or use a logarithmic scale.

\image html paraview05.gif " " 
\image latex paraview05.eps " " width=\textwidth


You can split a window by clicking on 
\image html split_button.gif " " 
\image latex split_button.eps " " width=0.07\textwidth
. Clicking on the plot window will make
it active. You can switch on/off the display of a pipeline member by
clicking on the eye
\image html eye_button.gif " " 
\image latex eye_button.eps " " width=0.02\textwidth
on the left. You can display different values and states in different windows:

\image html paraview051.gif " " 
\image latex paraview051.eps " " width=\textwidth


<hr>
<hr>

\section py_mult The oomph-convert and makePvd scripts for multiple files and animations


\subsection py_sample_mult An example session for data from a serial computation

Here is a quick demonstration of \c oomph-convert and \c makePvd
scripts in action

-# Add \c oomph-lib's bin directory to your path (in the example shown
   here, \c oomph-lib
   is installed in the directory \c /home/mheil/version185/oomph): 
   \n\n
   \code
   biowulf: 10:31:50$ PATH=$PATH:/home/mheil/version185/oomph/bin
   \endcode
   \n
-# Here is what's in the current directory at the moment: \c soln?. \c dat are
   the \c oomph-lib output files that illustrate the progress
   of the mesh adaptation during
   <a href="../../poisson/fish_poisson/html/index.html">the adaptive 
   solution of a Poisson equation in a fish-shaped domain.</a>
   \n\n
   \code
   biowulf: 11:05:10$ ll
   total 824
   -rw-r--r--    1 mheil    users        2292 May 21 09:19 soln0.dat
   -rw-r--r--    1 mheil    users      176776 May 21 09:19 soln1.dat
   -rw-r--r--    1 mheil    users      278117 May 21 09:19 soln2.dat
   -rw-r--r--    1 mheil    users      367408 May 21 09:19 soln3.dat
   \endcode
   \n
-# Run \c oomph-convert on all files (the -z option adds zeroes to the
   numbers -- this is only required if the files are to combined
   into an animation by paraview)
   \n\n
   \code
   biowulf: 11:16:13$ oomph-convert -z soln*.dat


   -- Processing soln0.dat
   * oomph-convert.py, ver. 20110615
   Parse input file for Tecplot zones........done
   * 0 lines ignored    
   Write nodal coordinates...................done
   Write cell connectivity..................done
   Write cell offsets.......................done
   Write cell types.........................done
   Write field 01/01.........................done
   * Conversion done in 0 seconds
   * Output file name: soln00000.vtu 
 

   -- Processing soln1.dat
   * oomph-convert.py, ver. 20110615
   Parse input file for Tecplot zones........done
   * 0 lines ignored    
   Write nodal coordinates...................done
   Write cell connectivity..................done
   Write cell offsets.......................done
   Write cell types.........................done
   Write field 01/01.........................done
   * Conversion done in 0 seconds
   * Output file name: soln00001.vtu 


   [further output suppressed]
  


   \endcode
   \n
-# We now have the corresponding \c * \c .vtu files
   \n\n
   \code
   biowulf: 12:37:05$ ll
   total 2568
   -rw-r--r--    1 mheil    users        5979 Jun 21 12:37 soln00000.vtu
   -rw-r--r--    1 mheil    users      377490 Jun 21 12:37 soln00001.vtu
   -rw-r--r--    1 mheil    users      592990 Jun 21 12:37 soln00002.vtu
   -rw-r--r--    1 mheil    users      789325 Jun 21 12:37 soln00003.vtu
   -rw-r--r--    1 mheil    users        2292 Jun 21 09:19 soln0.dat
   -rw-r--r--    1 mheil    users      176776 Jun 21 09:19 soln1.dat
   -rw-r--r--    1 mheil    users      278117 Jun 21 09:19 soln2.dat
   -rw-r--r--    1 mheil    users      367408 Jun 21 09:19 soln3.dat
   \endcode
   \n
   These \c * \c .vtu files can be displayed individually
   as discussed above.
   \n\n
-# To produce an animation of the results with paraview, create a \c * \c .pvd 
   file using \c makePvd
   \n\n   
   \code
   biowulf: 12:40:56$ makePvd soln mysoln.pvd
   --> File mysoln.pvd created
   \endcode
   \n\n
-# ...and visualise it:
   \n\n
   \code
   biowulf: 12:42:08$ paraview --data=mysoln.pvd  
   \endcode
.


<hr>
<hr>


\subsection screenshots_mult Screenshots from the paraview session

Here's a screenshot from the paraview session: once the \c * \c .pvd
file is loaded you can customise the plot style as discussed
in the previous example, and then use the
\image html play_buttons.gif " " 
\image latex play_buttons.eps " " width=0.1\textwidth
\c Play/Stop/... buttons to 
animate the progress of the mesh adaptation. 

\image html paraview_animation_select.gif " " 
\image latex paraview_animation_select.eps " " width=\textwidth

<hr>
<hr>

\subsection py_sample_mult_par An example session for data from a parallel computation
 	 
\c oomph-lib typically outputs results from parallel (distributed)
computations on a processor-by-processor basis, resulting in filenames 
of the form
\code
 	 soln_proc0_0.dat 	 	 \ 
 	 soln_proc1_0.dat 	 	 | 
 	 	 : 	 	 	 | Data for timestep 0 
 	 soln_proc[NPROC-1]_0.dat 	 / 
 	 	 
 
 	 soln_proc0_1.dat 	 	 \ 
 	 soln_proc1_1.dat 	 	 | 
 	 	 :                    	 | Data for timestep 1 
 	 soln_proc[NPROC-1]_1.dat, 	 /

                 :

\endcode
where NPROC is the number of processors. An animation of such data 
obviously requires the output from different processors (but for the 
the same timestep) to be combined. Provided, the filenames have 
the pattern 
\code
	 [stem]proc[processor_number]_[timestep_number].dat
\endcode
(note the "proc" and "_", both of which are required), the pvd file 
can be generated by first processing the files with \c oomph-convert,
\code
        oomph-convert -z [stem]proc* 
\endcode
followed by 
\code
         makePvd [NPROC] [stem] myplot.pvd 
\endcode
So, for the files listed above, to produce a pvd file that 
contains data from a computation with four processors the commands
\code
biowulf: 12:40:56$ oomph-convert soln_proc*
\endcode
followed by
\code
biowulf: 12:40:59$ makePvd 4 soln_ soln.pvd 
\endcode
would create the file soln.pvd from which paraview 
can create an animation of the solution. 

<hr>
<hr>

\section filters Data analysis with filters

In order to analyse the data, we can apply filters. Some filters are
accessible directly via the navigation bar; a full list is available
in the \c Filters menu:

\image html paraview06.gif " " 
\image latex paraview06.eps " " width=\textwidth

Here are few examples of filters available:

-# \c Calculator: Evaluates a user-defined expression e.g \f$
   \frac{1}{2} V1^2 \f$ and creates a new data array, called here
   "Energy", containing the result of this expression:
   \n
\image html paraview07.gif " " 
\image latex paraview07.eps " " width=\textwidth
   \n
-# \c Contour: Extracts the points, curves or surfaces where a scalar
   field is equal to a user-defined value e.g \f$ V1 = -0.5 \f$:
   \n
\image html paraview08.gif " " 
\image latex paraview08.eps " " width=\textwidth
   \n
-# \c Clip: Intersects the geometry with a half space. (Warning:
   with some versions of Paraview, zooming on the clipped surface can
   cause the X server to crash.)
   \n
\image html paraview09.gif " " 
\image latex paraview09.eps " " width=\textwidth
\image html paraview091.gif " " 
\image latex paraview091.eps " " width=\textwidth
   \n
-# \c Slice: Intersects the geometry with a plane. (Warning:
   with some versions of Paraview, zooming on the clipped surface can
   cause the X server to crash.)
   \n
\image html paraview10.gif " " 
\image latex paraview10.eps " " width=\textwidth
\image html paraview11.gif " " 
\image latex paraview11.eps " " width=\textwidth
   \n
-# \c Threshold: Extracts cells that lie within a specified range of
values
   \n
\image html paraview12.gif " " 
\image latex paraview12.eps " " width=\textwidth


\section howto How to ...

\subsection outlines Select and extract elements

Click on the button:
\n
\image html selec_button.gif " " 
\image latex selec_button.eps " " width=0.03\textwidth
\c Select \c Cells \c On to select elements on the surface (2D selection)
\n
\image html selec2_button.gif " " 
\image latex selec2_button.eps " " width=0.03\textwidth
\c Select \c Points \c On to select points on the surface (2D selection)
\n
\image html selec3_button.gif " " 
\image latex selec3_button.eps " " width=0.03\textwidth
\c Select \c Cells \c Through to select elements through the region
selected (3D selection)
\n
\image html selec4_button.gif " " 
\image latex selec4_button.eps " " width=0.03\textwidth
\c Select \c Points \c Through to select points through the region
selected (3D selection)
\n
When your selection is  highlighted, go to \c Filters->Data \c
Analysis->Extract \c Selection and 
\image html apply_button.gif " " 
\image latex apply_button.eps " " width=0.07\textwidth
\c Apply the filter to extract the selected elements. You can now
modify or apply filters on the extracted data only.

Here is an example of extraction of the surface elements of the curved
pipe data:

\image html how_to_01.gif " " 
\image latex how_to_01.eps " " width=\textwidth


<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

