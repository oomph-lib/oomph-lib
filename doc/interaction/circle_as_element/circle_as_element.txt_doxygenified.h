/**

\mainpage Demo problem: Upgrading a GeomObject to a GeneralisedElement -- how to parametrise unknown domain boundaries

In many previous examples we demonstrated how to use \c GeomObjects
to parametrise curvilinear, moving domain boundaries. We
used this functionality extensively to perform simulations of 
problems in moving domains. The techniques illustrated in these examples 
are adequate for problems in which the motion of the boundary is 
prescribed; see e.g. the  
<A HREF="../../../navier_stokes/collapsible_channel/html/index.html">
simulation of fluid flow in a 2D channel in which a part of one
channel wall performs a prescribed oscillation</A>. 

  Here we shall demonstrate how to use \c GeomObjects to
parametrise curvilinear, moving domain boundaries in problem in which
the position of the boundary is unknown and has to be determined as
part of the overall solution. We will also address the following
questions:
- What is a \c GeomObject's geometric \c Data?
- What is a \c GeneralisedElement's external and internal \c Data?
.


<HR>
<HR>

\section geom_data Part 1: A closer look at the GeomObject class: Geometric Data

In an  <A HREF="../../../poisson/fish_poisson2/html/index.html">
earlier example</A>, we introduced the \c GeomObject as 
a base class for all objects that provide an analytical parametrisation of a
geometric object by specifying its  Eulerian coordinates (i.e.
the position vector to "material points" in the object,
\f$ {\bf r}\f$) as a function of some intrinsic (Lagrangian)
coordinates, \f$ {\bf \zeta} \f$, so that \f$ {\bf r} = {\bf r}({\bf
\zeta})\f$. Every specific \c GeomObject must implement the pure
virtual function 

\code
void GeomObject::position(const Vector<double>& zeta, Vector<double>& r);
\endcode

where the number of Lagrangian and Eulerian coordinates involved in 
the parametrisation (i.e. the sizes of the vectors \c zeta and \c r) are
stored as unsigned integers in the protected member data, 
\c GeomObject::NLagrangian and \c GeomObject::Ndim, respectively.


Most specific geometric objects store some additional parameters
as private data members. For instance, a \c GeomObject that 
represents a 2D circle, parametrised by the 1D Lagrangian
coordinate \f$ \zeta \f$ as
\f[ {\bf r}(\zeta)= \left( \begin{array}{c} X_c + R \cos(\zeta) \\
 Y_c + R \sin(\zeta)
\end{array} \right)
\f]
may be implemented as follows:

\dontinclude circle.h
\skipline start_of_circle 
\until end of SimpleCircle class

 
The shape and position of this \c GeomObject is defined by
the parameters \f$ X_c, Y_c \f$ and \f$ R \f$ which are stored
as double precision numbers in the object's private member data. 


In many applications (e.g. in free-boundary problems) it is
necessary to regard some (or all) of the parameters that 
specify the object's shape as (an) unknown(s) in the overall problem.
Within \c oomph-lib all unknowns are values in \c Data objects.
We refer to any \c Data objects whose values affect the
shape of a \c GeomObject as the \c GeomObject's geometric \c Data.
The \c GeomObject class defines interfaces for two functions that provide
access to such \c Data. The function

\code
virtual unsigned GeomObject::ngeom_data() const;
\endcode

should return the number of geometric \c Data objects that affect 
a \c GeomObject's shape and position.  The function

\code
virtual Data* GeomObject::geom_data_pt(const unsigned& j);
\endcode

should return a pointer to the \c GeomObject's \c j - th geometric 
\c Data object. Both functions are implemented as 
"broken" virtual functions in the \c GeomObject base class in order to
facilitate the creation of new \c GeomObjects. The functions should be 
overloaded in \c GeomObjects that actually contain geometric 
\c Data. If this is not done, the broken versions in the base 
class issue a suitable explanation before throwing an error.



Here is a re-implementation of the above \c GeomObject in a form that
allows the the parameters \f$ X_c, Y_c \f$ and \f$ R \f$ to be
regarded as unknowns, by storing them as values 
in a suitable \c Data object:



\skipline //==========
\until {
 
The object has two constructors which we will discuss separately.
The arguments to the first constructor specify the x- and
y-coordinates of the ring's centre (\f$ X_c\f$ and \f$ Y_c\f$),
and its radius, \f$ R \f$, as double-precision numbers. The constructor creates
storage for these values in the \c Geom_data_pt
vector which is stored in the object's protected member data.
First, the \c Geom_data_pt vector (empty by
default) is resized to provide storage for
a (pointer to a) single Data item. Next we create a Data object that
provides storage for three values and store the pointer to this \c Data
object in the first (and only) component of the \c Geom_data_pt vector.
Finally, we set the values to those specified by the arguments and
pin them to reflect the fact that, by default, the values
are constants rather than unknowns in the problem.
[\b Note: Clearly, there is some ambiguity as to how to distribute
the values of the geometric parameters among the geometric \c Data.
Here we chose to store the three values in a single \c Data object,
but we could also have stored each value in its own \c Data object.]

\skipline Constructor
\until }

With this constructor, the geometric \c Data that controls the shape
of the object is created internally -- the "user" only specifies
the values of the parameters. Their values are accessible via
the access functions 
 
\skipline /// Access function to x-coordinate of centre 
\until value_pt(2);}

but their internal representation as (pinned) values of a \c Data 
object remains hidden. Access to the geometric \c Data is, 
however, possible via the functions

\skipline How many items
\until Data* geom_data_pt(

which overload the broken versions in the \c GeomObject base class.
Both functions access the protected vector of pointers to the object's
geometric \c Data:

\skipline protected:
\until Geom_data_pt;

The second constructor is appropriate for cases in which 
the \c Data object that specifies the object's shape has already been
created elsewhere. In this case, we simply pass the pointer
to the Data object to the constructor and store it in the
first entry of the \c Geom_data_pt vector:
\dontinclude circle.h
\skipline   Alternative constructor
\until end of alternative constructor



The boolean flag \c Must_clean_up is stored as private member data in
the class. It is used in the destructor to decide
if the geometric \c Data should be deleted. If the \c Data was created
internally by the object (i.e. if the first constructor was used), the \c Data
must also be deleted by the object; if the \c Data was created externally
(i.e. if the second constructor was used), the \c Data must not be 
deleted as it may still be required in other parts of the code.

\skipline Destructor
\until } // end of destructor
 
The function \c GeomObject::position(...) 
simply extracts the relevant parameters (\f$ X_c, Y_c, R \f$) from the
\c Geom_data_pt vector and computes the position vector as a function
of the 1D Lagrangian coordinate \f$ \zeta \f$:

\skipline Position Vector at Lagrangian coord
\until end of position(...)
 
 
<HR>
<HR>

\section example Part 2: Upgrading a GeomObject to a GeneralisedElement

We will now consider a (toy!) problem in which the ability to
treat a parameter that specifies the shape and/or the position 
of a geometric object as an 
unknown (i.e. as a value in a \c Data object) is essential:
A circular ring of radius \f$ R \f$ is mounted on a linearly elastic
spring of stiffness \f$ k_{stiff} \f$. We wish to find the ring's
equilibrium position when it is subjected to a vertical load, 
\f$ f_{load} \f$. 

\image html circle_as_element_sketch.gif "Sketch of the example problem. " 
\image latex circle_as_element_sketch.eps "Sketch of the example problem. " width=0.75\textwidth


To solve this problem with \c oomph-lib [Yes, the solution is 
\f$ Y_c = f_{load}/k_{stiff}\f$ but let's pretend we 
don't know this...], we employ C++ inheritance
to "upgrade" the \c GeneralCircle object, discussed above, to a \c
GeneralisedElement in which the vertical position of the ring, \f$ Y_c
\f$, is one of the unknowns whose value is determined by the
equilibrium equation (in residual form),
\f[
 f_{load} - Y_c \ k_{stiff} = 0. \ \ \ \ \ \ \ \ \ \ \ (1)
\f]
In the present problem, the load on the ring is a control parameter,
but in different applications it may be an unknown whose value
has to be determined as part of the solution. (For instance, 
in fluid-structure interaction problems, the load
on an elastic structure depends on the unknowns in the adjacent fluid
elements.) Therefore, we also represent \f$ f_{load} \f$ by a \c Data
object. 

To solve the above problem with \c oomph-lib, the residual equation (1)
 must be implemented in the \c get_residuals(...) function of a suitable 
\c GeneralisedElement. Recall that the elemental residual
vector computed by a \c GeneralisedElement can depend on two types
of \c Data:
- <B>Internal <EM>Data</EM></B> stores values that are "internal" to the 
  element, i.e. values that the element is "in charge of". 
  Such \c Data is accessible via the pointer-based access function
  \c GeneralisedElement::internal_data_pt(...)
- <B>External <EM>Data</EM></B> stores values that are "external" to the 
  element. The values stored in such \c Data \e affect the element's 
  residual vector but its values are not determined \e by the element. 
  Such \c Data is accessible via the pointer-based access function
  \c GeneralisedElement::external_data_pt(...)
.
In the present context, it is most natural to regard the 
\c GeneralCircle's geometric \c Data as internal to the element and
the load \c Data as external. 


Here is the class definition for the \c ElasticallySupportedRingElement
which combines the \c GeneralisedElement and \c GeneralCircle classes by 
multiple inheritance. Its role as a \c GeomObject allows the object to be used
in the parametrisation of domain boundaries; its role as a 
\c GeneralisedElement allows the value of \f$ Y_c \f$ to be determined
as part of the solution. 
 

\dontinclude circle_as_generalised_element.h 
\skipline start_of_general_circle
\until public GeneralCircle


The arguments to the constructor specify the initial geometric
parameters, \f$ X_c, Y_c \f$ and \f$ R \f$. We pass them to the 
constructor of the \c GeneralCircle object and assign a default value for
the spring stiffness  \f$ k_{stiff}\f$ (stored as a private data
member of the class; see below):
\dontinclude circle_as_generalised_element.h 
\skipline Constructor
\until K_stiff(1.0)

Next, we add the \c GeneralCircle's geometric \c Data (created automatically
by the constructor of the \c GeneralCircle object; see above) 
into the element's storage for its internal \c Data. This ensures
that the geometric \c Data is included in the element's
equation numbering procedures. Within the 
context of the \c GeneralCircle, all geometric parameters were 
regarded as constants and their values were pinned. Here, the 
vertical position [stored in the second entry in the geometric \c
Data] is unknown, therefore we un-pin it.

\skipline {
\until unpin(1);

Once \c Data has been added to a \c GeneralisedElement's internal 
\c Data, it is deleted by the destructor of the
\c GeneralisedElement when the \c GeneralisedElement goes out of
scope. The \c Data must therefore not be deleted again
by the destructor of the \c GeneralCircle class, and we change the
cleanup responsibilities accordingly:

\dontinclude circle_as_generalised_element.h 
\skipline Change cleanup
\until }


Since the \c GeneralisedElement's destructor will delete
the internal \c Data, the destructor can remain empty:
\skipline Destructor
\until }
 
The \c Data whose one-and-only value represents the load
must be set by the "user", using the function \c set_load_pt(...)
As discussed above, we store the pointer to the load \c Data object in
the \c GeneralisedElement's external \c Data and record its index
within that storage scheme.

\dontinclude circle_as_generalised_element.h 
\skipline  Set pointer to Data 
\until end of set_load_pt(...)

[Note the sanity check which asserts that the load \c Data object 
only contains a single value; see \ref notes for a further
discussion of this aspect.]

 
The \c load() function provides access to the load specified by the 
load \c Data. It returns zero if no load was set -- a sensible
default.
 
\skipline acting on the 
\until end of load()


Next, we provide an access functions to the spring stiffness parameter,
\skipline spring stiffness
\until double

and functions that allow the vertical displacement of the ring to 
be pinned and un-pinned:
 
\skipline Pin
\until end of unpin_yc()
 
Finally, we implement the pure virtual functions 
\c GeneralisedElement::get_jacobian(...) and 
\c GeneralisedElement::get_residuals(...) which are required to
make the element functional. These functions must compute the
element's contributions to the \c Problem's global 
residual vector and Jacobian matrix. As usual, we implement them as wrappers
to a single function that computes the residual and (optionally)
the Jacobian matrix:

\skipline /// Compute element residual vector (wrapper)
\until } // end of get_jacobian(...)

The "real work" is done in the protected member function
\c get_residuals_generic(...), where we distinguish two cases.
-# The load is prescribed, i.e. the (single) value in the load \c Data
   object is pinned.
   In this case, the element's \c Data contains only one 
   unknown -- the vertical
   displacement \f$ Y_c \f$, stored as value 1 in the internal \c Data. 
   The element's residual vector has a single entry, given by
   \f[
   {\tt r}_0 = f_{load} - Y_c \ k_{stiff}
   \f]
   and the element's Jacobian matrix is a 1x1 matrix whose single
   entry is given by
   \f[  
   {\tt J}_{00} = \frac{\partial {\tt r}_0}{\partial Y_c} =  - k_{stiff} 
   \f]
-# If the load is unknown (i.e. an unknown in the overall problem)
   the element's \c Data contains two unknowns. Recall that an element 
   only makes a
   contribution to the residuals associated with unknowns that 
   it is "in charge of" -- the external \c Data is assumed to be 
   "determined" by another element. 
   Therefore, the elements residual vectors (i.e. its contribution
   to the \c Problem's global residual vector) is given by
   \f[
   \left( \begin{array}{c}
   {\tt r}_0 \\
   {\tt r}_1
   \end{array} \right)
   =
   \left( \begin{array}{c}
   f_{load} - Y_c \ k_{stiff}\\
   0
   \end{array} \right).
   \f]
   The element's 2x2 Jacobian matrix contains the derivatives
   of the residuals with respect to the element's unknowns, \f$ Y_c
   \f$ and \f$ f_{load}\f$, 
   \f[  
   \left( \begin{array}{cc}
   {\tt J}_{00} &    {\tt J}_{01} \\
   {\tt J}_{10} &    {\tt J}_{11} 
   \end{array} \right)
   =
   \left( \begin{array}{cc}
   \partial {\tt r}_0 / \partial Y_c&
   \partial {\tt r}_0 / \partial f_{load}   \\
   \partial {\tt r}_1 / \partial Y_c  &
   \partial {\tt r}_1 / \partial f_{load}  
   \end{array} \right)
   =
   \left( \begin{array}{cc}
   - k_{stiff} & 1 \\
   0 & 0
   \end{array} \right)
   \f]
.
Note that we have (correctly!) assumed that the (fully automatic)
local equation numbering procedure, implemented in 
\c GeneralisedElement::assign_local_eqn_numbers(), 
regards the internal degree of freedom as local unknown "0" and 
the external one (if it exists) as unknown "1". However, to be on 
the safe side, we determine the local equation 
numbers via the access functions \c internal_local_eqn(i,j)
and \c external_local_eqn(i,j), which return the local equation
numbers of the j-th value, stored in the element's i-th internal
(or external) \c Data object:

\skipline Compute element 
\until // end of get_residuals_generic(...)


<HR>
<HR>

\section driver Part 3: The driver code

Finally, we provide an example that shows our 
\c ElasticallySupportedRingElement in action. The animation below
shows the result of the computational simulation of the toy problem
described at the beginning of the previous section: An elastically
supported ring, subjected to a vertical load. The animation shows
the position of the ring for various values of the load.
 
\image html circle_as_element.gif "Vertical displacement of a Circle, mounted on an elastic support. " 
\image latex circle_as_element.eps "Vertical displacement of a Circle, mounted on an elastic support. " width=0.75\textwidth

 

<HR>

\subsection main The driver code

The driver builds the \c Problem object and solves the problem
for various load levels:

\dontinclude geom_object_element.cc
\skipline start_of_driver
\until end of driver


<HR>

\subsection class The problem class definition

Here is the problem class definition which requires little comment: 

\dontinclude geom_object_element.cc
\skipline start_of_problem
\until };

<HR>

\subsection constructor The problem constructor

In the problem constructor, we create the (single)
\c ElasticallySupportedRingElement and and assign a value to its 
"spring stiffness" parameter. 
\skipline =============start
\until k_stiff()

Next, we create the problem's mesh: We start by creating an (empty)
\c Mesh object and then add the (pointer to the) newly created
\c ElasticallySupportedRingElement (in its incarnation as
a \c GeneralisedElement) to it:
 
\skipline Build mesh
\until add_element

As discussed in the previous section, the load must be specified
as a \c Data object that stores a single value. In our example,
the load is prescribed (i.e. not an unknown in the problem)
therefore the value must be pinned:

\skipline Create the load
\until set_load_pt

Finally, we set up the equation numbering scheme, set 
the output directory and open a trace file in which we will document
the load/displacement characteristics. 


\skipline Setup equation
\until // end of constructor


<HR>

\subsection post The post-processing function

The post-processing routine is straightforward and is listed only to 
highlight the need to cast the pointer to the \c GeneralisedElement
(as returned by \c Mesh::element_pt(...)) to the specific
element used here -- obviously, a \c GeneralisedElement does not
have member functions that return the load on the ring or its vertical
displacement.


\skipline start_of_doc_solution
\until // end of doc_solution


<HR>
<HR>

\section notes Comments and Exercises
 
\subsection ex Exercises
-# In the current implementation of the \c ElasticallySupportedRingElement,
   the load \c Data must only store a single value: Generalise 
   this to the case where
   the \c Data object can contain an arbitrary number of values. 
   You could add an additional argument to the 
   \c set_load_pt(...) function to specify which value in the \c Data object
   represents the load. Think carefully in which other parts of the code 
   this index will be required. \n\n
.
  
  

<HR>
<HR>


\section sources Source files for this tutorial
- The source files for this tutorial are located in the directory:\n\n
<CENTER>
<A HREF="
../../../../
demo_drivers/interaction/free_boundary_poisson/
">
demo_drivers/interaction/free_boundary_poisson/
</A>
</CENTER>\n
- The driver code is: \n\n
<CENTER>
<A HREF="
../../../../
demo_drivers/interaction/free_boundary_poisson/geom_object_element.cc
">
demo_drivers/interaction/free_boundary_poisson/geom_object_element.cc
</A>
</CENTER>
.

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

