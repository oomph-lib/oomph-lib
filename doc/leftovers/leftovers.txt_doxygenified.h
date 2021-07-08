/**

\mainpage Leftovers

\htmlonly

The following table provides a detailed side-by-side comparison 
between the procedural and object-oriented 
implementation of finite-element solution of problem P2.

\anchor table

 <CENTER> <B> Algorithm 1 </B> </CENTER>

 <B> Phase 1: Set up the problem </B> 

 <B> Phase 1a: Problem specification</B> 

<TABLE BORDER="1">
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- Choose the number of elements, \f$ N_E \f$, and the number of 
  nodes per element, \f$ n. \f$ This defines the total number of 
  nodes, \f$N\f$, and the local shape functions  
  \f$ \psi_j(s) \ \ (j=1,...,n)\f$  for all elements. 
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- Choose the  Element type, e.g. \c QPoissonElement<3,2>, a
  nine-node quadrilateral Poisson element. 
- Choose a suitable  Mesh, e.g. \c FishMesh.
</TD>
</TR>
</TABLE>


<B> Phase 1b: Mesh generation</B> 

<TABLE BORDER="1">
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- Discretise the domain by specifying the positions 
  \f$ X_{ij} \ (j=1,..,N; i=1,2)\f$,
  of the \f$ N \f$ nodes. 
- Generate the lookup scheme \f$ {\cal J}(j,e)\f$ that 
   establishes the relation
   between global and local node numbers. 
- Identify which nodes are located on which domain boundaries. 
</TD>
<TD WIDTH="50%"  VALIGN="TOP">
- The  Mesh is typically built in the  Problem constructor; 
  the  Problem stores a pointer to the  Mesh:
\code 
Problem::mesh_pt()=new FishMesh<QPoissonElement<3,2> >; 
\endcode
[Note how the  Element type is passed to the  Mesh constructor as
a template parameter. The \c FishMesh has a fixed number of
 Elements; for many other  Meshes, the number of  Elements can
be passed as an argument to the  Mesh constructor.]
- Once the  Mesh is built, we can find out how many distinct boundaries
  it has by calling
\code 
unsigned n_bound=Problem::mesh_pt()->nboundaries();
\endcode
- The \f$i\f$-th node on boundary \f$i_{bound}\f$ is 
accessible via
\code 
unsigned i;
unsigned i_bound;
Node* node_pt=Problem::mesh_pt()->boundary_node_pt(i_bound,i);
\endcode
</TD>
</TR>
</TABLE>


<B> Phase 1c: "Pin" nodes with essential (Dirichlet) 
     boundary conditions </B>  

<TABLE BORDER="1">
<TR>
<TD WIDTH="50%" VALIGN="TOP">
  - Loop over all global nodes j that are located on Dirichlet boundaries:
    - Assign a negative equation number to reflect the node's "pinned" 
      status:
      \f[ {\cal E}(j) = -1 \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
  - Pinning of  Nodes on (Dirichlet) boundaries is typically done
    in the  Problem constructor, using the  Mesh's access 
    functions to the boundary  Nodes.
  - As an example, the following code fragment pins the  Nodes
    on all boundaries
    \code 
    // How many boundaries are there in the mesh?
    unsigned num_bound = Problem::mesh_pt()->nboundary();

    // Loop over the boundaries
    for(unsigned ibound=0;ibound<num_bound;ibound++)
     {
      // How many nodes are there on this boundary?
      unsigned num_nod= Problem::mesh_pt()->nboundary_node(ibound);

      // Loop over nodes on the boundary
      for (unsigned inod=0;inod<num_nod;inod++)
       {
        // Pin the node
        Problem::mesh_pt()->boundary_node_pt(ibound,inod)->pin(0); 
       }
     }
    \endcode 
    [See \ref systems_of_PDEs for an explanation of the argument
     to the \c pin() function.]
</TD>
</TR>
</TABLE>


 <B> Phase 1d: Apply boundary conditions and provide initial guesses 
     for all unknowns</B> 

 <TABLE BORDER="1">
 <TR>
 <TD WIDTH="50%" VALIGN="TOP">
  - Loop over all global nodes \f$j=1,...,N\f$:
    - Provide an initial guess for the unknown nodal value (e.g.
      \f$ U_j=0 \f$) and set the nodal value to the (known)
      boundary value for nodes on the boundary, \f$ U_j = g(X_{ij}). \f$
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- The  Node constructor initialises all nodal values to zero. 
  These default values must be overwritten for nodes on the
  boundary. The nodal values in the interior of the domain 
  can be over-written, if a "better" initial guess for the solution 
  is available. Typically both steps are performed in the Problem
  constructor.
</TD>
</TR>
</TABLE>


 <B> Phase 1e: Set up the global equation numbering scheme</B> 

 <TABLE BORDER="1">
 <TR>
 <TD WIDTH="50%" VALIGN="TOP">
- Initialise the total number of unknowns, \f$M=0.\f$
- Loop over all global nodes \f$j=1,...,N\f$:
  - If global node j is not pinned
    (i.e. if \f$ {\cal E}(j) \ne -1 \f$) 
     - Increment the number of unknowns
        \f[ M=M+1 \f]
     - Assign the global equation number
        \f[ {\cal E}(j) = M \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
</TD>
- The equation numbering procedure is generic and fully 
  implemented. The equation numbering 
  is usually performed in the  Problem constructor by calling
  \code
  Problem::assign_eqn_numbers();
  \endcode 
</TR>
</TABLE>


 <B> Phase 1f: Set up the local equation numbering scheme</B> 

 <TABLE BORDER="1">
 <TR>
 <TD WIDTH="50%" VALIGN="TOP">
- Loop over the elements  \f$ e=1,...,N_e \f$
  <TABLE BORDER="0">
  <TR>
  <TD bgcolor="cornsilk">
  - Initialise the counter for the number of degrees of 
    freedom in this element, \f$ j_{dof}=0 \f$.
  - Loop over the element's local nodes \f$ j_{local}=1,...,n\f$
    - Determine the global node number 
      \f$ j_{global} = {\cal J}(j_{local},e) \f$
    - Determine the global equation number \f${\cal E}(j_{global})\f$
    - If  \f${\cal E}(j_{global}) \ne -1\f$:
      - Increment the number of degrees of 
        freedom in this element \f$ j_{dof}=j_{dof}+1\f$
      - Add the entry to the lookup scheme that relates local 
        and global equation numbers
        \f$ \widehat{\cal E}(j_{dof},e) = {\cal E}(j_{global})\f$
      - Store the local equation number associated with the
        current local node: \f$ {\cal L}(j_{local},e)=j_{dof}.\f$
      .
    - Else:
      - Set the local equation number associated with the
        current local node to -1 to indicate that it is pinned:
        \f$ {\cal L}(j_{local},e)=-1\f$
    .
  - Assign the number of degrees of freedom in the element
    \f$ {\cal N}_{dof}(e)=j_{dof} \f$
  </TD>
  </TR>
  </TABLE>
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- The local  equation numbering is performed automatically by
  \c Problem::assign_eqn_numbers() which executes the
   Elements' member function 
   \c GeneralisedElement::assign_local_eqn_numbers().
  \n \n
  This (virtual) function is defined in \c GeneralisedElement 
  and must be implemented for every specific  Element. For instance,
  <TABLE BORDER="0">
  <TR>
  <TD bgcolor="cornsilk">
  \c
  QPoissonElement<3,2>::assign_local_eqn_numbers()
  </TD>
  </TR>
  </TABLE>
  implements the equation numbering procedure, highlighted in the
  left column, for nine-node quadrilateral Poisson elements. 
  [See \ref systems_of_PDEs for a more detailed discussion of this
   aspect.]
</TD>
</TR>
</TABLE>


 <B> End of setup:</B> 

<TABLE BORDER="1">
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- The setup phase is now complete and we can determine the 
  (current) FE representation of all quantities in element \f$ e\f$:\n\n
  - The global coordinates \f$ x_i \f$ are given by
    \f[
    x_i = \sum_{j=1}^{n}  X_{i{\cal J}(j,e)} \ \psi_j(s)
    \f]
  - Their derivatives with respect to the local coordinates are
    \f[
     \widehat{\cal J}_{ik} =  \frac{\partial x_i}{\partial s_j} =
     \sum_{k=1}^n X_{i{\cal J}(k,e)}
     \ \frac{\partial \psi_k(s)}{\partial s_j}
    \f]
  - The Jacobian of the mapping between the local and global coordinates
    is given by
    \f[
      \widehat{\cal J} =  \det \widehat{\cal J}_{ik}
    \f]
  - The function \f$ u\f$ is obtained by interpolation between the
    the nodal values:
    \f[
    u = \sum_{j=1}^{n}  U_{{\cal J}(j,e)} \ \psi_j(s_1,s_2) 
    \f]
  - Its derivatives with respect to the global coordinates are
    \f[ u = \frac{\partial u}{\partial x_i} =
    {\widehat J}^{-1}_{ij}
    \left(\sum_{k=1}^n U_{{\cal J}(k,e)}
    \ \frac{\partial \psi_k(s_1,s_2)}{\partial s_j}
    \right) 
    \f]
  - Finally, the derivatives of the local shape functions, 
    \f$ \psi_i \f$, with respect to the global coordinates 
    \f$ x_j \f$ are given by
    \f[  \frac{\partial \psi_i}{\partial x_j} = {\widehat J}^{-1}_{ik} \ 
         \frac{\partial \psi_i(s_1,s_2)}{\partial s_k}
          \mbox{\ \ \ for $j=1,..,n$} 
    \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- The FE representation of the unknown function, the global coordinates, 
  and their derivatives are available from  Element
  member functions. Here are a few representative member functions for our 
  nine-node quadrilateral Poisson element:\n\n
  - The FE representation of the global (Eulerian) coordinate \f$x_i\f$ at 
    the local coordinate
    \f$ (s_1,s_2)\f$:
    \code
    Vector<double> s(2);
    double x=QPoissonElement<3,2>::interpolated_x(s,0);
    double y=QPoissonElement<3,2>::interpolated_x(s,1);
    \endcode
    [\b Note: In the C++ implementation, all indices start at 0, so the
    global coordinates are \f$(x_0,x_1)\f$.]
  - The Jacobian of the mapping between the global (Eulerian) and
    local coordinates, \f$x_i(s_j)\f$, at the local coordinate
    \f$ (s_0,s_`)\f$:
    \code
    Vector<double> s(2);
    double jacobian=QPoissonElement<3,2>::J_eulerian(s);
    \endcode
  - The FE representation of the solution at the local coordinate
    \f$ (s_1,s_2)\f$:
    \code
    Vector<double> s(2);
    double u=QPoissonElement<3,2>::interpolated_u(s);
    \endcode
  - The shape functions and their first derivatives with respect to
    the local coordinates at the local coordinate
   \f$ (s_0,s_1)\f$:
   \code
    Shape psi(9);     // We have 9 shape functions 
    DShape dpsi(9,2); // The 9 shape functions 
                      // have 2 derivatives each
                      // (w.r.t to s_0 and s_1)
    Vector<double> s(2);
    QPoissonElement<3,2>::dshape_local(s,psi,dpsi);
    \endcode
  - The shape functions and their first derivatives with respect to
    the global (Eulerian) coordinates at the local coordinate
    \f$ (s_0,s_1)\f$:
    \code
    Shape psi(9);       // We have 9 shape functions 
    DShape dpsidx(9,2); // The 9 shape functions 
                        // have 2 derivatives each 
                        // (w.r.t to x_0 and x_1)
    Vector<double> s(2);
    QPoissonElement<3,2>::dshape_eulerian(s,psi,dpsidx);
   \endcode
 </TD>
 </TR>
 </TABLE>

 <B> Phase 2: Solution  </B> 
 
 <TABLE BORDER="1" >
 <TR>
 <TD WIDTH="50%" VALIGN="TOP">
 - The solution process has three phases:
   - Phase 2a: Setup the linear system
   - Phase 2b: Solve the linear system
   - Phase 2c: Update the solution.

   These phases are detailed below.
 </TD>
 <TD WIDTH="50%" VALIGN="TOP">
 - The three phases listed on the left, 
   represent the central steps in the \ref Newton iteration.
   \c Problem::newton_solve() provides a "black box" implementation 
   of Newton's method and represents \c oomph-lib's default nonlinear solver. 
   Once a Problem has been set up 
   as discussed above, the solution only requires the specification
   of a linear solver for the Newton iteration.
   It's as easy as this:
   \code 
   // Build/set a linear solver
   Problem::linear_solver_pt() = new SomeSolver();

   // Solve the Problem
   Problem::newton_solve();

   // Done -- go away and document your solution...
         :
         :
   \endcode
   The Newton solver uses a default tolerance of \f$10^{-8} \f$.
   This default value (stored as a protected data member of the
   Problem class) can be over-written in any of the Problem's
   member functions:
   \code
   Problem::Newton_solver_tolerance=1.0e-14;
   \endcode
   Details of the implementation are given below:
 </TD>
 </TR>
 </TABLE>


 <B> Phase 2a: Setup the linear system </B> 

 <TABLE BORDER="1" >
 <TR>
 <TD WIDTH="50%" VALIGN="TOP">

- Initialise the global residual vector, \f$ r_j=0 \f$ for \f$j=1,...,M \f$
  and the global Jacobian matrix \f$ J_{jk}=0 \f$ for \f$j,k=1,...,M \f$
- Loop over the elements  \f$ e=1,...,N_E\f$
  \n  \n   
  <TABLE BORDER="0">
  <TR>
  <TD bgcolor="cornsilk">
  <CENTER> <B> Compute the element's residual vector and 
  Jacobian matrix </B>  </CENTER>
  \n
  - Determine the number of degrees of freedom in this element,
    \f$ n_{dof}={\cal N}_{dof}(e) \f$.
  - Initialise the element residual vector, 
     \f$ r_{j}^{(e)}=0 \f$ for \f$j=1,...,n_{dof} \f$
     and the element Jacobian matrix \f$ J_{jk}^{(e)}=0 \f$ 
     for \f$j,k=1,...,n_{dof} \f$
  - Loop over the Gauss points \f$ i_{int}=1,...,N_{int}\f$
    - Determine the local coordinate of the integration point
      \f$ s_j = S_{j\, i_{int}} \mbox{\ \ \ for $j=1,2$.} \f$
    - Compute
      - the global coordinates  \f$ x_i = x_i(s_j) \f$
      - the source function    \f$ f = f(x_i) =f(x_i(s_j)) \f$
      - the derivative \f$ \partial u/\partial x_i \f$ and 
      - the shape functions  \f$ \psi_j \f$ and their derivatives
        \f$ \partial \psi_j/\partial x_k \f$
      - the Jacobian of the mapping between local and global
        coordinates, \f$ \widehat{\cal J} = 
        \det \partial x_i/\partial s_j\f$
      .
    - Loop over the local nodes \f$ j_{local}=1,...,n \f$
      - Determine the global node number 
        \f$ j_{global}={\cal J}(j_{local},e)\f$
      - Determine the global equation number 
        \f$ {\cal E}( j_{global}) \f$
      - If  \f$ {\cal E}( j_{global}) \ne -1 \f$
        - Determine the local equation number from the element's 
          lookup scheme \f$ i_{dof}= {\cal L}(j_{local},e)\f$.
        -  Add the contribution to the element's residual vector
             \f[ r_{i_{dof}}^{(e)} = r_{i_{dof}}^{(e)} +
             \left( \sum_{l=1}^2 \frac{\partial 
             u}{\partial x_l} \ 
             \frac{\partial \psi_{j_{local}}}{\partial x_l} + 
             f \  \psi_{j_{local}} \right) 
             \widehat{\cal J} \ W_{i_{int}} \f]
          - Loop over the local nodes \f$ k_{local}=1,...,n \f$
            - Determine the global node number 
              \f$ k_{global}={\cal J}(k_{local},e)\f$
            - Determine the global equation number 
              \f$ {\cal E}( k_{global}) \f$
            - If  \f$ {\cal E}( k_{global}) \ne -1 \f$
              - Determine the local equation number from the element's 
                lookup scheme \f$ j_{dof}= {\cal L}(j_{local},e)\f$.
              - Add the contribution to the element's Jacobian matrix
                \f[ J_{i_{dof}j_{dof}}^{(e)} = J_{i_{dof}j_{dof}}^{(e)} + 
                \left( \sum_{l=1}^2  
                \frac{\partial \psi_{k_{local}}}{\partial x_l} \ 
                \frac{\partial \psi_{j_{local}}}{\partial x_l}
                \right)
                \widehat{\cal J} \ W_{i_{int}}  \f]
  </TD>
  </TR>
  </TABLE>
  \n  \n  \n \n \n
  <CENTER> <B> Add the element's contribution to the global residual
  vector and Jacobian matrix </B>  </CENTER>
  \n
  - Loop over the local degrees of freedom \f$ i_{dof}=1,...,n_{dof} \f$
    - Add the element's contribution to the global residual vector
      \f[ r_{\widehat{\cal E}(i_{dof},e)} = 
      r_{\widehat{\cal E}(i_{dof},e)} + r_{i_{dof}}^{(e)} 
      \f]
      - Loop over the local degrees of freedom \f$ j_{dof}=1,...,n_{dof}\f$
        - Add the element's contribution to the global Jacobian matrix
          \f[ J_{\widehat{\cal E}(i_{dof},e)\widehat{\cal E}(j_{dof},e)}=  
          J_{\widehat{\cal E}(i_{dof},e)\widehat{\cal E}(j_{dof},e)}+ 
          J_{i_{dof} j_{dof}}^{(e)} \f]
        .
      .
    .
  .
  </TD>
  <TD WIDTH="50%" VALIGN="TOP">
   - The assembly process is completely generic:
     - Find out how many degrees of freedom we have from
       \code
       unsigned ndof_global=Problem::ndofs();
       \endcode
       This allows us to create and initialise the global Jacobian matrix
       and the global residual vector.
       \code
       Vector<double> residual(ndof_global);
       DenseDoubleMatrix jacobian(ndof_global,ndof_global);
       \endcode
       [\b Note: This is just an illustration! Do \b not store
       the global Jacobian matrix as a \c DenseMatrix --
       the \c oomph-lib implementation of a dense matrix!]. 
     - Find out how many elements there are (in the  Mesh)
       from
       \code
       Problem::mesh_pt()->nelement();
       \endcode
     - Now we can loop over all  Elements.  Elements are accessed
       via the  Mesh's access function \c Mesh::element_pt(), so
       we obtain a pointer to the \f$ e \f$-th element in the
        Problem's  Mesh, by calling
       \code
       FiniteElement* element_pt=Problem::mesh_pt()->element_pt(e);
       \endcode
     - How many degrees of freedom does this  Element have?
       \code
       unsigned ndof=element_pt->ndof();
       \endcode
     - We can now allocate the  Element's Jacobian matrix, 
       \c el_jacobian, and its residual vector, \c el_residual.
       \code
       Vector<double> el_residual(ndof);
       DenseDoubleMatrix el_jacobian(ndof);
       \endcode
     - Get their entries from the  Element by calling
       \code
       element_pt->get_jacobian(el_jacobian,el_residual);
       \endcode
       Note that the implementation of this function
       depends on the (mathematical) problem to be solved but its interface
       is completely generic. The function is therefore defined as
       a virtual function in the abstract base class \c GeneralisedElement,
       and needs to be implemented in a derived class. 
       For instance, 
       <TABLE BORDER="0">
       <TR>
       <TD bgcolor="cornsilk">
       \n
       \c QPoissonElement<NNODE_1D,DIM>::get_jacobian(el_jacobian,el_residual);
       \n \n
       </TD>
       </TR>
       </TABLE>
       implements the computation of the element residual and 
       Jacobian matrix
       for the \c DIM -dimensional Poisson equation, discretised with 
       isoparametric elements from the \c QElement family, so that
       - \c DIM=1 corresponds to 1D line elements
       - \c DIM=2 corresponds to 2D quadrilateral elements
       - \c DIM=3 corresponds to 3D brick  elements.
       .
       The second template parameter, \c NNODE_1D indicates the
       number of nodes along the element edges so that
       - \c NNODE_1D=2 corresponds to elements with linear 
       (or {bi/tri}linear) shape functions
       - \c NNODE_1D=3 corresponds to elements with quadratic 
       (or {bi/tri}quadratic) shape functions
       - etc.
       .
     \n \n 

     - Copying of the  Element's Jacobian matrix and residual vector 
       into their global counterparts is again completely generic
       and involves calls to the  Element's \c eqn_number() function
       which returns the global equation number corresponding
       to a local equation number. For instance, the assembly
       loop for the global residual vector \c residual can be
       implemented as follows
       \code
       // Loop over the number of dofs in the element
       for (unsigned i_local=0;i_local<ndof;i_local++)
        {
         // Get the corresponding global equation number
         unsigned i_global=element_pt->eqn_number(i_local);

         // Add the entry to the global residual vector
         residual[i_global]+=el_residual[i_local];
        }
        \endcode
   </TD>
   </TR>
   </TABLE>


 <B> Phase 2b: Solve the linear system </B> 

 <TABLE BORDER="1">
 <TR>
 <TD WIDTH="50%" VALIGN="TOP">
- Solve the \f$ M \times M\f$ linear system
  \f[ \sum_{j=1}^{M}  J_{jk} \ y_k = - r_j
        \mbox{\ \  \ for $j=1,...,M.$} \f]
  for \f$ y_k \ \ (k=1,...,M).\f$ 

 </TD>
 <TD WIDTH="50%" VALIGN="TOP">
 - The solution of the fully-assembled linear system
   can now be obtained, using the linear solver specified.
   \n \n
   [\b Note: Frontal solvers such as HSL_MA42 and HSL_MP42
   interlace the assembly and LU decomposition of the linear system;
   in such cases there is no clear distinction between phases
   2a and 2b.]
 </TD>
 </TR>
 </TABLE>


 <B> Phase 2c: Update the initial guess </B> 

 <TABLE BORDER="1">
 <TR>
 <TD WIDTH="50%" VALIGN="TOP">

-  Loop over all global nodes \f$ j\f$:
        - Determine the equation number: \f$ {\cal E}(j) \f$
        - If \f$ {\cal E}(j) \ne -1 \f$:
   \f[ U_{{\cal E}(j)} = U_{{\cal E}(j)}  +  y_{{\cal E}(j)} \f] 

 </TD>
 <TD WIDTH="50%" VALIGN="TOP">
 - The Newton solver \c Problem::newton_solve(...) updates the nodal 
   values directly via pointer-based access to the unknowns.
 </TD>
 </TR>
 </TABLE>

 <B> Phase 3: Postprocessing (document the solution)</B> 
 <TABLE BORDER="1">
 <TR>
 <TD WIDTH="50%" VALIGN="TOP">
- The finite-element solution in element \f$e\f$ is 
  given by 
  \f[ u = \sum_{j=1}^{n}  U_{{\cal J}(j,e)} \ \psi_j(s) \f]
 </TD>
 <TD WIDTH="50%" VALIGN="TOP">
 - All finite elements and meshes have an \c output() function,
   so the postprocessing merely requires a call to the  Mesh's
   \c output function:
   \code
   // Write the solution in all elements to file
   ostream output_file;
   Problem::mesh_pt()->output(output_file);
   \endcode
 </TD>
 </TR>
 </TABLE>

<HR>
<HR>

\section extensions Extensions to more complicated problems
We will now briefly discuss further extensions of
the method. ***expand***

<HR>

\subsection non_dirichlet_bc Non-Dirichlet boundary conditions
The examples considered so far have all been problems with
Dirichlet boundary conditions, i.e. problems in which the
value of the unknown function was prescribed on the entire
domain boundary. To illustrate how non-Dirichlet
boundary conditions are incorporated into a finite element discretisation,
we consider the following variant of problem P2:


<CENTER>
<TABLE> 
<TR>
<TD>
<CENTER>
<B> Problem P4 </B>
</CENTER>
Solve
\f[ {\cal R}(x_i; u(x_i)) = \sum_{k=1}^2
\frac{\partial^2 u(x_i)}{\partial x_k^2} -
f(x_i) =0  \mbox{ \ \ \ \ for $x_i\in D$} \f]
subject to 
\f[ u|_{\partial D_{Dirichlet}} = g(\zeta)
\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (1)\f]
and 
\f[ \left. \frac{\partial u}{\partial n}\right|_{\partial
D_{Neuman}} = h(\zeta),
\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \   (2) \f]
where \f$ f(x_i), g(\zeta)\f$ and \f$ h(\zeta)\f$ are given.
The domain boundary \f$ \partial D\f$ is parametrised by
the arclength \f$ \zeta \f$, so that the position vector
to the boundary is given by
\f[
{\bf r}^{\partial D}(\zeta) = \left( x_1^{\partial D}(\zeta),
                                     x_2^{\partial D}(\zeta) \right)^T
\f]
The boundary consists of two parts so that
\f$\partial D = \partial D_{Dirichlet} \cup \partial D_{Neuman} \f$
and \f$ \partial D_{Dirichlet} \cap \partial D_{Neuman} = 0. \f$ \n
\f$ \partial/\partial n \f$ denotes the derivative in
direction of the outer unit normal to the domain.
\image html dirichlet_and_neumann_bc.gif "Sketch of the domain and its " 
\image latex dirichlet_and_neumann_bc.eps "Sketch of the domain and its " width=0.55\textwidth

</TD> 
</TR>
</TABLE>
</CENTER>

In the context of the physical problem of steady heat conduction,
where \f$ u(x_i) \f$ represents the temperature, 
the Dirichlet boundary condition \f$ u|_{\partial D_{Dirichlet}} =
g\f$ prescribes the temperature along a certain part of the
boundary, while the Neumann condition 
\f$\partial u/\partial n|_{\partial D_{Neuman}} = h\f$ prescribes
the heat flux into the body along the remainder of its
surface.

As before, we define a weak solution of this problem 
as a(ny) function \f$ u_w(x_i) \f$ that satisfies the
Dirichlet boundary conditions  (1), 
\f[ u_w|_{\partial D_{Dirichlet}} = g(\zeta),\f]
and for which the weighted residual
\f[ r = \int_{D}  {\cal R}(x_i; u_w(x_i)) \ \phi^{(test)}(x_i) 
\ \mbox{d}x_1 \mbox{d}x_2 \ \ \ \ \ \ \ \ (3) \f]
vanishes for \e any  "test function" 
\f$ \phi^{(test)}(x_i) \f$ that satisfies homogeneous
boundary conditions on the \f$\partial D_{Dirichlet}\f$,
\f[ \phi^{(test)}|_{\partial D_{Dirichlet}} = 0. \f]

The definition is consistent with 
P2<SUB>weak</SUB>: the two are identical if
\f$ \partial D = \partial D_{Dirichlet} \f$,
or, equivalently, if \f$ \partial D_{Neumann}=0 \f$. 
The role of the test function in the weak
formulation of our
problem is to test whether the candidate solution violates any of the
conditions imposed on it. There are two types of test: 
-# Does the candidate solution satisfy the PDE in the interior of 
   the domain? 
-# Does the candidate solution satisfy the boundary conditions?
.
We assume the Dirichlet conditions to be satisfied a priori, so 
it is not necessary to perform any tests on  
\f$ \partial D_{Dirichlet} \f$. The easiest way to exclude
\f$ \partial D_{Dirichlet} \f$ from the "test" is to set
\f$ \phi^{(test)}|_{\partial D_{Dirichlet}} = 0 \f$ as we have
done. If \f$ \partial D_{Neumann} \ne 0 \f$, we need to 
test whether our candidate solution satisfies 
the Neumann boundary condition (2). 
Certain types of Neumann boundary conditions can be
naturally incorporated into the weighted residual.

Performing the integration by parts in our present problem
converts  (3) into
\f[ \int_D \left[ \sum_{k=1}^2 \frac{\partial u(x_i)}{\partial x_k} \ 
 \frac{\partial \phi^{(test)}(x_i)}{\partial x_k}  +
 f(x_i) \  \phi^{(test)}(x_i) \right]  \mbox{d}x_1 \mbox{d}x_2 = \f]
\f[
= \int_{\partial D_{Dirichlet}}  \phi^{(test)} \
 \frac{\partial u}{\partial n}  \ \mbox{d}\zeta \ + \
 \int_{\partial D_{Neumann}}  
 \ \phi^{(test)} \ \frac{\partial u}{\partial n} \ \mbox{d}\zeta, \ \ \ \ \ 
(4) \f]
where \f$\zeta\f$ is the arclength along the domain boundary 
\f$ \partial D\f$. We have split the integral along \f$ \partial D\f$
into two separate integrals over 
\f$ \partial D_{Dirichlet}\f$ and \f$ \partial D_{Neumann}\f$.
Since \f$ \phi^{(test)}|_{\partial D_{Dirichlet}} = 0, \f$
the integral over \f$\partial D_{Dirichlet}\f$ vanishes
automatically. The second integral has to be retained because
\f$ \phi^{(test)} \f$ can take any value on  
\f$\partial D_{Neumann}\f$. However, \f$ \partial u/\partial n \f$
in its integrand is specified by the Neumann boundary 
condition (2). We can, therefore, incorporate
(2) into the weighted residual by writing
\f[ \int_D \left[ \sum_{k=1}^2 \frac{\partial u(x_i)}{\partial x_k} \ 
 \frac{\partial \phi^{(test)}(x_i)}{\partial x_k}  +
 f(x_i) \  \phi^{(test)}(x_i) \right]  \mbox{d}x_1 \mbox{d}x_2 =
 \int_{\partial D_{Neumann}}  
 \phi^{(test)}(\zeta) \ h(\zeta) \ d\zeta. \ \ \ \ \ 
(5) \f]
We now introduce the usual finite-element expansion for \f$ u(x_i)\f$
in terms of the \f$ N \f$ nodal values and the global
shape functions \f$ \psi_j(x_i) \f$, 
\f[ u(x_i) = \sum_{j=1}^{N} \ U_j \ \psi_j(x_i), \f]
where the Dirichlet
boundary conditions determine the nodal values of all nodes
that are located on \f$\partial
D_{Dirichlet}\f$. The nodal values of
the nodes in the interior of the domain and on 
\f$\partial D_{Neumann}\f$ are unknown and must
be determined from the solution of the discrete equations
\f[ f_{{\cal E}(l)} = \int_D \left[ \sum_{j=1}^N U_j \sum_{k=1}^2 
             \frac{\partial \psi_j(x_i)}{\partial x_k} \ 
 \frac{\partial \psi_l(x_i)}{\partial x_k}  +
 f(x_i) \  \psi_l(x_i) \right]  \mbox{d}x_1 \mbox{d}x_2\ \f]
 \f[ - \
 \int_{\partial D_{Neumann}}  
 \psi_l(\zeta) \ h(\zeta) \ \mbox{d}\zeta =0 \\ \mbox{\ \ \ for all $l$ for
  which ${\cal E}(l)>0$,}
\ \ \ \ \ \ 
(6)  \f]
where we have used the global equation number \f${\cal E}(l)\f$ 
to identify nodes whose values are not pinned.
Let us now assume that we have discretised the domain by using
isoparametric elements, so that
within element \f$ e \f$ the solution is represented as
\f[ u(s_1,s_2) = \sum_{j_{local}=1}^n U_{{\cal J}(j_{local},e)} 
\ \psi_{j_{local}}(s_1,s_2) \f]
and 
\f[ x_i(s_1,s_2) = \sum_{j_{local}=1}^n X_{i{\cal J}(j_{local},e)} 
\ \psi_{j_{local}}(s_1,s_2), \ \ \ \ \ \ \
(7) \f]
where \f$ n \f$ is the number of nodes in the element. 
We employ the element-by-element assembly procedure described 
in Algorithm 1 to evaluate the domain integral
in (6).
During this procedure, we visit each element in the mesh, loop over 
the element's \f$ n \f$ nodes 
\f$ l_{local}=1,...,n \f$, check whether the node's global equation number
\f$ i_{global} = {{\cal E}({\cal J}(l_{local},e))}\f$ is positive, and, if so,
add the "bulk" contribution
\f[
\int_{-1}^{1}   \int_{-1}^{1} 
\left[ \sum_{j_{local}=1}^n U_{{\cal J}(j_{local},e)} \sum_{k=1}^2 
\frac{\partial \psi_{j_{local}}}{\partial x_k} \ 
\frac{\partial \psi_{l_{local}}}{\partial x_k}  +
f \  \psi_{l_{local}} \right]  \widehat{\cal J} \ \mbox{d}s_1 \mbox{d}s_2 
\f]
to the \f$ i_{global}\f$-th entry of the global residual vector. 

It remains to assemble the contributions from the line integral around
the boundary of the domain, which can also be done element-by-element.
To illustrate the
procedure, we consider the example shown below: A (topologically 
rectangular) domain has been 
decomposed into 12 four-node elements. Neumann boundary 
conditions are to be applied along the domain's "southern " and 
"eastern" boundaries. The "northern" and "western" boundaries
are subject to Dirichlet boundary conditions.
\image html dirichlet_and_neumann_mesh.gif "Sketch of a discretised domain with Dirichlet and Neumann boundary " 
\image latex dirichlet_and_neumann_mesh.eps "Sketch of a discretised domain with Dirichlet and Neumann boundary " width=0.75\textwidth
A discrete representation of the domain's southern boundary 
is given by the southern edges of elements 1,2 and 3. 
Along these edges the local coordinate \f$ s_1\f$ varies 
from -1 to +1 while \f$ s_2 \equiv -1\f$. We obtain
a piecewise linear approximation for the continuous
domain boundary \f$ {\bf r}^{\partial D}\f$ by substituting
\f$ s_1 =  S \in [-1,1]\f$ and  \f$ s_2 = -1 \f$ into the
mapping (7) between
the local and global coordinates,
\f[
x_i^{\partial D} \approx x_i(s_1=S,s_2=-1) 
\f]
The integral
of any function  \f$ {\cal F}\big(x_1,x_2\big) \f$
along the elements' southern edges can be performed
in terms of element's local coordinates by writing
\f[
\int_{\mbox{Southern edge of element $e$}} {\cal F}\big(x_1(\zeta) ,x_2(\zeta)\big) \ d\zeta = 
\int_{-1}^{1} \left. {\cal F}\big(x_1(s_1,s_2),x_2(s_1,s_2)\big) \right|_{
s_2=-1}
\ \widehat{\cal J}^{(1D)}(s_1) \ ds_1,
\f] 
where
\f[  \widehat{\cal J}^{(1D)}(s_1) = \left. \sqrt{
\left( \frac{\partial x_1(s_1,s_2)}{\partial s_1} \right)^2 + 
\left( \frac{\partial x_2(s_1,s_2)}{\partial s_1} \right)^2 }
\right|_{s_2=-1} 
\f]
is the Jacobian of the mapping between the local coordinate
that parametrises the edge of the element, \f$s_1\f$,  and the arclength,
\f$ \zeta\f$.
Similarly, the domain's eastern boundary is represented
by the eastern edges of elements 3, 6, 9 and 12 whose position
is given by \f$ x_i(s_1=1,s_2=S) \f$ where  \f$ S \in [-1,1]\f$,
and
\f[
\int_{\mbox{Eastern edge of element $e$}} {\cal
F}\big(x_1(\zeta) ,x_2(\zeta)\big) \ d\zeta = 
\int_{-1}^{1} \left. {\cal F}\big(x_1(s_1,s_2),x_2(s_1,s_2)\big) \right|_{
s_1=1}
\ \widehat{\cal J}^{(1D)}(s_2) \ ds_2,
\f]
where
\f[  \widehat{\cal J}^{(1D)}(s_2) = \left. \sqrt{
\left( \frac{\partial x_1(s_1,s_2)}{\partial s_2} \right)^2 + 
\left( \frac{\partial x_2(s_1,s_2)}{\partial s_2} \right)^2 }
\right|_{s_1=1} 
\f]
In general, the line integrals are evaluated
numerically, e.g. using the 1D Gauss rules discussed in 
section \ref local_coords.

Inspection of equation
(6) shows that the line 
integral only makes a contribution to the discrete residuals whose
nodes are located on the Neumann boundaries. [Recall that the 
global shape function \f$ \psi_l \f$ 
is equal to one at global node \f$ l \f$ and vanishes at all
other nodes. Furthermore,  \f$ \psi_l \f$ is nonzero
only in the patch of elements adjacent to node \f$ l \f$ and
vanishes on the boundary of that patch.]
The line integral in (6)
can, therefore, be evaluated by augmenting the assembly 
of the elemental' residual vector as follows:

- Compute the "standard" element residual vector that arises
  from the area integral, as discussed in algorithm 1.
- If the element's southern edge is located on a Neumann boundary:
  - Set \f$s_2=-1\f$.
  - Loop over the 1D Gauss points \f$ i_{int}=1,...,N_{int} \f$.
    - Determine the local coordinate of the integration point
      \f$ s_1 = S_{i_{int}} \f$ and the associated 
      weight \f$W_{i_{int}}.\f$
    - Compute
      - the global coordinates  \f$ x_i = x_i(s_1,s_2) \f$
      - the prescribed-flux function \f$ h=h(x_1,x_2)\f$
      - the shape functions  \f$ \psi_j(s_1,s_2) \f$
      - the Jacobian of the 1D mapping
        \f[
        \widehat{\cal J}^{(1D)} = \left. \sqrt{
        \left( \frac{\partial x_1(s_1,s_2)}{\partial s_1} \right)^2 + 
        \left( \frac{\partial x_2(s_1,s_2)}{\partial s_1} \right)^2 }
        \right|_{s_2=-1} 
        \f]
      .
    - Loop over those local nodes that are located on the
      element's southern edge (e.g. \f$ j_{local}=1,2\f$ for a four-node 
      quadrilateral element).
      - Determine the global node number 
        \f$ j_{global} = {\cal J}(j_{local},e)\f$.
      - Determine the global equation number \f${\cal E}(j_{global})\f$.
      - If \f${\cal E}(j_{global})\ne -1\f$:
        - Determine the local equation number 
        \f$i_{dof}={\cal L}(j_{local},e)\f$.
        - Add the contribution to the element's residual vector
          \f[
          r^{(e)}_{i_{dof}} = r^{(e)}_{i_{dof}} + 
          h \ \psi_{j_{local}} \ \widehat{\cal J}^{(1D)} \ W_{i_{int}}
          \f]
        .
      .
    .
  .
- Perform equivalent operations for the eastern, northern
  and western edges, if required.

Note that it is not necessary to change the setup of the Jacobian matrix
because the boundary terms are independent of the unknown variables. 


This procedure can easily be generalised to higher-dimensional
elements such as 3D Poisson elements for which 
Neumann boundary conditions are applied on the elements'
2D faces.



\subsubsection neumann_implementation Implementation
The procedure outlined in the previous section could easily
be implemented within the general framework of \c oomph-lib. 
One possible implementation would be the following:
- Add a vector of bools, \c Vector<bool> \c Neumann_bc, say,
  to the private data of \c PoissonEquations and use
  the value of \c Neumann_bc[i] to indicate whether
  the element's \c i-th face is subject to 
  a Neumann boundary condition (assuming that \c
  Neumann_bc[i]=false by default.)
- Provide storage for function pointers to the (global) functions
  that specify the prescribed fluxes on the element's
  faces. 
- When assembling the element's residual vector
  in \c PoissonEquations::get_residuals(), execute the 
  assembly loop for the boundary integrals only for those edges
  for which  \c Neumann_bc[i]=true.
- To solve a problem with Neumann boundary conditions we would
  discretise the domain as before, loop over
  all elements that have nodes on \f$ \partial D_{Neumann}\f$, 
 setting the appropriate 
  \c Neumann_bc flag to \c true and specifying
  a function pointer to the prescribed-flux function.
.
The implementation is straightforward and actually
presents an instructive exercise. [Suggestion: Rather than 
incorporating the additional functionality directly
into the existing Poisson elements, use inheritance to create a new
class, \c PoissonEquationsWithFluxBC, say, which contains
the additional data and provides the necessary access functions. 
The derived class must
overload \c get_residuals(...) with a function that initially calls \c
PoissonEquations:get_residuals(...) and then adds the 
contribution from the line integrals]. 

Unfortunately, the implementation suggested above has
several serious disadvantages
-# The test of \c Neumann_bc[i] must be performed 
   in every element, even though only a relatively small number of
   elements are likely to be located on the domain boundaries.
-# The code that implements the computation of the integrals
   along the four possible element edges is very similar
   (the main difference being which local coordinate is held
   constant along the appropriate element edge, and which 
   local nodes make a contribution to the element residual vector)
   but this similarly is not exploited.
-# Finally, and most significantly, many steps in the above
   procedure, such as the identification
   of the nodes that are located on specific element edges, are 
   completely generic and will also occur in the application of
   Neumann-type boundary conditions in other equations (see
   section ***). Incorporating the procedure directly into 
   the Poisson elements makes it impossible to "reuse" the
   relevant algorithms elsewhere.

Our implementation employs a slightly different
approach: rather than incorporating the assembly of the
boundary integral into the Poisson element itself, we treat
the boundary integral over the element's face as a distinct
lower-dimensional \c FiniteElement. We separate
the specification of the element geometry from the implementation
of the specific boundary conditions, and use templating and inheritance
to create Neumann-type elements for entire families 
of elements. Here will discuss the \c PoissonFluxElement class,
a class of \c FiniteElements that impose flux boundary conditions
on the "faces" of associated "bulk" Poisson elements. The
\c PoissonFluxElement class is
templated by the element type of the associated "bulk" element
so that a \c PoissonFluxElement<QPoissonElement<2,3> > 
is the three-node line element (from the \c QElement family) that 
imposes a prescribed-flux boundary condition along one of the
1D edges of a nine-node quadrilateral Poisson element. 
Similarly, a \c PoissonFluxElement<QPoissonElement<3,2> > 
is the four-node quadrilateral element that 
imposes a prescribed-flux boundary condition on a
2D face of an eight-node brick-shaped Poisson element. 


Our implementation involves two new classes:
- The \c FaceElement class implements all the functionality that could 
  conceivably be required by any \c FiniteElement that imposes
  Neumann-type boundary conditions on the faces of higher-dimensional 
  "bulk" elements. (The \c PoissonFluxElements requires none
  of this functionality; and we defer a complete explanation of 
   the structure of this class.)
- The \c FaceGeometry class, a policy class that defines
  the geometry of the lower-dimensional \c FiniteElements
  that form the faces for a given "bulk" element. 
.
We shall illustrate the use of these classes in a "top-down"
example in which we solve Poisson's equation in a rectangular
domain, with Neumann boundary conditions 
along the domain's southern boundary. The full driver code
for this example
(<A HREF="../../to_be_written/html/index.html">prescribed_flux_poisson_driver.cc</A>)
can be found in the 
<A HREF="../../to_be_written/html/index.html">demo directory</A>.
The code fragment below illustrates the application of Neumann 
boundary conditions for this problem and presents the only
change to <A HREF="../../to_be_written/html/index.html">toy_poisson_driver.cc</A>, 
in which the problem is solved with Dirichlet
boundary conditions everywhere.


The code implements the problem by defining a suitable 
problem class, \c FluxPoissonProblem, which we template 
by the element type so that the problem
can be solved with a variety of quadrilateral Poisson elements,
\e e.g. the four-node quad element \c QPoissonElement<2,2>, 
the nine-node quad element  \c QPoissonElement<2,3>, etc.
Here is the constructor for the \c FluxPoissonProblem class:
\code

//========================================================================
/// \short Constructor for Poisson problem with prescribed flux boundary
/// conditions. 
//========================================================================
template<class ELEMENT>
FluxPoissonProblem<ELEMENT>::FluxPoissonProblem()
{ 
 
 // Build and assign mesh (3 x 4 elements in unit square)
 Problem::mesh_pt()=new SimpleRectangularQuadMesh<ELEMENT>(3,4,1.0,1.0);

\endcode

As always, we start by creating a mesh, specifying the type of the
"bulk" element as the template argument in the mesh constructor. 
When the mesh generation process is complete, the "bulk" elements
have been created and are accessible via the various pointer-based 
access functions defined in the \c Mesh class. Now we must
create the prescribed-flux elements and
attach them to those "bulk" elements that have edges which coincide
with the southern domain boundary (mesh boundary 1 in this example).
We use \c Mesh member functions to determine
- the "bulk" elements that are located on that boundary,
- the local coordinate of these elements that remains fixed along the
boundary, and
- the value of the "fixed" local coordinate; in fact, this can only be
the maximum or minimum possible value.
.
\code 

// How many elements are adjacent to boundary 1?
unsigned b=1;
unsigned n_element=mesh_pt()->nboundary_element(b);

// Loop over the bulk elements on mesh boundary b:
for(unsigned e=0;e<n_element;e++)
 {

  // Get pointer to the e-th bulk element on mesh boundary b. (Need to 
  // recast to the specific bulk element type, because 
  // boundary_element_pt() returns a pointer to a FiniteElement. 
  ELEMENT* bulk_elem_pt=dynamic_cast<ELEMENT*>(
   mesh_pt()->boundary_element_pt(b,e));
  
  // Which local coordinate (in the bulk element) is constant
  // along the boundary?
  unsigned s_fixed_index=mesh_pt()->s_fixed_index_at_boundary(b,e);

  // Is the constant coordinate fixed at its max. or min. value?
  int s_limit=mesh_pt()->s_limit_at_boundary(b,e);

  \endcode

 Now we construct the \c PoissonFluxElements (templated by
 the type of the bulk element). Their constructor takes a pointer to
 the bulk element, and the information about the fixed local 
 coordinate that we have just obtained from the \c Mesh:

 \code

  // Build the corresponding prescribed-flux element
  PoissonFluxElement<ELEMENT>* flux_element_pt = new 
  PoissonFluxElement<ELEMENT>(bulk_elem_pt,s_fixed_index,s_limit);

  \endcode
 
  The constructor creates a prescribed-flux element with the
  appropriate element geometry and sets its node pointers
  to the nodes on the appropriate face of the
  "bulk" element. We set the function pointer to the prescribed flux
  function -- here assumed to be the function \c get_flux(...) 
  in the namespace \c GlobalFunctions.

  \code

  // Set the pointer to the prescribed-flux function
  flux_element_pt->flux_fct_pt()=&GlobalFunctions::get_flux;

  \endcode

  Finally, we add the newly created element to the mesh:
 
  \code

  //Add the prescribed-flux element to the mesh
   mesh_pt()->add_element_pt(flux_element_pt);

 }
\endcode

We still need to apply the Dirichlet boundary conditions on the
the remaining boundary nodes. This follows the usual pattern
of looping over the boundary nodes and pinning their values, though
we omit the nodes on boundary 1, because their values
are unknown.

\code 
// Set the boundary conditions for this problem: All nodes are
// free by default -- just pin the ones that have Dirichlet conditions
// here. 
unsigned num_bound = mesh_pt()->nboundary();
for (unsigned b=0;b<num_bound;b++)
 {
  //Leave nodes on boundary 1 free
  if (b!=1)
   {
    unsigned num_nod= mesh_pt()->nboundary_node(b);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      mesh_pt()->boundary_node_pt(b,inod)->pin(0); 
     }
   }
 }

 [...]


} // End of constructor
\endcode

[We have omitted the (standard) code for the assignment of the 
linear solver, the assignment of equation numbers, etc.]

The fine details of the procedure are embodied in
class definition of the \c PoissonFluxElements:

\code

//======================================================================
/// \short A class for elements that allow the imposition of a 
/// prescribed flux on the boundaries of Poisson elements.
/// The element geometry is obtained from the  FaceGeometry<ELEMENT> 
/// policy class.
//======================================================================
template <class ELEMENT>
class PoissonFluxElement : public virtual FaceGeometry<ELEMENT>, 
public virtual FaceElement
{

public:


 /// \short Constructor, takes the pointer to the "bulk" element, the 
 /// index of the fixed local coordinate and its value, represented
 /// by an integer (+/- 1), indicating that the face is located
 /// at the max. or min. value of the "fixed" local coordinate
 /// in the bulk element.
 PoissonFluxElement(FiniteElement* bulk_el_pt, 
                    const unsigned& fixed_local_coord_index, 
                    const int& s_limit);


[...]


}; // End of class definition for PoissonFluxElement

\endcode

\c PoissonFluxElements are 
derived from the \c FaceGeometry class, 
a policy class used to establish the geometry
of the face element as a function of the "bulk" element 
specified by the template parameter. This inheritance structure
is, in fact,  equivalent to the way in which \c QPoissonElements
are derived from geometric elements of the \c QElement family.
The \c FaceGeometry class should be 
defined for every conceivable type of "bulk" element but this
is usually a straightforward process. For
"bulk" elements from the \c QPoissonElement family, the
prescribed-flux elements have the same number of nodes along
their 1D edges as their "bulk" elements but their spatial
dimension is reduced by one. Hence, the complete class definition is:

\code
//=======================================================================
/// Face geometry for the QPoissonElement elements: The spatial 
/// dimension of the face elements is one lower than that of the
/// bulk element but they have the same number of points
/// along their 1D edges.
//=======================================================================
template<unsigned DIM, unsigned NNODE_1D>
class FaceGeometry<QPoissonElement<DIM,NNODE_1D> >: 
 public virtual QElement<DIM-1,NNODE_1D>
{

  public:
 
 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional QElement
 FaceGeometry() : QElement<DIM-1,NNODE_1D>() {}


}; // End of FaceGeometry definition for QPoissonElements

\endcode

hierher: probably need to handle the 0D case separately.

This construction ensures that the number of nodes, etc., of
\c PoissonFluxElements is consistent with that of their
higher-dimensional "bulk" elements. We now need to ensure that the 
\c PoissonFluxElement's node pointers refer to the (already existing)
nodes on the appropriate face of the "bulk" element. 
The process is purely geometric, and
is performed by \c FiniteElement::build_face_element(...) function
of the specific "bulk" geometric elements.

\code

//======================================================================
/// Constructor, takes the pointer to the "bulk" element, the 
/// index of the fixed local coordinate and its value represented
/// by an integer (+/- 1), indicating that the face is located
/// at the max. or min. value of the "fixed" local coordinate
/// in the bulk element.
//======================================================================
template<class ELEMENT>
PoissonFluxElement<ELEMENT>::PoissonFluxElement(
                                  FiniteElement *bulk_el_pt, 
                                  const unsigned& fixed_local_coord_index, 
                                  const int& s_limit) : 
  FaceGeometry<ELEMENT>(), FaceElement()
  { 

   // Let the bulk element build the FaceElement, 
   // i.e. setup the pointers 
   // to its nodes (by referring to the appropriate nodes in the bulk 
   // element), etc.
   bulk_el_pt->build_face_element(fixed_local_coord_index,s_limit,this);
 
   [...]

 }

\endcode

The constructor initially calls the constructor of the appropriate
\c FaceGeometry and thus establishes the element geometry. This is 
followed by a call to the \c FaceElement constructor
which allocates storage for any additional data that needs to be
stored on (generic) face elements. Finally, we call the
"bulk" element's \c build_face_element() function to setup the
node pointers for the present (\c this) element.
[\b Note: \c build_face_element() also assigns the additional data that is
stored in the \c FaceElement base class. As mentioned above,
none of this information is required in the present context,
so we postpone the detailed discussion of this data to 
section ***FREE SURFACE NAVIER-STOKES ELEMENTS***.]

<HR>



\subsection systems_of_PDEs Systems of PDEs: Multiple nodal values

So far, we have only considered scalar ODEs and PDEs but the finite
element method can easily be extended to problems with 
vector-valued unknowns. As a simple example, we first consider a
system of 1D  Poisson equations
\f[ 
\frac{\mbox{d}^2 u_k(x)}{\mbox{d} x} =
f_k(x)  \mbox{ \ \ \ \ for $x\in [0,1]$  
\ \ \ subject to \ \ \ $u_k(0) = g_{0k}$ \ \ \ 
and \ \ \ $u_k(1) = g_{1k}$, }
\mbox{\hspace{5cm}} \f]
where the coefficients \f$g_{0k} \f$ and
\f$g_{1k}\f$ and the functions \f$ f_k(x) \f$ are given, and 
\f$ k=1,...,K\f$.
Since the \f$ K \f$ ODEs have exactly the same form, we obtain
\f$ K \f$ (identical) weak problems
<CENTER>
<TABLE>
<TR>
<TD>
<CENTER>
<B> Problem P3 </B>\n\n
Find the functions  \f$ u_k(x) \in H(D) \f$ (where 
\f$ k=1,...,K\f$) that satisfy the essential boundary conditions
\f[ u_k(0) = g_{0k} \ \ \ \mbox{and} \ \ \   u_k(1) = g_{1k}
\mbox{ \ \ \ \ for $k=1,...,K$,} \f]
and for which
\f[ \int_0^1 \left( \frac{\mbox{d} u_k(x)}{\mbox{d} x} \ 
\frac{\mbox{d} \phi^{(test)}(x)}{\mbox{d} x} \
+ f_k(x) \  \phi^{(test)}(x) \right)\  dx = 0 
\mbox{ \ \ \ \ for $k=1,...,K$,}
\f]
for \e all test functions  \f$ \phi^{(test)}(x) \in H_0(D). \f$
</CENTER>
</TD>
</TR>
</TABLE>
</CENTER>
\n\n
We expand all unknown functions in terms of the same (global) 
shape functions
\f[
u^{(FE)}_k(x) = \sum_{j=1}^{N} U_{kj} \ \psi_j(x).
\f]
Galerkin's method then yields the discrete residuals
  \f[ r_{kj} =  \int_0^1 \left\{
  \sum_{l=1}^{N} U_{kl} \frac{\mbox{d} 
  \psi_l(x)}{\mbox{d} x}
  \frac{\mbox{d} \psi_j(x)}{\mbox{d} x}\ 
  + f_k(x) \  \psi_j(x) \right\} dx \mbox{\ \ \ \ \ 
  \ for $i=2,...,N-1,$ \ \ and $k=1,...,K,$  } \f]
and, as before, we obtain one algebraic equation \f$r_{kj} \f$
for each of the \f$ (N-2)\times K \f$ unknown coefficients
\f$U_{kj}\f$ whose values are not fixed by the boundary conditions.
The finite element solution of these equations can easily be
derived by modifying algorithm @RA[P1_Newton6]@ so that
- every Node stores \f$ K \f$ values (rather than a single scalar).
- the local and global equation numbering schemes are
  modified to generate a linear numbering scheme for the 
  the doubly-indexed unknowns and equations.
.
For this purpose we represent the (global) equation number 
of nodal value \f$ k \f$ at node \f$ j \f$ by \f$ {\cal E}(j,k)\f$
and denote the number of nodal values that are stored at node \f$ j \f$ by 
\f$ {\cal N}(j)\f$.  As before, we reflect the fact that nodal 
values are pinned by setting their global equation number to -1. 
The implementation only requires a few straightforward
changes to Algorithm @RA[P1_Newton6]@ 
\n \n
<TABLE BORDER="1">
<TR>
<TD WIDTH="50%" VALIGN="TOP">
<CENTER> <B> Revised Phase 1c: "Pin" nodal values with essential (Dirichlet) 
     boundary conditions</B></CENTER>  \n
- Loop over all global nodes \f$ j \f$ that are located 
  on Dirichlet boundaries:
  - <TABLE BORDER=0>
    <TR>
    <TD bgcolor="cornsilk">
    Determine the number of values stored at this node, 
    \f$ {\cal N}(j) \f$.
    </TD>
    </TR>
    </TABLE>
    - <TABLE BORDER=0>
      <TR>
      <TD bgcolor="cornsilk">
      Loop over the nodal values, \f$k=1,...,{\cal N}(j)\f$: 
      </TD>
      </TR>
      </TABLE>
    - Assign a negative equation number to reflect the node's "pinned" 
      status:
      \f[ {\cal E}(j,k) = -1 \f]
</TD>
</TR>
</TABLE>
\n
and the  modifications (relative to the original version in Algorithm 
@RA[P1_Newton6]@) have been highlighted. 
The changes to Phase 1d (the provision of an initial guess for
the unknowns and the application of the boundary conditions)
should be obvious. The implementation of 
the global equation numbering scheme requires only one additional 
loop over the nodal values: 
<TABLE BORDER="1">
<TR>
<TD WIDTH="50%" VALIGN="TOP">
<CENTER> <B> Revised Phase 1e: Set up the global equation numbers
with multiple nodal values</B></CENTER>  \n
- Initialise the total number of unknowns, \f$M=0.\f$
- Loop over all global nodes \f$j=1,...,N\f$:
  - <TABLE BORDER=0>
    <TR>
    <TD bgcolor="cornsilk">
    Determine the number of values stored at this node, 
    \f$ {\cal N}(j) \f$.
    </TD>
    </TR>
    </TABLE>
    - <TABLE BORDER=0>
      <TR>
      <TD bgcolor="cornsilk">
      Loop over the nodal values, \f$k=1,...,{\cal N}(j)\f$: 
      </TD>
      </TR>
      </TABLE>
      - If nodal value k is not pinned
        (i.e. if \f$ {\cal E}(j,k) \ne -1 \f$) 
        - Increment the number of unknowns
          \f[ M=M+1 \f]
        - Assign the global equation number
          \f[ {\cal E}(j,k) = M \f]
</TD>
</TR>
</TABLE>
\n
This is the form in 
which the global equation numbering procedure is  
implemented in \c Problem::assign_eqn_numbers(); [Actually, 
this function implements a few additional extensions that we will
discuss in sections \ref discontinuous_interpolations and \ref
solids_lagrange.]

While the global equation numbering scheme is generic and
therefore fully implemented in \c oomph-lib, the 
local equation numbering scheme
and the corresponding assembly of the elements' 
residual vectors and Jacobians needs to implemented differently
for each specific element. Here is an illustration of the
revisions to the relevant phases of Algorithm  @RA[P1_Newton6]@:

 <TABLE BORDER="1">
 <TR>
 <TD>
 <B> Revised Phase 1f: Set up the local equation numbering scheme</B> 
- Loop over the elements  \f$ e=1,...,N_e \f$
  - Initialise the counter for the number of degrees of 
    freedom in this element, \f$ j_{dof}=0 \f$.
  - Loop over the element's local nodes \f$ j_{local}=1,...,n\f$
    - Determine the global node number 
      \f$ j_{global} = {\cal J}(j_{local},e) \f$
    - <TABLE BORDER=0>
      <TR>
      <TD bgcolor="cornsilk">
      Determine the number of values stored at this node, 
      \f$ {\cal N}(j_{local}) \f$.
      </TD>
      </TR>
      </TABLE>
      - <TABLE BORDER=0>
        <TR>
        <TD bgcolor="cornsilk">
        Loop over the nodal values, \f$k=1,...,{\cal N}(j_{local})\f$: 
        </TD>
        </TR>
        </TABLE>
        - Determine the global equation number \f${\cal E}(j_{global},k)\f$
        - If  \f${\cal E}(j_{global},k) \ne -1\f$:
          - Increment the number of degrees of 
            freedom in this element \f$ j_{dof}=j_{dof}+1\f$
          - Add the entry to the lookup scheme that relates
            local and global equation numbers, 
            \f$ \widehat{\cal E}(j_{dof},e) = {\cal E}(j_{global})\f$
          - Store the local equation number associated with the
            current nodal value at the current local node
            \f$ {\cal L}(j_{local},k,e) = j_{dof}\f$
          .
        - Else
          - Set the local equation number associated with the
            current nodal value at the current node to -1 to
            indicate that it is pinned
            \f$ {\cal L}(j_{local},k,e) = -1 \f$
          .
        .
      .
    .
  - Assign the number of degrees of freedom in the element
    \f$ {\cal N}_{dof}(e)=j_{dof} \f$
</TD>
</TR>
</TABLE>

 <TABLE BORDER="1">
 <TR>
 <TD>
     <CENTER> <B> Revision of Phase 2a: Computation of the 
     element's residual vector and Jacobian matrix </B>  </CENTER>
    \n
  - Determine the number of degrees of freedom in this element,
    \f$ n_{dof}={\cal N}_{dof}(e) \f$.
  - Initialise the element residual vector, 
     \f$ r_{j}^{(e)}=0 \f$ for \f$j=1,...,n_{dof} \f$
     and the element Jacobian matrix \f$ J_{jk}^{(e)}=0 \f$ 
     for \f$j,k=1,...,n_{dof} \f$
  - Loop over the Gauss points \f$ i_{int}=1,...,N_{int}\f$
    - Determine the local coordinate of the integration point
      \f$ s = S_{i_{int}} \f$ and the associated weight \f$W_{i_{int}}.\f$
    - Compute
      - the global coordinate  \f$ x = x(s) \f$
      - the source functions   \f$ f_k = f_k(x) =f_k(x(s))
                               \mbox{\ \ \ $(k=1,...,K)$} \f$
      - the derivatives \f$ \mbox{d} u_k/\mbox{d} x 
        \mbox{\ \ \ $(k=1,...,K)$}\f$ and 
      - the shape functions  \f$ \psi_j \f$ and their derivatives
        \f$ \mbox{d} \psi_j/\mbox{d} x \f$
      - the Jacobian of the mapping between local and global
        coordinates, \f$ \widehat{\cal J} = \mbox{d} x/\mbox{d} s\f$
      .
    - Loop over the local nodes \f$ j_{local}=1,...,n \f$
      - Determine the global node number 
        \f$ j_{global}={\cal J}(j_{local},e)\f$
      - <TABLE BORDER=0>
        <TR>
        <TD bgcolor="cornsilk">
        Determine the number of values stored at this node, 
        \f$ {\cal N}(j_{local}) \f$.
        </TD>
        </TR>
        </TABLE>
      - <TABLE BORDER=0>
        <TR>
        <TD bgcolor="cornsilk">
        Loop over the nodal values, \f$k=1,...,{\cal N}(j_{local})\f$: 
        </TD>
        </TR>
        </TABLE>
        - Determine the global equation number 
          \f$ {\cal E}( j_{global},k) \f$
        - If  \f$ {\cal E}( j_{global},k) \ne -1 \f$
          - Determine the local equation number from the element's
            lookup scheme, \f$ i_{dof} = {\cal L}(j_{local},k,e) \f$.
          - Add the contribution to the element's residual vector
            \f[ r_{i_{dof}}^{(e)} = r_{i_{dof}}^{(e)} +
            \left( \frac{\mbox{d} 
            u_k}{\mbox{d} x} \ 
            \frac{\mbox{d} \psi_{j_{local}}}{\mbox{d} x} + 
            f_k \  \psi_{j_{local}} \right) 
            \widehat{\cal J}\  W_{i_{int}} \f]
            - Loop over the local nodes \f$ k_{local}=1,...,n \f$
              - Determine the global node number 
                \f$ k_{global}={\cal J}(k_{local},e)\f$
              - <TABLE BORDER=0>
                <TR>
                <TD bgcolor="cornsilk">
                Determine the number of values stored at this node, 
                \f$ {\cal N}(k_{local}) \f$.
                </TD>
                </TR>
                </TABLE>
              - <TABLE BORDER=0>
                <TR>
                <TD bgcolor="cornsilk">
                Loop over the nodal values, 
                \f$m=1,...,{\cal N}(k_{local})\f$: 
                </TD>
                </TR>
                </TABLE>
                - Determine the global equation number 
                \f$ {\cal E}( k_{global},m) \f$
                 - If  \f$ {\cal E}( k_{global},m) \ne -1 \f$
                   - Determine the local equation number from the element's
                     lookup scheme, \f$ j_{dof} = {\cal L}(k_{local},m,e) \f$.
                   - Add the contribution to the element's Jacobian matrix
                     \f[ J_{i_{dof}j_{dof}}^{(e)} = 
                         J_{i_{dof}j_{dof}}^{(e)} + 
                         \left( 
                         \frac{\mbox{d} \psi_{k_{local}}}{\mbox{d} x} \ 
                         \frac{\mbox{d} \psi_{j_{local}}}{\mbox{d} x}
                         \right)
                         \widehat{\cal J} \ W_{i_{int}}  \f]
                  .
                .
              .
            .
          .
        .
      .
    .
    \n  \n   
</TD>
</TR>
</TABLE>
\n
Again, the modifications are straightforward. At each node
we loop over the multiple nodal values; furthermore, we have added an
additional index to the local equation numbering scheme, so that
\f$ {\cal L}(j_{local},k,e)\f$ stores the local equation number
of the \f$k\f$-th value at local node  \f$j_{local}\f$ in element \f$e\f$.

The example considered here is, of course, very simple. In fact, the
\f$ K \f$ ODEs in Problem P3 are uncoupled and could (should!) 
have been solved independently. However, the 
algorithms developed above can be applied to any systems 
of PDEs in which the \f$ K \f$ differential operators have the same form. 
A more practically relevant example is given by the 3D Navier-Lame 
equations of linear elasticity:
\f[
\sum_{k=1}^3  \left[ (\lambda+\mu) 
\frac{\partial^2 u_k}{\partial x_k\partial x_i} + 
\mu \frac{\partial^2 u_i}{\partial x_k^2} \right]
+ F_i =0 \mbox{\ \ \ \ $(i=1,...,3)$}
\f]
The local equation numbering scheme for these equations is 
identical to that for Problem P3 shown above. Only the computation
of the element's residual vector and Jacobian matrix must  be
changed. This is left as an exercise. 


<B>Implementation</B> \n\n

To allow the object-oriented implementation of the above
algorithms, each Node stores its nodal values in a vector. 
The number of values that need to be stored at the Node
is provided as an argument to its constructor.
The argument to the
\c Node::value_pt(...) function identifies
individual nodal values:
\code
unsigned i_value;
double nodal_value = *(Node::value_pt(i_value));
\endcode
The Node object has a member function that returns the number
of values stored,
\code
unsigned nvalues = Node::nvalue();
\endcode
and the global equation numbers of the individual nodal values 
are accessible via a Node member function
\code
unsigned i_value;
unsigned i_eqn = Node::eqn_number(i_value);
\endcode
Finally, the nodal values' pinned status can be obtained from
\code
unsigned i_value;
bool pinned = Node::is_pinned(i_value);
\endcode
and nodal values are pinned/unpinned with
\code
unsigned i_value,j_value;
Node::pin(i_value);
Node::unpin(j_value);
\endcode
This finally provides the explanation for the mysterious argument
that we had to provide to these functions in all previous
examples: In a scalar problem, each Node only stores a single
value whose index is "0".

  One question still needs to be addressed, however. Nodes are
usually created during the Mesh generation process. How does the
Mesh constructor "know" how many values a given Node should store? 
Also, what should the spatial dimension of the newly created Node be?
(Note that the spatial dimension of a Node is not necessarily the same
as that of the Mesh: For instance, a 1D Mesh can be used to 
discretise a spatially one-dimensional Poisson equation, in which
case each node must store a single (Eulerian) coordinate.
However, the same Mesh could also be used to discretise a 
one-dimensional structural object, such as a beam, in a 3D Eulerian 
space, in which case each Node needs to store three (Eulerian)
coordinates. 

The example considered above shows that the number of 
values that are "required" at a given Node is determined by the 
system of equations that is being solved. This information becomes 
available when a specific system of equations is implemented in a specific
FiniteElement. 
The FiniteElement class therefore provides the pure virtual member function
\code
FiniteElement::required_nvalue(...)
\endcode

which returns the number of values required at a given local Node
in the element. This function should be implemented for every specific
FiniteElement. For a Poisson element (or any other element that implements
the discretisation of a scalar PDE), the function should return 1; for an 
element that discretises the 3D Navier-Lame equations, the function
should return 3; etc. Note that our implementation allows the possibility
that different Nodes in an element require storage for different
numbers of values, see \ref mixed_interpolation.

The required spatial dimension of the Node is, again, determined by each
specific element and should be specified by calling the
protected member function
\code
FiniteElement::set_nodal_dimension(...)
\endcode

For instance, for an element that solves the 2D Poisson equation, this
function must be called in the constructor with argument 2. 

Using the element-specific information provided by these 
virtual FiniteElement member functions, the construction of a Node 
of the appropriate type can be performed by the function
\code
FiniteElement::construct_node(...)
\endcode

This function performs the following tasks:
- it creates a new Node, determining the arguments to the Node
  constructor (spatial dimension, number of values, etc) 
  from the FiniteElement member functions just discussed.
- it stores a pointer to the newly created Node in the element's
  own lookup scheme.
- it returns a pointer to the newly created Node so that it can be
  added to the Mesh's storage scheme for Nodes. 
.

Here is a code fragment that illustrates the creation of Nodes during
the mesh generation process. (Recall that specific Meshes are usually 
templated by the element type, so that a Mesh that was originally 
developed for the solution of a Poisson equation can also be used to 
solve the Navier-Stokes equations, say).

\code
template<class ELEMENT>
class SomeMesh : public Mesh
{

 public:

  /// Mesh constructor
  SomeMesh()
   {

     // Create the first element (Note that element constructors 
     // never have any arguments, therefore this works
     // for any type of element!)
     ELEMENT* el_pt=new ELEMENT;

     // Add the (pointer to the) element to the Mesh's Element_pt
     // vector
     Element_pt.push_back(new ELEMENT);

     // How many nodes does the element have?
     unsigned nnod=element_pt->nnode();

     // Loop over the nodes (we're creating the first element
     // so all nodes are new)
     for (unsigned j=0;j<nnod;j++)
       {
        // Create the element's j-th local node,
        // using the parameters (spatial dimension, number
        // of values, etc) specified by  FiniteElement::required_ndim(...),
        // FiniteElement::required_nvalue(...), etc.
        // Also store a pointer to the newly created Node
        // in the element's own lookup scheme for nodes
        Node* nod_pt=element_pt->construct_node(j);

        // Add the newly created node to the Mesh's own storage scheme
        Node_pt.push_back(nod_pt);
       }


      // Now do the subsequent elements

      [...]


      // Assign the nodal coordinates

      [...]

   }

};
\endcode

 

<HR>

\subsection mixed_interpolation Mixed interpolation

In the examples considered so far, each equation in the system of
PDEs had the same structure and it was 
natural to use the same interpolation for all unknowns
and the same test functions for each PDE in the system. However, 
there are many systems of PDEs whose constituent equations 
have different structures. As an example, we consider the 
steady 2D Navier-Stokes equations, which comprise the two momentum equations
\f[
Re \sum_{k=1}^2  u_k  \frac{\partial u_i}{\partial x_k} +
\frac{\partial p}{\partial x_i} -  
\sum_{k=1}^2 \frac{\partial^2 u_i}{\partial x_k^2} 
 =0 \mbox{\ \ \ \ $(i=1,...,2)$}
\f]
and the continuity equation
\f[
\sum_{k=1}^2  \frac{\partial u_k}{\partial x_k} =0.
\f]
The structure of the momentum equations differs from 
that of the continuity equation. Furthermore, the equations contain 
second derivatives of the two velocity components, 
\f$ u_i \ (i=1,...2),\f$ but only first derivatives of the 
pressure, \f$p\f$. 
This raises two important questions:
- How should we define the weak form for this coupled problem? Should
  we use the same class of test functions for both types of equations
  or should we use a different class for each one?
- How should we represent the two types of unknowns? Should we use
  the same basis functions for the pressure and the velocity or should
  we use a different type of interpolation for each one?
. 
For a general system of PDEs, these are difficult questions which can only
be answered by careful analysis of the governing
equations. Luckily, in the case of the Navier-Stokes equations,
the analysis is well-established and we refer
to the standard literature [e.g. P.M. Gresho and R.L. Sani:
"Incompressible Flow and the Finite Element Method" Vol. 1 and 2.
Wiley (2000)] which provides several possible answers to these
questions. "Equal-order interpolations", in which the
same shape functions are used to represent the pressure
and the velocity components, are generally found to result
in unstable elements. Most stable discretisations are
based on "mixed interpolations" in which the velocity components are
represented by one set of basis functions,
\f[ 
u_i(x_k) = \sum_{j=1}^N U_{ij} \ \psi_j(x_k),\ \ \ \ \ \  (8)
\f]
and another set of basis functions, \f$ \psi^{(P)}_j(x_k) \ 
(j=1,...,N_P)\f$,  is employed to approximate the pressure
\f[ 
p(x_k) = \sum_{j=1}^{N_P} P_j \ \psi^{(P)}_j(x_k). \ \ \ \ \  (9)
\f]
The number of nodes involved in the representation of
the velocity, \f$N\f$, can (and usually does) differ from the
number of nodes involved in the representation of the pressure,
\f$N_P\f$.

Stability constraints (the infamous LBB *REFERENCE*condition) impose
certain restrictions on the choice of basis functions. 
One well-known stable discretisation employs piecewise quadratic interpolations
for the velocities and a piecewise linear interpolation for the
pressure. The resulting elements are known as "Taylor-Hood" elements. 
When implemented as an isoparametric, 2D quadrilateral element, the 
element geometry and the velocity components are represented 
by the bi-quadratic local shape functions \f$ \psi_j(x_k)
\ \ (j=1,...,9)\f$ which interpolate between the nodal values
and coordinates of the element's \f$ n=9 \f$ local nodes. 
The pressure is represented in terms of the bilinear local shape functions, 
\f$ \psi^{(P)}_j(x_k) \ \ (j=1,...,4)\f$, which interpolate between 
pressure values at the element's  \f$ n_p=4 \f$ corner nodes, as 
shown in this sketch:

\image html taylor_hood.gif "Velocity (blue) and pressure (red) nodes in a quadrilateral Taylor-Hood " 
\image latex taylor_hood.eps "Velocity (blue) and pressure (red) nodes in a quadrilateral Taylor-Hood " width=0.35\textwidth

In this case, the number of nodal values is not the same for each
node: the corner nodes store two velocity components and one pressure;
all other nodes only store the two velocity components.

 The weak form of the Navier-Stokes equations is obtained by
using the pressure shape functions \f$ \psi^{(P)}_j(s_k)\f$ as the
test functions for the continuity equation and the velocity shape
functions to test the momentum equations. Upon performing the
usual integration by parts and using the divergence
theorem (Exercise!), the (global) residuals
are then given by
\f[
f_{il} = \int_D  \left[
\left( Re \sum_{k=1}^2  u_k  \frac{\partial u_i}{\partial x_k} \right)
\psi_l - p \frac{\partial \psi_l}{\partial x_i} +  
\sum_{k=1}^2 \frac{\partial u_i}{\partial x_k} 
             \frac{\partial \psi_l}{\partial x_k} \right] dV = 0,
\ \ \ \ \ (10)
\f]
where we have ignored any boundary terms (hence, implying 
the velocity is prescribed along the entire domain boundary), 
and 
\f[
f_{l}^{(P)} = \int_D  \left( \sum_{k=1}^2
 \frac{\partial u_k}{\partial x_k} \right) \psi^{(P)}_l  dV = 0.
\ \ \ \ \ (11)
\f]
The system of nonlinear algebraic equations (10)
provides one equation, \f$ f_{il} \f$, for each velocity
component \f$ (i=1,2) \f$ at those nodes \f$(l=1,...,N)\f$ where
the velocity is not prescribed by a Dirichlet boundary condition.
The weak continuity equations (11) provides 
the remaining  equations \f$ f_{l}^{(P)} =0 \f$ for the
pressure values  \f$ P_j \ (j=1,...,N_p)\f$ that are not prescribed.
[Usually there are no boundary conditions for the pressure -- the
pressure can only be set indirectly through applied-traction boundary 
conditions in which the normal component of the applied traction
contains a contribution from the pressure and the viscous normal stresses; only
if the velocity is prescribed along the entire boundary is it legal,
and, in fact, necessary,  to pin a single pressure value.]

<B>Implementation</B> \n\n

The implementation of an element with continuous 
mixed interpolation is straightforward:
- The function \c FiniteElement:required_nvalue(...) needs to specify
  the required number of values for each node.
- The element can "recycle" the geometric shape functions to implement
  the isoparametric representation of the velocity components. 
- The element needs to define another set of basis functions 
  which is used to represent the pressure.  
- The local equation numbering scheme needs to be augmented to 
  include the pressure degrees of freedom. 
.

<HR>


\subsection discontinuous_interpolation Discontinuous interpolation -- "Internal Data"

The choice of basis functions in Taylor-Hood-type Navier-Stokes elements
appears to be perfectly natural. The bi-quadratic interpolation of the velocity
components between the element's nodal values enforces the inter-element 
continuity of the velocity, which ensures that the integrals in
(10) "make sense". The LBB constraint requires the 
use of a lower-order representation for the pressure; the "obvious"
choice is bilinear 
interpolation between the pressure values at the four corner nodes,
which also ensures inter-element continuity of the
pressure. An
inspection of the weak form of the Navier-Stokes equations, 
(10) and (11),
shows that this is not strictly required because the integration by parts has
removed the derivatives on the pressure. Hence the weak form "makes
sense" even if the pressure is globally discontinuous. This suggests the
use of different types of basis functions for the representation of
the pressure. An example is given by the Crouzeix-Raviart Navier-Stokes 
elements in which the pressure within each element is represented 
by the three basis functions
\f[
\psi_1^{(P)}(s_1,s_2) = 1,\ \ \ \ \ \  
\psi_2^{(P)}(s_1,s_2) = s_1,\ \ \ \ \ \  
\psi_3^{(P)}(s_1,s_2) = s_2,
\f]
so that the pressure in element \f$ e \f$ is given by the complete
linear polynomial,
\f[
p(s_1,s_2) = \sum_{j_p=1}^3 P_{j_p e}  \ \psi_j^{(P)}(s_1,s_2).
\f]
In this representation, the discrete pressure values,
\f$P_{j_p e}\f$, are not associated with any Nodes -- instead,
they represent "internal" data for each individual element.

<B>Implementation</B> \n\n

To facilitate the implementation of discontinuous interpolations
we introduce two modifications to our data structure:
-# We split the \c Node class into two separate classes:
   - The base class \c Data which contains only the "values"
     and their global equation numbers.
   - The derived class \c Node which inherits from \c Data and
     adds the additional functionality that is required by actual Nodes
     (storage for the Eulerian coordinates and the associated 
     access functions).
-# Each FiniteElement provides storage for "Internal Data" in
   a vector of (pointers to) \c Data objects. By default, the 
   vector of "Internal
   Data" is a vector of zero length. Specific FiniteElements
   can resize this vector in their constructor to provide a sufficient
   amount of storage for their Internal Data (if they require any).
   Internal Data is included in the
   global equation numbering procedure by providing 
   the function \c FiniteElement::assign_internal_eqn_numbers()
   which allocates distinct global equation numbers to all values
   that are stored in the element's "Internal Data". the
   function is executed for each element when the global equation 
   numbers are assigned by \c Problem::assign_eqn_numbers(). 
.

<HR>

\subsection higher_order_pdes Higher-order PDEs -- Hermite interpolants

All the examples considered so far involved second-order PDEs
whose weak form could be integrated by parts. Therefore we only needed
to evaluate integrals that involve first derivatives
of the basis functions. Such integrals can be evaluated
for basis functions that are globally continuous but have
e.g. discontinuous first derivatives at isolated points (at the
element boundaries). The piecewise
linear/quadratic/... (global) finite-element shape functions 
introduced in section \ref fe_basis_fct
satisfied these criteria. Since these representation of the function
with these shape functions implements a (Lagrange) interpolation 
between the nodal values, such finite elements are often referred
to as Lagrange finite elements. 

Let us now consider the fourth-order problem 
\f[ \frac{\mbox{d}^4 u(x)}{\mbox{d} x^4} = f(x)  \mbox{ \ \ \
for $x\in[0,1]$\ \ \ \ \ \ \ \ \ (1)}  \f]
subject to the essential boundary conditions
\f[ u(0)=g_0, \ \  u(1)=g_1, \ \ 
\left. \frac{\mbox{d}u}{\mbox{d}x}\right|_{x=0}=h_0, \ \  
\left.\frac{\mbox{d}u}{\mbox{d}x}\right|_{x=1}=h_1,\f]
where \f$f(x)\f$ and the constants \f$g_0,g_1,h_0\f$ and \f$h_1\f$ are
given. 

We define the weak form of (1) by requiring that the (weak) solution
satisfies the essential boundary conditions and that the
weighted residual
\f[
r = \int_0^1  \left( \frac{\mbox{d}^4 u(x)}{\mbox{d} x^4} - f(x) \right)
 \phi^{(test)}(x) \ \mbox{d}x \ \ \ \ \ \ \ \ \ (2)
\f]
vanishes for all test functions  \f$ \phi^{(test)}(x) \f$ which satisfy
the homogeneous boundary conditions
\f[ \phi^{(test)}(0)=0, \ \  \phi^{(test)}(1)=0, \ \ 
\left. \frac{d\phi^{(test)}}{dx}\right|_{x=0}=0, \ \  
\left.\frac{d\phi^{(test)}}{dx}\right|_{x=1}=0.
\ \ \ \ \ \ \ \ \ (3) \f]
We can integrate the weak form (2) twice by parts to obtain
\f[
r = \int_0^1  \left(
    \frac{\mbox{d}^2 u(x)}{\mbox{d} x^2} 
    \frac{\mbox{d}^2 \phi^{(test)}(x)}{\mbox{d} x^2}
     - f(x) \phi^{(test)}(x) \right) \ dx
   \ \  + \ \ \left[ \left. \frac{\mbox{d}^3 u(x)}{\mbox{d} x^3} 
              \phi^{(test)}(x) \right] \right|_{0}^1
   \ \  - \ \ \left[ \left. \frac{\mbox{d}^2 u(x)}{\mbox{d} x^2}
       \frac{\mbox{d} \phi^{(test)}(x)}{\mbox{d} x} \right]
      \right|_{0}^1.
\f]
As before, the boundary terms vanish because \f$ \phi^{(test)}(x)\f$  
satisfies the homogeneous boundary conditions (3), and we obtain
\f[
r = \int_0^1  \left(
    \frac{\mbox{d}^2 u(x)}{\mbox{d} x^2} 
    \frac{\mbox{d}^2 \phi^{(test)}(x)}{\mbox{d} x^2}
     - f(x) \phi^{(test)}(x) \right) \ \mbox{d}x \ \ \ \ \ \ \ \ \ (4).
\f]
We now follow the usual (Galerkin) procedure of
expanding the solution \f$ u(x) \f$ and the test function 
\f$ \phi^{(test)}(x)\f$ in terms of
the same set of basis functions \f$ \psi_j(x)\f$. The integral (4) shows
that the basis functions now
need to have square-integrable second derivatives, i.e. 
we have to require that \f$ \psi_j(x) \in H^2_0(D)\f$ rather
than  \f$ \psi_j(x) \in H^1_0(D)\f$ as in the case of the
second-order problem P1. 
The Lagrange interpolants that are generated by piecewise
linear or quadratic finite element shape functions are not members
of  \f$ H^2_0(D)\f$ because their first derivatives are discontinuous between
elements. 

Basis functions of sufficient differentiability can be
constructed by Hermite interpolation. For this  purpose we introduce
two different types of basis functions:
- Basis functions of type 1 provide an interpolation between the
  nodal values and thus satisfy the conditions
  \f[ \psi_{j1}(X_k) = \delta_{jk}, \f]
  implying that they are equal to one at exactly \e one node while they 
  vanish at all others. Furthermore, these basis functions satisfy
  \f[ \left. \frac{\mbox{d}\psi_{j1}}{\mbox{d}x} \right|_{X_k} = 0 
  \mbox{\ \ for $j,k=1,...,N$}, \f]
  so that their derivatives vanish at \e all nodal points.
- Basis functions of type 2 provide an independent interpolation
  for the derivatives of the function and thus satisfy the conditions
  \f[ \left. \frac{\mbox{d}\psi_{j2}}{\mbox{d}x} \right|_{X_k} = \delta_{jk} , \f]
  implying that their derivative is equal to one at exactly 
  \e one nodal point and zero at all others. Furthermore, these
  basis functions satisfy
  \f[ \psi_{j2}(X_k) = 0 \mbox{\ \ for $j,k=1,...,N$}, \f]
  so that they vanish at \e all nodal points.

The global Hermite shape functions are illustrated in this sketch:

\image html 1Dmesh_with_hermite_shape_fcts.gif "1D Hermite shape " 
\image latex 1Dmesh_with_hermite_shape_fcts.eps "1D Hermite shape " width=0.75\textwidth

The global Hermite shape functions have local equivalents which are
again expressed in terms of a local coordinate \f$s \in [0,1]\f$
(note that the coordinate range differs from that for the Lagrange 
elements) as

\todo ***ANDREW SHOULD WE CHANGE THE RANGE TO [-1,1] ?? ***

\f[ \psi_{11}(s) =  2s^3 - 3s^2 + 1 \f]
\f[ \psi_{12}(s) = s^3 - 2s^2 + s \f]
\f[ \psi_{21}(s) = 3s^2 - 2s^3 \f]
\f[ \psi_{22}(s) = s^3 - s^2. \f]

The shape functions for [-1,1] are
\f[ \psi_{11}(s) =  \frac{1}{4}\left(s^3 - 3s + 2\right) \f]
\f[ \psi_{12}(s) = \frac{1}{4}\left(s^3 - s^2 - s + 1 \right) \f]
\f[ \psi_{21}(s) = \frac{1}{4}\left(2 + 3s - s^3\right) \f]
\f[ \psi_{22}(s) = \frac{1}{4}\left(s^3 + s^2 - s - 1\right). \f]

Using these shape functions we can represent the unknown function \f$ u \f$ 
in element \f$e\f$ as
\f[ u = \sum_{j=1}^2  \sum_{k=1}^2 U_{{\cal J}(j,e)k}\ \psi_{jk}(s), \f]
where each node now stores two "types" of degree of freedom: \f$U_{j1} \f$
represents the nodal value at node \f$ j \f$ and  \f$U_{j2} \f$ represents
(a measure of) the derivative of the function at that node.  
In isoparametric elements, 
the same shape functions are used to represent the mapping 
between \f$x\f$ and \f$ s\f$ so we have
\f[ x(s) = \sum_{j=1}^2  \sum_{k=1}^2 X_{{\cal J}(j,e)k}\ \psi_{jk}(s). \f]
Here  \f$X_{j1} \f$ represents the \f$x\f$ coordinate of nodal point
\f$j\f$. We refer to the coefficients \f$X_{j2} \f$ as "generalised
coordinates". During the mesh-generation process the values of
these generalised coordinates must be chosen such that the
mapping \f$ x(s)\f$ is one-to-one over the entire element.

<HR>


\subsection hyperbol_parabol_timestep Parabolic and hyperbolic problems -- timestepping

The PDEs discussed so far were all of elliptic type. However,
the finite element method can also be used for the (semi-)discretisation
of parabolic and hyperbolic problems. We will demonstrate
this by considering the two classical examples of
these equations: i) the unsteady heat equation and ii) the wave equation.
\subsubsection unsteady_heat The unsteady heat equation
The 2D unsteady heat equation is given by 
\f[
\sum_{k=1}^2 \frac{\partial^2 u(x_i,t)}{\partial x_k^2} 
= \frac{\partial u(x_i,t)}{\partial t} + f(x_i,t)
\mbox{\ \ \ \ \ in $D$}
\f]
subject to the boundary condition
\f[
u\left(x_i^{\partial D}(\zeta),t\right) = g(\zeta,t) 
\ \ \ \ \ \ \ \ (12)
\f]
and the initial condition
\f[
u(x_i,t=t_0) = u^{(IC)}(x_i), \ \ \ \ \ \ \ \ (13)
\f]
for the given functions \f$ f(x_i,t), \ g(\zeta,t) \f$
and \f$ u^{(IC)}(x_i)\f$.
Superficially, the only difference between this parabolic problem
and the elliptic Poisson problem P2 is the presence of the 
time-derivative on the RHS of the PDE and the fact that the
boundary conditions and the source function are time-dependent. 
Motivated by this observation, we employ the same spatial discretisation
as in the Poisson problem by writing
\f[
u(x_i,t) = \sum_{j=1}^N U_j(t) \ \psi_j(x_i)
\f]
where the coefficients  \f$ U_j(t) \f$ are now
(continuous) functions of time \f$ t \f$. Following the
usual procedures, we obtain the residuals
\f[
r_l(t) = \int_D \left[ \sum_{j=1}^N \left( U_j(t) 
\sum_{k=1}^2
\frac{\partial \psi_j(x_i)}{\partial x_k} 
\frac{\partial \psi_l(x_i)}{\partial x_k} +
\frac{\mbox{d} U_j(t)}{\mbox{d} t} \psi_j(x_i) \ \psi_l(x_i) \right) +
f(x_i,t) \ \psi_l(x_i) \right] \mbox{d}x_1 \mbox{d}x_2 = 0 
\ \ \ \ \ \ (14)
\f]
for all coefficients \f$ U_l(t) \f$ whose values are not determined
by the boundary conditions. The residuals now form
a system of ODEs, rather than a system of algebraic equations --
the spatial discretisation has only achieved a semi-discretisation.
We obtain the initial conditions for the unknown coefficients 
\f$ U_l(t) \f$ by evaluating (13) at the nodal
positions
\f[
U_j(t=t_0) = u^{(IC)}(X_{ij})
\ \ \ \ \ \ (15)
\f]
The solution of the initial-value problem defined by
the system of first-order ODEs (14) and the 
initial conditions (15) can now be obtained by
any ODE solver. Within \c oomph-lib we provide the implicit 
timesteppers from the BDF family as the default
timesteppers for parabolic problems. In these schemes,
the time-derivatives are approximated by
\f[
\frac{\mbox{d} U_j}{\mbox{d} t} \approx 
\sum_{\tau=0}^{N_{history}}  \left.^\tau U_j\right.  \left.^\tau W\right.
\ \ \ \ \ \ \  (16)
\f]
where the values \f$ \left.^\tau U_j\right. \f$  represent
the history of the nodal value \f$  U_j \f$, i.e. its values
at current and previous (discrete) timesteps. 
If the time-integration is performed with a constant timestep
 \f$ \Delta t \f$, we have
\f[
 \left.^\tau U_j\right. =  U_j(t-\tau \Delta t)
\f]
so \f$ \left.^0 U_j\right. \f$
represents the (unknown) nodal value at the present time.
The weights \f$ \left.^\tau W\right. \f$ are constants that depend
only on the timestep \f$ \Delta t\f$.  And, yes, 
the \f$ \tau \f$ is a left superscript.

The different members of the BDF family are distinguished by the
number of previous values used in the approximation (16).
For the BDF1 scheme we have \f$ N_{history}=1\f$ with the
weights
\f[  
\left.^0 W\right. = \frac{1}{\Delta t} \mbox{\ \ \ \ and \ \ \ }
\left.^1 W\right. = - \frac{1}{\Delta t}
\f]
Inserting the weights into (16), we recover the 
classical first-order-accurate backward Euler scheme:
\f[
\frac{\mbox{d} U_j}{\mbox{d} t} \approx 
\frac{ \left.^0 U_j\right. -  \left.^1 U_j\right. }{\Delta t}
\f]
Similarly, in the BDF2 scheme we have \f$ N_{history}=2\f$ with the
weights
\f[  
\left.^0 W\right. = \frac{3}{2\Delta t}, \ 
\left.^1 W\right. = - \frac{2}{\Delta t}\mbox{\ \ \ \ and \ \ \ }
\left.^2 W\right. = \frac{1}{2\Delta t}
\f]
which yields the second-order backward Euler scheme:
\f[
\frac{\mbox{d} U_j}{\mbox{d} t} \approx 
\frac{ 3\left.^0 U_j\right. -  4 \left.^1 U_j\right.  +  
1\left.^2 U_j\right. }{2\Delta t}
\f]
Inserting the approximations (16) for the time-derivatives
into  (14) gives a system of algebraic 
equations 
 \f[
 r_l = \int_D \left\{ \sum_{j=1}^N \left.^{0}U_j \right.  \left[
 \sum_{k=1}^2 \left(
 \frac{\partial \psi_j(x_i)}{\partial x_k} 
 \frac{\partial \psi_l(x_i)}{\partial x_k} \right) +
 \left.^{0} W \right. 
 \psi_j(x_i) \ \psi_l(x_i) \right] +
 \sum_{\tau=1}^{N_{history}} 
 \left.^{\tau} W \right.  \left.^{\tau} U_j\right. 
 f(x_i,t) \ \psi_l(x_i) \right\} dx_1 dx_2 = 0 
 \ \ \ \ \ \ (17)
 \f]
that can be solved to determine the unknown values of the coefficients 
\f$ \left.^{0} U_j\right. =  U_j \f$ at the present time.

<CENTER> <B>Algorithm *** Timestepping</B></CENTER>
-# Set time \f$ t \f$ to its initial value, \f$ t=t_0. \f$
-# Choose a timestep \f$\Delta t \f$.
-# Compute the weights \f$ \left.^{\tau} W \right. \f$
   for \f$ \tau=0,...,N_{history}\f$.
-# Assign initial values for all history 
   values so that the initial conditions are satisfied
   at  \f$ t=t_0 \f$.
-# Move to the next timestep:
   - Increase the value of the continuous time: 
     \f[
     t=t+\Delta t
     \f]
   - Update the history values. For a BDF scheme, this
     involves shifting the history values back by one
     time level:
     - Loop over the time levels backwards, \f$ \tau=N_{history},...,1 \f$
       \f[
       \left.^{\tau} U_j\right. = \left.^{\tau-1} U_j\right.
       \f]
     .
   - Update the nodal values of all nodes on the boundary according to the
     time-dependent boundary conditions (12).
   - Solve ((17)) (evaluating the
     source function at the present time \f$ t \f$) for the
     unknown nodal values \f$ \left.^{0} U_j\right. =  U_j \f$.
   - Document the solution.
-# Go to *** for the next timestep.

A key step in this algorithm is the provision of
history values for the nodal values at \f$ t=t_0 \f$. Recall that the
state of the system at this time is determined
by the initial condition (13) which determines
the values of \f$ \left.^{0} U_j\right. \f$. 
The value of \f$ \left.^{1} U_j\right. \f$ is
not determined by  (13) but could be obtained
by solving (***) 

Before advancing to the next (=first) timestep, we shift
this history values back therefore the initial value for
 \f$ \left.^{0} U_j\right. \f$ is copied into 
\f$ \left.^{1} U_j\right. \f$ and  \f$ \left.^{0} U_j\right. \f$
is regarded as an unknown -- we can (and usually do) retain its initial
value as the initial guess for the solution at \f$ t = t_0+\Delta t\f$.
If we timestepping is performed with a BDF1 scheme in which
\f$ \left.^{1} U_j\right. \f$ is the only history value,
the scheme is therefore "self-starting" in the sense that
the provision of the initial conditions (13)
completely specify initial values for all required history
values. 

\subsubsection timestep_implementation Implementation
The implementation of timestepping procedures in \c oomph-lib
is fairly straightforward

<HR>

\subsection solids_lagrange Solid mechanics problems -- Lagrangian coordinates
All problems discussed so far were formulated in an Eulerian frame,
i.e. we assumed all unknowns to be functions of the spatially-fixed,
Eulerian coordinates, \f$x_i \ (i=1,...,3)\f$, and of time, 
\f$ t \f$. Many problems 
in solid mechanics are formulated more conveniently in terms of 
body-attached, Lagrangian coordinates. In this section we will
briefly review the essential concepts of nonlinear continuum mechanics
and present the principle of virtual displacements which forms the
basis for all large-displacement solid mechanics computations in \c oomph-lib.

Throughout this section we will use the summation convention 
so that repeated indices are assumed to be summed over the range
of the three spatial coordinates and all free indices range 
from 1 to 3. We will retain the summation signs
for all other sums, such as sums over the nodes etc.

\subsubsection geometry The geometry
The figure below introduces the essential geometrical concepts: 
We employ a set of Lagrangian coordinates, \f$ \xi^i\f$, to
parametrise the (Eulerian) position vector to material points in the
 body's undeformed position:
\f[
{\bf r} = {\bf r}(\xi^i).
\f]
The specific choice of the Lagrangian coordinates is 
irrelevant for the subsequent development. For analytical studies, it is
advantageous to choose a body-fitted Lagrangian coordinate
system (as shown in the sketch) because this allows boundary
conditions to be applied on iso-levels of the coordinates;
in computational studies, it is usually preferable to use a 
coordinate system in which the governing equations are as compact
as possible. A cartesian coordinate system is best suited for this 
purpose.

We denote the tangent vectors to the coordinate lines \f$ \xi^i = const.\f$
in the undeformed configuration by
\f[
{\bf g}_i = \frac{\partial {\bf r}}{\partial \xi^i}
\f]
and define the components of the covariant metric tensor via 
the inner products
\f[
g_{ij} = {\bf g}_i \cdot {\bf g}_j.
\f]
This tensor defines the "metric" because the (square of the) infinitesimal 
length, \f$ ds \f$, of a line element with coordinate increments 
\f$ d\xi^i\f$ is given by
\f[
(ds)^2 = g_{ij} \  d\xi^i d\xi^j.
\f]
The volume of the infinitesimal parallelepiped formed by
the coordinate increments  \f$ d\xi^i\f$ is given by
\f[
dv = \sqrt{g} \ d\xi^1 d\xi^2 d\xi^3,
\f]
where 
\f[
g = \det g_{ij}. 
\f]
 
\image html Lagrangian_coord.gif "2D sketch of the Lagrangian coordinates. The Lagrangian coordinates parametrise the position vector to material points in the body. As the body deforms, the Lagrangian coordinates remain attached to the same material points and the coordinate lines become distorted. The change in the length of infinitesimal material line elements provides an objective measure of the body's " 
\image latex Lagrangian_coord.eps "2D sketch of the Lagrangian coordinates. The Lagrangian coordinates parametrise the position vector to material points in the body. As the body deforms, the Lagrangian coordinates remain attached to the same material points and the coordinate lines become distorted. The change in the length of infinitesimal material line elements provides an objective measure of the body's " width=0.75\textwidth
As the body deforms, the Lagrangian coordinates remain "attached" to the
same material points. The body's deformation can therefore be
described bt the vector field that specifies the position vectors 
to material particles in the deformed configuration,
\f[
{\bf R} = {\bf R}(\xi^i).
\f]
As in the undeformed coordinate system, we  form the 
tangent vectors to the deformed coordinate lines \f$ \xi^i = const.\f$
and denote them by
\f[
{\bf G}_i = \frac{\partial {\bf R}}{\partial \xi^i}.
\f]
The inner product of these vectors defines the metric tensor in 
the deformed configuration
\f[
G_{ij} = {\bf G}_i \cdot {\bf G}_j
        = \frac{\partial {\bf R}}{\partial \xi^i} \cdot  
          \frac{\partial {\bf R}}{\partial \xi^j},
\f]
and we have equivalent relations for the lengths of line elements,
\f[
(dS)^2 = G_{ij} \  d\xi^i d\xi^j
\f]
and the volume of infinitesimal parallelepipeds,
\f[
dV = \sqrt{G} \ d\xi^1 d\xi^2 d\xi^3,
\f]
where
\f[
G = \det G_{ij}. 
\f] 
Since the metric tensors \f$ G_{ij }\f$ and \f$ g_{ij }\f$ provide 
a measure of the length of 
material line elements in the deformed and undeformed configurations,
respectively, their difference
\f[
\gamma_{ij} =  \frac{1}{2} (G_{ij} - g_{ij})  
\f]
provides an objective measure of the strain and is known as 
the Green strain tensor. 

\subsubsection equilibrium Equilibrium and the Principle of Virtual Displacements

Let us now assume that the body is subjected to
- an applied surface traction  \f$ {\bf t}   \f$ --
  a force per unit deformed surface area, applied on (part of) the
  body's deformed surface area \f$ A_{tract} \subset \partial V\f$,
- a body force \f$ {\bf f} \f$ --
  a force per unit volume of the undeformed (!) body [This can
  easily be expressed in terms of a force per unit deformed volume, 
 \f$ {\bf F }\f$, by invoking conservation of mass, which shows that
 \f$ {\bf f } = \sqrt{G/g} \ {\bf F }\f$.],
. 
while its displacements are prescribed on the remaining part of
the boundary, \f$ A_{displ}  \subset \partial V \f$ (where 
\f$ A_{tract} \cap A_{displ}=0 \f$ and  \f$A_{tract} \cup A_{displ}
= \partial V\f$ ), so that
\f[
{\bf R}(\xi^k) = {\bf R}^{(BC)}(\xi^k) 
\mbox{ \ \ \ \ on \ \ \ $ A_{displ},$ }
  \ \ \ \ \ \ \ \ (19)
\f]
where \f$ {\bf R}^{(BC)}(\xi^k) \f$ is given.
\image html Solid_boundary_conds.gif "Sketch of the boundary conditions: The surface of the body is either subject to a prescribed traction or its displacement is " 
\image latex Solid_boundary_conds.eps "Sketch of the boundary conditions: The surface of the body is either subject to a prescribed traction or its displacement is " width=0.75\textwidth
The deformation is governed by the principle of virtual displacements
\f[
\int \left\{ \sigma^{ij}  \ \delta \gamma_{ij} -  \left( {\bf f} - 
 \rho \frac{\partial^2 {\bf R}}{\partial t^2} \right) \cdot 
 \delta {\bf R} \right\} \ dv -
  \oint_{A_{tract}} {\bf t} \cdot \delta {\bf R} \ dA =0,
  \ \ \ \ \ \ \ \ (20)
 \f]
 where \f$ \sigma^{ij} = \sigma^{ji}\f$ is the (symmetric) 
 second Piola-Kirchhoff stress tensor, \f$ \rho \f$ is the density, and
 \f$ \delta (.) \f$ represents the variation of \f$ (.)\f$.
 See, e.g., Green, A.E. \& Zerna, W. "Theoretical Elasticity", 
 Dover (1992); or Wempner, G. "Mechanics of Solids with Applications 
 to Thin Bodies", Kluwer (1982) for more details.

 Upon choosing the particles' position vector in the deformed 
 configuration, \f${\bf R}(\xi^j)\f$, as the unknown, the variation 
 of the strain tensor becomes
 \f[
 \delta \gamma_{ij} = \frac{1}{2} \left( 
         \frac{\partial {\bf R}}{\partial \xi^i} \cdot  
         \delta \ \frac{\partial {\bf R}}{\partial \xi^j} +
         \frac{\partial {\bf R}}{\partial \xi^j} \cdot  
         \delta \ \frac{\partial {\bf R}}{\partial \xi^i}
       \right)
 \f]
 and we have
 \f[
  \delta {\bf R}  = {\bf 0} \mbox{ \ \ \ \ on \ \ \ $ A_{displ}.$ }
  \ \ \ \ \ \ \ \ (21)
 \f]
 The 2nd Piola Kirchhoff stress tensor is symmetric, therefore
 we can write the variation of the strain energy in
 (20) as
 \f[
 \int  \sigma^{ij}  \ \delta \gamma_{ij} \ dv = 
 \int \sigma^{ij}  \ \frac{\partial {\bf R}}{\partial \xi^i} \cdot  
          \delta \ \frac{\partial {\bf R}}{\partial \xi^j} \ dv
  \f]
 and obtain the displacement form of the principle of 
 virtual displacements:
 \f[
 \int \left\{ \sigma^{ij} \ \frac{\partial {\bf R}}{\partial \xi^i} \cdot  
          \delta \ \frac{\partial {\bf R}}{\partial \xi^j} 
      - \left( {\bf f} - 
  \rho \frac{\partial^2 {\bf R}}{\partial t^2} \right) \cdot 
 \delta {\bf R} \right\} \ dv 
 - \oint_{A_{tract}} {\bf t} \cdot \delta {\bf R} \ dA =0.
 \f]
This needs to be augmented by a constitutive equation which
determines the stress as a function of the body's deformation,
(and possibly the history of its deformation). Here we will only
consider elastic behaviour, where the stress is only a function of the
strain.

\subsubsection elastic_constitutive Constitutive Equations for Purely Elastic Behaviour
For purely
elastic behaviour, the stress is only a function of the
instantaneous, local strain and the constitutive equation
has the form
\f[
\sigma^{ij} =  \sigma^{ij} (\gamma_{kl}).
\f]

The functional form of the constitutive equation is different
for compressible/incompressible/near-incompressible behaviour:
-# \b Compressible \b Behaviour: \n If the material is compressible,
   the stress can be computed directly from the deformed and undeformed
   metric tensors,
   \f[
  \sigma^{ij} = \sigma^{ij}(\gamma_{kl}) =  
  \sigma^{ij}\bigg( \frac{1}{2} (G_{kl} - g_{kl})\bigg).
  \f]
-# \b Incompressible \b Behaviour: \n If the material is incompressible,
   its deformation is constrained by the condition that
   \f[
   \det G_{ij} - \det g_{ij}= 0 
   \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (22) 
   \f]
   which ensures that the volume of infinitesimal material
   elements remains constant during the deformation. This
   condition is typically enforced by a Lagrange multiplier which
   plays the role of a pressure. In such cases, the 
   stress tensor has the form
   \f[
   \sigma^{ij} = -p G^{ij} + 
   \overline{\sigma}^{ij}\big(\gamma_{kl}\big),
   \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (23)  
   \f]
   where only the deviatoric part of the stress tensor, 
   \f$ \overline{\sigma}^{ij}, \f$ depends directly on the
   strain. The pressure \f$ p \f$ needs to be determined
   independently by enforcing the incompressibility constraint
   (22).
   Given the deformed and undeformed metric tensors,
   the computation of the stress tensor \f$ \sigma^{ij} \f$ 
   for an incompressible
   material therefore requires the computation of the following quantities:
   - The deviatoric stress \f$ \overline{\sigma}^{ij} \f$
   - The contravariant deformed metric tensor  \f$ G^{ij} \f$
   - The determinant of the deformed 
     metric tensor \f$ \det G_{ij} \f$ which 
     is required in equation (22)
     whose solution determines the pressure.
   .
  \n
-# \b Nearly \b Incompressible \b Behaviour: \n If the material is nearly
   incompressible, it is advantageous to split the stress into
   its deviatoric and hydrostatic parts by writing the
   constitutive law in the form
   \f[
   \sigma^{ij} = -p G^{ij} + 
   \overline{\sigma}^{ij}\big(\gamma_{kl}\big),
   \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (24)
   \f]
   where the deviatoric part of the stress tensor, 
   \f$ \overline{\sigma}^{ij}, \f$ depends on the
   strain. This form of the constitutive
   law is identical to that of the incompressible
   case and it involves a pressure \f$ p \f$ which needs to be
   determined from an additional equation. In the 
   incompressible case, this equation was given by the incompressibility
   constraint (22). Here, we need to augment 
   the constitutive law for the deviatoric stress by
   an additional equation for the pressure. Generally this takes the
   form
   \f[
   p = \kappa \ d
   \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (25)
   \f]
   where \f$ \kappa \f$ is the "bulk modulus", a material property
   that needs to be specified by the constitutive law.   
   \f$ d \f$ is the (generalised) dilatation, i.e. the relative change
   in the volume of an infinitesimal material element (or some
   suitable generalised quantity that is related to it). As the 
   material approaches incompressibility, \f$ \kappa \to \infty\f$, so that
   infinitely large pressures would be required to achieve
   any change in volume. To facilitate the implementation of
   (25) as the equation for the pressure, we re-write 
   it in the form
   \f[
   p \ \frac{1}{\kappa} -  d\big(g_{ij},G_{ij}\big) = 0
   \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (26) 
   \f]
   which only involves quantities that remain finite 
   as we approach true incompressibility.
   \n 
   Given the deformed and undeformed metric tensors,
   the computation of the stress tensor \f$ \sigma^{ij} \f$ 
   for a nearly incompressible
   material therefore requires the computation of the following quantities:
   - The deviatoric stress \f$ \overline{\sigma}^{ij} \f$
   - The contravariant deformed metric tensor  \f$ G^{ij} \f$
   - The generalised dilatation \f$ d \f$
   - The inverse of the bulk modulus \f$ \kappa \f$
   .
  \n

The abstract base class \c ConstitutiveLaw provides
interfaces for the computation of the stress in all three forms.

\subsubsection non-dim_solid Non-dimensionalisation
The principle of virtual displacements (20)
is written in dimensional form.
We generally prefer to work with non-dimensional quantities and
will now discuss the non-dimensionalisation used in the
actual implementation of the equations in \c oomph-lib. For this purpose we 
first re-write equation (20) as
\f[
\int_{v} \left\{ \sigma^{*ij} \ \delta \gamma_{ij} -  \left( {\bf f}^* - 
 \rho \frac{\partial^2 {\bf R}^*}{\partial t^{*2}} \right) \cdot 
 \delta {\bf R}^* \right\} \ dv^* 
 - \oint_{A_{tract}} {\bf t}^* \cdot \delta {\bf R}^* \ dA^* =0,
  \ \ \ \ \ \ \ \ (27)
 \f]
where have used asterisks to label those
dimensional quantities that will have non-dimensional equivalents.
(Some quantities, such as the strain, are already
dimensionless, while others, such as the density will not have
any non-dimensional counterparts. We do not introduce modifiers
for these).

We now non-dimensionalise all lengths with a problem-specific
length-scale, \f$ {\cal L}\f$, given e.g. the length of the solid
body, so that
\f[
\xi^{*i} = {\cal L} \ \xi^{i}, \ \ \  {\bf R}^* = {\cal L} \ {\bf R}, 
\ \ \ dA^* = {\cal L}^2 \ dA, \mbox{ \ \ and \ \ } 
dv^* = {\cal L}^3 \ dv.\ \ \ 
\f]
We use a characteristic stiffness, \f$ {\cal S}\f$, (e.g. the
material's Young's modulus \f$ E \f$) to non-dimensionalise
the stress and the loads as
\f[
\sigma^{*ij} = {\cal S} \ \sigma^{ij}, \ \ \  {\bf t}^* = {\cal S}\
{\bf t}, \mbox{ \ \ and \ \ } {\bf f}^* = {\cal S} / {\cal L} \
{\bf f},
\f]
and we non-dimensionalise time with a problem-specific
timescale \f$ {\cal T}\f$ (e.g. the period of some external forcing),
so that
\f[
t^* = {\cal T} t.
\f]
This transforms (27) into
\f[
\int_{v} \left\{ \sigma^{ij} \ \delta \gamma_{ij} - \left( {\bf f} - 
 \Lambda^2 \frac{\partial^2 {\bf R}}{\partial t^{2}} \right) \cdot 
 \delta {\bf R} \right\} \ dv 
 -\oint_{A_{tract}} {\bf t} \cdot \delta {\bf R} \ dA =0,
  \ \ \ \ \ \ \ \ (28)
 \f]
 where 
\f[
\Lambda = \frac{1}{\cal T} \sqrt{\frac{\rho {\cal L}}{\cal S}}
\f]
is the ratio of system's "intrinsic" timescale, 
\f$T_{intrinsic} = \sqrt{\rho {\cal L}/{\cal S}},\f$ to the
time, \f${\cal T}\f$, used in the non-dimensionalisation of the
equations. If a given problem has no externally imposed timescale
(e.g. in the free vibrations of an elastic body) \f$T_{intrinsic}\f$
(or some suitable problem-dependent multiple thereof) provides a natural 
timescale for the non-dimensionalisation. Therefore we use 
\f$ \Lambda=1 \f$ as the default value in all solid mechanics
equations. If preferred, computations can, of course, be
be performed with dimensional units, provided all quantities
are expressed in consistent units (e.g. from the SI system).
In this case \f$ \Lambda^2 \f$ represents the dimensional
density of the material.

We adopt a similar approach for non-dimensionalisation of the 
constitutive equations. Typically, the constitutive parameters
(e.g. Young's modulus and Poisson's ratio for a Hookean
constitutive equation) are passed to the \c ConstitutiveLaw
as arguments to its constructor. Where possible, we select one of these
parameters as the reference stress \f$ {\cal S} \f$ and give it
a default value of 1.0. Hence, if a Hookean constitutive
law is instantiated with just one argument (the Poisson ratio \f$
\nu \f$), the stress is assumed to have been scaled on Young's
modulus. If two arguments are provided, the second argument
should be interpreted as the ratio of
the material's Young's modulus to the reference stress 
\f$ {\cal S} \f$ used in the non-dimensionalisation of the equations.

\subsubsection 2D_solid 2D problems: Plane strain.
Many solid mechanics problems can be regarded as two-dimensional
in the sense that the quantities of interest only depend on two
spatial coordinates. In such problems it is important to consider what
constraints the system is subjected to in the third direction -- clearly
all real solid bodies are three-dimensional!
In plane strain problems, the displacements of material points are assumed
to be parallel to the 2D plane spanned by the coordinates \f$x_1
\f$ and \f$ x_2\f$, so that any displacements 
normal to this plane are suppressed. In this case, significant 
transverse stresses can develop.
Conversely, in plane stress problems, it is assumed that no stresses
develop in the transverse direction; in this case we need
to allow material particles to be displaced transversely.
Since we have formulated our problem in terms of positions 
(i.e. in terms of the displacement of material points), our formulation
naturally produces a <B>plane strain</B> problem if we reduce
the equations to two dimensions:
We assume that the transverse displacement vanishes, 
while the remaining (in-plane) displacements are only 
functions of the in-plane coordinates, \f$ x_1\f$ and \f$ x_2\f$.

It is important to remember that the 2D version of all equations needs to 
produce plane-strain behaviour when new, strain-energy-based 
constitutive equation classes are formulated:
When implementing such strain-energy functions, recall that any
invariants of metric tensors etc. are the invariants of the full
3D quantities which are not necessarily the same as those
of the corresponding 2D quantities.



\subsubsection isotropic_growth Isotropic growth.
Many biological tissues undergo growth processes. For instance, the cells 
that make up a solid tumour divide regularly and as a result the
total mass of the tumour increases. If the growth 
occurs non-uniformly so that certain parts of the tumour grow faster
than others, regions that grow more slowly restrain the expansion
of their neighbours. This process can induce significant growth stresses.
The scenario is similar (but not identical) to that of thermal 
growth in which a non-uniform temperature distribution in 
a body creates thermal stresses. (The important difference between
these two cases is that in the latter process mass is conserved --
thermal expansion only leads to an increase in the volume occupied
by the material, whereas biological growth via cell division 
increases the mass of the tumour). 

It is easy to incorporate such growth processes into our 
theoretical framework. The general idea is sketched in the 
following figure:




<TABLE>
<TR>
<TD COLSPAN=2>
<CENTER>
\image html Isotropic_growth_0.gif " " 
\image latex Isotropic_growth_0.eps " " width=0.35\textwidth
\n <B> State 0:</B> Undeformed, ungrown and 
stress-free reference configuration.
</CENTER>
</TD>
</TR>
</TABLE>

<TABLE>
<TR>
<TD>
\image html Isotropic_growth_U1.gif " " 
\image latex Isotropic_growth_U1.eps " " width=0.35\textwidth
\n <B> (Hypothetical) state U1:</B> The individual infinitesimal 
material elements
have expanded isotropically and the elements are in 
a stress-free state. The isotropic growth
is spatially uniform -- all elements have expanded by the same amount.
</TD>
<TD>
\image html Isotropic_growth_N1.gif " " 
\image latex Isotropic_growth_N1.eps " " width=0.35\textwidth
\n<B> (Hypothetical) state N1:</B>  All infinitesimal material elements
have expanded isotropically and the elements are in 
a stress-free state. The isotropic growth
is spatially non-uniform, so individual elements
have grown by different amounts.
</TD>
</TR>
</TABLE>
<TABLE>
<TR>
<TD>
\image html Isotropic_growth_U2.gif " " 
\image latex Isotropic_growth_U2.eps " " width=0.35\textwidth
\n<B> State U2:</B> Since the isotropically-grown infinitesimal
 material elements have grown by the same amount, they can
 be (re-)assembled to form a continuous, stress-free body. 
</TD>
<TD>
\image html Isotropic_growth_N2.gif " " 
\image latex Isotropic_growth_N2.eps " " width=0.35\textwidth
\n<B> State N2:</B> Since the individual material elements have
grown by different amounts they can (in general) not be (re-)assembled
to form a continuous body without undergoing some deformation.
The deformation of the material elements (relative to their
stress-free shape in the hypothetical, stress-free state N1) 
induces internal stresses -- the so-called growth-stresses.
</TD>
</TR>
</TABLE>
<TABLE>
<TR>
<TD>
\image html Isotropic_growth_UE.gif " " 
\image latex Isotropic_growth_UE.eps " " width=0.35\textwidth
\n<B> State UE:</B> When subjected to external tractions 
and to body forces, the uniformly-grown material elements deform. Their
deformation (relative to their stress-free shape in state
U1 (or, equivalently, U2), generates internal stresses which
- balance the applied loads,
- keep all infinitesimal material elements in local equilibrium.
</TD>
<TD>
\image html Isotropic_growth_NE.gif " " 
\image latex Isotropic_growth_NE.eps " " width=0.35\textwidth
\n<B> State NE:</B>  When subjected to external tractions 
and to body forces, the material elements deform further. Their
deformation (relative to their stress-free shape in state
N1), generates internal stresses which
- balance the applied loads,
- keep the material elements in local equilibrium.
</TD>
</TR>
</TABLE>




We start our analysis with the stress-free (and 
"ungrown") reference configuration "0" at the top of the diagram,
and initially follow the deformation shown in the left half
of the sketch. In a first step, each infinitesimal material
element in  the body is subjected to the same isotropic growth
which changes its mass from  \f$\mbox{d}M\f$ to
\f$ \Gamma \ \mbox{d}M \f$. Assuming that the growth process does not
change the density of the material, \f$ \Gamma\f$ also specifies
the volumetric growth of the material. All material
elements grow by the same amount, therefore the individual elements
can be (re-)assembled to form a continuous body (state
U2). In this state, the body is 
stress-free and, compared to the reference configuration "0",
all lengths have increased by a factor of \f$\Gamma^{1/d}\f$ 
where \f$ d \f$ is the body's spatial dimension (i.e. \f$d=2\f$ in 
the sketch above). The covariant basis vectors in this 
uniformly-grown, stress-free configuration are therefore given by
\f[
\widetilde{\bf g}_i = \Gamma^{1/d} \ {\bf g}_i,
\f] 
the metric tensor is given by
\f[
\widetilde{g_{ij}} = \Gamma^{2/d} \ g_{ij},
\f] 
and the volume of an infinitesimal material element
has increased from 
\f[
dv = \sqrt{g}  \ d\xi^1 ... \ d\xi^D
\f]
to 
\f[
\widetilde{dv} = \Gamma \sqrt{g} \ d\xi^1 ... \ d\xi^D 
\ \ \ \ \ \ \ \ \ \ \ \ (29)
\f]
We now subject the stress-free, uniformly-grown body to
external loads and determine its deformation, using the principle
of virtual displacements. Since the uniformly grown state "U2"
is stress-free, we may regard it as the reference state for
the subsequent deformation. The strain tensor that describes
the deformation from the stress-free (and uniformly-grown) state "U2"
to the final equilibrium configuration "UE" must therefore
be defined as
\f[
\gamma_{ij} = \frac{1}{2} \left(G_{ij} - \widetilde{g_{ij}}\right)
 = \frac{1}{2} \left(G_{ij} - \Gamma^{2/D} \ g_{ij} \right)
\ \ \ \ \ \ \ \ \ \ \ \ (30)
\f]
Equations (29) and (30) allow
us to express the principle of virtual displacements
in terms of
- the metric of the undeformed (and un-grown) reference
  state "0", 
- the volumetric growth \f$ \Gamma \f$, and
- quantities associated with the final deformed configuration
  "UE":
.
\f[
\int_{v} \left\{ \sigma^{ij} \ \delta \gamma_{ij} - \left( {\bf f} - 
 \Lambda^2 \frac{\partial^2 {\bf R}}{\partial t^{2}} \right) \cdot 
 \delta {\bf R} \right\} \ \Gamma \ dv 
 -\oint_{A_{tract}} {\bf t} \cdot \delta {\bf R} \ dA =0,
  \ \ \ \ \ \ \ \ (31)
\f]
Note that this equation does not contain any references to quantities
in the intermediate states "U1" and "U2". 

We will now consider the case of spatially non-uniform growth,
illustrated in the right half of the sketch. If the isotropic
growth is spatially non-uniform, \f$ \Gamma = \Gamma(\xi^i)\f$,
the growth will try to expand all infinitesimal material 
elements isotropically -- but each one by a different amount as
illustrated by the hypothetical state N1 in which each material
element has expanded to its stress-free state in which its
metric tensor is given by
\f[
\widetilde{g_{ij}} = \Gamma(\xi^k) \ g_{ij}.
\f] 
Material elements will only be stress-free if the strain 
\f[
\gamma_{ij} = \frac{1}{2} \left(G_{ij} - \widetilde{g_{ij}}\right)
 = \frac{1}{2} \left(G_{ij} - \Gamma^{2/D}(\xi^k) \ g_{ij} \right)
\ \ \ \ \ \ \ \ \ \ \ \ (32)
\f]
relative to their isotropically grown shape in state N1 is zero.
In general, the displacements induced by such an isotropic
expansion will be
incompatible and it would be impossible to (re-)assemble the 
individually grown material elements to a continuous body
unless the material elements undergo some deformation. 
The elements' deformation relative to their
stress-free shape in N1 will generate internal "growth-stresses"
(stage N2). When subjected to external loads and body forces
the body will undergo further deformations until the
the stress (generated by the particles' total deformation relative to
their stress-free state N1) balances the applied loads.

It is important to realise that, as in the case of spatially
uniform growth, the strain defined by (30) is 
an intrinsic quantity that provides a measure of each particles' 
\e local deformation relative to its stress-free shape in N1. 
The intermediate (and in the current case hypothetical), isotropically 
grown state N1 does not appear in the analysis -- it only serves
to define the stress-free shape for each infinitesimal
material element. Equation (31) therefore
remains valid.


\subsubsection CartesianLagrangian Specialisation to a cartesian basis and finite element discretisation
If the problem does not have any
symmetries (e.g. axisymmetry) whose exploitation would reduce
the spatial dimension of the problem, the most compact form
of the equations is obtained by resolving all vectors into
a fixed cartesian basis so that the undeformed position vector
is given by
\f[
{\bf r}(\xi^j) = r_i(\xi^j) \ {\bf e}_i,
 \ \ \ \ \ \ \ \ (33)
\f]
Similarly, we write 
\f[
{\bf f} = f_i \ {\bf e}_i,
\f]
\f[
{\bf t} = t_i \ {\bf e}_i,
\f]
and 
\f[
{\bf R}^{(BC)} = R^{(BC)}_i \ {\bf e}_i.
\f]

We use the Eulerian coordinates in the undeformed 
configuration as the Lagrangian coordinates so that
\f[
r_i(\xi^j) = \xi^i. \ \ \ \ \ \ \ \ (34)
\f]
With this choice, the tangent vectors to the undeformed coordinate lines
are the cartesian basis vectors
\f[
{\bf g}_i =  {\bf e}_i, 
\f]
and the undeformed metric tensor is the Kronecker delta
(the unit matrix)
\f[
g_{ij} = \delta_{ij}.
\f]
We now expand the (unknown) deformed position vector in the same basis,
\f[
{\bf R}(\xi^j) = R_i(\xi^j) \ {\bf e}_i,
\f]
and derive a finite element approximation for this vector field
from the principle of virtual displacements. For this purpose we decompose
the undeformed body into a number of finite elements, using the
standard mesh generation process described previously.
This establishes the Eulerian position \f$ X_{ij}^{(0)} \
(j=1,...,N) \f$ of the \f$ N \f$
nodes in the body's undeformed configuration.
Since the Eulerian coordinates of material points in the undeformed 
configuration coincide with their Lagrangian coordinates (see 
(34)), a finite-element representation of the 
Lagrangian coordinates is established by writing
\f[
\xi^i = \sum_{j=1}^{N} \Xi_{ij} \ \psi_j,
\f]
where \f$ \Xi_{ij} = X_{ij}^{(0)} \f$ is the \f$i\f$-th 
Lagrangian coordinate of global node \f$j\f$, and the \f$\psi_j\f$
are the global finite-element shape functions. In practice, the
\f$\psi_j\f$ are, of course, represented by local shape functions, 
\f$\psi_j(s_k)\f$, so that the Lagrangian coordinate at local 
coordinate \f$ s_k \f$ in element \f$ e \f$ is given by
\f[
\xi^i(s_k) = \sum_{j=1}^{n} \Xi_{i{\cal J}(j,e)} \ \psi_j(s_k).
\f]
We employ the same basis functions to represent the components of 
the unknown vector field \f$ {\bf R}(\xi^j)\f$, by writing
\f[
R_i(\xi^k) = \sum_{j=1}^{n} X_{ij} \ \psi_j(\xi^k),
\ \ \ \ \ \ \ \   (35)
\f]
and treat the (Eulerian) nodal positions \f$ X_{ij} \f$ as 
the unknowns. With this discretisation, the 
variations in  \f$ {\bf R}(\xi^j) \f$ correspond to variations
in the nodal positions \f$ X_{ij} \f$ so that
\f[
\delta {\bf R} =   \sum_{j=1}^{N} \delta X_{ij} \ \psi_j \ {\bf e}_i
\f]
\f[
\delta \frac{\partial {\bf R}}{\partial \xi^k} = 
\sum_{j=1}^{N} \delta X_{ij} \ 
\frac{\partial \psi_j}{\partial \xi^k} \ {\bf e}_i,
\f]
\f[
\mbox{ etc. }
\f]
The principle of virtual displacement  (31) 
therefore becomes
 \f[\sum_{m=1}^{N}
 \left\{
  \int \left[ \sigma^{ij} \ 
  \left(\sum_{l=1}^N X_{kl} \frac{\partial \psi_l}{\partial \xi^i} \right)
       \frac{\partial \psi_m}{\partial \xi^j} \
      - \left( f_k - 
  \Lambda^2 \left(\sum_{l=1}^N
    \frac{\partial^2 X_{kl}}{\partial t^2} \psi_l \right) \right) 
 \psi_m \right] \ \Gamma \ dv 
 - \oint_{A_{tract}} t_k \ \psi_m \ dA \right\} 
 \delta X_{km}  =0. \ \ \ \  (36)
 \f]
[Note that summation convention enforces the summation over the index 
\f$ k \f$.]
The displacement boundary condition (19)
determines the positions of all nodes that are located on the
boundary \f$A_{displ}\f$, 
\f[
 X_{ij} = R^{(BC)}_i(\Xi_{lj}) \mbox{ \ \ \ if node $j$ is 
 located on $A_{displ},$}
\f]
and equation (21) 
requires their variations to vanish,
\f[
 \delta X_{ij} = 0 \mbox{ \ \ \ if node $j$ is located on $A_{displ}.$}
\f]
The variations of all other nodal positions are arbitrary (and independent
of each other), therefore the terms in the curly brackets in 
(36) need to vanish individually. This provides one 
(discrete) equation for each unknown \f$X_{km}\f$
 \f[ f_{km} = 
  \int \left[ \sigma^{ij} \ 
  \left(\sum_{k=1}^N X_{kl} \frac{\partial \psi_l}{\partial \xi^i} \right)
       \frac{\partial \psi_m}{\partial \xi^j} \
      - \left( f_k - 
  \Lambda^2 \left(\sum_{k=1}^N
    \frac{\partial^2 X_{kl}}{\partial t^2} \psi_l \right) \right) 
 \psi_m \right] \ \Gamma \ dv \
 - \oint_{A_{tract}} t_k \ \psi_m \ dA 
  =0. \ \ \ \ \  \ \ \ \ (37)
 \f]
These equations can again be assembled in an element-by-element
fashion.

\subsubsection solid_implementation Implementation
We will now discuss how the discrete equations (37) are
implemented in \c oomph-lib. To facilitate the analysis
of multi-physics problems, we introduce generalisations of the \c Node, 
\c FiniteElement and \c Mesh classes which provide
separate storage (and access functions) for all solid mechanics
Data. The resulting \c SolidFiniteElements can 
be used as stand-alone elements for the simulation of 
pure solid mechanics problems. More importantly, however, the design
makes it easy to employ multiple inheritance to create 
more complex elements that solve the equations of solid mechanics
together with 
any other field equations. For instance, if we combine
a \c FiniteElement that solves the unsteady heat equation
with an \c SolidFiniteElement that describes the elastic 
deformations, we obtain an element that can be used to simulate 
unsteady heat conduction in an elastic body that performs 
large-amplitude oscillations, say.

<CENTER> <B>The SolidNode class:</B></CENTER> 
The class \c SolidNode is derived from \c Node and 
implements the additional functionality required for solid
mechanics problems. The key new feature is that each Node
needs to store its (fixed) Lagrangian coordinate \f$ \Xi_{ij}\f$,
while its Eulerian position  \f$ X_{ij}\f$ needs to be regarded
as an unknown. This requires the following changes to
member functions of the \c Node class:

- The function \c Node::x(...) returns the
  Eulerian coordinates of the \c Node. Internally, this function
  accesses the nodal coordinates via pointers to double
  precision numbers. This data structure is adequate for
  problems in which the nodal positions are either known a priori 
  or determined by an algebraic node update function. However, the
  scheme does not allow us to regard the nodal positions as
  unknowns -- as required by our expansion
  (35). In \c SolidNodes
  the nodal positions are therefore stored as \c Data (whose 
  values can either be unknowns or pinned values, and which can have a time
  history). We overload the \c Node::x(...) function by 
  \code
  SolidNode::x(...)
  \endcode
  which returns the (now variable) Eulerian nodal position of the node.
  The function
  \code
  SolidNode::position_eqn_number(...)
  \endcode
  gives access to the global equation number for each
  nodal coordinate and, following our usual convention, the function
  returns -1 if a coordinate is pinned.
- We introduce a new member function 
  \code
  SolidNode::xi(...)
  \endcode
  which returns the (fixed) Lagrangian coordinates of the node.
- Both functions have the usual extensions to generalised
  coordinates
  \code
  SolidNode::x_gen(...)
  \endcode
  and 
  \code
  SolidNode::xi_gen(...)
  \endcode
  which are required for Hermite elements and any other elements
  that use generalised nodal coordinates to interpolate the
  element's geometry.
- Similar to the \c Data member function \c Data::pin(...) which can
  be used to pin specific nodal values, we provide the function
  \code
  SolidNode::position_pin(...)
  \endcode
  which allows pinning of selected nodal coordinates.
- Dynamic problems require the evaluation of time-derivatives
  of the nodal positions, such as \f$ \partial^2 X_{ij}/\partial t^2\f$,
  see (37). These time-derivatives are evaluated by the
  \c TimeStepper of the positional \c Data.  By default,
  we use the same \c TimeStepper for the nodal \c Data and for the nodal 
  positions. In multi-physics problems this might not be 
  appropriate. Consider, for instance, solving an unsteady 
  heat equation in a dynamically deforming, elastic body. In this problem
  the 2nd time-derivatives of the nodal position might be evaluated
  by a \c Newmark scheme, acting on the history values of the nodal
  positions, whereas the time-derivatives of the temperature
  might be determined by a \c BDF scheme, operating on the 
  history values of the nodal \c Data.
  In such cases, the default assignment for the two timesteppers 
  can be overwritten with the access functions
  \code
  SolidNode::position_time_stepper_pt()
  \endcode
  and 
  \code
  SolidNode::time_stepper_pt()
  \endcode
  where the latter is inherited from \c Data::time_stepper_pt().
- Our implementation is based on the displacement form of the
  principle of virtual displacements in which the position vector 
  \f$ {\bf R}(\xi^i) \f$ in the deformed configuration is regarded
  as the unknown vector field. Equation  (35)
  defines the representation of this vector
  field within each finite element in terms of the nodal coordinates.
  Some constitutive equations require the representation
  of additional (non-positional) variables.
  For instance, for incompressible (or nearly incompressible) materials,
  the stress \f$ \sigma^{ij} \f$ contains a contribution from the
  (scalar) pressure field \f$ p \f$; see the discussion of the
  constitutive equations above. If we choose a continuous
  representation for the pressure in which its value is
  interpolated between nodal values (as in Taylor-Hood-type
  elements), the nodal pressure values could in principle be stored
  in the "normal" nodal Data. However, to retain a clear 
  distinction between
  solid-mechanics and non-solid-mechanics quantities we store
  such "additional nodal solid values" in their own \c Data object 
  to which we provide access with the function
  \code
  SolidNode::additional_solid_data_pt()
  \endcode
  As with all other nodal quantities, the number of additional 
  solid data values, needs to be specified as an argument to the
  \c SolidNode constructor when the \c SolidNode is created.
  The access function
  \code
  SolidNode::nadditional_solid_value()
  \endcode
  returns the number of additional solid values that are stored
  at an \c SolidNode; the function
  \code
  SolidNode::additional_solid_eqn_number(...)
  \endcode
  returns the global equation number of the SolidNode's 
  additional solid values, and there are further member functions
  that allow pinning/unpinning of these values, etc.
- Finally, \c SolidNodes overload the function 
  \c Data::assign_eqn_numbers()
  with 
  \code
  SolidNode::assign_eqn_numbers()
  \endcode
  which creates global equation numbers for all (non-pinned) 
  positional values and any additional (nodal) solid values, and
  then deals with the "normal" nodal Data by calling
  \c Data::assign_eqn_numbers().

<CENTER><B>The SolidFiniteElement class:</B></CENTER>
The class \c SolidFiniteElement is derived from \c FiniteElement and 
implements the additional functionality required for solid
mechanics problems. Again, most of the additional (or revised)
functionality is related to the presence of the two coordinate
systems which requires the following changes to \c FiniteElement
member functions:
- The nodes of \c SolidFiniteElements are \c SolidNodes, therefore
  we overload the function \c FiniteElements::construct_node(...)
  to ensure that an \c SolidNode with the appropriate 
  amount of storage is built. As in the case of \c FiniteElements, 
  the required spatial dimension of the elements' constituent 
  nodes, their number of nodal values etc. are specified via
  the \c FiniteElement's (pure) virtual member functions 
  \c FiniteElement::required_ndim(...),
  \c FiniteElement::required_nvalue(...), etc, which need to be
  implemented for all specific elements that are derived from
  the \c SolidFiniteElement base class.  As discussed above, the 
  constructor of the \c SolidNodes requires additional parameters, such
  as the number of Lagrangian coordinates. These need to be specified 
  by implementing \c SolidFiniteElement::required_nlagrangian(...) 
  and similar other functions. As in the 
  case of \c FiniteElements, many of these functions are already
  implemented as virtual (rather than pure virtual) functions which provide
  sensible default values. For instance, by default, the number
  of additional nodal solid values is set to zero, etc. Such
  functions need to be overloaded in specific derived elements
  if the default values are not appropriate.
- The interpolation of the Eulerian coordinates, implemented
  in \c FiniteElement::interpolated_x(...), can remain unchanged
  because of the way in which we have overloaded \c Node::x(...) with
  SolidNode::x(...). \c SolidFiniteElements provide additional
  functions, such as 
  \code
  SolidFiniteElement::interpolated_xi(...)
  \endcode
  which determines the interpolated Lagrangian coordinates at
  a specified local coordinate within the element, or
  \code
  SolidFiniteElement::xi(...)
  \endcode
  which provides a wrapper for the nodal values of the
  Lagrangian coordinates. 
- The displacement form of the principle of virtual displacements 
  (37) contains derivatives of the shape functions
  with respect to the Lagrangian (rather than the Eulerian)
  coordinates. Their computation is implemented in
  \code
  SolidFiniteElement::dshape_lagrangian(...)
  \endcode
- We have already discussed the need to represent
  additional, non-positional solid variables (such as pressures)
  within an \c SolidFiniteElement. The storage 
  of the "additional solid values" at the \c SolidNodes
  allows a globally continuous representation of such variables
  (similar to the continuous pressure representation in Taylor-Hood
  elements). Discontinuous representations
  (as in Crouzeix-Raviart elements) require the 
  storage of additional internal solid variables. In principle, these
  could again be stored in the \c FiniteElement's internal Data
  but for the reasons outlined at the beginning of this section, 
  we store \e all solid variables in separate Data items.
  The function 
  \code
  SolidFiniteElement::internal_solid_data_pt(...)
  \endcode
  provides access to the additional internal solid \c Data.
  The total number of \c Data items stored as 
  internal solid \c Data  can be obtained from 
  \code
  SolidFiniteElement::ninternal_solid_data() 
  \endcode
- We have now created storage and access functions for all the 
  additional \c Data that is required for solid mechanics problems: 
  - Data that represents the nodal positions,
  - Data that represents additional, non-positional nodal Data
    such as continuously interpolated pressures, and
  - Data that represents additional internal Data, used, e.g.,
    to represent discontinuously interpolated pressures.
  .
  We need to ensure that
  these Data items are included in the various equation numbering
  schemes. This only requires overloading of 
  \c FiniteElement::assign_internal_eqn_numbers(...)
  with a version that also sets up the global equation numbers
  for the additional internal solid values. 
  Furthermore, we provide the function
  \c SolidFiniteElement::assign_solid_local_eqn_numbers() 
  which sets up the local equation numbering scheme for all
  solid Data associated with an element. This function will be
  called from the overloaded 
  \c SolidMesh::assign_local_eqn_numbers()
  function to be discussed in the next section.
- We're done! \c SolidFiniteElements now form a suitable basis for
  all elements whose deformation is determined
  by the equations of solid mechanics (or some variant thereof). 
  To implement a specific
  solid mechanics element, we need to represent its geometry,
  its state of stress, etc., in terms of the \c
  SolidFiniteElement's nodal and additional internal solid
  \c Data. This requires the specification of the shape functions
  and the functions that compute the element's Jacobian matrix
  and its residual vector -- the latter implementing the 
  element's contribution to the global residual vector defined by the 
  discretised principle of virtual displacements, (37).
  As for "normal" \c
  FiniteElements it is sensible to construct specific \c
  SolidFiniteElements in a hierarchy which separates between
  the implementation of the governing equations and the 
  representation of the element geometry. For instance,
  the \c SolidQElement family represents the generalisation of the 
  \c QElement family to \c SolidFiniteElements, while \c
  PVDEquations implement the principle of virtual
  displacements (37). The two are combined
  by multiple inheritance to form the \c QPVDElement class.
- The computation of the
  element's Jacobian matrix requires the evaluation of
  the derivatives of the discrete residuals (37) with respect
  to the unknown nodal positions, and with respect to any
  additional unknown solid mechanics variables (e.g. pressures). The
  derivatives with respect to the nodal positions result in 
  fairly complex algebraic expressions and in our experience it
  is more efficient to evaluate these entries by finite
  differencing. This operation is completely generic
  and is implemented in 
  \code
  SolidFiniteElement::get_jacobian(...)
  \endcode 
  It is therefore not necessary to implement this function for
  \c SolidFiniteElements that only involve positional degrees of
  freedom (unless the setup of the element's Jacobian matrix 
  by finite differencing
  is deemed to be too inefficient -- however, as stated above, we very much
  doubt that its analytical computation would be faster). 
  An example is given by the \c SolidEquations class in which
  only the function \c SolidFiniteElement::get_residual(...)
  is implemented -- the element's Jacobian matrix for 
  \c SolidFiniteElements that are derived from these equations
  is computed automatically by the finite differencing operations
  implemented in \c SolidFiniteElement::get_jacobian(...)
  \n \n
  \c SolidFiniteElements that involve non-positional degrees
  of freedom can still call this function to compute the
  entries that correspond to derivatives with respect to
  the nodal positions. Any other entries (such as the derivatives
  of the residuals with respect to the solid pressure degrees
  of freedom) need to be computed separately. The
  \c SolidEquationsWithPressure class provides an example for
  this approach.


<CENTER><B> The SolidMesh class: </B></CENTER>
Finally, we introduce the \c SolidMesh class as a generalisation
of the \c Mesh class. The key features of this class are:
- It overloads the \c Mesh::node_pt(...) function with
  \code
  SolidMesh::node_pt(...)
  \endcode
  which returns a pointer to an \c SolidNode, rather than a
  "normal" \c Node. Equivalent access functions are implemented for
  all other \c Mesh member functions that return pointers to
  \c Nodes.
- It overloads the \c Mesh::assign_local_eqn_numbers()
  function with
  \code
  SolidMesh::assign_local_eqn_numbers(...)
  \endcode
  This function calls the elements' "normal" 
  \c GeneralisedElement::assign_local_eqn_numbers() function,
  followed by a call to  
  \c SolidFiniteElement::assign_solid_local_eqn_numbers(),
  which implements the local equation numbering for the 
  solid mechanics degrees of freedom. 
- Finally, we provide the function
  \code
  SolidMesh::set_lagrangian_nodal_coordinates()
  \endcode
  which assigns the current Eulerian coordinates of all \c Nodes
  to their Lagrangian coordinates, thus turning the current
  configuration into the reference configuration). This function
  greatly facilitates the construction of \c SolidMeshes via
  inheritance from existing \c Meshes. If, for instance, \c SomeMesh
  is an existing, fully functional \c Mesh, the corresponding
  \c SolidMesh can be constructed with a few lines
  of code, as in this example:
  \code
  //=================================================================
  // SolidMesh version of SomeMesh
  //=================================================================
  template<class ELEMENT>
  class SomeSolidMesh : public virtual SomeMesh<ELEMENT>,
                          public virtual SolidMesh
   {
     public:

     // Constructor: Call the constructor to the underlying Mesh
     // then assign the Lagrangian coordinates -- for the
     // PARANOID user, check that the element specified
     // in the template argument is derived from an 
     // SolidFiniteElement
     SomeSolidMesh() : SomeMesh()
      {
  #ifdef PARANOID
     // Check that the element type is derived from
     // the SolidFiniteElement class
     SolidFiniteElement* el_pt=dynamic_cast<SolidFiniteElement*>
      (finite_element_pt(0));
     if (el_pt==0)
      {
       cout << "Element needs to be derived from SolidFiniteElement " << endl;
       assert(false);
      }
  #endif
  
     // Make the current configuration the undeformed one by
     // setting the nodal Lagrangian coordinates to their current
     // Eulerian ones
     set_lagrangian_nodal_coordinates();
    }
 };
\endcode



\subsubsection Solid_IC Timestepping and the generation of initial conditions
In time-dependent problems, the boundary value problem defined by
the variational principle (28) needs to be
augmented by suitable initial conditions which specify the
state of the system at time \f$ t=t_0. \f$ The initial
conditions specify the initial shape of the solid body,
\f[
{\bf R}(\xi^i,t=t_0) =  {\bf R}^{(IC)}(\xi^i),
\ \ \ \ \ \ \ \ \ \ \ \ \ (38)
\f]
and its initial velocity,
\f[
\left. \frac{\partial {\bf R}(\xi^i,t)}{\partial t} \right|_{t=t_0} = 
{\bf V}^{(IC)}(\xi^i),
\ \ \ \ \ \ \ \ \ \ \ \ (39)
\f]
where \f$ {\bf R}^{(IC)}(\xi^i) \f$ and \f$ {\bf V}^{(IC)}(\xi^i) \f$
are given. The accelerations \f$ \partial^2 {\bf R}/\partial t^2
\f$ at \f$t=t_0\f$ follow from the solution of
(28) and can therefore not be enforced,
unless we wish to initialise the time-stepping procedure with a known
exact solution \f$ {\bf R}^{(exact)}(\xi^i,t) \f$. (Only!) in this case
are we allowed to assign an initial value for the acceleration via
\f[
\left. \frac{\partial^2 {\bf R}(\xi^i,t)}{\partial t^2} \right|_{t=t_0} = 
{\bf A}^{(IC)}(\xi^i) \equiv
\left. \frac{\partial^2 {\bf R}^{(exact)}(\xi^i,t)}{\partial t^2} 
\right|_{t=t_0}.
\ \ \ \ \ \ \ \ \ \ \ \ (40)
\f]
We will assume that time-stepping is performed with the \c Newmark 
method which is our default timestepper for hyperbolic
problems. In this case the time-derivatives of the nodal positions in
(28) are replaced by an approximation
which involves the current and three "history values" of the nodal 
positions. To start the time-integration, we need to assign 
suitable values to these quantities to ensure that the intial 
state of the system is represented correctly.

To assign the intial values for the nodal positions, we
(temporarily) remove all boundary conditions for the nodal 
positions and determine their initial values by solving equation
(38) in its weak form,
\f[
f_{il}^{(0)} = \int \left( \sum_{j=1}^N X_{ij}^{(0)} \psi_j(\xi^k) 
- R^{(IC)}_i(\xi^k)
\right) \psi_l(\xi^k) \ dv = 0
\mbox{ \ \ \ \ for \ \ \ $l=1,...,N; \ \ i=1,..,3$}
\ \ \ \ \ \ \ \ \ \ \ \ \ (41)
\f]
where  \f$ R^{(IC)}_i \f$ is the \f$i\f$-th component of 
\f${\bf R}^{(IC)}.\f$ Equation (41) provides
\f$ 3 \times N\f$ equations for the \f$ 3 \times N \f$ components
of the initial nodal positions, \f$ X_{ij}^{(0)} \f$ (where
\f$i=1,...3; \ j=1,...,N)\f$. To determine the initial nodal
velocities, we repeat the same procedure with the prescribed
velocities and solve
\f[
f_{il}^{(1)} = \int \left( \sum_{j=1}^N X_{ij}^{(1)} \psi_j(\xi^k) 
- V^{(IC)}_i(\xi^k)
\right) \psi_l(\xi^k) \ dv = 0
\mbox{ \ \ \ \ for \ \ \ $l=1,...,N;  \ \ i=1,..,3$}
\ \ \ \ \ \ \ \ \ \ \ \ \ (42)
\f]
for the initial nodal velocities, \f$ X_{ij}^{(1)} \f$ (where
\f$i=1,...3; \ j=1,...,N)\f$. Finally, assuminig that we have 
an exact solution for the accelerations, we solve
\f[
f_{il}^{(2)} = \int \left( \sum_{j=1}^N X_{ij}^{(2)} \psi_j(\xi^k) 
- A^{(IC)}_i(\xi^k)
\right) \psi_l(\xi^k) \ dv = 0
\mbox{ \ \ \ \ for \ \ \ $l=1,...,N; \ \ i=1,..,3$}
\ \ \ \ \ \ \ \ \ \ \ \ \ (43)
\f]
for the initial nodal accelerations, \f$ X_{ij}^{(2)} \f$ (where
\f$i=1,...3; \ j=1,...,N)\f$. Having determined the 
nodal positions and their first and second time-derivatives
at \f$ t=t_0 \f$, we can use the functions
\c Newmark::assign_initial_data_values_stage1(...) and
\c Newmark::assign_initial_data_values_stage2(...) to compute the 
positional history values which ensure that the Newmark 
approximations for the
initial velocity and acceleration are correct.
This procedure is fully implemented in the function
\code
SolidMesh::Solid_IC_problem.set_newmark_initial_condition_directly(...)
\endcode
whose arguments are:
- The pointer to the problem being solved. This is needed because
  the solution of equations (41) -
  (43) requires a temporary change to the
  boundary conditions and to the equation numbering scheme. Once
  the history values are assigned, the original
  boundary conditions are restored and the original equation numbers 
  are re-generated by executing
  the Problem's \c Problem::assign_eqn_numbers().
- A pointer to the \c SolidMesh on which the initial conditions
  are assigned.
- A pointer to the \c TimeStepper (which has to be of the \c Newmark
  family).
- A pointer to the "Elasic initial condition" object (discussed below).
- The initial timestep. 

Here is a brief outline of the implementation:
All \c SolidFiniteElements store a pointer to an
\c SolidInitialCondition object (hierher: change Cond to
Condition in code). By default this pointer is set to NULL,
indicating that \c FiniteElement::get_residual(...)
should compute the residuals of the "normal" governing equations.
\c SolidFiniteElements whose initial conditions are to be set with 
the above function, need re-direct the computation of the residual to
\code 
SolidFiniteElement::get_residuals_for_ic(...)
\endcode
\n
whenever the pointer to the \c SolidInitialCondition is non-NULL, as 
illustrated in this code fragment:
\n \n
\code
//=======================================================================
// Compute the residuals for the elasticity equations.
//=======================================================================
template <unsigned DIM>
 void SolidEquations<DIM>::get_residuals(Vector<double> &residuals)
 {
 
  // Simply set up initial condition?
  if (Solid_ic_pt!=0)
   {
    get_residuals_for_ic(residuals);
    return;
   }

  // Set up residuals for principle of virtual displacements

  [...]


  }
 \endcode

The \c SolidInitialCondition object stores a (pointer to a) \c GeomObject
and a flag that indicates which time-derivative (0th, 1st
or 2nd) of the \c GeomObject's position vector is to be assigned to
the nodal coordinates. Based on the value of this flag,
the function \c SolidFiniteElement::get_residuals_for_ic(...),
is able to compute the residuals corresponding to equations
(41), (42) or 
(43). 

This all sounds very complicated (and it is!) but luckily
all the hard work has already been done and the relevant
procedures are fully implemented. Hence,
the actual assignment of the initial conditions is as simple as
this:

\code

// Somebody has created a Problem somewhere
Problem* problem_pt = new SomeProblem;

// Somebody has created an SolidMesh somewhere
SolidMesh* solid_mesh_pt = new SomeSolidMesh;

// Somebody has created a Newmark timestepper somewhere
TimeStepper* time_stepper_pt=new Newmark<1>; 

[...]

// Create the GeomObject whose time-dependent deformation 
// specifies the initial conditions for our solid body:
// The position vector and its time-derivatives
// are available from GeomObject::dpositiondt(...).
GeomObject* geom_obj_pt=new SomeGeomObject;

//Setup object that specifies the initial conditions:
SolidInitialCondition* ic_pt = new SolidInitialCondition(geom_obj_pt);

// Choose the initial timestep
double dt=0.01; 

// Assign the initial conditions:
SolidMesh::Solid_IC_problem.set_newmark_initial_condition_directly(
             problem_pt,solid_mesh_pt,time_stepper_pt,IC_pt,dt);

// Done!

[...]

 \endcode

If we do not know an exact solution to our problem (and in most
cases we obviously won't...), we can only use the procedure described above
to determine the initial nodal positions and velocities. 
We then need to solve the equations (37) for the
remaining "history value". Since the equations (37)
are linear in the accelerations, this is a linear problem
whose Jacobian matrix is proportional to the mass matrix
\f[
M_{ij} = \int {\cal M} \ \psi_i \ \psi_j \ dv. 
\ \ \ \ \ \ \ (44)
\f]
The procedure which determines the initial "history values" 
from the given initial positions and velocities while ensuring
consistency with the governing equation at \f$ t = t_0\f$
is implemented in 
\c SolidMesh::Solid_IC_problem.set_newmark_initial_condition_consistently(...)
which takes the same arguments as the function that assigns 
the acceleration directly, but also requires a function pointer
to the "multiplier" \f$  {\cal M} \f$.  If there
is no growth, the multiplier is given by the timescale ratio
\f$\Lambda^2\f$; if the body is subjected to uniform isotropic
growth, \f$\Gamma \ne 1\f$, the multiplier is equal to 
\f$\Gamma \Lambda^2\f$. 
If the wrong multiplier is specified (or if it is omitted, in
which the default value of 1.0 is used) the residuals (37)
will be nonzero (or at least larger than the tolerance specified in \c
SolidICProblem). In this case a warning is issued and the
code execution terminates. This behaviour can be suppressed
by increasing the tolerance suitably, but you do this at your own risk!


<B>Important:</B> 
The above procedures can only handle the assignment of 
initial conditions in problems that are formulated in
terms of displacements and do \b not involve any additional variables
such as solid pressures. We do not believe that
it is possible to implement the assignment of initial conditions
for such problems without additional knowledge about the precise 
form of the constitutive equations. If the presence of "internal
solid Data" is detected during the assignment of the initial
guess, code execution terminates with an appropriate warning message.




\subsubsection solid_traction Solid traction elements
 write after Poisson flux elements...

\subsubsection beam_and_shell_elements Beam and shell elements














Here is another table that contrasts the various objects with their
procedural representation: 

<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- "The Problem" and "The Mesh".
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- The  Problem stores a pointer to a  Mesh object which is accessible
  via a pointer-based access function,
  \code
  Mesh* mesh_pt=Problem::mesh_pt();
  \endcode
- The  Mesh provides pointer-based access
  to its constituent  Elements and  Nodes,
  \code
  unsigned j;
  Node* node_pt=Mesh::node_pt(j);
  \endcode
  and
  \code
  unsigned e;
  FiniteElement* element_pt=Mesh::finite_element_pt(e);
  \endcode
  or
  \code
  unsigned e;
  GeneralisedElement* element_pt=Mesh::element_pt(e);
  \endcode
- The Mesh provides pointer-based access
  to Nodes that are located on its boundaries:
  \code
  unsigned i_boundary;
  unsigned j;
  Node* boundary_node_pt=Mesh::boundary_node_pt(i_boundary,j);
  \endcode
</TD>
</TR>
</TABLE>


<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
  - Global node \f$j\f$
</TD>
<TD WIDTH="50%" VALIGN="TOP">
  - Global node numbers identify Nodes within the Mesh
    (Nodes themselves do not have numbers!). Access to the 
    Nodes is via pointers:
    \code
    unsigned j;
    Node* node_pt=Problem::mesh_pt()->node_pt(j);
    \endcode
</TD>
</TABLE>


<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- Element \f$e\f$
</TD>
<TD WIDTH="50%" VALIGN="TOP">
-  Element numbers identify  Elements within the  Mesh
  ( Elements themselves do not have numbers!). Access to the 
   Elements is via pointers:
\code
unsigned e;
GeneralisedElement el_pt=Problem::mesh_pt()->element_pt(e);
\endcode
or
\code
unsigned e;
FiniteElement el_pt=Problem::mesh_pt()->finite_element_pt(e);
\endcode
</TD>
</TR>
</TABLE>


<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- The \f$ i\f$-th coordinate of node \f$j\f$
  \f[ X_{ij} \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- The nodal coordinates are member data of the  Node object and are
  accessed via member functions
\code
double x=Node::x(0);
double y=Node::x(1);
\endcode
</TD>
</TR>
</TABLE>


<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- The total number of nodes
  \f[ N \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- The total number of  Nodes in the  Mesh are 
  accessible via a  Mesh member function
\code
unsigned n_nodes=Mesh::nnode();
\endcode
</TD>
</TR>
</TABLE>


<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- The total number of elements
  \f[ N_E \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- The total number of  Elements in the  Mesh are 
  accessible via a  Mesh member function
\code
unsigned n_elem=Mesh::nelement();
\endcode
</TD>
</TR>
</TABLE>


<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- The total number of unknowns in the problem
  \f[ M \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- The total number of unknowns (degrees of freedom)
  in the  Problem are accessible accessible via a
   Problem member function
\code
unsigned n_dof=Problem::ndof();
\endcode
</TD>
</TR>
</TABLE>


<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- The total number of nodes in an element
  \f[ n \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- The total number of  Nodes in an  Element are
  accessible via a  FiniteElement member function
  \code
  unsigned n_node=FiniteElement::nnode();
  \endcode
</TD>
</TR>
</TABLE>


<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- Global node number \f$i_{global}\f$ of local node 
  \f$i_{local}\f$ in element \f$ e \f$
  \f[
  i_{global}= {\cal J}(i_{local},e)
  \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- \c FiniteElements provide pointer-based access to their  Nodes
   \code
   unsigned i_local;
   Node* node_pt=FiniteElement::node_pt(i_local);
   \endcode
</TD>
</TR>
</TABLE>


<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- The equation number of node \f$j\f$
  \f[ {\cal E}(j) \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- The global equation number of a  Node is
  accessible via a  Node member function
\code
unsigned i_eqn=Node::eqn_number(0);
\endcode
[See \ref systems_of_PDEs for an explanation of the argument
 to this function.]
</TD>
</TR>
</TABLE>


<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- The fact that node \f$ j \f$ is pinned  
  is indicated by
  \f[ {\cal E}(j) = -1 \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- The fact that a  Node is pinned  
  is indicated by
  \code
  bool pinned=Node::is_pinned(0);
  \endcode
  or, equivalently,
  \code
  bool pinned=(Node::eqn_number(0)<0);
  \endcode
  [See \ref systems_of_PDEs for an explanation of the argument
   to these functions.]
</TD>
</TR>
</TABLE>


<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- Nodal value at node \f$ j \f$
  \f[ U_j \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- Access to nodal values is via pointers
  \code
  double u=*(Node::value_pt(0));
  \endcode
  [See \ref systems_of_PDEs for an explanation of the argument
   to this function.]
</TD>
</TR>
</TABLE>


<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- The global equation number \f$ i_{global}\f$ 
  of the local unknown \f$i_{local}\f$
  in element  \f$e\f$
  \f[  i_{global} = \widehat{\cal E}(i_{local},e) \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- The global equation number of local unknown \f$ i_{local} \f$ 
  in an  Element is
  accessible via a \c GeneralisedElement member function
\code
unsigned i_local;
unsigned i_global=GeneralisedElement::eqn_number(i_local);
\endcode
</TD>
</TR>
</TABLE>


<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- All shape functions in element \f$ e \f$
  \f[ \psi_i(s_1,s_2) \ \ (i=1,...,n) \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- Shape functions are 
  implemented in concrete  Elements,
  \code
  Vector<double> s(2); // 2 local coordinates
  Shape psi(9);        // 9 shape functions
  QPoissonElement<3,2>::shape(s,psi);
  \endcode
  using the standard interfaces defined in the \c FiniteElement
  base class.
</TD>
</TR>
</TABLE>


<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- Derivatives of all shape functions in element \f$ e \f$
  with respect to the local coordinates \f$ (s_1,s_2) \f$, 
  \f[ \partial \psi_i/\partial s_j \ \ (i=1,...,n; \ j=1,2) \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- Shape functions and their derivatives are implemented in
  concrete  Elements,
  \code
  Vector<double> s(2); // 2 local coordinates
  Shape psi(9);        // 9 shape functions
  DShape dpsi(9,2);    // 9 shape functions have
                       // 2 derivatives each
  QPoissonElement<3,2>::dshape_local(s,psi,dpsi);
  \endcode
  using the standard interfaces defined in the \c FiniteElement
  base class.
</TD>
</TR>
</TABLE>


<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- Derivatives of all shape functions in element \f$ e \f$
  with respect to the global (Eulerian) coordinates \f$ (x_1,x_2) \f$, 
  \f[ \partial \psi_i/\partial x_j \ \ (i=1,...,n; \ j=1,2) \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- Shape functions and their derivatives are implemented in concrete
   Elements,
  \code
  Vector<double> s(2); // 2 local coordinates
  Shape psi(9);        // 9 shape functions
  DShape dpsidx(9,2);  // 9 shape functions have
                       // 2 derivatives each
  QPoissonElement<3,2>::dshape_eulerian(s,psi,dpsidx);
  \endcode
  using the standard interfaces defined in the \c FiniteElement
  base class.
</TD>
</TR>
</TABLE>


<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- Number of degrees of freedom  in element \f$ e \f$
  \f[ {\cal N}_{dof}(e) \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- The number of degrees of freedom in an element are available  from a 
  \c GeneralisedElement member function
  \code
  unsigned ndof=GeneralisedElement::ndof();
  \endcode
</TD>
</TR>
</TABLE>


<TABLE BORDER="1" >
<TR>
<TD WIDTH="50%" VALIGN="TOP">
- The local residual vector and Jacobian matrix in element \f$ e \f$
  \f[ J^{(e)}_{ij} \mbox{ \ \ and \ \ }  r^{(e)}_{i} \ \ \ 
     (i,j=1,...,{\cal N}_{dof}(e)) \f]
</TD>
<TD WIDTH="50%" VALIGN="TOP">
- The computation of the element Jacobian matrix and 
  residual vector are implemented in concrete  Elements,
  using the standard interfaces defined in the \c GeneralisedElement
  base class:
  \code

  // Which element are we dealing with?
  unsigned e;
  QPoissonElement<3,2> * el_pt=Problem::mesh_pt()->element_pt(e);

  // How many dofs are there in the element?
  unsigned ndof=el_pt->ndof();

  // Allocate residual vector and Jacobian matrix
  Vector<double> residual(ndof);
  DenseDoubleMatrix jacobian(ndof,ndof);

  // Get the residual vector and Jacobian matrix from the element
  el_pt->get_jacobian(jacobian,residual);

  \endcode
</TD>
</TR>
</TABLE>
<HR>
<HR>

\endhtmlonly

\latexonly
{\bf
\begin{center}
Sorry -- latex version of this document is switched off. As the title
suggests, it contains only a collection of leftover bits of
documentation, is therefore incomplete, and tends to break
when the pdf version is made.
\end{center}
}
\endlatexonly

<hr>
<hr>
\section pdf PDF file
A <a href="../latex/refman.pdf">pdf version</a> of this document is available.
**/

