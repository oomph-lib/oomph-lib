//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================



//====start_of_physical_parameters=====================
/// Namespace for phyical parameters
//======================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re=250.0;

 /// Womersley = Reynolds times Strouhal
 double ReSt=250.0;
 
 /// Non-dimensional wall thickness. As in Heil (2004) paper.
 double H=1.0e-2;
 
 /// 2nd Piola Kirchhoff pre-stress. As in Heil (2004) paper.
 double Sigma0=1.0e3;

 /// Pointer to Data object that stores external pressure
 Data* P_ext_data_pt=0;

 /// Min. pressure. Only used in steady runs during parameter
 /// incrementation. Use 1.5 for values of Re to
 /// span the range in Heil (2004) paper. 
 double Pmin=1.5;

 /// Max. pressure. Only used in steady runs during parameter
 /// incrementation. Use 2.0 for Re=250; 3.7 for Re=0 to
 /// span the range in Heil (2004) paper. 
 double Pmax=2.0;

 /// Jump in pressure after a restart --  used to give a steady
 /// solution a kick before starting a time-dependent run
 double P_step=0.0;

 /// Current prescribed vertical position of control point 
 /// (only used for displacement control)
 double Yprescr = 1.0;

 /// Min. of prescribed vertical position of conrol point 
 /// (only used during parameter study with displacement control).
 /// 0.6 corresponds to the value in Heil (2004) paper for static runs.
 double Yprescr_min=0.6; 

 /// Load function: Apply a constant external pressure to the wall.
 /// Note:  This is the load without the fluid contribution!
 /// Fluid load gets added on by FSIWallElement.
 void load(const Vector<double>& xi, const Vector<double>& x,
           const Vector<double>& N, Vector<double>& load)
 { 
  for(unsigned i=0;i<2;i++) 
   {
    load[i] = -P_ext_data_pt->value(0)*N[i];
   }
 } //end of load


 /// Fluid structure interaction parameter: Ratio of stresses used for
 /// non-dimensionalisation of fluid to solid stresses. As in Heil (2004) 
 /// paper 
 double Q=1.0e-2;


} // end of namespace



/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////




//==========start_of_BL_Squash =========================================
/// Namespace to define the mapping [0,1] -> [0,1] that re-distributes
/// nodal points across the channel width.
//======================================================================
namespace BL_Squash
{
 
 /// Boundary layer width
 double Delta=0.1;

 /// Fraction of points in boundary layer
 double Fract_in_BL=0.5; 

 /// Mapping [0,1] -> [0,1] that re-distributes
 /// nodal points across the channel width
 double squash_fct(const double& s)
 {
  // Default return
  double y=s;
  if (s<0.5*Fract_in_BL)
   {
    y=Delta*2.0*s/Fract_in_BL;
   }
  else if (s>1.0-0.5*Fract_in_BL)
   {
    y=2.0*Delta/Fract_in_BL*s+1.0-2.0*Delta/Fract_in_BL;
   }
  else
   {
    y=(1.0-2.0*Delta)/(1.0-Fract_in_BL)*s+
     (Delta-0.5*Fract_in_BL)/(1.0-Fract_in_BL);
   }

  return y;
 }
}// end of BL_Squash






/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////


//====start_of_underformed_wall============================================
/// Undeformed wall is a steady, straight 1D line in 2D space 
///  \f[ x = X_0 + \zeta \f]
///  \f[ y = H \f]
//=========================================================================
class UndeformedWall : public GeomObject
{

public:

 /// Constructor: arguments are the starting point and the height
 /// above y=0.
 UndeformedWall(const double& x0, const double& h): GeomObject(1,2)
  {
   X0=x0;
   H=h;
  }
 

 /// Position vector at Lagrangian coordinate zeta 
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   // Position Vector
   r[0] = zeta[0]+X0;
   r[1] = H;
  }


 /// Parametrised position on object: r(zeta). Evaluated at
 /// previous timestep. t=0: current time; t>0: previous
 /// timestep. Calls steady version.
 void position(const unsigned& t, const Vector<double>& zeta,
               Vector<double>& r) const
  {
   // Use the steady version
   position(zeta,r);

  } // end of position


 /// Posn vector and its  1st & 2nd derivatives
 /// w.r.t. to coordinates:
 /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i). 
 /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ = 
 /// ddrdzeta(alpha,beta,i). Evaluated at current time.
 virtual void d2position(const Vector<double>& zeta,
                         Vector<double>& r,
                         DenseMatrix<double> &drdzeta,
                         RankThreeTensor<double> &ddrdzeta) const
  {
   // Position vector
   r[0] = zeta[0]+X0;
   r[1] = H;

   // Tangent vector
   drdzeta(0,0)=1.0;
   drdzeta(0,1)=0.0;

   // Derivative of tangent vector
   ddrdzeta(0,0,0)=0.0;
   ddrdzeta(0,0,1)=0.0;

  } // end of d2position

private :

 /// x position of the undeformed beam's left end. 
 double X0;

 /// Height of the undeformed wall above y=0.
 double H;

}; //end_of_undeformed_wall






//====Namespace_for_flags================================
/// Namespace for flags
//======================================================
namespace Flags
{

 /// Resolution factor (multiplier for number of elements across channel)
 unsigned Resolution_factor=1;

 /// Use displacement control (1) or not (0)
 unsigned Use_displ_control=1;

 /// Steady run (1) or unsteady run (0)
 unsigned Steady_flag=1;

 /// Number of steps in parameter study
 unsigned Nsteps=5;

 /// String to identify the run type in trace file
 string Run_identifier_string="";

 /// Name of restart file
 string Restart_file_name="";

 /// Doc flags
 void doc_flags()
 {

  std::cout << "\nFlags:\n"
            <<   "======\n"; 

  std::cout << "-- Resolution factor: " << Resolution_factor << std::endl;
   
  if (Steady_flag)
   { 
    std::cout << "-- Steady run " << std::endl;
    if (Use_displ_control)
     {
      std::cout << "-- Using displacement control " << std::endl;
     }
    else
     {
      std::cout << "-- Not using displacement control " << std::endl;
     }
   }
  else
   {
    std::cout << "-- Unsteady run " << std::endl;
    if (Use_displ_control)
     {
      std::cout << "-- Not using displacement control (command line flag\n"
                << "   overwritten because it's an unsteady run!) " 
                << std::endl;
     }
   }

  std::cout << "-- Reynolds number: " 
            << Global_Physical_Variables::Re  << std::endl;

  std::cout << "-- FSI parameter Q: " 
            << Global_Physical_Variables::Q  << std::endl;


  if (Restart_file_name!="")
   {
    std::cout << "-- Performing restart from: " << Restart_file_name 
              << std::endl;
    std::cout << "-- Jump in pressure: " << Global_Physical_Variables::P_step
              << std::endl;
   }
  else
   {
    std::cout << "-- No restart " << std::endl;
   }
  std::cout << std::endl;
 }

}



//====start_of_problem_class==========================================
/// Problem class
//====================================================================
template <class ELEMENT>
class FSICollapsibleChannelProblem : public virtual Problem
{

public :

/// Constructor: The arguments are the number of elements and
/// the lengths of the domain.
 FSICollapsibleChannelProblem(const unsigned& nup, 
                              const unsigned& ncollapsible,
                              const unsigned& ndown,
                              const unsigned& ny,
                              const double& lup,
                              const double& lcollapsible, 
                              const double& ldown,
                              const double& ly,
                              const bool& displ_control,
                              const bool& steady_flag);
 
 /// Destructor
 ~FSICollapsibleChannelProblem()
  { 
   // Mesh gets killed in general problem destructor
  }

 /// Steady run; virtual so it can be overloaded in derived problem
 /// classes
 virtual void steady_run();

 /// Unsteady run; virtual so it can be overloaded in derived problem
 /// classes. Specify timestep or use default of 0.1
 virtual void unsteady_run(const double& dt=0.1);

 /// Access function for the specific bulk (fluid) mesh
 AlgebraicCollapsibleChannelMesh<ELEMENT>* bulk_mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<
    AlgebraicCollapsibleChannelMesh<ELEMENT>*>
    (Bulk_mesh_pt);
  }

 /// Access function for the wall mesh
 OneDLagrangianMesh<FSIHermiteBeamElement>* wall_mesh_pt() 
  {
   return Wall_mesh_pt;

  } // end of access to wall mesh


 /// Actions before solve. Reset counter for number of Newton iterations
 void actions_before_newton_solve()
  {
   Newton_iter=0;
  }

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}


 /// Update before checking Newton convergence: Update the
 /// nodal positions in the fluid mesh in response to possible 
 /// changes in the wall shape.
 void actions_before_newton_convergence_check()
  {
   // Update mesh
   Bulk_mesh_pt->node_update();

   // Increment counter
   Newton_iter++;
  }


 /// Doc the steady solution
 virtual void doc_solution_steady(DocInfo& doc_info, ofstream& trace_file, 
                                  const double& cpu, const unsigned &niter);

 /// Doc the unsteady solution
 virtual void doc_solution_unsteady(DocInfo& doc_info, ofstream& trace_file, 
                                    const double& cpu, const unsigned &niter);

 /// Apply initial conditions
 void set_initial_condition();


  protected:
 
 /// Dump problem to disk to allow for restart.
 void dump_it(ofstream& dump_file);

 /// Read problem for restart from specified restart file.
 void restart(ifstream& restart_file);

 /// Number of elements in the x direction in the upstream part of the channel
 unsigned Nup;

 /// Number of elements in the x direction in the collapsible part of 
 /// the channel
 unsigned Ncollapsible;

 /// Number of elements in the x direction in the downstream part of the channel
 unsigned Ndown;

 /// Number of elements across the channel
 unsigned Ny;

 /// x-length in the upstream part of the channel
 double Lup;

 /// x-length in the collapsible part of the channel
 double Lcollapsible;

 /// x-length in the downstream part of the channel
 double Ldown;

 /// Transverse length
 double Ly;

 /// Pointer to the "bulk" mesh
 AlgebraicCollapsibleChannelMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the mesh that contains the displacement control element
 Mesh* Displ_control_mesh_pt; 
 
 /// Use displacement control?
 bool Displ_control;

 /// Pointer to the "wall" mesh
 OneDLagrangianMesh<FSIHermiteBeamElement>* Wall_mesh_pt; 

 /// Pointer to the left control node
 Node* Left_node_pt;
 
 /// Pointer to right control node
 Node* Right_node_pt;
 
 /// Pointer to control node on the wall
 Node* Wall_node_pt;

 /// Flag for steady run
 bool Steady_flag;

 /// Pointer to GeomObject at which displacement control is applied
 /// (or at which wall displacement is monitored in unsteady runs)
 GeomObject* Ctrl_geom_obj_pt;

 /// Vector of local coordinates of displacement control point
 /// in Ctrl_geom_obj_pt
 Vector<double> S_displ_ctrl;

 /// Pointer to geometric object (one Lagrangian, two Eulerian 
 /// coordinates) that will be built from the wall mesh
 MeshAsGeomObject* Wall_geom_object_pt;

 /// Counter for Newton iterations
 unsigned Newton_iter;

 /// DocInfo object
 DocInfo Doc_info;

};//end of problem class




//=====start_of_constructor======================================
/// Constructor for the collapsible channel problem
//===============================================================
template <class ELEMENT>
FSICollapsibleChannelProblem<ELEMENT>::FSICollapsibleChannelProblem(
 const unsigned& nup, 
 const unsigned& ncollapsible,
 const unsigned& ndown,
 const unsigned& ny,
 const double& lup,
 const double& lcollapsible, 
 const double& ldown,
 const double& ly,
 const bool& displ_control,
 const bool& steady_flag)
{


 // Initialise number of Newton iterations
 Newton_iter=0;

 // Store problem parameters
 Nup=nup;
 Ncollapsible=ncollapsible;
 Ndown=ndown;
 Ny=ny;
 Lup=lup;
 Lcollapsible=lcollapsible;
 Ldown=ldown;
 Ly=ly;
 Steady_flag=steady_flag;

 // Displacement control only makes sense for steady problems
 if (Steady_flag)
  {
   Displ_control=displ_control;
  }
 else
  {
   Displ_control=false;
   if (displ_control)
    {
     std::cout << "Switched off displacement control for unsteady run!"
               << std::endl;
    }
  }
 

 // Overwrite maximum allowed residual to accomodate bad initial guesses
 Problem::Max_residuals=1000.0;

 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. Note: This is appropriate even for
 // the steady problem as we'll explicitly call the *steady*
 // Newton solver which disables the timesteppers
 // during the solve.
 add_time_stepper_pt(new BDF<2>);

 // Create a dummy Steady timestepper that stores two history values
 Steady<2>* wall_time_stepper_pt = new Steady<2>;

 // Add the wall timestepper to the Problem's collection of timesteppers.
 add_time_stepper_pt(wall_time_stepper_pt);

 // Geometric object that represents the undeformed wall: 
 // A straight line at height y=ly; starting at x=lup.
 UndeformedWall* undeformed_wall_pt=new UndeformedWall(lup,ly);

 //Create the "wall" mesh with FSI Hermite elements
 Wall_mesh_pt = new OneDLagrangianMesh<FSIHermiteBeamElement>
  (Ncollapsible,Lcollapsible,undeformed_wall_pt,wall_time_stepper_pt);

 // Build a geometric object (one Lagrangian, two Eulerian coordinates)
 // from the wall mesh
 Wall_geom_object_pt=
  new MeshAsGeomObject(Wall_mesh_pt); 

 // Get pointer to/local coordinate in wall geom object that contains 
 // control node -- adjusted for different values of Q, so that
 // the point is located near the point of strongest collapse.
 Vector<double> zeta_displ_ctrl(1);
 zeta_displ_ctrl[0]=3.5;
 if (std::abs(Global_Physical_Variables::Q-1.0e-3)<1.0e-10)
  {
   zeta_displ_ctrl[0]=3.0;
  }
 //if (std::abs(Global_Physical_Variables::Q-1.0e-4)<1.0e-10) 
 if (Global_Physical_Variables::Q<=1.0e-4) 
  {
   zeta_displ_ctrl[0]=2.5;
  }
 std::cout << "Wall control point located at zeta=" <<zeta_displ_ctrl[0] 
           << std::endl;
 S_displ_ctrl.resize(1);

 // Locate control point (pointer to GeomObject and local coordinate in it)
 Wall_geom_object_pt->locate_zeta(zeta_displ_ctrl,
                                  Ctrl_geom_obj_pt,
                                  S_displ_ctrl);


 // Normal load incrementation or unsteady run
 //===========================================
 Displ_control_mesh_pt=new Mesh;

 // Choose element in which displacement control is applied:
 SolidFiniteElement* controlled_element_pt=
  dynamic_cast<SolidFiniteElement*>(Ctrl_geom_obj_pt);
 
 // Fix the displacement in the vertical (1) direction...
 unsigned controlled_direction=1;
 
 // Pointer to displacement control element
 DisplacementControlElement* displ_control_el_pt;
 
 // Build displacement control element
 displ_control_el_pt=
  new DisplacementControlElement(controlled_element_pt,
                                 S_displ_ctrl,
                                 controlled_direction,
                                 &Global_Physical_Variables::Yprescr);
 
 // The constructor of the  DisplacementControlElement has created
 // a new Data object whose one-and-only value contains the
 // adjustable load: Use this Data object in the load function:
 Global_Physical_Variables::P_ext_data_pt=displ_control_el_pt->
  displacement_control_load_pt();
 
 // Add the displacement-control element to its own mesh
 Displ_control_mesh_pt->add_element_pt(displ_control_el_pt); 
 

 if (!Displ_control) 
  {
   // Create Data object whose one-and-only value contains the
   // (in principle) adjustable load
   Global_Physical_Variables::P_ext_data_pt=new Data(1);
   
   //Pin the external pressure because it isn't actually adjustable.
   Global_Physical_Variables::P_ext_data_pt->pin(0);
  }

 //Build bulk (fluid) mesh
 Bulk_mesh_pt = 
  new AlgebraicCollapsibleChannelMesh<ELEMENT>
  (nup, ncollapsible, ndown, ny,
   lup, lcollapsible, ldown, ly,
   Wall_geom_object_pt,
   time_stepper_pt());


 // Add the sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Wall_mesh_pt);
 add_sub_mesh(Displ_control_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();
   

 // Complete build of fluid mesh
 //----------------------------- 
 
 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 unsigned n_element=Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   // Set the Womersley number
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;

   // Switch off mesh velocity in steady runs
   if (Flags::Steady_flag)
    {
     el_pt->disable_ALE();
    }
   else
    {
     // Is element in rigid part?
     bool is_in_rigid_part=true;
      
     // Number of nodes
     unsigned nnod=el_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       double x=el_pt->node_pt(j)->x(0);
       if ((x>=Lup)&&(x<=(Lup+Lcollapsible)))
        {
         is_in_rigid_part=false;
         break;
        }
      }
     if (is_in_rigid_part)
      {
       el_pt->disable_ALE();
      }
    }

  } // end loop over elements



 // Apply boundary conditions for fluid
 //------------------------------------

 //Pin the velocity on the boundaries
 //x and y-velocities pinned along boundary 0 (bottom boundary) :
 unsigned ibound=0; 
 unsigned num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   for(unsigned i=0;i<2;i++)
    {
     bulk_mesh_pt()->boundary_node_pt(ibound, inod)->pin(i);
    }
  }
  
 //x and y-velocities pinned along boundaries 2, 3, 4 (top boundaries) :
 for(ibound=2;ibound<5;ibound++)
  { 
   num_nod= bulk_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     for(unsigned i=0;i<2;i++)
      {
       bulk_mesh_pt()->boundary_node_pt(ibound, inod)->pin(i);
      }
    }
  }

 //y-velocity pinned along boundary 1 (right boundary):
 ibound=1; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   bulk_mesh_pt()->boundary_node_pt(ibound, inod)->pin(1);
  }


 //Both velocities pinned along boundary 5 (left boundary):
 ibound=5; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   bulk_mesh_pt()->boundary_node_pt(ibound, inod)->pin(0);
   bulk_mesh_pt()->boundary_node_pt(ibound, inod)->pin(1);
  }
//end of pin_velocity


 // Complete build of wall elements
 //--------------------------------
  
 //Loop over the elements to set physical parameters etc.
 n_element = wall_mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast to the specific element type
   FSIHermiteBeamElement *elem_pt = 
    dynamic_cast<FSIHermiteBeamElement*>(wall_mesh_pt()->element_pt(e));
    
   // Set physical parameters for each element:
   elem_pt->sigma0_pt() = &Global_Physical_Variables::Sigma0;
   elem_pt->h_pt() = &Global_Physical_Variables::H;
    
   // Set the load vector for each element
   elem_pt->load_vector_fct_pt() = &Global_Physical_Variables::load;

   // Function that specifies the load ratios
   elem_pt->q_pt() = &Global_Physical_Variables::Q;

   // Set the undeformed shape for each element
   elem_pt->undeformed_beam_pt() = undeformed_wall_pt;

   // The normal on the wall elements as computed by the FSIHermiteElements
   // points away from the fluid rather than into the fluid (as assumed
   // by default)
   elem_pt->set_normal_pointing_out_of_fluid();

   // Displacement control? If so, the load on *all* elements
   // is affected by an unknown -- the external pressure, stored
   // as the one-and-only value in a Data object: Add it to the
   // elements' external Data.
   if (Displ_control)
    {
     //The external pressure is external data for all elements
     elem_pt->add_external_data(Global_Physical_Variables::P_ext_data_pt);
    }


  } // end of loop over elements



 // Boundary conditions for wall mesh
 //----------------------------------

 // Set the boundary conditions: Each end of the beam is fixed in space
 // Loop over the boundaries (ends of the beam)
 for(unsigned b=0;b<2;b++)
  {
   // Pin displacements in both x and y directions
   wall_mesh_pt()->boundary_node_pt(b,0)->pin_position(0); 
   wall_mesh_pt()->boundary_node_pt(b,0)->pin_position(1);
  }
  


 //Choose control nodes
 //---------------------
  
 // Left boundary
 ibound=5; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 unsigned control_nod=num_nod/2;
 Left_node_pt= bulk_mesh_pt()->boundary_node_pt(ibound, control_nod);
  
 // Right boundary
 ibound=1; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 control_nod=num_nod/2;
 Right_node_pt= bulk_mesh_pt()->boundary_node_pt(ibound, control_nod);
  
 
 // Set the pointer to the control node on the wall
 num_nod= wall_mesh_pt()->nnode();
 Wall_node_pt=wall_mesh_pt()->node_pt(Ncollapsible/2);



 // Setup FSI
 //----------

 // The velocity of the fluid nodes on the wall (fluid mesh boundary 3)
 // is set by the wall motion -- hence the no-slip condition must be
 // re-applied whenever a node update is performed for these nodes. 
 // Such tasks may be performed automatically by the auxiliary node update 
 // function specified by a function pointer:
 ibound=3; 
 num_nod= bulk_mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   static_cast<AlgebraicNode*>(
    bulk_mesh_pt()->boundary_node_pt(ibound, inod))->
    set_auxiliary_node_update_fct_pt(
     FSI_functions::apply_no_slip_on_moving_wall);
  }

 // Work out which fluid dofs affect the residuals of the wall elements:
 // We pass the boundary between the fluid and solid meshes and 
 // pointers to the meshes. The interaction boundary is boundary 3 of the 
 // 2D fluid mesh.
 FSI_functions::setup_fluid_load_info_for_solid_elements<ELEMENT,2>
  (this,3,Bulk_mesh_pt,Wall_mesh_pt);
  
 // Setup equation numbering scheme
 cout <<"Total number of equations: " << assign_eqn_numbers() << std::endl; 
 
}//end of constructor



//====start_of_doc_solution_steady============================================
/// Doc the solution for a steady run
//============================================================================
template <class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>:: 
doc_solution_steady(
		    DocInfo &doc_info,
		    ofstream& trace_file,
		    const double& cpu, const unsigned &niter)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output fluid solution 
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 bulk_mesh_pt()->output(some_file,npts);
 some_file.close();

 // Document the wall shape
 sprintf(filename,"%s/beam%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 wall_mesh_pt()->output(some_file,npts);
 some_file.close();
 
// Write restart file
 sprintf(filename,"%s/restart%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 some_file.precision(16);                          
 dump_it(some_file);
 some_file.close();

 // Write trace file 
 trace_file << Global_Physical_Variables::P_ext_data_pt->value(0)  << " ";
 trace_file << Global_Physical_Variables::Yprescr  << " ";

 // Write trace file 
 trace_file << Left_node_pt->value(0) << " "
            << Right_node_pt->value(0) << " "
            << cpu << " "
            << Newton_iter << " " 
            << std::endl; 

 
} // end_of_doc_solution_steady








//====start_of_doc_solution_unsteady==========================================
/// Doc the solution for an unstady run
//============================================================================
template <class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>:: 
doc_solution_unsteady(
		      DocInfo& doc_info,
		      ofstream& trace_file,
		      const double& cpu,
		      const unsigned &niter)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output fluid solution 
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 bulk_mesh_pt()->output(some_file,npts);
 some_file.close();

 // Document the wall shape
 sprintf(filename,"%s/beam%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 wall_mesh_pt()->output(some_file,npts);
 some_file.close();
 
// Write restart file
 sprintf(filename,"%s/restart%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 dump_it(some_file);
 some_file.close();

 // Write trace file 
 trace_file << time_pt()->time() << " ";

 // Get/doc y-coordinate of control point
 Vector<double> r(2);
 Ctrl_geom_obj_pt->position(S_displ_ctrl,r);
 trace_file << r[1]  << " ";

 // Write trace file 
 trace_file << Left_node_pt->value(0) << " "
            << Right_node_pt->value(0) << " "
            << cpu << " "
            << Newton_iter << " " 
            << std::endl; 

 
} // end_of_doc_solution_steady




//=====start_of_dump_it===================================================
/// Dump the solution to disk to allow for restart
//========================================================================
template <class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>::dump_it(ofstream& dump_file)
{

 // Number of submeshes must agree when dumping/restarting so 
 // temporarily add displacement control mesh back in before dumping...
 if (!Displ_control)
  {
   flush_sub_meshes();
   add_sub_mesh(Bulk_mesh_pt);
   add_sub_mesh(Wall_mesh_pt);
   add_sub_mesh(Displ_control_mesh_pt);
   rebuild_global_mesh();
   assign_eqn_numbers();
  }

 // Write current external pressure
 dump_file <<  Global_Physical_Variables::P_ext_data_pt->value(0) 
           << " # external pressure" << std::endl;

 // Call generic dump()
 Problem::dump(dump_file); 

 // ...strip displacement control mesh back out after dumping if
 // we don't actually need it
 if (!Displ_control)
  {
   flush_sub_meshes();
   add_sub_mesh(Bulk_mesh_pt);
   add_sub_mesh(Wall_mesh_pt);
   rebuild_global_mesh();
   assign_eqn_numbers();
  }


} // end of dump_it



//=================start_of_restart=======================================
/// Read solution from disk for restart
//========================================================================
template <class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>::restart(ifstream& restart_file)
{



// Read external pressure
  
// Read line up to termination sign 
 string input_string;
 getline(restart_file,input_string,'#');
 restart_file.ignore(80,'\n');
 
 if (Global_Physical_Variables::P_step!=0.0)
  {
   std::cout 
    << "Increasing external pressure from "
    << double(atof(input_string.c_str())) << " to "
    << double(atof(input_string.c_str()))+Global_Physical_Variables::P_step
    << std::endl;
  }
 else
  {
   std::cout << "Running with unchanged external pressure of " 
             << double(atof(input_string.c_str())) << std::endl;
  }
    
 // Set external pressure
 Global_Physical_Variables::P_ext_data_pt->
  set_value(0,double(atof(input_string.c_str()))+
            Global_Physical_Variables::P_step);

 // Read the generic problem data from restart file
 Problem::read(restart_file); 

 //Now update the position of the nodes to be consistent with
 //the possible precision loss caused by reading in the data from disk
 this->Bulk_mesh_pt->node_update();

 // Strip out displacement control mesh if we don't need it
 if (!Displ_control)
  {
   flush_sub_meshes();
   add_sub_mesh(Bulk_mesh_pt);
   add_sub_mesh(Wall_mesh_pt);
   rebuild_global_mesh();
   assign_eqn_numbers();
  }


} // end of restart



//====start_of_apply_initial_condition========================================
/// Apply initial conditions
//============================================================================
template <class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>::set_initial_condition()
{ 

 // Check that timestepper is from the BDF family
 if (!Steady_flag)
  {
   if (time_stepper_pt()->type()!="BDF")
    {
     std::ostringstream error_stream;
     error_stream << "Timestepper has to be from the BDF family!\n"
                  << "You have specified a timestepper from the "
                  << time_stepper_pt()->type() << " family" << std::endl;
     
     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
    }
  }


 // Pointer to restart file
 ifstream* restart_file_pt=0;
 
 // Restart?
 //---------

 if (Flags::Restart_file_name!="")
  {
   // Open restart file
   restart_file_pt= new ifstream(Flags::Restart_file_name.c_str(),
                                 ios_base::in);
   if (restart_file_pt!=0)
    {
     cout << "Have opened " << Flags::Restart_file_name << 
      " for restart. " << std::endl;
     restart(*restart_file_pt);
     return;
    }
   else
    {
     std::ostringstream error_stream;
     error_stream 
      << "ERROR while trying to open " << Flags::Restart_file_name << 
      " for restart." << std::endl;

     throw OomphLibError(
      error_stream.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);
    }
  }


 // No restart
 else
  {
   // Update the mesh
   bulk_mesh_pt()->node_update();
   
   // Loop over the nodes to set initial guess everywhere
   unsigned num_nod = bulk_mesh_pt()->nnode();
   for (unsigned n=0;n<num_nod;n++)
    {
     // Get nodal coordinates
     Vector<double> x(2);
     x[0]=bulk_mesh_pt()->node_pt(n)->x(0);
     x[1]=bulk_mesh_pt()->node_pt(n)->x(1);
     
     // Assign initial condition: Steady Poiseuille flow
     bulk_mesh_pt()->node_pt(n)->set_value(0,6.0*(x[1]/Ly)*(1.0-(x[1]/Ly)));
     bulk_mesh_pt()->node_pt(n)->set_value(1,0.0);
    } 
   
   // Assign initial values for an impulsive start
   bulk_mesh_pt()->assign_initial_values_impulsive();
  }


} // end of set_initial_condition




//====steady_run==============================================================
/// Steady run
//============================================================================
template <class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>::steady_run()
{ 

 // Set initial value for external pressure (on the wall stiffness scale). 
 // This can be overwritten in set_initial_condition.
 Global_Physical_Variables::P_ext_data_pt->
  set_value(0,Global_Physical_Variables::Pmin);
 
 // Apply initial condition
 set_initial_condition();
 
 //Set output directory
 Doc_info.set_directory("RESLT");
 
 // Open a trace file 
 ofstream trace_file;
 char filename[100];   
 sprintf(filename,"%s/trace.dat",Doc_info.directory().c_str());
 trace_file.open(filename);

 // Write trace file header
 trace_file << "VARIABLES=\"p<sub>ext</sub>\","
            << "\"y<sub>ctrl</sub>\",";
 trace_file << "\"u_1\","
            << "\"u_2\","
            << "\"CPU time for solve\","
            << "\"Number of Newton iterations\","
            << std::endl;

 trace_file << "ZONE T=\"";
 trace_file << "Re=" << Global_Physical_Variables::Re << ", ";
 trace_file << "Q=" << Global_Physical_Variables::Q << ", ";
 trace_file << "resolution factor: " << Flags::Resolution_factor << ". ";
  trace_file << Flags::Run_identifier_string << "\" ";
 trace_file << std::endl;

 // Output the initial solution (dummy for CPU time)
 doc_solution_steady(Doc_info,trace_file,0.0,0);
 
 // Increment step number
 Doc_info.number()++;
 
 
 // Increment for external pressure (on the wall stiffness scale)
 double delta_p=(Global_Physical_Variables::Pmax-
                 Global_Physical_Variables::Pmin)/double(Flags::Nsteps-1);
 
 // Initial and final values for control position
 Global_Physical_Variables::Yprescr=1.0;
 
 // Steady run
 double delta_y=
  (Global_Physical_Variables::Yprescr_min-Global_Physical_Variables::Yprescr)/
  double(Flags::Nsteps-1);
 
 
 // Parameter study
 //----------------
 for (unsigned istep=0;istep<Flags::Nsteps;istep++)
  {
   
   // Displacement control?
   if (Displ_control)
    {
     std::cout << "Solving for control displ = " 
               << Global_Physical_Variables::Yprescr 
               << std::endl;
    }
   else
    {
     std::cout << "Solving for p_ext = " 
               << Global_Physical_Variables::P_ext_data_pt->value(0) 
               << std::endl;
    }
      
   // Solve the problem
   //------------------
   clock_t t_start = clock();
   
   // Explit call to the steady Newton solve.
   steady_newton_solve();

   clock_t t_end= clock();
   double cpu=double(t_end-t_start)/CLOCKS_PER_SEC;
   
      
   // Outpt the solution
   //-------------------
   doc_solution_steady(Doc_info,trace_file,cpu,0);
   
   // Step number
   Doc_info.number()++;
   
   // Adjust control parameter
   if (Displ_control)
    {
     // Increment control position
     Global_Physical_Variables::Yprescr+=delta_y;
    }
   else
    {
     // Increment external pressure 
     double old_p=Global_Physical_Variables::P_ext_data_pt->value(0);
     Global_Physical_Variables::P_ext_data_pt->set_value(0,old_p+delta_p);
    }
   
  }
 
 // Close trace file.
 trace_file.close();

}






//====unsteady_run============================================================
/// Unsteady run. Specify timestep or use default of 0.1.
//============================================================================
template <class ELEMENT>
void FSICollapsibleChannelProblem<ELEMENT>::unsteady_run(const double& dt)
{ 
  
  // Set initial value for external pressure (on the wall stiffness scale). 
  // Will be overwritten by restart data if a restart file (and pressure
  // jump) are specified
  Global_Physical_Variables::P_ext_data_pt->
   set_value(0,Global_Physical_Variables::Pmax);

  // Initialise timestep -- also sets the weights for all timesteppers
  // in the problem.
  initialise_dt(dt);

  std::cout << "Pressure before set initial: " 
            << Global_Physical_Variables::P_ext_data_pt->value(0)
            << std::endl;

  // Apply initial condition
  set_initial_condition();

  std::cout << "Pressure after set initial: " 
            << Global_Physical_Variables::P_ext_data_pt->value(0)
            << std::endl;

  //Set output directory
  Doc_info.set_directory("RESLT");

  // Open a trace file 
  ofstream trace_file;
  char filename[100];   
  sprintf(filename,"%s/trace.dat",Doc_info.directory().c_str());
  trace_file.open(filename);


 // Write trace file header
 trace_file << "VARIABLES=\"time\","
            << "\"y<sub>ctrl</sub>\",";
 trace_file << "\"u_1\","
            << "\"u_2\","
            << "\"CPU time for solve\","
            << "\"Number of Newton iterations\""
            << std::endl;

  trace_file << "ZONE T=\"";
  trace_file << "Re=" << Global_Physical_Variables::Re << ", ";
  trace_file << "Q=" << Global_Physical_Variables::Q << ", ";
  trace_file << "resolution factor: " << Flags::Resolution_factor << ". ";
  trace_file << Flags::Run_identifier_string << "\" ";
  trace_file << std::endl;
  
  // Output the initial solution (dummy for CPU time)
  doc_solution_unsteady(Doc_info,trace_file,0.0,0);

  // Increment step number
  Doc_info.number()++;

  // Timestepping loop
  //------------------
  for (unsigned istep=0;istep<Flags::Nsteps;istep++)
   {
    
    // Solve the problem
    //------------------
    clock_t t_start = clock();
    
    // Explit call to the unsteady Newton solve.
    unsteady_newton_solve(dt);
        
    clock_t t_end= clock();
    double cpu=double(t_end-t_start)/CLOCKS_PER_SEC;
    
 
    // Output the solution
    //--------------------
    doc_solution_unsteady(Doc_info,trace_file,cpu,0);
 
    // Step number
    Doc_info.number()++;
    
   }
  
  // Close trace file.
  trace_file.close();
  
}


