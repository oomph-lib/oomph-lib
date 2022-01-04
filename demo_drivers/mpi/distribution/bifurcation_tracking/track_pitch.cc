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
//Driver function for a simple test Navier--Stokes problem
//Basic C++ input/output libraries


//IMPORTANT NOTE: The convergence of the block-decomposed solver
//degrades under SuperLU_DIST because it uses only static pivoting
//MUMPS appears to be much more stable, or the full-augmented system
//should be used.

#include <iostream>
#include <fstream>
#include <cstdio>

//FSI++ include files

//All the generic stuff -- geometric elements, base classes, etc
#include "generic.h"
//Fluid-based elements (i.e. Navier--Stokes)
#include "axisym_navier_stokes.h"
//A standard library for 2D meshes
#include "meshes/rectangular_quadmesh.h"

//Use the std namespace
using namespace std;

using namespace oomph;

#ifdef OOMPH_HAS_HYPRE
//=============================================================================
/// helper method for the block diagonal F block preconditioner to allow 
/// hypre to be used for as a subsidiary block preconditioner
//=============================================================================
namespace Hypre_Subsidiary_Preconditioner_Helper
{
 Preconditioner* set_hypre_preconditioner()
 {
  return new HyprePreconditioner;
 }
}
#endif



//======================================================================
//A Class for the mesh to be used in the rotating Disc equations
//======================================================================
template <class ELEMENT>
class CylinderMesh : public RectangularQuadMesh<ELEMENT>
{
 //Private data for the mesh
private:
 //Data to hold the number of elements in each section of the mesh
 unsigned Nx1, Nx2, Ny1, Ny2;
 //Data that holds the fractions along the disc in which to introduce
 //the second zones
 double Xfraction, Yfraction;
 //A Vector of pointers to node that line along the centreline of the disc
 Vector<Node *> Centreline_node_pt;

 //Public functions
 public:
 //Constructor. This "builds" the mesh. The arguments are the number
 //of elements in each zone, the radius of the disc and the
 //timestepper that will be used. It may seem odd to pass a timestepper
 //to a mesh. This is because the storage for any previous timesteps
 //must be assigned when the mesh is constructed.
 CylinderMesh(const unsigned &nx1, const unsigned &nx2,
                const unsigned &ny1, const unsigned &ny2,
                const double &lx) :
  RectangularQuadMesh<ELEMENT>(nx1+nx2,ny1+ny2,
                               0.0,lx,0.0,1.0,false,false)
  {
   //Set mesh variables, number of elements in each section
   Nx1 = nx1; Nx2 = nx2; Ny1 = ny1; Ny2 = ny2;

   //Minimum and maximum values of x and fraction of length 
   Xfraction=0.5;
   //RectangularQuadMesh<ELEMENT>::Ymax = 1.0; 
   Yfraction=0.5;

   //Call the generic mesh construction routine
   //(Defined in RectangularQuadMesh<ELEMENT>
   this->build_mesh();
   
   //Set up the Vector Centreline_node_pt containing the nodes on the 
   //channel centreline. Note that this must be called before the 
   //elements are reordered because it relies on the element being in
   //the order defined in build_mesh()
   assign_centreline_nodes();
     
   //Now reorder the elements in the best form for the frontal solver
   //This arranges in elements in vertical slices starting from x=0 and moving
   //along the channel to x=lx
   RectangularQuadMesh<ELEMENT>::element_reorder();
  }

 //A member function that defines the spacing of the nodes in the
 //r (radial) direction
 //In this mesh, there are two different regions in the r-direction
 //with uniform spacing in each. 
 //Note that this can easily be changed to generate non-uniform meshes
 //if desired
 double x_spacing_function(unsigned xelement, unsigned xnode,
                           unsigned yelement, unsigned ynode)
  {
   double Xmin = RectangularQuadMesh<ELEMENT>::Xmin;
   double Xmax = RectangularQuadMesh<ELEMENT>::Xmax;
   unsigned Np = RectangularQuadMesh<ELEMENT>::Np;
   //Set up some spacing parameters
   //Region one starts at Xmin 
   //Region two starts at Xmin + Xfraction(Xmax-Xmin)
   double x1init = Xmin, x2init = Xmin + Xfraction*(Xmax-Xmin);
   //Calculate the spacing between the nodes in each region
   //Assuming uniform spacing
   //Region one has a length Xfraction*(Xmax-Xmin)
   double x1step = Xfraction*(Xmax-Xmin)/((Np-1)*Nx1);
   //Region two has a length (1.0-Xfraction)*(Xmax-Xmin)
   double x2step = (1.0-Xfraction)*(Xmax-Xmin)/((Np-1)*Nx2);
   
   //Now set up the particular spacing 
   //(it's different in the two different regions)
   if(xelement < Nx1) {return (x1init + x1step*((Np-1)*xelement + xnode));}
   else {return (x2init + x2step*((Np-1)*(xelement-Nx1) + xnode));}
  }

 //A member function that defines the spacing of the nodes in the z 
 //(axial) direction.In this mesh there are two different regions
 //in the y-direction, with uniform spacing in each
 //Note that this can easily be changed to generate non-uniform meshes
 //if desired. This may be desirable if sharp boundary, or shear layers,
 //develop. Alternatively we could learn how to use Matthias's adaptive 
 //refinement
 double y_spacing_function(unsigned xelement, unsigned xnode,
                           unsigned yelement, unsigned ynode)
  {
   double Ymin = RectangularQuadMesh<ELEMENT>::Ymin;
   double Ymax = RectangularQuadMesh<ELEMENT>::Ymax;
   unsigned Np = RectangularQuadMesh<ELEMENT>::Np;

   //Set up some spacing parameters
   //The lower region starts at Ymin
   //The upper region starts at Ymin + Yfraction*(Ymax-Ymin)
   double y1init = Ymin, y2init = Ymin + Yfraction*(Ymax-Ymin);
   //Calculate the space between each node in each region,
   //Assumming uniform spacing
   //Region one has a length Yfraction(Ymax-Ymin)
   double y1step = Yfraction*(Ymax-Ymin)/((Np-1)*Ny1);
   //Region two has a length (1.0-Yfraction)*(Ymax-Ymin)
   double y2step = (1.0-Yfraction)*(Ymax-Ymin)/((Np-1)*Ny2);
 
   //Now return the actual node position, it's different in the two
   //regions, of course
   if(yelement < Ny1) {return (y1init + y1step*((Np-1)*yelement + ynode));}
   else {return (y2init + y2step*((Np-1)*(yelement-Ny1) + ynode));}
  }

 //A function to pass all the nodes on the centreline into the 
 //centreline_node_pt Vector
 void assign_centreline_nodes()
  {
   //Find the number of nodes in each element
   unsigned Nnode=this->finite_element_pt(0)->nnode_1d();
   unsigned Nx = RectangularQuadMesh<ELEMENT>::Nx;
   //Loop over all elements on the centreline
   for(unsigned e=0;e<Nx;e++)
    {
     //Loop over all but the last node on each
     for(unsigned l=0;l<(Nnode-1);l++)
      {
       //Add each node to the centreline Vector
       Centreline_node_pt.
        push_back(this->finite_element_pt(Ny1*Nx+e)->node_pt(l));
      }
    }
      //Add the final node to the centreline Vector
   Centreline_node_pt.push_back(this->finite_element_pt(Ny1*Nx + Nx-1)
                                ->node_pt(Nnode-1));
  }

 //Access function for the centreline node Vector
 Node* centreline_node_pt(const unsigned &i) const
  {return Centreline_node_pt[i];}

 //Rescale the mesh to a new length
 double rescale_length(const double &length)
  {
   //Save the present length
   double old_length = this->Xmax;
   //Set the new length
   this->Xmax = length;
   //Find the ratio
   double ratio = length/old_length;

   //Loop over all nodes and scale x position by ratio
   unsigned long n_node = this->nnode();
   for(unsigned n=0;n<n_node;n++) {this->Node_pt[n]->x(0) *= ratio;}
   return ratio;
  }
};


//A Problem class that solves the Navier--Stokes equations
//in an axisymmetric geometry
//N.B. we usually template problems by element-type and a timestepper in
//any problems where different elements or timesteppers MAY be used
template<class ELEMENT>
class RotatingProblem : public Problem
{
private:
 //The Reynolds number will be private member data
 double Re;

 //Length of the domain
 double Length;

 //Position of the centre node
 unsigned Central_node;

public:
 //Constructor:
 //Nr: Number of elements in the r (radial) direction
 //Nz: Number of elements in the z (axial) direction
 RotatingProblem(const unsigned &Nr1, const unsigned &Nr2,
                 const unsigned &Nz1, const unsigned &Nz2);

 //Set boundary conditions on the wall. 
 void set_boundary_conditions();

 //Solve the problem
 void solve_system();

 //Finish full specification of the elements
 void finish_problem_setup();

 //Overload access function for the mesh
 CylinderMesh<ELEMENT>* mesh_pt() 
  {return dynamic_cast<CylinderMesh<ELEMENT>*>(Problem::mesh_pt());}
 
 //--pure virtual functions in Problem class that MUST be overloaded

 //No actions to be taken after each solve step 
 void actions_after_newton_solve() {}

 //No actions to be taken before each solve step
 void actions_before_newton_solve() {} //{change_length(); set_boundary_conditions();}

 //No actions after change in Re
 void actions_after_change_in_bifurcation_parameter() {}

 //No actions to be taken before each solve step
 void actions_before_newton_convergence_check() 
  {change_length(); set_boundary_conditions();
#ifdef OOMPH_HAS_MPI
   this->synchronise_all_dofs();
#endif
  }

 //Function to change the length of the domain
 void change_length() 
  {
   double ratio = mesh_pt()->rescale_length(Length);
   //If the ratio is actually bigger than 1.0
   //if(std::abs(ratio - 1.0) > 1.0e-5)
   {
     //rescale the velocities
     unsigned long n_node = mesh_pt()->nnode();
     for(unsigned n=0;n<n_node;++n)
      {
       //If the node is not pinned then rescale
       for(unsigned i=0;i<3;i+=2)
        {
         if(!mesh_pt()->node_pt(n)->is_pinned(i))
          {
           (*mesh_pt()->node_pt(n)->value_pt(i)) *= ratio;
          }
        }
      }
   }
  }
};

//Constructor 
template<class ELEMENT>
RotatingProblem<ELEMENT>::RotatingProblem
(const unsigned &Nr1, const unsigned &Nr2,
 const unsigned &Nz1, const unsigned &Nz2) :
 Re(0.0), Length(12.0) //Initialise value of Re to zero
{

 //Set the Central Node
 Central_node = 2*Nr1;

 //Now create the mesh, a generic square mesh, the boundaries of the mesh
 //are labelled 0 to 3 starting from the "bottom" and proceeding in an
 //anti-clockwise direction. The parameters of the constructor are passed 
 //directly to the mesh constructor, and the width of the channel is 
 //always 1.0. The mesh is defined in stdmesh.h
 Problem::mesh_pt() = 
  new CylinderMesh<ELEMENT>(Nr1,Nr2,Nz1,Nz2,Length);

 //dynamic_cast<SuperLUSolver*>(this->linear_solver_pt())
 //->set_solver_type(SuperLUSolver::Serial);

//  MumpsSolver* solver_pt=new MumpsSolver;
//  this->linear_solver_pt() = solver_pt;
//  solver_pt->enable_suppress_warning_about_MPI_COMM_WORLD();

//  //Build iterative linear solver
//  /*GMRES<CRDoubleMatrix>* iterative_linear_solver_pt = new
//   GMRES<CRDoubleMatrix>;
 
//   // Set maximum number of iterations
//   iterative_linear_solver_pt->max_iter() = 500;
 
// //  // Set tolerance
//   iterative_linear_solver_pt->tolerance() = 1.0e-8;   
 
//   AxisymmetricNavierStokesLSCPreconditioner* prec_pt = 
//    new AxisymmetricNavierStokesLSCPreconditioner;
//   //Set the mesh
//   prec_pt->set_navier_stokes_mesh(this->mesh_pt());

// //  //Set the preconditioner
//   iterative_linear_solver_pt->preconditioner_pt() = prec_pt;
 
// //  //Set the linear solver
// this->linear_solver_pt() = iterative_linear_solver_pt;*/

/*#ifdef OOMPH_HAS_HYPRE
//Trap because HYPRE can't handle the case when OOMPH_HAS_MPI, but we
//run in serial
#ifndef OOMPH_HAS_MPI

 //Set up the internal preconditioners
 Preconditioner* P_matrix_preconditioner_pt = new HyprePreconditioner;
 
 // Set parameters for use as preconditioner on Poisson-type problem
 Hypre_default_settings::set_defaults_for_2D_poisson_problem(
  static_cast<HyprePreconditioner*>(P_matrix_preconditioner_pt));
 
 // Use Hypre for the Schur complement block
 prec_pt->set_p_preconditioner(P_matrix_preconditioner_pt);
 
 // Shut up!
 static_cast<HyprePreconditioner*>(P_matrix_preconditioner_pt)->
   disable_doc_time();
#endif
#endif
*/ 
//Preconditioner* F_matrix_preconditioner_pt = 
//  new BlockDiagonalPreconditioner<CRDoubleMatrix>;
/*#ifdef OOMPH_HAS_HYPRE
//Trap because HYPRE can't handle the case when 
//OOMPH_HAS_MPI, but we run in serial
#ifndef OOMPH_HAS_MPI
 dynamic_cast<BlockDiagonalPreconditioner<CRDoubleMatrix>* >
  (F_matrix_preconditioner_pt)->set_subsidiary_preconditioner_function
  (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_preconditioner);
#endif
#endif*/
 

/* 
   Preconditioner* F_matrix_preconditioner_pt = new HyprePreconditioner;
 
 // Shut up!
 static_cast<HyprePreconditioner*>(F_matrix_preconditioner_pt)->
  disable_doc_time();
 
 // Set parameters for use as preconditioner in for momentum 
 // block in Navier-Stokes problem
 Hypre_default_settings::set_defaults_for_navier_stokes_momentum_block(
 static_cast<HyprePreconditioner*>(F_matrix_preconditioner_pt));*/
 
 // Use Hypre for momentum block 
 //prec_pt->set_f_preconditioner(F_matrix_preconditioner_pt);

  
  //Set the linear solver (frontal solver) (function defined in Problem class)
  /*linear_solver_pt() = new HSL_MA42;
  

 static_cast<HSL_MA42*>(linear_solver_pt())->lenbuf_factor0() = 3.0;
 static_cast<HSL_MA42*>(linear_solver_pt())->lenbuf_factor1() = 3.0;
 static_cast<HSL_MA42*>(linear_solver_pt())->lenbuf_factor2() = 3.0;
 static_cast<HSL_MA42*>(linear_solver_pt())->front_factor() = 3.0;*/

 //Specify the eigensolver shift
 static_cast<ARPACK*>(eigen_solver_pt())->set_shift(5.0);

 //Complete the build of the elements 
 finish_problem_setup();
 
 //Set the boundary conditions
 //Pin the u- velocity on all boundaries 
 //Pin the v- velocity on all boundaries apart from the outside edge
 //and the w-velocity on all boundaries but the axis and outside edge
 for(unsigned i=0;i<4;i++)
  {
   //Find the number of nodes on the boundary
   unsigned Nboundary_node = mesh_pt()->nboundary_node(i);
   //Loop over the nodes on the boundary
   for(unsigned n=0;n<Nboundary_node;n++)
    {
     //Pin the u velocity on all boundaries
     mesh_pt()->boundary_node_pt(i,n)->pin(0);

     //Pin the v velocity on all cylinder boundaries
     //apart from the outlet
     if(i!=1)
      {
       mesh_pt()->boundary_node_pt(i,n)->pin(2);
      }
     
     //Pin the swirl velocity on the top and bottom
     if((i==0) || (i==2))
      {
        mesh_pt()->boundary_node_pt(i,n)->pin(1);
       }
    }
  }
 
 //Pin a single pressure value which should really be on the centreline
 //The fact that it isn't will introduce a small imperfection into the 
 //system, but it shouldn't matter because pinning the value to zero means
 //that the unbalanced term should also be zero, which will not contribute
 //to the symmetry dot product.
 dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0))->fix_pressure(0,0.0);
 
 //Setup all the equation numbering and look-up schemes 
 //(function defined in Problem class)
 oomph_info << assign_eqn_numbers() << std::endl; 
}

//Finish the specification of the elements in the problem
template<class ELEMENT>
void RotatingProblem<ELEMENT>::finish_problem_setup()
{
 //Loop over all the (fluid) elements 
 unsigned long Nfluid = mesh_pt()->nelement();
 for(unsigned long e=0;e<Nfluid;e++)
  {
   //Cast to the particular element type, this is necessary because
   //the base elements don't have the member functions that we're about
   //to call!
   ELEMENT *temp_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
   
   //Set the Reynolds number for each element 
   //(yes we could have different Reynolds number in each element!!)
   temp_pt->re_pt() = &Re;
   //Need to have non-zero timescale
   temp_pt->re_st_pt() = &Re;
   //The mesh is NOT moving
   temp_pt->disable_ALE();

  }
}

//Set the boundary conditions
template<class ELEMENT>
void RotatingProblem<ELEMENT>::set_boundary_conditions()
{
 //NOTE: The default value of all parameters is zero, so we need only 
 //set the values that are non-zero on the boundaries, i.e. the swirl
 
 //Top boundary
 {
  unsigned i=2;
  //Find the number of nodes on the boundary
  unsigned Nboundary_node = mesh_pt()->nboundary_node(i);
  //Loop over the nodes on the boundary
  for(unsigned n=0;n<Nboundary_node;n++)
   {
    //Get the radial values
    double r = mesh_pt()->boundary_node_pt(i,n)->x(0);
    mesh_pt()->boundary_node_pt(i,n)->set_value(2,-1.0*r);
   }
 }
 
 //Bottom boundary
 {
  unsigned i=0;
  //Find the number of nodes on the boundary
  unsigned Nboundary_node = mesh_pt()->nboundary_node(i);
  //Loop over the nodes on the boundary
  for(unsigned n=0;n<Nboundary_node;n++)
   {
    //Get the radial values
    double r = mesh_pt()->boundary_node_pt(i,n)->x(0);
    mesh_pt()->boundary_node_pt(i,n)->set_value(2,1.0*r);
   }
 }
}

//Solve the system for a number of different values of the Reynolds number
template<class ELEMENT>
void RotatingProblem<ELEMENT>::solve_system()
{
 //Let's calculate the derivatives analytically wrt the Reynolds number
 this->set_analytic_dparameter(&Re);
 this->set_analytic_hessian_products();

 //Define a string that we can set to be the name of the output file
 //char filename[100];

 //Set the boundary conditions (only need to do this once)
 //If the boundary conditions depend up on time, or Re, we need to
 //reset them every time something changes. This is most easily
 //achieved using the actions_before_newton_solve() {} function.
 set_boundary_conditions();
 
 //Reset Reynolds number
 Re = 0.0; 

 //Continuation parameters

 //Control parameter for contiuation
 double *cont_param = &Re;
 //Initial continuation step
 double ds = 0.1;
 reset_arc_length_parameters();
 //Bifurcation_detection = true;
 Desired_newton_iterations_ds = 4;
 Max_residuals = 1000.0;

 //Open trace file
 std::ostringstream trace_filename;
 trace_filename << "trace" << this->communicator_pt()->my_rank() << ".dat";
 ofstream trace(trace_filename.str().c_str());
 ofstream file;
 //Loop over paramters
 for(unsigned l=0;l<2;l++)
  {
   Re += 12.50;
   try
    {
     steady_newton_solve();
    }
   catch (OomphLibError &error)
    {
     error.disable_error_message();
     cout << "Caught solver error -- continuing regardless \n.";
    }

   //Output to the trace file
   if(this->communicator_pt()->my_rank()==0)
    {
     trace << Re << " " << Length << " " 
           << mesh_pt()->centreline_node_pt(0)->value(1) << " "
           << mesh_pt()->centreline_node_pt(Central_node)->value(1) << " "
           << mesh_pt()->centreline_node_pt(Central_node)->x(0)
           << std::endl;
    }
   else
    {
     trace << Re << " " << Length << " " 
           << mesh_pt()->centreline_node_pt(Central_node)->value(1) << " "
           << mesh_pt()->centreline_node_pt(Central_node)->x(0)
           << std::endl;
    }

   //Output data at each step
   //Create the filename, including the array index
   /*sprintf(filename,"step%g_%g.dat",Length,Re);
   //Actually, write the data
   file.open(filename);
   mesh_pt()->output(file,5);
   file.close();*/
  }

 trace.close();

 //Analytic vs FD test
/* {
  unsigned n_dof = ndof();
  LinearAlgebraDistribution dist(this->communicator_pt(),n_dof,false);
  DoubleVector peturb(&dist,0.0);

  peturb[0] = 1.0;
  Vector<DoubleVector> input(1), output(1);;
  input[0].build(&dist,0.0); output[0].build(&dist,0.0);
  input[0][0] = 1.0;

  this->get_hessian_vector_products(peturb,input,output);
  exit(1);
  }*/


 
 //Specify the symmetry the hard way
 unsigned n_dof_local = this->dof_distribution_pt()->nrow_local();
 DoubleVector symm(this->dof_distribution_pt(),0.0);
 Vector<double> backup(n_dof_local);
   
 {
  for(unsigned n=0;n<n_dof_local;n++) {backup[n] = dof(n);}
  
  //Now sort out the problem
  unsigned n_node = mesh_pt()->nnode();
  for(unsigned n=0;n<n_node;n++)
   {
    Node* nod_pt = mesh_pt()->node_pt(n);
    double y = nod_pt->x(1);
    //Loop over the nodal freedoms
    //Anti-symmetric in first
    if(!nod_pt->is_pinned(0)) {nod_pt->set_value(0,0.01*(y-0.5));}
    //Symmetric in second and third
    if(!nod_pt->is_pinned(1)) {nod_pt->set_value(1,0.01*(y-0.5)*(y-0.5));}
    if(!nod_pt->is_pinned(2)) {nod_pt->set_value(2,0.01*(y-0.5)*(y-0.5));}
   }
  
  //Now set all the pressures to zero
  unsigned n_fluid = mesh_pt()->nelement();
  for(unsigned e=0;e<n_fluid;e++)
   {
    ELEMENT* elem_pt = 
     dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
//Ignore halo elements
#ifdef OOMPH_HAS_MPI
    if(!elem_pt->is_halo())
     {
#endif
      unsigned n_pres = 3;
      for(unsigned p=0;p<n_pres;p++)
       {
        if(!elem_pt->internal_data_pt(0)->is_pinned(p))
         {
          elem_pt->internal_data_pt(0)->set_value(p,0.0);
         }
       }
      //Let's set the middle value for an antisymmetric pressure
      double y = elem_pt->node_pt(4)->x(1);
      if(!elem_pt->internal_data_pt(0)->is_pinned(0))
       {
        elem_pt->internal_data_pt(0)->set_value(0,0.01*(y-0.5));
       }
#ifdef OOMPH_HAS_MPI
     }
#endif
   }
 }
 //Now copy the values in
 for(unsigned n=0;n<n_dof_local;n++)
  {
   symm[n] = dof(n);
   dof(n) = backup[n];
   }
   
 //Let's try to find the pitchfork
 activate_pitchfork_tracking(&Re,symm);
 try
  {
   steady_newton_solve();
  }
 catch (OomphLibError &error)
  {
   error.disable_error_message();
   cout << "Caught solver error -- continuing regardless \n.";
  }
 

 
 //Pitchfork found at
 oomph_info << "Pitchfork found at "
           << Re << " " << std::endl;

 if(communicator_pt()->my_rank()==0)
  {
   oomph_info << "The slack parameter is " << dof(2*n_dof_local+1) << std::endl;
  }



 //Now stop and use the first eigenfunction as the appropriate symmetry
 //vector
 Vector<DoubleVector> efn;
 this->assembly_handler_pt()->get_eigenfunction(efn);
 activate_pitchfork_tracking(&Re,efn[0]);
 try
  {
   steady_newton_solve();   
  }
 catch (OomphLibError &error)
  {
   error.disable_error_message();      
   cout << "Caught solver error -- continuing regardless \n.";
  }
 
 oomph_info << "Pitchfork found at "
            << Re << " " << std::endl;
 if(communicator_pt()->my_rank()==0)
  {
   oomph_info << "The slack parameter is " << dof(2*n_dof_local+1) 
              << std::endl;
  }

 
 //Now continue in the length of the domain
 cont_param = &Length;
 ds = 0.01;
 Desired_proportion_of_arc_length=0.5;
 //Parameter_derivative = 0.5;
 //Max_newton_iterations = 6;
 //Desired_newton_iterations_ds=3;
 trace_filename.str("");
 trace_filename << "trace_pitch" << this->communicator_pt()->my_rank() 
                << ".dat";
 trace.open(trace_filename.str().c_str());
 //Output the current value
 if(this->communicator_pt()->my_rank()==0)
  {
   trace << Re << " " << Length << " " 
         << dof(2*n_dof_local+1) << " " 
         << mesh_pt()->centreline_node_pt(0)->value(1) << " "
         << mesh_pt()->centreline_node_pt(Central_node)->value(1) << " "
         << mesh_pt()->centreline_node_pt(Central_node)->x(0)
         << std::endl;
  }
 else
  {
   trace << Re << " " << Length << std::endl;
  }
   

 //Loop over a single continuation step
 for(unsigned l=0;l<1;l++)
  {
   std::cout << "Solving with ds " << ds << std::endl;
   double ds_next = arc_length_step_solve(cont_param,ds);
   ds = ds_next;
   
   //Output the current value
   if(this->communicator_pt()->my_rank()==0)
    {
     trace << Re << " " << Length << " " 
           << dof(2*n_dof_local+1) << " " 
           << mesh_pt()->centreline_node_pt(0)->value(1) << " "
           << mesh_pt()->centreline_node_pt(Central_node)->value(1) << " "
           << mesh_pt()->centreline_node_pt(Central_node)->x(0)
           << std::endl;
    }
   else
    {
     trace << Re << " " << Length << std::endl;
    }
  }
 
 
 //Close the trace file
 trace.close();
}

//Main driver loop
int main(int argc, char *argv[])
{
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif


 //Construct the problem, you can use either element types
 RotatingProblem<AxisymmetricQCrouzeixRaviartElement> problem(20,20,10,10);

 //Set up a dummy partition

#ifdef OOMPH_HAS_MPI
 //Set up a dummy partition
 unsigned n_element = problem.mesh_pt()->nelement();
 Vector<unsigned> element_partition(n_element);
 for(unsigned e=0;e<n_element/2;e++) {element_partition[e]=0;}
 for(unsigned e=n_element/2;e<n_element;e++) {element_partition[e]=1;}

 //DocInfo mesh_doc_info;
 //bool report_stats=true;
 //mesh_doc_info.set_directory("RESLT_MESH");
 problem.distribute(element_partition);//,mesh_doc_info,report_stats);
 //problem.check_halo_schemes(mesh_doc_info);
 //problem.distribute();
#endif


 //Solve the problem :)
 problem.solve_system();


#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif
}








