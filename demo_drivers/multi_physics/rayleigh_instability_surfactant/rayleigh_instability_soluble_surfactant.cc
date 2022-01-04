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
// Driver for axisymmetric single-layer fluid problem. Plateau-Rayleigh 
// instability (unstable if H>2*pi*R -> forming drops)
// in the presence of a soluble surfactant. The problem is described in
// Numerical analaysis of the Rayleigh instability in capillary tubes:
// The influence of surfactant solubility by
// D. M. Campana & F. A. Saita, Phys. Fluids.
// vol 18, 022104 (2006)
// The parameters are taken from their Table 3

// Generic oomph-lib header
#include "generic.h"

// Axisymmetric Navier-Stokes headers
#include "axisym_navier_stokes.h"

// Interface headers
#include "fluid_interface.h"
//Header for the pesudo-solid mesh updates
#include "solid.h"
#include "constitutive.h"

// The standard rectangular mesh
#include "meshes/rectangular_quadmesh.h"

// Equations for bulk surfactant transport
#include "axisymmetric_advection_navier_stokes_elements.h"

//Use the oomph and std namespaces 
using namespace oomph;
using namespace std;

//==start_of_namespace===================================================
/// Namespace for physical parameters
/// The parameter values are chosen to be those used in Figures 7, 12
/// in Campana & Saita
//=======================================================================
namespace Global_Physical_Variables
{

  //Film thickness parameter
 double Film_Thickness = 0.18;

 /// Reynolds number
 double Re = 1.0;
 
 /// Womersley number
 double ReSt = Re; // (St = 1)
 
 /// Product of Reynolds number and inverse of Froude number
 double ReInvFr = 0.0; // (Fr = 0)

 /// Capillary number
 double Ca = pow(Film_Thickness,3.0);

 /// External pressure
 double P_ext = 0.0;

 /// Direction of gravity
 Vector<double> G(3);

 /// Wavelength of the domain
 double Alpha = 0.8537;
 
 /// Free surface cosine deformation parameter
 double Epsilon = 1.0e-3;

 /// Surface Elasticity number
 double Beta = 0.1;

 /// Alpha for absorption kinetics
 double Alpha_absorption = 1.0;

 /// K parameter that
 /// describes solubility number
 double K = 0.01;

 // Bulk Peclet number
 double Pe = 10000.0;
 // Bulk Peclet Strouhal number
 double PeSt = Pe;

 //We shall adopt the same scaling
 //as Campana & Saita, so divide
 //the bulk equation through by
 //Peclet number to give 1 as the
 //nominal peclet number
 double Pe_reference_scale=1.0;
 
 //Chose our diffusion in the bulk equations
 //to be 1/Pe, after dividing through by Pe
 double Diff = 1.0/Pe;
 
 /// Surface Peclet number
 double Peclet_S = 10000.0;

 /// Sufrace Peclet number multiplied by Strouhal number
 double Peclet_St_S = 1.0; 
 
 /// Pvd file -- a wrapper for all the different
 /// vtu output files plus information about continuous time
 /// to facilitate animations in paraview
 ofstream Pvd_file;

 /// Pseudo-solid Poisson ratio
 double Nu = 0.1;
 
} // End of namespace



//================================================================
/// Interface class to handle the mass transport between bulk
/// and surface as well as the surfactant transport along the
//interface
//================================================================
template<class ELEMENT>
class ElasticAxisymmetricSolubleSurfactantTransportInterfaceElement :
public ElasticAxisymmetricSurfactantTransportInterfaceElement<ELEMENT>
{
private:
 /// Pointer to adsorption number
 double *Alpha_pt;
 /// Pointer to solubility number
 double *K_pt;

 /// Storage for the index at
 /// which the bulk concentration is stored
 unsigned C_bulk_index;

protected:

  /// Get the bulk surfactant concentration
 double interpolated_C_bulk(const Vector<double> &s)
  {
     //Find number of nodes
   unsigned n_node = this->nnode();

   //Get the nodal index at which the unknown is stored
   const unsigned c_index = C_bulk_index;

   //Local shape function
   Shape psi(n_node);

   //Find values of shape function
   this->shape(s,psi);

   //Initialise value of C
   double C = 0.0;

   //Loop over the local nodes and sum
   for(unsigned l=0;l<n_node;l++) 
     {
       C += this->nodal_value(l,c_index)*psi(l);
     }

   return(C);
  }
 
 /// The time derivative of the bulk surface concentration
  double dc_bulk_dt_surface(const unsigned &l) const
  {
   // Get the data's timestepper
    TimeStepper* time_stepper_pt= this->node_pt(l)->time_stepper_pt();
   
   //Initialise dudt
    double dcdt=0.0;
   //Loop over the timesteps, if there is a non Steady timestepper
    if (time_stepper_pt->type()!="Steady")
    {
     //Find the index at which the variable is stored
     const unsigned c_index = C_bulk_index;

     // Number of timsteps (past & present)
     const unsigned n_time = time_stepper_pt->ntstorage();
     
     for(unsigned t=0;t<n_time;t++)
      {
       dcdt += time_stepper_pt->weight(1,t)*this->nodal_value(t,l,c_index);
      }
    }
   return dcdt;
  }

  /// Overload the Helper function to calculate the residuals and 
  /// jacobian entries. This particular function ensures that the
  /// additional entries are calculated inside the integration loop
  void add_additional_residual_contributions_interface(
   Vector<double> &residuals, DenseMatrix<double> &jacobian,
   const unsigned &flag,const Shape &psif, const DShape &dpsifds,
   const DShape &dpsifdS, const DShape &dpsifdS_div,
   const Vector<double> &s,
   const Vector<double> &interpolated_x, const Vector<double> &interpolated_n, 
   const double &W,const double &J)
  {
   //Call the standard transport equations
   ElasticAxisymmetricSurfactantTransportInterfaceElement<ELEMENT>::
    add_additional_residual_contributions_interface(
     residuals, jacobian, flag,psif,dpsifds,dpsifdS,dpsifdS_div,
     s, interpolated_x, interpolated_n, W, J);

   //Add the additional mass transfer terms
   const double k = this->k();
   const double alpha = this->alpha();
   //Find out how many nodes there are
   unsigned n_node = this->nnode();
   
   //Storage for the local equation numbers and unknowns
   int local_eqn = 0, local_unknown = 0;
   
   //Surface advection-diffusion equation
   
   //Find the index at which the bulk concentration is stored
   unsigned c_bulk_index = this->C_bulk_index;
   
   //Now calculate the bulk and surface concentrations at this point
   //Assuming the same shape functions are used (which they are)
   double interpolated_C = 0.0;
   double interpolated_C_bulk = 0.0;
   
   //Loop over the shape functions
   for(unsigned l=0;l<n_node;l++)
    {
     const double psi = psif(l);
     const double C_ = this->nodal_value(l,this->C_index[l]);
     
     interpolated_C += C_*psi;
     interpolated_C_bulk += this->nodal_value(l,c_bulk_index)*psi;
    }

   //Pre compute the flux between the surface and bulk
   double flux = alpha*(interpolated_C_bulk - interpolated_C);
  
  //Now we add the flux term to the appropriate residuals
  for(unsigned l=0;l<n_node;l++)
   {
    //BULK FLUX
    
    //Read out the appropriate local equation
    local_eqn = this->nodal_local_eqn(l,c_bulk_index);
    
    //If not a boundary condition
    if(local_eqn >= 0)
     {
      //Add flux to the bulk
      residuals[local_eqn] -= k*flux*psif(l)*W*J;

      if(flag)
       {
        for(unsigned l2=0;l2<n_node;l2++)
         {
          local_unknown = this->nodal_local_eqn(l2,this->C_index[l2]);
          if(local_unknown >= 0)
           {
            jacobian(local_eqn,local_unknown) += k*alpha*psif(l2)*psif(l)*W*J;
           }
          
          local_unknown = this->nodal_local_eqn(l2,c_bulk_index);
          if(local_unknown >= 0)
           {
            jacobian(local_eqn,local_unknown) -= k*alpha*psif(l2)*psif(l)*W*J;
           }
         }
       }
     } //End of contribution to bulk
    

    //SURFACE FLUX
    //Read out the appropriate local equation
    local_eqn = this->nodal_local_eqn(l,this->C_index[l]);
    
    //If not a boundary condition
    if(local_eqn >= 0)
     {
      //Add flux to the surface
      residuals[local_eqn] -= flux*psif(l)*W*J;

      if(flag)
       {
        for(unsigned l2=0;l2<n_node;l2++)
         {
          local_unknown = this->nodal_local_eqn(l2,this->C_index[l2]);      
          if(local_unknown >= 0)
           {
            jacobian(local_eqn,local_unknown) += alpha*psif(l2)*psif(l)*W*J;
           }
          
          local_unknown = this->nodal_local_eqn(l2,c_bulk_index);
          if(local_unknown >= 0)
           {
            jacobian(local_eqn,local_unknown) -= alpha*psif(l2)*psif(l)*W*J;
           }
         }
       } //End of contribution to surface

     }
   }

  }

public:
 /// Constructor that passes the bulk element and face index down
 /// to the underlying
  ElasticAxisymmetricSolubleSurfactantTransportInterfaceElement(
   FiniteElement* const &element_pt, const int &face_index) : 
   ElasticAxisymmetricSurfactantTransportInterfaceElement<ELEMENT>
   (element_pt,face_index)
  {
   //Initialise the values
   Alpha_pt = &this->Default_Physical_Constant_Value;
   K_pt = &this->Default_Physical_Constant_Value;

   //Hack because I know the bulk index in this case
   C_bulk_index = 3;
  }

 
 /// Return the adsorption number
 double alpha() {return *Alpha_pt;}
 
 /// Return the solubility nubmer
 double k() {return *K_pt;}
 
 /// Access function for pointer to adsorption number
 double* &alpha_pt() {return Alpha_pt;}
 
 /// Access function for pointer to solubility number
 double* &k_pt() {return K_pt;}

 /// Overload the output function
 void output(std::ostream &outfile, const unsigned &n_plot)
  {
   outfile.precision(16);
   
   //Set output Vector
   Vector<double> s(1);
   
   //Tecplot header info 
   outfile << "ZONE I=" << n_plot << std::endl;
   
   const unsigned n_node = this->nnode();
   const unsigned dim = this->dim();
   
   Shape psi(n_node);
   DShape dpsi(n_node,dim);
   
   //Loop over plot points
   for(unsigned l=0;l<n_plot;l++)
    {
	s[0] = -1.0 + l*2.0/(n_plot-1);
        
        this->dshape_local(s,psi,dpsi);
        Vector<double> interpolated_tangent(2,0.0);
        for(unsigned l=0;l<n_node;l++)
         {
          const double dpsi_ = dpsi(l,0);
          for(unsigned i=0;i<2;i++)
           {
            interpolated_tangent[i] += this->nodal_position(l,i)*dpsi_;
           }
         }

        //Tangent
        double t_vel = (this->interpolated_u(s,0)*interpolated_tangent[0] + this->interpolated_u(s,1)*interpolated_tangent[1])/
         sqrt(interpolated_tangent[0]*interpolated_tangent[0] + interpolated_tangent[1]*interpolated_tangent[1]);

        
        //Output the x,y,u,v 
	for(unsigned i=0;i<2;i++) outfile << this->interpolated_x(s,i) << " ";
	for(unsigned i=0;i<2;i++) outfile << this->interpolated_u(s,i) << " ";      
        //Output a dummy pressure
	outfile << 0.0 << " ";
        //Output the concentrations
	outfile << this->interpolated_C(s) << " ";
        outfile << interpolated_C_bulk(s) << " ";
        //Output the interfacial tension
        outfile << this->sigma(s) << " " << t_vel << std::endl;
    }
   outfile << std::endl;
  }
 
};

//==start_of_problem_class===============================================
/// Single axisymmetric fluid interface problem including the
/// transport of an soluble surfactant.
//=======================================================================
template<class ELEMENT, class TIMESTEPPER>
class InterfaceProblem : public Problem
{
 
public:
 
 /// Constructor: Pass the number of elements in radial and axial directions 
 /// and the length of the domain in the z direction)
 InterfaceProblem(const unsigned &n_r, const unsigned &n_z, 
                  const double &l_z);
 
 /// Destructor (empty)
 ~InterfaceProblem() {}

 /// Set initial conditions: Set all nodal velocities to zero and
 /// initialise the previous velocities to correspond to an impulsive
 /// start
 void set_initial_condition()
  {
   // Determine number of nodes in mesh
   const unsigned n_node = Bulk_mesh_pt->nnode();

   // Loop over all nodes in mesh
   for(unsigned n=0;n<n_node;n++)
    {
     // Loop over the three velocity components
     for(unsigned i=0;i<3;i++)
      {
       // Set velocity component i of node n to zero
       Bulk_mesh_pt->node_pt(n)->set_value(i,0.0);
      }
     //Set the bulk concentration to be 1
     Bulk_mesh_pt->node_pt(n)->set_value(3,1.0);
    }

   // Initialise the previous velocity values for timestepping
   // corresponding to an impulsive start
   assign_initial_values_impulsive();

  } // End of set_initial_condition


 /// The global temporal error norm, based on the movement of the nodes
 /// in the radial direction only (because that's the only direction
 /// in which they move!)
 double global_temporal_error_norm()
  {
   //Temp
   double global_error = 0.0;
   
   //Find out how many nodes there are in the problem
   const unsigned n_node = Bulk_mesh_pt->nnode();
   
   //Loop over the nodes and calculate the errors in the positions
   for(unsigned n=0;n<n_node;n++)
    {
     //Set the dimensions to be restricted to the radial direction only
     const unsigned n_dim = 1; 
     //Set the position error to zero
     double node_position_error = 0.0;
     //Loop over the dimensions
     for(unsigned i=0;i<n_dim;i++)
      {
       //Get position error
       double error = 
        Bulk_mesh_pt->node_pt(n)->position_time_stepper_pt()->
        temporal_error_in_position(Bulk_mesh_pt->node_pt(n),i);
     
       //Add the square of the individual error to the position error
       node_position_error += error*error;
      }
     
     //Divide the position error by the number of dimensions
     node_position_error /= n_dim;
     //Now add to the global error
     global_error += node_position_error;
    }
 
   //Now the global error must be divided by the number of nodes
   global_error /= n_node;
   
   //Return the square root of the errr
   return sqrt(global_error);
  }


 /// Access function for the specific mesh
 ElasticRectangularQuadMesh<ELEMENT>*  Bulk_mesh_pt; 

 /// Mesh for the free surface (interface) elements
 Mesh* Interface_mesh_pt;

 /// Pointer to the constitutive law
 ConstitutiveLaw* Constitutive_law_pt;

 /// Node for documentatin
 Node* Document_node_pt;
 
 /// Doc the solution
 void doc_solution(DocInfo &doc_info);
 
 /// Do unsteady run up to maximum time t_max with given timestep dt
 void unsteady_run(const double &t_max, const double &dt); 

 /// Compute the total mass of the insoluble surfactant
 double compute_total_mass()
  {
   //Initialise to zero
   double mass = 0.0;
   
   // Determine number of 1D interface elements in mesh
   const unsigned n_interface_element = Interface_mesh_pt->nelement();
   
   // Loop over the interface elements
   for(unsigned e=0;e<n_interface_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     ElasticAxisymmetricSurfactantTransportInterfaceElement<ELEMENT>* el_pt = 
      dynamic_cast<ElasticAxisymmetricSurfactantTransportInterfaceElement<ELEMENT>*>
      (Interface_mesh_pt->element_pt(e));
     //Add contribution from each element
     mass += el_pt->integrate_c();
    }
   return mass;
  } // End of compute_total_mass
    

private:

 /// Deform the mesh/free surface to a prescribed function
 void deform_free_surface(const double &epsilon)
  {
   //Loop over all nodes in the mesh
   const unsigned n_node = Bulk_mesh_pt->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     //Find out the r value
     double r = Bulk_mesh_pt->node_pt(n)->x(0);
     // Determine z coordinate of spine
     double z_value = Bulk_mesh_pt->node_pt(n)->x(1);
     
     //Set the thickess of the filme
     double thickness = 
      Global_Physical_Variables::Film_Thickness*(1.0 +  
                                                 + epsilon*cos(Global_Physical_Variables::Alpha*z_value));

     //Now scale the position accordingly
     Bulk_mesh_pt->node_pt(n)->x(0) = 1.0 - (1.0-r)*thickness;
    }

   //Reset the lagrangian coordinates
   Bulk_mesh_pt->set_lagrangian_nodal_coordinates();
  } // End of deform_free_surface

  /// Trace file
ofstream Trace_file;
 
}; // End of problem class



//==start_of_constructor==================================================
/// Constructor for single fluid interface problem
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
InterfaceProblem<ELEMENT,TIMESTEPPER>::
InterfaceProblem(const unsigned &n_r, 
                 const unsigned &n_z,
                 const double &l_z) 

{
 // Allocate the timestepper (this constructs the time object as well)
 add_time_stepper_pt(new TIMESTEPPER(true));

 // Build and assign mesh (the "false" boolean flag tells the mesh
 // constructor that the domain is not periodic in r)
 Bulk_mesh_pt = new ElasticRectangularQuadMesh<ELEMENT>(n_r,n_z,1.0,l_z,time_stepper_pt());

 //Create "surface mesh" that will only contain the interface elements
 Interface_mesh_pt = new Mesh;
 {
  // How many bulk elements are adjacent to boundary b?
   // Boundary 1 is the outer boundary
   // The boundaries are labelled
   //                   2
   //                 3   1
   //                   0

  unsigned n_element = Bulk_mesh_pt->nboundary_element(3);
  
  // Loop over the bulk elements adjacent to boundary b?
  for(unsigned e=0;e<n_element;e++)
   {
    // Get pointer to the bulk element that is adjacent to boundary b
    ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
     Bulk_mesh_pt->boundary_element_pt(3,e));
    
    // Find the index of the face of element e along boundary b
    int face_index = Bulk_mesh_pt->face_index_at_boundary(3,e);
    
    // Build the corresponding free surface element
    ElasticAxisymmetricSolubleSurfactantTransportInterfaceElement<ELEMENT>* interface_element_pt = new 
      ElasticAxisymmetricSolubleSurfactantTransportInterfaceElement<ELEMENT>(bulk_elem_pt,face_index);
    
    //Add the prescribed-flux element to the surface mesh
     Interface_mesh_pt->add_element_pt(interface_element_pt);
    
   } //end of loop over bulk elements adjacent to boundary b
 }


 //Use the first node as our document
 Document_node_pt = Bulk_mesh_pt->node_pt(0);
 
 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Interface_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();


 // --------------------------------------------
 // Set the boundary conditions for this problem
 // --------------------------------------------

 //Pin all azimuthal velocities
 //and all vertical positions so the
 //nodes can only move horizontally
 {
  unsigned n_node = this->Bulk_mesh_pt->nnode();
  for(unsigned n=0;n<n_node;++n)
   {
    this->Bulk_mesh_pt->node_pt(n)->pin(2);
    this->Bulk_mesh_pt->node_pt(n)->pin_position(1);
   }
 }
 
 // All nodes are free by default -- just pin the ones that have
 // Dirichlet conditions here
 unsigned ibound = 3;
 unsigned n_node = Bulk_mesh_pt->nboundary_node(ibound);
 for(unsigned n=0;n<n_node;n++)
   {
     Bulk_mesh_pt->boundary_node_pt(ibound,n)->set_value(4,1.0);
   }

 // Determine number of mesh boundaries
 const unsigned n_boundary = Bulk_mesh_pt->nboundary();
 
 // Loop over mesh boundaries
 for(unsigned b=0;b<n_boundary;b++)
  {
   // Determine number of nodes on boundary b
   const unsigned n_node = Bulk_mesh_pt->nboundary_node(b);

   // Loop over nodes on boundary b
   for(unsigned n=0;n<n_node;n++)
    {
     // Pin azimuthal velocity on bounds
     Bulk_mesh_pt->boundary_node_pt(b,n)->pin(2);

     // Pin velocity on the outer wall
     if(b==1)
      {
	Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0);
	Bulk_mesh_pt->boundary_node_pt(b,n)->pin(1);
        Bulk_mesh_pt->boundary_node_pt(b,n)->pin_position(0);
      }
     // Pin axial velocity on bottom and top walls (no penetration)
     if(b==0 || b==2)
      {
       Bulk_mesh_pt->boundary_node_pt(b,n)->pin(1);
      }
    } // End of loop over nodes on boundary b
  } // End of loop over mesh boundaries
 
 // ----------------------------------------------------------------
 // Complete the problem setup to make the elements fully functional
 // ----------------------------------------------------------------
 Constitutive_law_pt = new GeneralisedHookean(&Global_Physical_Variables::Nu);
 
 // Determine number of bulk elements in mesh
 const unsigned n_bulk = Bulk_mesh_pt->nelement();

 // Loop over the bulk elements
 for(unsigned e=0;e<n_bulk;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   // Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   // Set the Womersley number
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;

   // Set the product of the Reynolds number and the inverse of the
   // Froude number
   el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;

   // Set the direction of gravity
   el_pt->g_pt() = &Global_Physical_Variables::G;

   // Set the constitutive law
   el_pt->constitutive_law_pt() = Constitutive_law_pt;
   
   //Set the peclet numbers
   el_pt->pe_pt() = &Global_Physical_Variables::Pe_reference_scale;

   el_pt->pe_st_pt() = &Global_Physical_Variables::Pe_reference_scale;
   //Set the diffusion scale
   el_pt->d_pt() = &Global_Physical_Variables::Diff;
  } // End of loop over bulk elements

 // Create a Data object whose single value stores the external pressure
 Data* external_pressure_data_pt = new Data(1);
 
 // Pin and set the external pressure to some arbitrary value
 double p_ext = Global_Physical_Variables::P_ext;

 external_pressure_data_pt->pin(0);
 external_pressure_data_pt->set_value(0,p_ext);

 // Determine number of 1D interface elements in mesh
 const unsigned n_interface_element = Interface_mesh_pt->nelement();

 // Loop over the interface elements
 for(unsigned e=0;e<n_interface_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ElasticAxisymmetricSolubleSurfactantTransportInterfaceElement<ELEMENT>* el_pt = 
    dynamic_cast<ElasticAxisymmetricSolubleSurfactantTransportInterfaceElement<ELEMENT>*>
    (Interface_mesh_pt->element_pt(e));

   // Set the Capillary number
   el_pt->ca_pt() = &Global_Physical_Variables::Ca;

   // Pass the Data item that contains the single external pressure value
   el_pt->set_external_pressure_data(external_pressure_data_pt);

   // Set the surface elasticity number
   el_pt->beta_pt() = &Global_Physical_Variables::Beta;

   // Set the sorption
   el_pt->alpha_pt() = &Global_Physical_Variables::Alpha_absorption;

   // Set the K parameter
   el_pt->k_pt() = &Global_Physical_Variables::K;

   // Set the surface peclect number
   el_pt->peclet_s_pt() = &Global_Physical_Variables::Peclet_S;

   // Set the surface peclect number multiplied by strouhal number
   el_pt->peclet_strouhal_s_pt() = &Global_Physical_Variables::Peclet_St_S;

  } // End of loop over interface elements

 // Setup equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << std::endl;

} // End of constructor


   
//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::
doc_solution(DocInfo &doc_info)
{ 

 // Output the time
 double t= time_pt()->time();
 cout << "Time is now " << t << std::endl;

 // Document in trace file
 Trace_file << time_pt()->time() << " "
            << Document_node_pt->x(0) 
            << " " << this->compute_total_mass() << std::endl;

 ofstream some_file;
 char filename[100];

 // Set number of plot points (in each coordinate direction)
 const unsigned npts = 5;

 // Open solution output file
// sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
//         doc_info.number());
// some_file.open(filename);

 // Output solution to file
// Bulk_mesh_pt->output(some_file,npts);
// some_file.close();
 //Put interface in separate file
 sprintf(filename,"%s/int%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Interface_mesh_pt->output(some_file,npts);
 some_file.close();

 // Write file as a tecplot text object...
// some_file << "TEXT X=2.5,Y=93.6,F=HELV,HU=POINT,C=BLUE,H=26,T=\"time = " 
 //         << time_pt()->time() << "\"";
 // ...and draw a horizontal line whose length is proportional
 // to the elapsed time
 //some_file << "GEOMETRY X=2.5,Y=98,T=LINE,C=BLUE,LT=0.4" << std::endl;
 //some_file << "1" << std::endl;
 //some_file << "2" << std::endl;
 //some_file << " 0 0" << std::endl;
 //some_file << time_pt()->time()*20.0 << " 0" << std::endl;

 // Close solution output file
// some_file.close();

 // Output solution to file in paraview format
 //sprintf(filename,"%s/soln%i.vtu",doc_info.directory().c_str(),
 //        doc_info.number());
 //some_file.open(filename);
 //Bulk_mesh_pt->output_paraview(some_file,npts);
 //some_file.close();
 
 // Write pvd information 
 //string file_name="soln"+StringConversion::to_string(doc_info.number())
 // +".vtu";
 //ParaviewHelper::write_pvd_information(Global_Physical_Variables::Pvd_file,
 //                                      file_name,t);

} // End of doc_solution

 

//==start_of_unsteady_run=================================================
/// Perform run up to specified time t_max with given timestep dt
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::
unsteady_run(const double &t_max, const double &dt)
{

 // Set value of epsilon
 double epsilon = Global_Physical_Variables::Epsilon;

 // Deform the mesh/free surface
 deform_free_surface(epsilon);

 // Initialise DocInfo object
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");
 
 // Initialise counter for solutions
 doc_info.number()=0;
 
 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);

 // Initialise trace file
 Trace_file << "time" << ", "
            << "edge spine height" << ", "
            << "mass " << ", " << std::endl;

 // Initialise timestep
 initialise_dt(dt);

 // Set initial conditions
 set_initial_condition();

 // Determine number of timesteps
 const unsigned n_timestep = unsigned(t_max/dt);

 // Open pvd file -- a wrapper for all the different
 // vtu output files plus information about continuous time
 // to facilitate animations in paraview
 //sprintf(filename,"%s/soln.pvd",doc_info.directory().c_str());
 //Global_Physical_Variables::Pvd_file.open(filename);
 //ParaviewHelper::write_pvd_header(Global_Physical_Variables::Pvd_file);

 // Doc initial solution
 doc_solution(doc_info);

 // Increment counter for solutions 
 doc_info.number()++;

 //double dt_desired = dt;

 // Timestepping loop
 for(unsigned t=1;t<=n_timestep;t++)
  {
   // Output current timestep to screen
   cout << "\nTimestep " << t << " of " << n_timestep << std::endl;
   
   // Take one fixed timestep
   unsteady_newton_solve(dt);

   //double dt_actual = 
   // adaptive_unsteady_newton_solve(dt_desired,1.0e-6);
   //dt_desired = dt_actual;


   // Doc solution
   doc_solution(doc_info);

   // Increment counter for solutions 
   doc_info.number()++;

   //Reset the lagrangian coordinates
   Bulk_mesh_pt->set_lagrangian_nodal_coordinates();

  } // End of timestepping loop

 // write footer and close pvd file
 //ParaviewHelper::write_pvd_footer(Global_Physical_Variables::Pvd_file);
 //Global_Physical_Variables::Pvd_file.close();

} // End of unsteady_run


/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////


//==start_of_main=========================================================
/// Driver code for single fluid axisymmetric horizontal interface problem 
//========================================================================
int main(int argc, char* argv[]) 
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 /// Maximum time
 double t_max = 1000.0;

 /// Duration of timestep
 const double dt = 0.1;

 // If we are doing validation run, use smaller number of timesteps
 if(CommandLineArgs::Argc>1) 
  { 
   t_max = 0.5; 
  }

 // Number of elements in radial (r) direction
 const unsigned n_r = 10;
   
 // Number of elements in axial (z) direction
 const unsigned n_z = 80;

 // Height of domain
 const double l_z = MathematicalConstants::Pi/Global_Physical_Variables::Alpha;
 
 // Set direction of gravity (vertically downwards)
 Global_Physical_Variables::G[0] = 0.0;
 Global_Physical_Variables::G[1] = 0.0;
 Global_Physical_Variables::G[2] = 0.0;

 // Set up the spine test problem with AxisymmetricQCrouzeixRaviartElements,
 // using the BDF<2> timestepper
 InterfaceProblem<PseudoSolidNodeUpdateElement<AxisymmetricQAdvectionCrouzeixRaviartElement, QPVDElement<2,3> > ,BDF<2> >
  problem(n_r,n_z,l_z);
 
 // Run the unsteady simulation
 problem.unsteady_run(t_max,dt);
} // End of main

