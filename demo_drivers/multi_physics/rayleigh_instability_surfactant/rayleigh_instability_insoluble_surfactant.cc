//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
// in the presence of an insoluble surfactant. The problem is described in
// A 2-D model of Rayleigh instability in capillary tubes --- surfactant 
// effects by D. Campana, J. Di Paolo & F. A. Saita, Int. J. Multiphase Flow,
// vol 30, pp 431--454, (2004).
// The initial version of this code was developed by Aman Rajvardhan.

// Generic oomph-lib header
#include "generic.h"

// Axisymmetric Navier-Stokes headers
#include "axisym_navier_stokes.h"

// Interface headers
#include "fluid_interface.h"

// The mesh, including horizontal spines
#include "meshes/horizontal_single_layer_spine_mesh.h"

//Use the oomph and std namespaces 
using namespace oomph;
using namespace std;

//==start_of_namespace===================================================
/// Namespace for physical parameters
/// The parameter values are chosen to be those used in Figures 8, 9 
/// in Campana et al.
//=======================================================================
namespace Global_Physical_Variables
{

  //Film thickness parameter
 double Film_Thickness = 0.2;

 /// Reynolds number
 double Re = 40.0;
 
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
 double Alpha = 1.047;
 
 /// Free surface cosine deformation parameter
 double Epsilon = 1.0e-3;

 /// Surface Elasticity number (weak case)
 double Beta = 3.6e-3;

 /// Surface Peclet number
 double Peclet_S = 4032.0;

 /// Sufrace Peclet number multiplied by Strouhal number
 double Peclet_St_S = 1.0; 
 
 /// Pvd file -- a wrapper for all the different
 /// vtu output files plus information about continuous time
 /// to facilitate animations in paraview
 ofstream Pvd_file;

} // End of namespace


 
namespace oomph
{

//==================================================================
///Spine-based Marangoni surface tension elements that add
///a linear dependence on concentration
///of a surface chemical to the surface tension, 
///which decreases with increasing concentration.
///The non-dimensionalisation is the same as Campana et al (2004)
///but we may wish to revisit this.
//=================================================================
template<class ELEMENT>
class SpineAxisymmetricMarangoniSurfactantFluidInterfaceElement : 
  public SpineAxisymmetricFluidInterfaceElement<ELEMENT>
{
private:
 /// Pointer to an Elasticity number
  double *Beta_pt;

 /// Pointer to Surface Peclet number
  double *Peclet_S_pt;

 /// Pointer to the surface Peclect Strouhal number
  double *Peclet_Strouhal_S_pt;

 /// Index at which the surfactant concentration is stored at the
 /// nodes
 unsigned C_index;
 
 /// Default value of the physical constants
 static double Default_Physical_Constant_Value;

protected:

 ///Get the surfactant concentration
 double interpolated_C(const Vector<double> &s)
  {
     //Find number of nodes
   unsigned n_node = this->nnode();

   //Get the nodal index at which the unknown is stored
   const unsigned c_index = C_index;

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
 
 /// The time derivative of the surface concentration
  double dcdt_surface(const unsigned &l) const
  {
   // Get the data's timestepper
    TimeStepper* time_stepper_pt= this->node_pt(l)->time_stepper_pt();
   
   //Initialise dudt
    double dcdt=0.0;
   //Loop over the timesteps, if there is a non Steady timestepper
    if (time_stepper_pt->type()!="Steady")
    {
     //Find the index at which the variable is stored
     const unsigned c_index = C_index;

     // Number of timsteps (past & present)
     const unsigned n_time = time_stepper_pt->ntstorage();
     
     for(unsigned t=0;t<n_time;t++)
      {
       dcdt += time_stepper_pt->weight(1,t)*this->nodal_value(t,l,c_index);
      }
    }
   return dcdt;
  }

 /// The surface tension function is linear in the
 /// concentration with constant of proportionality equal
 /// to the elasticity  number.
 double sigma(const Vector<double> &s)
  {
   //Find the number of shape functions
   const unsigned n_node = this->nnode();
   //Now get the shape fuctions at the local coordinate
   Shape psi(n_node);
   this->shape(s,psi);
   
   //Now interpolate the surfactant concentration
   double C=0.0;
   for(unsigned l=0;l<n_node;l++)
    {
     C += this->nodal_value(l,C_index)*psi(l);
    }
   
   //Get the Elasticity number
   double Beta = this->beta();
   //Return the variable surface tension
   return (1.0 - Beta*(C-1.0));
  } // End of sigma

 /// Fill in the contribution to the residuals
 /// Calculate the contribution to the jacobian
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
    this->fill_in_generic_residual_contribution_interface(residuals,jacobian,1);
    {
    //Use finite differences to handle concentration variations in the
    //bulk equations
    const unsigned n_node = this->nnode();
    //Find the number of dofs in the element
    const unsigned n_dof = this->ndof();
    //Create newres vector
    Vector<double> newres(n_dof);
    
    //Integer storage for local unknown
    int local_unknown=0;
    
    //Use the default finite difference step
    const double fd_step = this->Default_fd_jacobian_step;
    
    //Loop over the nodes again
    for(unsigned n=0;n<n_node;n++)
     {
      //Get the number of values stored at the node
      unsigned c_index = this->C_index;
      
      //Get the local equation number
      local_unknown = this->nodal_local_eqn(n,c_index);
      //If it's not pinned
      if(local_unknown >= 0)
       {
        //Store a pointer to the nodal data value
        double *value_pt = this->node_pt(n)->value_pt(c_index);
        
        //Save the old value of the Nodal data
        double old_var = *value_pt;
        
        //Increment the value of the Nodal data
        *value_pt += fd_step;
        
        //Calculate the new residuals
        this->get_residuals(newres);
        
        //Do finite differences
        for(unsigned m=0;m<n_dof;m++)
         {
          double sum = (newres[m] - residuals[m])/fd_step;
          //Stick the entry into the Jacobian matrix
          jacobian(m,local_unknown) = sum;
         }
        
        //Reset the Nodal data
        *value_pt = old_var;
       }
     }
    }

   //Call the generic routine to handle the spine variables
    SpineElement<FaceGeometry<ELEMENT> >::fill_in_jacobian_from_geometric_data(jacobian);
  }

 
  /// Overload the Helper function to calculate the residuals and 
  /// jacobian entries. This particular function ensures that the
  /// additional entries are calculated inside the integration loop
  void add_additional_residual_contributions_interface(Vector<double> &residuals, DenseMatrix<double> &jacobian,
                                                       const unsigned &flag,const Shape &psif, const DShape &dpsifds,
                                                       const DShape &dpsifdS, const DShape &dpsifdS_div,
                                                       const Vector<double> &s,
                                                       const Vector<double> &interpolated_x, const Vector<double> &interpolated_n, 
                                                       const double &W,const double &J)
  {
   //Flag to control whether the Campana formulation (false) 
   // or our own (true) is used
   bool Integrated_curvature = true;
   
  //Find out how many nodes there are
   unsigned n_node = this->nnode();
   
   //Storage for the local equation numbers and unknowns
   int local_eqn = 0, local_unknown = 0;
   
   //Surface advection-diffusion equation
   
   //Find the index at which the concentration is stored
   unsigned c_index = this->C_index;
   Vector<unsigned> u_index = this->U_index_interface;

   //Read out the surface peclect number
   const double Pe_s = this->peclet_s();
   //const double PeSt_s = this->peclet_strouhal_s();
   
   //Read out the radial position
   const double r = interpolated_x[0];

   //Rescale the jacobian to be the "raw" version
   const double J_raw = J/r;
   
   //Now calculate the concentration at this point
   //Assuming the same shape functions are used (which they are)
   double interpolated_C = 0.0;
   double interpolated_dCds = 0.0;
   double dCdt = 0.0;
   //The tangent vector
   const unsigned ndim = this->node_pt(0)->ndim();
   Vector<double> interpolated_tangent(ndim,0.0);
   Vector<double> interpolated_u(ndim,0.0);
   Vector<double> mesh_velocity(ndim,0.0);
   Vector<double> interpolated_duds(ndim,0.0);
   if(ndim+1 != u_index.size())
    {
      throw OomphLibError("Dimension Incompatibility",
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }
   
   for(unsigned l=0;l<n_node;l++)
     {
       const double psi = psif(l);
       const double dpsi = dpsifds(l,0);
       interpolated_C += this->nodal_value(l,c_index)*psi;
       interpolated_dCds += this->nodal_value(l,c_index)*dpsi;
       dCdt += dcdt_surface(l)*psi;
       for(unsigned i=0;i<ndim;i++)
        {
         interpolated_tangent[i] += this->nodal_position(l,i)*dpsi;
         interpolated_u[i] += this->nodal_value(l,u_index[i])*psi;
         interpolated_duds[i] += this->nodal_value(l,u_index[i])*dpsi;
        }
       //Mesh Velocity
       for(unsigned j=0;j<ndim;j++)
        {
         mesh_velocity[j] += this->dnodal_position_dt(l,j)*psi;
        }
     }
   
   
   double u_tangent = 0.0, t_length = 0.0;
   for(unsigned i=0;i<ndim;i++) 
     {
       u_tangent += interpolated_u[i]*interpolated_tangent[i];
       t_length  += interpolated_tangent[i]*interpolated_tangent[i];
     }

   //Work out the second derivative of position
   Vector<double> d2xds(2,0.0);
   for(unsigned i=0;i<2;i++)
    {
     d2xds[i] = this->nodal_position(0,i) + 
      this->nodal_position(2,i) - 2.0*this->nodal_position(1,i);
    }
   
   //Let's do the first component of the curvature
   double k1 = 0.0;
   for(unsigned i=0;i<2;i++)
    {
     k1 += (d2xds[i]/(J_raw*J_raw) - interpolated_tangent[i]*(
             interpolated_tangent[0]*d2xds[0] 
             + interpolated_tangent[1]*d2xds[1])/(J_raw*J_raw*J_raw*J_raw))*interpolated_n[i];
    }
   
   //Second component of the curvature
   double k2 = - (interpolated_n[0] / r);
   
   //Now we add the flux term to the appropriate residuals
   for(unsigned l=0;l<n_node;l++)
    {
     //Read out the apprporiate local equation
     local_eqn = this->nodal_local_eqn(l,c_index);
     
     //If not a boundary condition
     if(local_eqn >= 0)
      {
       //Time derivative
       residuals[local_eqn] += dCdt*psif(l)*r*W*J_raw;

       double tmp = 0.0;
       //Compute our version in which second derivatives are not required
       if(Integrated_curvature)
        {
         for(unsigned i=0;i<2;i++)
          {
           tmp += 
            (interpolated_u[i] - mesh_velocity[i])*interpolated_tangent[i];
          }
         //First Advection term
         residuals[local_eqn] += tmp*interpolated_dCds*psif(l)*r*W/J_raw;
         //Additional term from axisymmetric formulation
         residuals[local_eqn] += interpolated_C*interpolated_u[0]*psif(l)*W*J_raw;
         //Second Advection term
         residuals[local_eqn] += interpolated_C*
          (interpolated_duds[0]*interpolated_tangent[0]
           + interpolated_duds[1]*interpolated_tangent[1])*r*W*psif(l)/J_raw;
          }
       //This is the Campana formulation
       else
        {
         //ALE term
         for(unsigned i=0;i<2;i++)
          {
           tmp += -mesh_velocity[i]*interpolated_tangent[i];
          }
         residuals[local_eqn] += tmp*interpolated_dCds*psif(l)*r*W/J_raw;
         // Curvature (normal velocity) term
         residuals[local_eqn] -= 
          interpolated_C*(k1 + k2)*psif(l)*r*W*J_raw*(
           interpolated_u[0]*interpolated_n[0]
           + interpolated_u[1]*interpolated_n[1]);
         //Integrated by parts tangential advection term
         residuals[local_eqn] -= interpolated_C*(
          interpolated_tangent[0]*interpolated_u[0] + 
          interpolated_tangent[1]*interpolated_u[1])*dpsifds(l,0)*r*W/J_raw;
        }
       
       //Diffusion term
       residuals[local_eqn] += (1.0/Pe_s)*interpolated_dCds*dpsifds(l,0)*r*W/J_raw;
       
       //We also need to worry about the jacobian terms
       if(flag)
        {
         //Loop over the nodes again
         for(unsigned l2=0;l2<n_node;l2++)
          {
           //Get the time stepper
           TimeStepper* time_stepper_pt=this->node_pt(l2)->time_stepper_pt();
           
           //Get the unknown c_index
           local_unknown =this->nodal_local_eqn(l2,c_index);
           
           if(local_unknown >=0)
            {
             jacobian(local_eqn,local_unknown) += 
              time_stepper_pt->weight(1,0)*psif(l2)*psif(l)*r*W*J_raw;
             
             if(Integrated_curvature)
              {
               jacobian(local_eqn,local_unknown) += ((interpolated_u[0] - mesh_velocity[0])*interpolated_tangent[0]
                                                     + (interpolated_u[1] - mesh_velocity[1])*interpolated_tangent[1])*dpsifds(l2,0)
                *psif(l)*r*W/J_raw;


               jacobian(local_eqn,local_unknown) += psif(l2)*interpolated_u[0]*psif(l)*W*J_raw;
               
               jacobian(local_eqn,local_unknown) += psif(l2)*(interpolated_tangent[0]*interpolated_duds[0]
                                                              + interpolated_tangent[1]*interpolated_duds[1])*psif(l)*r*W/J_raw;
              }
             else
              {
               jacobian(local_eqn,local_unknown) -=
                (mesh_velocity[0]*interpolated_tangent[0] + mesh_velocity[1]*interpolated_tangent[1])*dpsifds(l2,0)*psif(l)*r*W/J_raw;
               
               jacobian(local_eqn,local_unknown) -= psif(l2)*(k1 + k2)*psif(l)*r*W*J_raw*(interpolated_u[0]*interpolated_n[0]
                                                                                      + interpolated_u[1]*interpolated_n[1]);
               
               jacobian(local_eqn,local_unknown) -= psif(l2)*(interpolated_tangent[0]*interpolated_u[0] + interpolated_tangent[1]*
                                                              interpolated_u[1])*dpsifds(l,0)*r*W/J_raw;
              }
             
             jacobian(local_eqn,local_unknown) += (1.0/Pe_s)*dpsifds(l2,0)*dpsifds(l,0)*r*W/J_raw;
             
            }


           //Loop over the velocity components
           for(unsigned i2=0;i2<ndim;i2++)
            {
             
             //Get the unknown
             local_unknown = this->nodal_local_eqn(l2,u_index[i2]);
             
             
             //If not a boundary condition
             if(local_unknown >= 0)
              {
               
               // jacobian(local_eqn,local_unknown) += time_stepper_pt->weight(1,0)*psif(l2)*PeSt_s*psif(l)*r*W*J_raw;
               if(Integrated_curvature)
                {
                 jacobian(local_eqn,local_unknown) += 
                  interpolated_dCds*psif(l2)*interpolated_tangent[i2]*psif(l)*r*W/J_raw;

                 if(i2==0)
                  {
                   jacobian(local_eqn,local_unknown) += interpolated_C*psif(l2)*W*J_raw;
                  }
                }
               else
                {
                 jacobian(local_eqn,local_unknown) -= interpolated_C*(k1 + k2)*psif(l)*r*W*J_raw*psif(l2)*interpolated_n[i2];
                 
                 jacobian(local_eqn,local_unknown) -=  interpolated_C*interpolated_tangent[i2]*psif(l2)*dpsifds(l,0)*r*W/J_raw;
                }
               
               jacobian(local_eqn,local_unknown) += 
                interpolated_C*interpolated_tangent[i2]*dpsifds(l2,0)*psif(l)*r*W/J_raw;
              }
            }
          }
        }
      }
    } //End of loop over the nodes
   
  }  

 
 /// Add the element's contribution to its residuals vector,
 /// jacobian matrix and mass matrix
  void fill_in_contribution_to_jacobian_and_mass_matrix( Vector<double> &residuals, DenseMatrix<double> &jacobian, 
							 DenseMatrix<double> &mass_matrix)
  {
    //Add the contribution to the jacobian
    this->fill_in_contribution_to_jacobian(residuals,jacobian);
    //No mass matrix terms, but should probably do kinematic bit here
  }

public:
 ///Constructor that passes the bulk element and face index down
 ///to the underlying
  SpineAxisymmetricMarangoniSurfactantFluidInterfaceElement(FiniteElement* const &element_pt, const int &face_index) : SpineAxisymmetricFluidInterfaceElement<ELEMENT>
														       (element_pt,face_index)
  {
   //Initialise the values
   Beta_pt = &Default_Physical_Constant_Value;
   Peclet_S_pt = &Default_Physical_Constant_Value;
   Peclet_Strouhal_S_pt = &Default_Physical_Constant_Value;

   //Add the additional surfactant terms to these surface elements
   
   //Read out the number of nodes on the face
   //For some reason I need to specify the this pointer here(!)
   unsigned n_node_face = this->nnode();
   //Set the additional data values in the face
   //There is one additional values at each node --- the lagrange multiplier
   Vector<unsigned> additional_data_values(n_node_face);
   for(unsigned i=0;i<n_node_face;i++) additional_data_values[i] = 1;
   //Resize the data arrays accordingly 
   this->resize_nodes(additional_data_values);

   //The C_index is the new final value
   //Minor HACK HERE
   C_index = this->node_pt(0)->nvalue()-1;
}
 
 ///Return the Elasticity number
  double beta() {return *Beta_pt;}

 ///Return the surface peclect number
 double peclet_s() {return *Peclet_S_pt;}
 
 ///Return the surface peclect strouhal number
 double peclet_strouhal_s() {return *Peclet_Strouhal_S_pt;}
 
 ///Access function for pointer to the Elasticity number
 double* &beta_pt() {return Beta_pt;}
 
 ///Access function for pointer to the surface Peclet number
 double* &peclet_s_pt() {return Peclet_S_pt;}
 
 ///Access function for pointer to the surface Peclet x Strouhal number
 double* &peclet_strouhal_s_pt() {return Peclet_Strouhal_S_pt;}
 
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
        //Output the concentration
	outfile << interpolated_C(s) << " ";
        //Output the interfacial tension
        outfile << sigma(s) << " " << t_vel << std::endl;
      }
 outfile << std::endl;
  }

 
 /// Compute the concentration intergated over the area
 double integrate_c() const
  {
   //Find out how many nodes there are
   unsigned n_node = this->nnode();
   
   //Set up memeory for the shape functions
   Shape psif(n_node);
   DShape dpsifds(n_node,1);

   //Set the value of n_intpt
   unsigned n_intpt = this->integral_pt()->nweight();
   
   //Storage for the local coordinate
   Vector<double> s(1);

   //Storage for the total mass
   double mass = 0.0;
   
   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Get the local coordinate at the integration point
     s[0] = this->integral_pt()->knot(ipt,0);
     
     //Get the integral weight
     double W = this->integral_pt()->weight(ipt);
     
     //Call the derivatives of the shape function
     this->dshape_local_at_knot(ipt,psif,dpsifds);
     
     //Define and zero the tangent Vectors and local velocities
     Vector<double> interpolated_x(2,0.0);
     Vector<double> interpolated_t(2,0.0);
     double interpolated_c = 0.0;
     
     
     //Loop over the shape functions to compute concentration and tangent
     for(unsigned l=0;l<n_node;l++)
      {
       interpolated_c += this->nodal_value(l,C_index)*psif(l);
       //Calculate the tangent vector
       for(unsigned i=0;i<2;i++)
        {
         interpolated_x[i] += this->nodal_position(l,i)*psif(l);
         interpolated_t[i] += this->nodal_position(l,i)*dpsifds(l,0);
        }
      }
     
     //The first positional coordinate is the radial coordinate
     double r = interpolated_x[0];

   //Calculate the length of the tangent Vector
   double tlength = interpolated_t[0]*interpolated_t[0] + 
    interpolated_t[1]*interpolated_t[1];
   
   //Set the Jacobian of the line element
   double J = sqrt(tlength);

   mass += interpolated_c*r*W*J;
    }
   return mass;
  }

};


//Define the default physical value to be one
template<class ELEMENT>
double SpineAxisymmetricMarangoniSurfactantFluidInterfaceElement<ELEMENT>::Default_Physical_Constant_Value = 1.0;

}


namespace oomph
{

//======================================================================
/// Inherit from the standard Horizontal single-layer SpineMesh
/// and modify the spine_node_update() function so that it is appropriate
/// for an annular film, rather than a fluid fibre.
//======================================================================
template <class ELEMENT>
class MyHorizontalSingleLayerSpineMesh : 
  public HorizontalSingleLayerSpineMesh<ELEMENT>
{

public:

 /// Constructor: Pass number of elements in x-direction, number of
 /// elements in y-direction, radial extent, axial length , and pointer 
 /// to timestepper (defaults to Steady timestepper)
 MyHorizontalSingleLayerSpineMesh(const unsigned &nx, 
                                  const unsigned &ny,
                                  const double &lx,
                                  const double &ly,
                                  TimeStepper* time_stepper_pt=
                                  &Mesh::Default_TimeStepper) :
  HorizontalSingleLayerSpineMesh<ELEMENT>(nx,ny,lx,ly,time_stepper_pt) {}

 /// Node update function assumed spines rooted at the wall
 /// fixed to be at r=1 and directed inwards to r=0.
 virtual void spine_node_update(SpineNode* spine_node_pt)
  {
   //Get fraction along the spine
   double W = spine_node_pt->fraction();

   //Get spine height
   double H = spine_node_pt->h();
   
   //Set the value of y
   spine_node_pt->x(0) = 1.0-(1.0-W)*H;
  }

};

}


//==start_of_problem_class===============================================
/// Single axisymmetric fluid interface problem including the
/// transport of an insoluble surfactant.
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

 /// Spine heights/lengths are unknowns in the problem so their values get
 /// corrected during each Newton step. However, changing their value does
 /// not automatically change the nodal positions, so we need to update all
 /// of them here.
 void actions_before_newton_convergence_check()
  {
   Bulk_mesh_pt->node_update();
  }

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
 MyHorizontalSingleLayerSpineMesh<ELEMENT>*  Bulk_mesh_pt; 

 /// Mesh for the free surface (interface) elements
 Mesh* Interface_mesh_pt;

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
     SpineAxisymmetricMarangoniSurfactantFluidInterfaceElement<ELEMENT>* el_pt = 
      dynamic_cast<SpineAxisymmetricMarangoniSurfactantFluidInterfaceElement<ELEMENT>*>
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
   // Determine number of spines in mesh
   const unsigned n_spine = Bulk_mesh_pt->nspine();
   
   // Loop over spines in mesh
   for(unsigned i=0;i<n_spine;i++)
    {
     
     // Determine z coordinate of spine
     double z_value = Bulk_mesh_pt->boundary_node_pt(1,i)->x(1);

     // Set spine height
     Bulk_mesh_pt->spine_pt(i)->height() = 
      Global_Physical_Variables::Film_Thickness*(1.0 +  
                                                 + epsilon*cos(Global_Physical_Variables::Alpha*z_value));

    } // End of loop over spines
   
   // Update nodes in bulk mesh
   Bulk_mesh_pt->node_update();

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
 Bulk_mesh_pt = new MyHorizontalSingleLayerSpineMesh<ELEMENT>(n_r,n_z,1.0,l_z,time_stepper_pt());
 
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
    SpineAxisymmetricMarangoniSurfactantFluidInterfaceElement<ELEMENT>* interface_element_pt = new 
      SpineAxisymmetricMarangoniSurfactantFluidInterfaceElement<ELEMENT>(bulk_elem_pt,face_index);
    
    //Add the prescribed-flux element to the surface mesh
     Interface_mesh_pt->add_element_pt(interface_element_pt);
    
   } //end of loop over bulk elements adjacent to boundary b
 }


 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Interface_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();


 // --------------------------------------------
 // Set the boundary conditions for this problem
 // --------------------------------------------

 //Pin all azimuthal velocities
 {
  unsigned n_node = this->Bulk_mesh_pt->nnode();
  for(unsigned n=0;n<n_node;++n)
   {
    this->Bulk_mesh_pt->node_pt(n)->pin(2);
   }
 }
 
 // All nodes are free by default -- just pin the ones that have
 // Dirichlet conditions here
 unsigned ibound = 3;
 unsigned n_node = Bulk_mesh_pt->nboundary_node(ibound);
 for(unsigned n=0;n<n_node;n++)
   {
     Bulk_mesh_pt->boundary_node_pt(ibound,n)->set_value(3,1.0);
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
   SpineAxisymmetricMarangoniSurfactantFluidInterfaceElement<ELEMENT>* el_pt = 
    dynamic_cast<SpineAxisymmetricMarangoniSurfactantFluidInterfaceElement<ELEMENT>*>
    (Interface_mesh_pt->element_pt(e));

   // Set the Capillary number
   el_pt->ca_pt() = &Global_Physical_Variables::Ca;

   // Pass the Data item that contains the single external pressure value
   el_pt->set_external_pressure_data(external_pressure_data_pt);

   // Set the surface elasticity number
   el_pt->beta_pt() = &Global_Physical_Variables::Beta;

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
            << Bulk_mesh_pt->spine_pt(0)->height()
            << " " << this->compute_total_mass() << std::endl;

 ofstream some_file;
 char filename[100];

 // Set number of plot points (in each coordinate direction)
 const unsigned npts = 5;

 // Open solution output file
 //sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
 //        doc_info.number());
 //some_file.open(filename);

 // Output solution to file
 //Bulk_mesh_pt->output(some_file,npts);
 //some_file.close();
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
 string file_name="soln"+StringConversion::to_string(doc_info.number())
  +".vtu";
 ParaviewHelper::write_pvd_information(Global_Physical_Variables::Pvd_file,
                                       file_name,t);

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
            << "contact angle left" << ", "
            << "contact angle right" << ", " << std::endl;

 // Initialise timestep
 initialise_dt(dt);

 // Set initial conditions
 set_initial_condition();

 // Determine number of timesteps
 const unsigned n_timestep = unsigned(t_max/dt);

 // Open pvd file -- a wrapper for all the different
 // vtu output files plus information about continuous time
 // to facilitate animations in paraview
 sprintf(filename,"%s/soln.pvd",doc_info.directory().c_str());
 Global_Physical_Variables::Pvd_file.open(filename);
 ParaviewHelper::write_pvd_header(Global_Physical_Variables::Pvd_file);

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
  } // End of timestepping loop

 // write footer and close pvd file
 ParaviewHelper::write_pvd_footer(Global_Physical_Variables::Pvd_file);
 Global_Physical_Variables::Pvd_file.close();

} // End of unsteady_run


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


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
 InterfaceProblem<SpineElement<AxisymmetricQCrouzeixRaviartElement >,BDF<2> >
  problem(n_r,n_z,l_z);
 
 // Run the unsteady simulation
 problem.unsteady_run(t_max,dt);
} // End of main

