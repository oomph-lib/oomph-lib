//Three-dimensional free-surface test case including insoluble
//surfactant transport. This is a 3D implementation of the
//axisymmetric Rayleigh--Plateau problem solved in 
//rayleigh_instability_insoluble_surfactant.cc, which 
//should be a hard test
//because it's implemented in Cartesian coordinates. 
//The main aim is to test the implementation of the surface transport
//equations in 3D (2D surface).
 
// The oomphlib headers   
#include "generic.h"
#include "navier_stokes.h"
#include "fluid_interface.h"

// The basic mesh
#include "meshes/single_layer_cubic_spine_mesh.h"

using namespace std;
using namespace oomph;


//==start_of_namespace==================================================
/// Namepspace for global parameters, chosen from Campana et al. as
/// in the axisymmetric problem.
//======================================================================
namespace Global_Physical_Variables
{
 //Film thickness parameter
 double Film_Thickness = 0.2;

 /// Reynolds number
 double Re = 40.0;

 /// Womersley = Reynolds times Strouhal
 double ReSt = Re; 
 
 /// Product of Reynolds and Froude number
 double ReInvFr = 0.0;

 /// Capillary number
 double Ca =  pow(Film_Thickness,3.0);

 /// Direction of gravity
 Vector<double> G(3);

// Wavelength of the domain
 double Alpha = 1.047;
 
 /// Free surface cosine deformation parameter
 double Epsilon = 1.0e-3;

 /// \short Surface Elasticity number (weak case)
 double Beta = 3.6e-3;

 /// \short Surface Peclet number
 double Peclet_S = 4032.0;

 /// \shorT Sufrace Peclet number multiplied by Strouhal number
 double Peclet_St_S = 1.0; 
 
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
class SpineSurfaceMarangoniSurfactantFluidInterfaceElement : 
  public SpineSurfaceFluidInterfaceElement<ELEMENT>
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
   
   //Now interpolate the temperature and surfactant concentration
   double C=0.0;
   for(unsigned l=0;l<n_node;l++)
    {
     C += this->nodal_value(l,C_index)*psi(l);
    }
   
   //Get the Elasticity numbers
   double Beta = this->beta();
   //Return the variable surface tension
   return (1.0 - Beta*(C-1.0));
  }

 /// \short Fill in the contribution to the residuals
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

 
  /// \short Overload the Helper function to calculate the residuals and 
  /// jacobian entries. This particular function ensures that the
  /// additional entries are calculated inside the integration loop
  void add_additional_residual_contributions_interface(Vector<double> &residuals, DenseMatrix<double> &jacobian,
                                                       const unsigned &flag,const Shape &psif, const DShape &dpsifds,
                                                       const Vector<double> &interpolated_x, const Vector<double> &interpolated_n, 
                                                       const double &W,const double &J)
  {
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

   //Now calculate the concentration and derivatives at this point
   //Assuming the same shape functions are used (which they are)
   double interpolated_C = 0.0;
   Vector<double> interpolated_dCds(2,0.0);
   double dCdt = 0.0;
   //The tangent vectors and velocity
   const unsigned ndim = this->node_pt(0)->ndim();
   DenseMatrix<double> interpolated_tangent(2,ndim,0.0);
   Vector<double> interpolated_u(ndim,0.0);
   Vector<double> mesh_velocity(ndim,0.0);
   DenseMatrix<double> interpolated_duds(ndim,2,0.0);
   if(ndim != u_index.size())
    {
      throw OomphLibError("Dimension Incompatibility",
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
    }

   //Loop over the shape functions
   for(unsigned l=0;l<n_node;l++)
     {
      const double psi = psif(l);
      interpolated_C += this->nodal_value(l,c_index)*psi;
      dCdt += dcdt_surface(l)*psi;
      //Velocity and Mesh Velocity
      for(unsigned j=0;j<ndim;j++)
       {
        interpolated_u[j] += this->nodal_value(l,u_index[j])*psi;
        mesh_velocity[j] += this->dnodal_position_dt(l,j)*psi;
       }
      //Loop over surface coordinates
      for(unsigned j=0;j<2;j++)
       {
        const double dpsi = dpsifds(l,j);
        interpolated_dCds[j] += this->nodal_value(l,c_index)*dpsi;

        for(unsigned i=0;i<ndim;i++)
         {
          interpolated_tangent(j,i) += this->nodal_position(l,i)*dpsi;
          interpolated_duds(i,j) += this->nodal_value(l,u_index[i])*dpsi;
         }
       }
     }
   
       

   
   //Now we need to work out surface metric tensor components
   DenseMatrix<double> a(2,2,0.0);
   for(unsigned al=0;al<2;al++)
    {
     for(unsigned be=0;be<(al+1);be++)
      {
       //Take the dot products of the tangent vectors
       for(unsigned i=0;i<3;++i)
        {
         a(al,be) += interpolated_tangent(al,i)*interpolated_tangent(be,i);
        }
      }
    }
   
   const double adet = a(0,0)*a(1,1) - a(1,0)*a(0,1); 
   
   //Work out the contravariant form by inverting the matrix
   DenseMatrix<double> a_up(2,2,0.0);
   a_up(0,0) = a(1,1)/adet; a_up(0,1) = -a(0,1)/adet;
   a_up(1,0) = a_up(0,1); a_up(1,1) = a(0,0)/adet;
   
   
   //No need for normal vector because we already have J
   
   //Pre-compute advection term
   double interpolated_advection_term = 0.0;
   for(unsigned al=0;al<2;al++)
    {
     for(unsigned be=0;be<2;be++)
      {
       for(unsigned i=0;i<3;i++)
        {
         interpolated_advection_term +=
          a_up(al,be)*interpolated_tangent(be,i)*(interpolated_dCds[al]*(interpolated_u[i] - mesh_velocity[i])
                                                  + interpolated_C*interpolated_duds(i,al));
        }
      }
    }
    
   
   //Now we add the flux term to the appropriate residuals
   for(unsigned l=0;l<n_node;l++)
    {
     //Read out the apprporiate local equation
     local_eqn = this->nodal_local_eqn(l,c_index);
     
     //If not a boundary condition
     if(local_eqn >= 0)
      {
       //Time derivative term
       residuals[local_eqn] += dCdt*psif(l)*W*J;
       
       //First Advection term
       residuals[local_eqn] += interpolated_advection_term*psif(l)*W*J;
       
       //Diffusion term
       double diffusion_term = 0.0;
       for(unsigned al=0;al<2;al++)
        {
         for(unsigned be=0;be<2;be++)
          {
           diffusion_term += a_up(al,be)*interpolated_dCds[al]*dpsifds(l,be);
          }
        }
       residuals[local_eqn] += (1.0/Pe_s)*diffusion_term*W*J;
       
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
              time_stepper_pt->weight(1,0)*psif(l2)*psif(l)*W*J;

             double tmp = 0.0;
             for(unsigned al=0;al<2;al++)
              {
               for(unsigned be=0;be<2;be++)
                {
                 for(unsigned j=0;j<3;j++)
                  {
                   tmp += 
                    a_up(al,be)*interpolated_tangent(be,j)*(dpsifds(l2,al)*(interpolated_u[j] - mesh_velocity[j])
                                                            + psif(l2)*interpolated_duds(j,al));
                  }
                }
              }

              jacobian(local_eqn,local_unknown) += tmp*psif(l)*W*J;
             
             //Reset tmp term
             tmp = 0.0;
             for(unsigned al=0;al<2;al++)
              {
               for(unsigned be=0;be<2;be++)
                {
                 tmp += a_up(al,be)*dpsifds(l2,al)*dpsifds(l,be);
                }
              }
             
             jacobian(local_eqn,local_unknown) += (1.0/Pe_s)*tmp*W*J;
             
            }


           //Loop over the velocity components
           for(unsigned i2=0;i2<ndim;i2++)
            {
             
             //Get the unknown
             local_unknown = this->nodal_local_eqn(l2,u_index[i2]);
             
             
             //If not a boundary condition
             if(local_unknown >= 0)
              {
               //Compute the bits from the advection term
               double tmp = 0.0;
               for(unsigned al=0;al<2;al++)
                {
                 for(unsigned be=0;be<2;be++)
                  {
                   tmp +=
                    a_up(al,be)*interpolated_tangent(be,i2)*(interpolated_dCds[al]*psif(l2)
                                                             + interpolated_C*dpsifds(l2,al));
                  }
                }

               //Add to the jacobian
               jacobian(local_eqn,local_unknown) += tmp*psif(l)*W*J;
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
  SpineSurfaceMarangoniSurfactantFluidInterfaceElement(FiniteElement* const &element_pt, const int &face_index) :
   SpineSurfaceFluidInterfaceElement<ELEMENT>(element_pt,face_index)
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
   Vector<double> s(2);
   Vector<double> n(3);
   Vector<double> u(3);
  
   
   outfile << this->tecplot_zone_string(n_plot);
   
   // Loop over plot points
   unsigned num_plot_points=this->nplot_points(n_plot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     // Get local coordinates of plot point
     this->get_s_plot(iplot,n_plot,s);
     //Get the outer unit normal
     this->outer_unit_normal(s,n);

     double u_n = 0.0;
     for(unsigned i=0;i<3;i++)
      {
       u[i] = this->interpolated_u(s,i);
       u_n += u[i]*n[i];
      }
     
     
     
     //Output the x,y,u,v 
     for(unsigned i=0;i<3;i++) outfile << this->interpolated_x(s,i) << " ";
     for(unsigned i=0;i<3;i++) outfile << u[i] << " ";      
     //Output a dummy pressure
     outfile << 0.0 << " ";
     //Output the concentration
     outfile << interpolated_C(s) << " ";
     //Output the interfacial tension
     outfile << sigma(s) << " " << u[0] - u_n*n[0] << " " << u[1] - u_n*n[1] << " " << u[2] - u_n*n[2] << std::endl;
    }
   outfile << std::endl;
  }

 
 /// \short Compute the concentration intergated over the surface area
 double integrate_c() const
  {
   //Find out how many nodes there are
   unsigned n_node = this->nnode();
   
   //Set up memeory for the shape functions
   Shape psif(n_node);
   DShape dpsifds(n_node,2);

   //Set the value of n_intpt
   unsigned n_intpt = this->integral_pt()->nweight();
   
   //Storage for the local coordinate
   Vector<double> s(2);

   //Storage for the total mass
   double mass = 0.0;
   
   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Get the local coordinate at the integration point
     for(unsigned i=0;i<2;i++) {s[i] = this->integral_pt()->knot(ipt,i);}
     
     //Get the integral weight
     double W = this->integral_pt()->weight(ipt);
     
     //Call the derivatives of the shape function
     this->dshape_local_at_knot(ipt,psif,dpsifds);


     //Define and zero the tangent Vectors and local velocities
     double interpolated_g[2][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0}};
     double interpolated_c = 0.0;
     
     
     //Loop over the shape functions
     for(unsigned l=0;l<n_node;l++)
      {
       interpolated_c += this->nodal_value(l,C_index)*psif(l);
       //Loop over directional components (now three)
       for(unsigned i=0;i<3;i++)
        {
         //Calculate the local tangent vectors
         interpolated_g[0][i] += this->nodal_position(l,i)*dpsifds(l,0);
         interpolated_g[1][i] += this->nodal_position(l,i)*dpsifds(l,1);  
        }
      }


     //Calculate the local metric tensor
     //The dot product of the two tangent vectors
     double gmet[2][2];
     for(unsigned al=0;al<2;al++)
      {
       for(unsigned be=0;be<2;be++)
        {
         //Initialise to zero
         gmet[al][be] = 0.0;
         //Add the dot product contributions
         for(unsigned i=0;i<3;i++)
          {
           gmet[al][be] += interpolated_g[al][i]*interpolated_g[be][i];
          }
        }
      }

     // Define the normal vector (cross product of tangent vectors)
     Vector<double> interpolated_n(3); 
     interpolated_n[0] = 
      interpolated_g[0][1]*interpolated_g[1][2] - 
      interpolated_g[0][2]*interpolated_g[1][1];
     interpolated_n[1] = 
      interpolated_g[0][2]*interpolated_g[1][0] - 
      interpolated_g[0][0]*interpolated_g[1][2]; 
     interpolated_n[2] = interpolated_g[0][0]*interpolated_g[1][1] - 
      interpolated_g[0][1]*interpolated_g[1][0];
     
     // Calculate the length of the vector
     double slength =  interpolated_n[0]*interpolated_n[0] +
      interpolated_n[1]*interpolated_n[1] +
      interpolated_n[2]*interpolated_n[2];
     
     //Set the determinant of the local metric tensor, 
     //which is equal to the length of the normal vector
     double J = sqrt(slength);
     
     
     mass += interpolated_c*W*J;
    }
   return mass;
  }



};


//Define the default physical value to be one
template<class ELEMENT>
double SpineSurfaceMarangoniSurfactantFluidInterfaceElement<ELEMENT>::Default_Physical_Constant_Value = 1.0;

}



//=======================================================================
/// Function-type-object to perform comparison of elements in y-direction
//=======================================================================
class ElementCmp
{
public:

 /// Comparison. Are the values identical or not?
 bool operator()(GeneralisedElement* const &x, GeneralisedElement* const &y) 
  const
  {
   FiniteElement* cast_x = dynamic_cast<FiniteElement*>(x);
   FiniteElement* cast_y = dynamic_cast<FiniteElement*>(y);

   if((cast_x ==0) || (cast_y==0)) {return 0;}
   else
    {return (cast_x->node_pt(0)->x(1)< cast_y->node_pt(0)->x(1));}
             
  }
};


namespace oomph
{
//===================================================================
/// Deform the existing cubic spine mesh into a annular section
/// with spines directed radially inwards from the wall
//===================================================================

template<class ELEMENT>
class AnnularSpineMesh : public SingleLayerCubicSpineMesh<ELEMENT>
{

public:
  AnnularSpineMesh(const unsigned &n_r, const unsigned &n_y, 
                   const unsigned &n_theta, 
		   const double &r_min, const double &r_max, 
                   const double &l_y,
		   const double &theta_min, const double &theta_max, 
		   TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper):
   //This will make a cubic mesh with n_r in the x-direction
   //n_y in the y-direction and n_theta in the z-direction
   //The coordinates will run from 0 to r_max, 0 to l_y and 0 to theta_max
   SingleLayerCubicSpineMesh<ELEMENT>(n_theta,n_y,n_r,r_max,l_y,theta_max,
                                      time_stepper_pt)
  {  
   
    //Find out how many nodes there are
    unsigned n_node = this->nnode();

    //Loop over all nodes
    for (unsigned n=0;n<n_node;n++)
      {
	//pointer to node
	Node* nod_pt=this->node_pt(n);
	SpineNode* spine_node_pt=dynamic_cast<SpineNode*>(nod_pt);
	//Get x/y/z coordinates
	double x_old=nod_pt->x(0);
	double y_old=nod_pt->x(1);
	double z_old=nod_pt->x(2);

	//Mapping
	
	double r = r_min + (r_max-r_min)*z_old;
        double theta = (theta_min + (theta_max-theta_min)*x_old);
	double y = y_old;
	
	
	if(spine_node_pt->spine_pt()->ngeom_parameter()==0)
	  {
           spine_node_pt->h() =
            Global_Physical_Variables::Film_Thickness; 
	    spine_node_pt->spine_pt()->add_geom_parameter(theta);
	  }
	
	//cout << spine_node_pt->spine_pt()->ngeom_parameter()  << std::endl;
	//Set new nodal coordinates
	nod_pt->x(0)=r*cos(theta);
	nod_pt->x(2)=r*sin(theta);
	nod_pt->x(1)=y;
	//cout << "one" << theta << std::endl;	
      }

  }


    virtual void spine_node_update(SpineNode* spine_node_pt)
  {
   //Get fraction along the spine
   double W = spine_node_pt->fraction();
   //Get theta along the spine
   double theta = spine_node_pt->spine_pt()->geom_parameter(0);
   //Get spine height
   double H = spine_node_pt->h();
   //Set the value of y
   spine_node_pt->x(0) = (1.0-W*H)*cos(theta);
   spine_node_pt->x(2) = (1.0-W*H)*sin(theta);
  }
};
}

//======================================================================
/// Single fluid interface problem including transport of an
/// insoluble surfactant.
//======================================================================
template<class ELEMENT, class TIMESTEPPER>
class InterfaceProblem : public Problem
{
 
public:

 /// Constructor: Pass number of elements in x and y directions. Also lengths
 /// of the domain in x- and y-directions and the height of the layer

 InterfaceProblem(const unsigned &n_r, const unsigned &n_y, 
                  const unsigned &n_theta, const double &r_min,
                  const double &r_max, const double &l_y,  
                  const double &theta_max);
 
 /// Spine heights/lengths are unknowns in the problem so their
 /// values get corrected during each Newton step. However,
 /// changing their value does not automatically change the
 /// nodal positions, so we need to update all of them
 void actions_before_newton_convergence_check(){Bulk_mesh_pt->node_update();}

 /// Run an unsteady simulation with specified number of steps
 void unsteady_run(const unsigned& nstep); 

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
 /// Compute the total mass
 double compute_total_mass()
  {
   double mass = 0.0;
   
   // Determine number of 1D interface elements in mesh
   const unsigned n_interface_element = Surface_mesh_pt->nelement();
   
   // Loop over the interface elements
   for(unsigned e=0;e<n_interface_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     SpineSurfaceMarangoniSurfactantFluidInterfaceElement<ELEMENT>* el_pt = 
      dynamic_cast<SpineSurfaceMarangoniSurfactantFluidInterfaceElement<ELEMENT>*>
      (Surface_mesh_pt->element_pt(e));

     mass += el_pt->integrate_c();
    }
   return mass;
  }


 
 private:

 /// Trace file
 ofstream Trace_file;

 /// Axial lengths of domain
 double R_max;

 double L_y;

 /// Pointer to bulk mesh
 AnnularSpineMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the surface mes
 Mesh* Surface_mesh_pt;

 /// Pointer to a node for documentation purposes
 Node* Document_node_pt;

};


//====================================================================
/// Problem constructor
//====================================================================
template<class ELEMENT, class TIMESTEPPER>
InterfaceProblem<ELEMENT,TIMESTEPPER>::InterfaceProblem
(const unsigned &n_r, const unsigned &n_y,const unsigned &n_theta,
 const double &r_min,
 const double &r_max, const double &l_y, const double &theta_max)
 : R_max(r_max), L_y(l_y)
{  
 //this->linear_solver_pt() = new HSL_MA42;

 //static_cast<HSL_MA42*>(this->linear_solver_pt())->lenbuf_factor0() = 3.0;
 //static_cast<HSL_MA42*>(this->linear_solver_pt())->lenbuf_factor1() = 3.0;
 //static_cast<HSL_MA42*>(this->linear_solver_pt())->lenbuf_factor2() = 3.0;
 
 //Allocate the timestepper
 add_time_stepper_pt(new TIMESTEPPER); 

 //Now create the bulk mesh -- this should be your new annular mesh
 Bulk_mesh_pt = new AnnularSpineMesh<ELEMENT>
   (n_r,n_y,n_theta,r_min,r_max,l_y,0.0,theta_max,time_stepper_pt());

 //Set the documented node
 Document_node_pt = Bulk_mesh_pt->boundary_node_pt(5,0);


 //Create the surface mesh that will contain the interface elements
 //First create storage, but with no elements or nodes
 Surface_mesh_pt = new Mesh;
   
 // Loop over those elements adjacent to the free surface,
 // which we shall choose to be the upper surface
 for(unsigned e1=0;e1<n_y;e1++)
  {
   for(unsigned e2=0;e2<n_theta;e2++)
    {
     // Set a pointer to the bulk element we wish to our interface
     // element to
     FiniteElement* bulk_element_pt =
      Bulk_mesh_pt->finite_element_pt(n_theta*n_y*(n_r-1) + e2 + e1*n_theta);

     // Create the interface element (on face 3 of the bulk element)
     FiniteElement* interface_element_pt =
      new SpineSurfaceMarangoniSurfactantFluidInterfaceElement<ELEMENT>(bulk_element_pt,3);

   // Add the interface element to the surface mesh
   this->Surface_mesh_pt->add_element_pt(interface_element_pt);
    }
  }
 
 // Add the two sub-meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all sub-meshes into a single mesh
 build_global_mesh();
 
 //Pin all nodes on the bottom
 unsigned long n_boundary_node = Bulk_mesh_pt->nboundary_node(0);
 for(unsigned long n=0;n<n_boundary_node;n++) 
  {
   for(unsigned i=0;i<3;i++)
    {
     Bulk_mesh_pt->boundary_node_pt(0,n)->pin(i);
    }
  }

 //On the front and back (y=const) pin in y-direction
 for(unsigned b=1;b<4;b+=2)
  {
   n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
   for(unsigned long n=0;n<n_boundary_node;n++)
    {
     Bulk_mesh_pt->boundary_node_pt(b,n)->pin(1);
    }
  }

 //On sides pin in z-direction
 for(unsigned b=4;b<5;b+=2)
  {
   n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
   for(unsigned long n=0;n<n_boundary_node;n++)
    {
     Bulk_mesh_pt->boundary_node_pt(b,n)->pin(2);
    }
  }

 // pinning the top surface
 for(unsigned b=2;b<3;b++)
  {
   n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
   for(unsigned long n=0;n<n_boundary_node;n++)
    {
     Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0);
    }
  }

 //Loop over the upper surface and set the surface concentration
 {
  unsigned b=5;
  n_boundary_node = Bulk_mesh_pt->nboundary_node(b);
  for(unsigned long n=0;n<n_boundary_node;n++) 
   {
    Bulk_mesh_pt->boundary_node_pt(b,n)->set_value(3,1.0);;
   }
 }

 
 //Create a Data object whose single value stores the
 //external pressure
 Data* external_pressure_data_pt = new Data(1);
 
 // Set and pin the external pressure to some random value
 external_pressure_data_pt->set_value(0,1.31);
 external_pressure_data_pt->pin(0);

 //Complete the problem setup to make the elements fully functional

 //Loop over the elements in the layer
 unsigned n_bulk=Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_bulk;e++)
  {
   //Cast to a fluid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;
   el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr;
   el_pt->g_pt() = &Global_Physical_Variables::G;
  }

 //Loop over 2D interface elements and set capillary number and 
 //the external pressure
 unsigned long interface_element_pt_range = 
  Surface_mesh_pt->nelement();
 for(unsigned e=0;e<interface_element_pt_range;e++)
  {
   //Cast to a interface element
   SpineSurfaceMarangoniSurfactantFluidInterfaceElement<ELEMENT>* el_pt = 
    dynamic_cast<SpineSurfaceMarangoniSurfactantFluidInterfaceElement<ELEMENT>*>
    (Surface_mesh_pt->element_pt(e));

   //Set the Capillary number
   el_pt->ca_pt() = &Global_Physical_Variables::Ca;
   
   //Pass the Data item that contains the single external pressure value
   el_pt->set_external_pressure_data(external_pressure_data_pt);

   // Set the surface elasticity number
   el_pt->beta_pt() = &Global_Physical_Variables::Beta;

   // Set the surface peclect number
   el_pt->peclet_s_pt() = &Global_Physical_Variables::Peclet_S;

   // Set the surface peclect number multiplied by strouhal number
   el_pt->peclet_strouhal_s_pt() = &Global_Physical_Variables::Peclet_St_S;

  }

 //Do equation numbering
 cout << assign_eqn_numbers() << std::endl; 

 //Now sort the elements to minimise frontwidth
 std::sort(mesh_pt()->element_pt().begin(),
           mesh_pt()->element_pt().end(),
           ElementCmp());

}

   

//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 

 //Output the time
 cout << "Time is now " << time_pt()->time() << std::endl;

 // Doc
 Trace_file << time_pt()->time();
 Trace_file << " "  << Document_node_pt->x(0) << " "
            << this->compute_total_mass()
            << std::endl;


 // Output solution, bulk elements followed by surface elements
 /*sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();*/
 //Change the file name and output surface in separate file
 sprintf(filename,"%s/surface%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Surface_mesh_pt->output(some_file,npts);
 some_file.close();
 
}

 


  //=============================================================================
/// Unsteady run with specified number of steps
//=============================================================================
template<class ELEMENT, class TIMESTEPPER>
void InterfaceProblem<ELEMENT,TIMESTEPPER>::unsteady_run(const unsigned& nstep)
{

 // Increase maximum residual
 Problem::Max_residuals=500.0;

 //Distort the mesh
 const double epsilon=Global_Physical_Variables::Epsilon;
 unsigned Nspine = Bulk_mesh_pt->nspine();
 for(unsigned i=0;i<Nspine;i++)
  {
    double y_value = Bulk_mesh_pt->boundary_node_pt(0,i)->x(1);
    
    Bulk_mesh_pt->spine_pt(i)->height() =
     Global_Physical_Variables::Film_Thickness*(1.0 +  
                                                + epsilon*cos(Global_Physical_Variables::Alpha*y_value));
  }

 //Make sure to update 
 Bulk_mesh_pt->node_update();

 // Doc info object
 DocInfo doc_info;

 // Set output directory
 doc_info.set_directory("RESLT");
 doc_info.number()=0;

 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);

 Trace_file << "VARIABLES=\"time\",";
 Trace_file << "\"h<sub>left</sub>\",\"h<sub>right</sub>\"";
 Trace_file << "\n";

 //Set value of dt
 double  dt = 0.1;

 //Initialise all time values
 assign_initial_values_impulsive(dt);
  
 //Doc initial solution
 doc_solution(doc_info);

//Loop over the timesteps
 for(unsigned t=1;t<=nstep;t++)
  {
   cout << "TIMESTEP " << t << std::endl;
   
   //Take one fixed timestep
   unsteady_newton_solve(dt);

   // Doc solution
   doc_info.number()++;
   doc_solution(doc_info);
   }

}


//==start_of_main========================================================
/// Driver code for unsteady two-layer fluid problem. If there are
/// any command line arguments, we regard this as a validation run
/// and perform only two steps.

// In my version we will change nsteps in the programs
//======================================================================
int main(int argc, char *argv[]) 
{

 // Set physical parameters:

 //Set direction of gravity: Vertically downwards
 Global_Physical_Variables::G[0] = 0.0;
 Global_Physical_Variables::G[1] = 0.0;
 Global_Physical_Variables::G[2] = 0.0;

 //Set the film thickness
 const double h = Global_Physical_Variables::Film_Thickness;
 
 //Axial lngth of domain
 double r_max = 1.0;
 double r_min = r_max - h;

 const double pi = 4.0*atan(1);

 double alpha = Global_Physical_Variables::Alpha;

 double l_y = pi/alpha;
 
 double theta_max = 0.5*pi;
  
 // Number of elements in r and theta direction
 unsigned n_r = 4;
 
 unsigned n_theta = 4;
 
 // Number of elements in axial direction
 unsigned n_y = 20;
 

 {
  //Set up the spine test problem
  //The minimum values are all assumed to be zero.
  InterfaceProblem<SpineElement<QCrouzeixRaviartElement<3> >,BDF<2> >
   problem(n_r,n_y,n_theta,r_min,r_max,l_y,theta_max);
  
  // Number of steps: 
  unsigned nstep;
  if(argc > 1)
   {
    // Validation run: Just five steps
    nstep=5;
   }
  else
   {
    // Full run otherwise
    nstep=1000;
   }
  
  // Run the unsteady simulation
  problem.unsteady_run(nstep);
 }
} // End of main


