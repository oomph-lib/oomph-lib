//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
// Driver for Solid deformation -- driven by boundary motion which
// is imposed via Lagrange multipliers

//Oomph-lib includes
#include "generic.h"
#include "solid.h"
#include "constitutive.h"

//The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

namespace oomph
{


//=========================================================================
/// Warped line in 2D space 
//=========================================================================
class WarpedLine : public GeomObject
{

public:


 /// Constructor: Pass height, length and warping amplitude
 WarpedLine(const double& height, 
            const double& length,
            const double& ampl) : GeomObject(1,2)
  {
   Height=height;
   Length=length;
   Ampl=ampl;
  }
 
 
 /// Broken copy constructor
 WarpedLine(const WarpedLine& dummy) 
  { 
   BrokenCopy::broken_copy("WarpedLine");
  } 
 
 /// Broken assignment operator
 void operator=(const WarpedLine&) 
  {
   BrokenCopy::broken_assign("WarpedLine");
  }


 /// Empty Destructor
 ~WarpedLine(){}

 /// \short Position Vector at Lagrangian coordinate zeta 
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   // Position Vector
   r[0] = zeta[0]+5.0*Ampl*zeta[0]*(zeta[0]-Length)*(zeta[0]-0.7*Length);
   r[1] = Height+
    Ampl*0.5*(1.0-cos(2.0*MathematicalConstants::Pi*zeta[0]/Length));
  }
 

 /// \short Parametrised position on object: r(zeta). Evaluated at
 /// previous timestep. t=0: current time; t>0: previous
 /// timestep. Forward to steady version
 void position(const unsigned& t, const Vector<double>& zeta,
                       Vector<double>& r) const
  {
   position(zeta,r);
  }

 /// Access to amplitude
 double& ampl() {return Ampl;}
 
private:

 /// Height
 double Height;

 /// Wavelength of perturbation
 double Length;

 /// Amplitude of perturbation
 double Ampl;

};




///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


//======================================================================
/// A class for elements that allow the imposition of a displacment
/// constraint for "bulk" solid elements via a Lagrange multiplier.
/// The geometrical information can be read from the FaceGeometry<ELEMENT> 
/// class and and thus, we can be generic enough without the need to have
/// a separate equations class.
/// \b NOTE: Currently (and for the foreseeable future) this 
/// element only works with bulk elements that do not have
/// generalised degrees of freedom (so it won't work with
/// Hermite-type elements, say). The additional functionality 
/// to deal with such elements could easily be added (once a 
/// a suitable test case is written). For now we simply throw
/// errors if an attempt is made to use the element with an unsuitable
/// bulk element.
//======================================================================
template <class ELEMENT>
class ImposeDisplacementLagrangeMultiplierElement : 
  public virtual FaceGeometry<ELEMENT>, 
  public virtual SolidFaceElement
{
 
public:

 /// \short Constructor, which takes a "bulk" element and the 
 /// value of the index and its limit
 ImposeDisplacementLagrangeMultiplierElement(FiniteElement* const &element_pt, 
                      const int &face_index) : 
  FaceGeometry<ELEMENT>(), FaceElement(), Boundary_shape_geom_object_pt(0)
  { 
   
   
#ifdef PARANOID
   {
    //Check that the bulk element is not a refineable 3d element
    //If it's three-d
    if(element_pt->dim()==3)
     {
      //Is it refineable
      if(dynamic_cast<RefineableElement*>(element_pt))
       {
        //Issue a warning
        OomphLibWarning(
         "This flux element will not work correctly if nodes are hanging\n",
         "ImposeDisplacementLagrangeMultiplierElement::Constructor",
         OOMPH_EXCEPTION_LOCATION);
       }
     }
   }
   {
    // Check that the bulk element does not require generalised positional
    // degrees of freedom
    if(element_pt->nnodal_position_type()!=1)
     {      
      throw OomphLibError(
       "ImposeDisplacementLagrangeMultiplierElement cannot (currently) be used with elements that have generalised positional dofs",
       "ImposeDisplacementLagrangeMultiplierElement::ImposeDisplacementLagrangeMultiplierElement()",
       OOMPH_EXCEPTION_LOCATION);
     }
   }
#endif
 
   //Attach the geometrical information to the element. N.B. This function
   //also assigns nbulk_value from the required_nvalue of the bulk element
   element_pt->build_face_element(face_index,this);

   // We need dim additional values (for the components of the
   // Lagrange multiplier vector) at each of the element's nodes
   unsigned n_node=this->nnode();
   Vector<unsigned> nadditional_data_values(n_node);
   for (unsigned j=0;j<n_node;j++)
    {
     nadditional_data_values[j]=element_pt->dim();
    }
   resize_nodes(nadditional_data_values); 
   
  }
 
 /// \short Access to GeomObject that specifies the prescribed 
 /// boundary displacement; GeomObject is assumed to be
 /// parametrised by the same coordinate the is used as
 /// the boundary coordinate in the bulk solid mesh to which
 /// this element is attached.
 GeomObject*& boundary_shape_geom_object_pt()
  {
   return Boundary_shape_geom_object_pt;
  }


 /// Fill in the residuals
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic routine with the flag set to 0
   fill_in_generic_contribution_to_residuals_displ_lagr_multiplier(
    residuals,GeneralisedElement::Dummy_matrix,0);
  }

 
 /// Fill in contribution from Jacobian
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_contribution_to_residuals_displ_lagr_multiplier(
    residuals,jacobian,1);
  }

 
 /// \short Output function
 void output(std::ostream &outfile, const unsigned &n_plot)
  {
   // Elemental dimension
   unsigned dim_el=dim();

   //Find the number of positional types
   unsigned n_position_type = this->nnodal_position_type();
   
#ifdef PARANOID
   if(n_position_type!=1)
    {      
     throw OomphLibError(
      "ImposeDisplacementLagrangeMultiplierElement cannot (currently) be used with elements that have generalised positional dofs",
      "ImposeDisplacementLagrangeMultiplierElement::output()",
      OOMPH_EXCEPTION_LOCATION);
    }
#endif


   //Local coord
   Vector<double> s(dim_el);
      
   // # of nodes, 
   unsigned n_node=nnode();
   Shape psi(n_node,n_position_type);

   //Tecplot header info 
   outfile << "ZONE I=" << n_plot << std::endl;
   
   //Loop over plot points
   for(unsigned l=0;l<n_plot;l++)
    {
     s[0] = -1.0 + l*2.0/(n_plot-1);
     
     // Get shape function
     shape(s,psi);
     
     //Calculate the Eulerian coordinates and Lagrange multiplier
     Vector<double> x(dim_el+1,0.0);
     Vector<double> lambda(dim_el+1,0.0);
     Vector<double> zeta(dim_el,0.0);
     for(unsigned j=0;j<n_node;j++) 
      {
       // higher dimensional quantities
       for(unsigned i=0;i<dim_el+1;i++)
        {
         x[i]+=nodal_position(j,i)*psi(j,0); // hierher need to sort
                                             // this out properly
         lambda[i]+=node_pt(j)->value(Nbulk_value[j]+i)*psi(j,0);
        }
       //In-element quantities
       for(unsigned i=0;i<dim_el;i++)
        {
         //Loop over positional types
         for (unsigned k=0;k<n_position_type;k++)
          {
           zeta[i]+=zeta_nodal(j,k,i)*psi(j,k);
          }
        }
      }

     // Get prescribed wall shape
     Vector<double> r_prescribed(dim_el+1);
     Boundary_shape_geom_object_pt->position(zeta,r_prescribed);

     //Output stuff
     outfile << x[0] << " "
             << x[1] << " "
             << -lambda[0] << " "
             << -lambda[1] << " "
             << r_prescribed[0] << " "
             << r_prescribed[1] << " "
             <<  std::endl;
    }
  }


 /// \short Output function
 void output(std::ostream &outfile)
  {
   unsigned n_plot=5;
   output(outfile,n_plot);
  }

 /// \short Set function for the boundary number in bulk mesh -- needed
 /// to get access to the appropriate boundary coordinate
 void set_boundary_number_in_bulk_mesh(const unsigned& b) 
  { 
   Boundary_number_in_bulk_mesh=b;
#ifdef PARANOID
   Boundary_number_in_bulk_mesh_has_been_set=true;
#endif
  }


 /// \short Specify the values of the boundary coordinate, zeta,
 /// at local node n in this element. k is the type of the
 /// coordinate, i identifies the coordinate direction). 
 /// Note: Boundary coordinates will have been set up when
 /// creating the underlying mesh, and their values will have 
 /// been stored at the nodes.
 double zeta_nodal(const unsigned &n,  const unsigned &k, const unsigned &i)
  const
  {
   //Vector in which to hold the intrinsic coordinate
   Vector<double> zeta(dim());
 
   //Get the boundary coordinate at node n
   node_pt(n)->get_coordinates_on_boundary(
    Boundary_number_in_bulk_mesh,k,zeta);

   //Return the individual coordinate
   return zeta[i];
  }


 /// \short Calculate the interpolated value of zeta, the boundary coordinate
 /// on the mesh boundary to which this element is attached, 
 /// at the local coordinate s of the element
 void interpolated_zeta(const Vector<double> &s, Vector<double> &zeta) const
  {
   //Find the number of nodes
   unsigned n_node = nnode();

   //Find the number of positional types
   unsigned n_position_type = this->nnodal_position_type();
  
   //Storage for the shape functions
   Shape psi(n_node,n_position_type);

   //Get the values of the shape functions at the local coordinate s
   shape(s,psi);
   
   //Find the number of coordinates
   unsigned ncoord = dim();
   
   //Initialise the value of zeta to zero
   for(unsigned i=0;i<ncoord;i++) {zeta[i] = 0.0;}

   //Add the contributions from each nodal dof to the interpolated value
   //of zeta.
   for(unsigned l=0;l<n_node;l++)
    {
     for(unsigned k=0;k<n_position_type;k++)
      {
       for(unsigned i=0;i<ncoord;i++)
        {
         zeta[i] += this->zeta_nodal(l,k,i)*psi(l,k);
        }
      }  
    }
  }

 /// \short For a given value of zeta, the boundary coordinate
 /// (assigned in the mesh and stored at the nodes), find the local 
 /// coordinate in this element that corresponds to zeta. This is achieved
 /// in generality by using Newton's method to find the value s such that
 /// interpolated_zeta(s) is equal to the desired value of zeta.
 /// If zeta cannot be located in this element, geom_object_pt is set
 /// to NULL. If zeta is located in this element, we return its "this"
 /// pointer.
 void locate_zeta(const Vector<double> &zeta, GeomObject*& geom_object_pt,
                  Vector<double> &s)
  {
   using namespace Locate_zeta_helpers;

   //Find the number of coordinates
   unsigned ncoord = this->dim();//DIM-1; 

   //Assign storage for the vector and matrix used in Newton's method
   Vector<double> dx(ncoord,0.0);
   DenseDoubleMatrix jacobian(ncoord,ncoord,0.0);

   //Initialise s to the middle of the element
   for(unsigned i=0;i<ncoord;i++) 
    {
     s[i] = 0.5*(s_max()+s_min());
    }
   
   //Counter for the number of Newton steps
   unsigned count=0;

   //Control flag for the Newton loop
   bool keep_going=true;
   
   //Storage for the interpolated value of zeta
   Vector<double> inter_zeta(ncoord);
   
   //Get the value of zeta at the initial guess
   interpolated_zeta(s,inter_zeta);
   
   //Set up the residuals
   for(unsigned i=0;i<ncoord;i++) {dx[i] = zeta[i] - inter_zeta[i];}
   
   //Main Newton Loop
   do    // start of do while loop
    {
     //Increase loop counter
     count++;

     //Bail out if necessary
     if(count > Max_newton_iterations)
      {
       std::ostringstream error_message;
       error_message << "Newton solver not converged in " 
                     << count << " steps\n"
                     << "Try adjusting the Tolerances, or refining the mesh."
                     << std::endl;
       throw OomphLibError(error_message.str(),
                           "FaceAsGeomObject::locate_zeta()",
                           OOMPH_EXCEPTION_LOCATION);
      }
	     
     //If it's the first time round the loop, check the initial residuals
     if(count==1)
      {
       double maxres = 
        std::abs(*std::max_element(dx.begin(),dx.end(),AbsCmp<double>()));

       //If it's small enough exit
       if(maxres < Newton_tolerance) 
        {
         keep_going=false;
         continue;
        }
      }
     
     //Compute the entries of the Jacobian matrix
     unsigned n_node = this->nnode();
     unsigned n_position_type = this->nnodal_position_type();
     Shape psi(n_node,n_position_type);
     DShape dpsids(n_node,n_position_type,ncoord);

     //Get the local shape functions and their derivatives
     dshape_local(s,psi,dpsids);
     
     //Calculate the values of dzetads
     DenseMatrix<double> interpolated_dzetads(ncoord,ncoord,0.0);
     for(unsigned l=0;l<n_node;l++)
      {
       for(unsigned k=0;k<n_position_type;k++)
        {
         for(unsigned i=0;i<ncoord;i++)
          {
           for(unsigned j=0;j<ncoord;j++)
            {
             interpolated_dzetads(i,j) += 
              this->zeta_nodal(l,k,i)*dpsids(l,k,j);
            }
          }
        }
      }
     
     //The entries of the Jacobian matrix are merely dresiduals/ds
     //i.e. - dzeta/ds
     for(unsigned i=0;i<ncoord;i++)
      {
       for(unsigned j=0;j<ncoord;j++)
        {
         jacobian(i,j) = - interpolated_dzetads(i,j);
        }
      }
     
     //Now solve the damn thing
     jacobian.solve(dx);
     
     //Add the correction to the local coordinates
     for(unsigned i=0;i<ncoord;i++) {s[i] -= dx[i];}
     
     //Get the new residuals
     interpolated_zeta(s,inter_zeta);
     for(unsigned i=0;i<ncoord;i++) {dx[i] = zeta[i] - inter_zeta[i];}
     
     //Get the maximum residuals
     double maxres = 
      std::abs(*std::max_element(dx.begin(),dx.end(),AbsCmp<double>()));
     //If we have converged jump straight to the test at the end of the loop
     if(maxres < Newton_tolerance) 
      {
       keep_going=false; 
       continue;
      } 
    }
   while(keep_going);


 
   //Test that the solution is within the element
   for(unsigned i=0;i<ncoord;i++)
    {
     // We're outside -- return the null pointer for the geom object
     if((s[i] - s_max() >  Rounding_tolerance) || 
        (s_min() - s[i] >  Rounding_tolerance)) 
      {
       geom_object_pt=0; 
       return;
      }
    }
   
   //Otherwise the required point is located in "this" element:
   geom_object_pt = this;
   
   //If we're over the limit by less than the rounding error, adjust
   for(unsigned i=0;i<ncoord;i++)
    {
     if(s[i] > s_max()) {s[i] = s_max();}
     if(s[i] < s_min()) {s[i] = s_min();}
    }
  }
 


protected:

 /// \short Helper function to compute the residuals and, if flag==1, the
 /// Jacobian 
 void fill_in_generic_contribution_to_residuals_displ_lagr_multiplier(
  Vector<double> &residuals, 
  DenseMatrix<double> &jacobian, 
  const unsigned& flag)
  {
   //Find out how many positional dofs there are
   unsigned n_position_type = this->nnodal_position_type();
   
#ifdef PARANOID
   if(n_position_type!=1)
    {      
     throw OomphLibError(
      "ImposeDisplacementLagrangeMultiplierElement cannot (currently) be used with elements that have generalised positional dofs",
      "ImposeDisplacementLagrangeMultiplierElement::fill_in_generic_contribution_to_residuals_displ_lagr_multiplier()",
      OOMPH_EXCEPTION_LOCATION);
    }
#endif

   //Find out how many nodes there are
   unsigned n_node = nnode();
   
   // Dimension of element
   unsigned dim_el=dim();

   //Set up memory for the shape functions
   Shape psi(n_node);
   DShape dpsids(n_node,dim_el); 

   //Set the value of n_intpt
   unsigned n_intpt = integral_pt()->nweight();
 
   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Get the integral weight
     double w = integral_pt()->weight(ipt);
     
     //Only need to call the local derivatives
     dshape_local_at_knot(ipt,psi,dpsids);

     //Calculate the Eulerian coordinates and Lagrange multiplier
     Vector<double> x(dim_el+1,0.0);
     Vector<double> lambda(dim_el+1,0.0);
     Vector<double> zeta(dim_el,0.0);
     DenseMatrix<double> interpolated_a(dim_el,dim_el+1,0.0);   

     // Loop over nodes
     for(unsigned j=0;j<n_node;j++) 
      {
       Node* nod_pt=node_pt(j);

       //Assemble higher-dimensional quantities
       for(unsigned i=0;i<dim_el+1;i++)
        {
         x[i]+=nodal_position(j,i)*psi(j);
         lambda[i]+=nod_pt->value(Nbulk_value[j]+i)*psi(j);
         for(unsigned ii=0;ii<dim_el;ii++)
          {
           interpolated_a(ii,i) += 
            lagrangian_position(j,i)*dpsids(j,ii);
          }
        }        
       for(unsigned k=0;k<n_position_type;k++)
        {   
         //Assemble in-element quantities: boundary coordinate
         for(unsigned i=0;i<dim_el;i++)
          {
           zeta[i]+=zeta_nodal(j,k,i)*psi(j,k);
          }
        }
      }
     
     
     //Now find the local undeformed metric tensor from the tangent Vectors
     DenseMatrix<double> a(dim_el);
     for(unsigned i=0;i<dim_el;i++)
      {
       for(unsigned j=0;j<dim_el;j++)
        {
         //Initialise surface metric tensor to zero
         a(i,j) = 0.0;
         //Take the dot product
         for(unsigned k=0;k<dim_el+1;k++)
          { 
           a(i,j) += interpolated_a(i,k)*interpolated_a(j,k);
          }
        }
      }

     
     //Find the determinant of the metric tensor
     double adet =0.0;
     switch(dim_el+1)
      {

      case 2:
       adet = a(0,0);
       break;

      case 3:
       adet = a(0,0)*a(1,1) - a(0,1)*a(1,0);
       break;

      default:
       throw 
        OomphLibError(
         "Wrong dimension fill_in_generic_contribution_to_residuals_displ_lagr_multiplier",
         "ImposeDisplacementLagrangeMultiplierElement::fill_in_generic_contribution_to_residuals_displ_lagr_multiplier()",
         OOMPH_EXCEPTION_LOCATION);
      }
     
     // Get prescribed wall shape
     Vector<double> r_prescribed(dim_el+1);
     Boundary_shape_geom_object_pt->position(zeta,r_prescribed);
     
     //Premultiply the weights and the square-root of the determinant of 
     //the metric tensor
     double W = w*sqrt(adet);

     // Assemble residuals and jacobian
     
     //Loop over directions
     for(unsigned i=0;i<dim_el+1;i++)
      {     
       //Loop over the nodes
       for(unsigned j=0;j<n_node;j++)
        {          
         
         // Assemble residual for Lagrange multiplier:
        
         // Local eqn number. Recall that the
         // (additional) Lagrange multiplier values are stored
         // after those that were created by the bulk elements:
         int local_eqn=nodal_local_eqn(j,Nbulk_value[j]+i);
         if (local_eqn>=0)
          {
           residuals[local_eqn]+=(x[i]-r_prescribed[i])*psi(j)*W;

           // Do Jacobian too?
           if (flag==1)
            {
             // Loop over the nodes again for unknowns (only diagonal
             // contribution to direction!).
             for(unsigned jj=0;jj<n_node;jj++)
              {     
               int local_unknown=position_local_eqn(jj,0,i);
               if (local_unknown>=0)
                {
                 jacobian(local_eqn,local_unknown)+=psi(jj)*psi(j)*W;
                }
              }
            }
          }

         
         // Add Lagrange multiplier contribution to bulk equations

         // Local eqn number: Node, type, direction
         local_eqn=position_local_eqn(j,0,i);
         if (local_eqn>=0)
          {
           // Add to residual
           residuals[local_eqn]+=lambda[i]*psi(j)*W;

           // Do Jacobian too?
           if (flag==1)
            {
             // Loop over the nodes again for unknowns (only diagonal
             // contribution to direction!).
             for(unsigned jj=0;jj<n_node;jj++)
              {     
               int local_unknown=nodal_local_eqn(jj,Nbulk_value[jj]+i);
               if (local_unknown>=0)
                {
                 jacobian(local_eqn,local_unknown)+=psi(jj)*psi(j)*W;
                }
              }
            }
          }

        }
      }
   
  
  } //End of loop over the integration points

  }


private:
 
 /// The boundary number in the bulk mesh to which this element is attached
 unsigned Boundary_number_in_bulk_mesh;

#ifdef PARANOID

 /// \short Has the Boundary_number_in_bulk_mesh been set? Only included if
 /// compiled with PARANOID switched on.
 bool Boundary_number_in_bulk_mesh_has_been_set;

#endif

 /// \short GeomObject that specifies the prescribed 
 /// boundary displacement; GeomObject is assumed to be
 /// parametrised by the same coordinate the is used as
 /// the boundary coordinate in the bulk solid mesh to which
 /// this element is attached.
 GeomObject* Boundary_shape_geom_object_pt;
 
 
}; 


} // end of oomph namespace extension


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Physical_Variables
{

 /// Domain height
 double H=1.0;

 /// Domain width
 double L=1.0;

 /// Boundary lline object
 WarpedLine Boundary_geom_object(H,L,0.0);

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt;

 /// Elastic modulus
 double E=1.0;

 /// Poisson's ratio
 double Nu=0.3;
 
} //end namespace



//=============begin_problem============================================ 
/// Problem class for deformation of elastic block by prescribed
/// boundary motion.
//====================================================================== 
template<class ELEMENT>
class PrescribedBoundaryDisplacementProblem : public Problem
{

public:

 /// Constructor:
 PrescribedBoundaryDisplacementProblem();
 
 /// Update function (empty)
 void actions_after_newton_solve() {}

 /// \short Update function (empty)
 void actions_before_newton_solve() {}

 /// Access function for the solid mesh
 ElasticRefineableRectangularQuadMesh<ELEMENT>*& solid_mesh_pt() 
  {return Solid_mesh_pt;} 

 /// Access function to the mesh of Lagrange multiplier elements
 SolidMesh*& lagrange_multiplier_mesh_pt()
  {return Lagrange_multiplier_mesh_pt;} 

 /// Actions before adapt: Wipe the mesh of Lagrange multiplier elements
 void actions_before_adapt();

 /// Actions after adapt: Rebuild the mesh of Lagrange multiplier elements
 void actions_after_adapt();

 /// Doc the solution
 void doc_solution();

private:

 /// \short Pass pointer to prescribed boundary GeomObject
 /// to elements that enforce prescribed boundary motion
 /// by Lagrange multiplilers
 void set_prescribed_boundary_shape();

 /// \short Create elements that enforce prescribed boundary motion
 /// by Lagrange multiplilers
 void create_lagrange_multiplier_elements();

 /// Delete elements that enforce prescribed boundary motion
 /// by Lagrange multiplilers
 void delete_lagrange_multiplier_elements();

 /// Trace file
 ofstream Trace_file;
 
 /// Pointers to node whose position we're tracing
 Node* Trace_node_pt;

 /// Pointer to solid mesh
 ElasticRefineableRectangularQuadMesh<ELEMENT>* Solid_mesh_pt;

 /// Pointers to meshes of Lagrange multiplier elements
 SolidMesh* Lagrange_multiplier_mesh_pt;

 /// DocInfo object for output
 DocInfo Doc_info;
 
};


//===========start_of_constructor======================================= 
/// Constructor: 
//====================================================================== 
template<class ELEMENT>
PrescribedBoundaryDisplacementProblem<ELEMENT>::PrescribedBoundaryDisplacementProblem() 
{

 // Create the mesh

 // # of elements in x-direction
 unsigned n_x=5;

 // # of elements in y-direction
 unsigned n_y=5;

 // Domain length in x-direction
 double l_x= Global_Physical_Variables::L;

 // Domain length in y-direction
 double l_y=Global_Physical_Variables::H;


 //Now create the mesh 
 solid_mesh_pt() = new ElasticRefineableRectangularQuadMesh<ELEMENT>(
  n_x,n_y,l_x,l_y);

 // Set error estimator
 solid_mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 //Assign the physical properties to the elements before any refinement
 //Loop over the elements in the main mesh
 unsigned n_element =solid_mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(i));
   
   // Set the constitutive law
   el_pt->constitutive_law_pt() =
    Global_Physical_Variables::Constitutive_law_pt;
  }


 // Choose a control node: The last node in the solid mesh
 unsigned nnod=solid_mesh_pt()->nnode();
 Trace_node_pt=solid_mesh_pt()->node_pt(nnod-1);

 // Refine the mesh uniformly
 solid_mesh_pt()->refine_uniformly();

 // Construct the mesh of elements that enforce prescribed boundary motion
 /// by Lagrange multiplilers
 Lagrange_multiplier_mesh_pt=new SolidMesh;
 create_lagrange_multiplier_elements();
 
 // Pass pointer to GeomObject that describes the imposed boundary 
 // shape to the elements
 set_prescribed_boundary_shape();
 
 // Solid mesh is first sub-mesh
 add_sub_mesh(solid_mesh_pt());

 // Add Lagrange multiplier sub-mesh
 add_sub_mesh(lagrange_multiplier_mesh_pt());

 // Build combined "global" mesh
 build_global_mesh();
 
 // Loop over all boundaries apart from the top one (2) 
 for (unsigned b=0;b<4;b++)
  {
   if (b!=2)
    {
     unsigned n_side = solid_mesh_pt()->nboundary_node(b);
     
     //Loop over the nodes
     for(unsigned i=0;i<n_side;i++)
      {
       solid_mesh_pt()->boundary_node_pt(b,i)->pin_position(0);
       solid_mesh_pt()->boundary_node_pt(b,i)->pin_position(1);
      }
    }
  }

 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
  solid_mesh_pt()->element_pt());

 //Attach the boundary conditions to the mesh
 cout << assign_eqn_numbers() << std::endl; 

 // Set output directory
 Doc_info.set_directory("RESLT");

 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",Doc_info.directory().c_str());
 Trace_file.open(filename);
 

} //end of constructor


//=====================start_of_actions_before_adapt======================
/// Actions before adapt: Wipe the mesh of elements that impose
/// the prescribed boundary displacements
//========================================================================
template<class ELEMENT>
void PrescribedBoundaryDisplacementProblem<ELEMENT>::actions_before_adapt()
{
 // Kill the  elements and wipe surface mesh
 delete_lagrange_multiplier_elements();
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();

}// end of actions_before_adapt



//=====================start_of_actions_after_adapt=======================
///  Actions after adapt: Rebuild the mesh of elements that impose
/// the prescribed boundary displacements
//========================================================================
template<class ELEMENT>
void PrescribedBoundaryDisplacementProblem<ELEMENT>::actions_after_adapt()
{
 // Create elements from all elements that are 
 // adjacent to boundary 2 and add them to surface meshes
 create_lagrange_multiplier_elements();
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();
 
 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
  solid_mesh_pt()->element_pt());

 // Set pointer to prescribed boundary shape
 set_prescribed_boundary_shape();
 
}// end of actions_after_adapt



//========start_of_set_prescribed_boundary_shape==========================
/// Set boundary information for the Lagrange multiplier FaceElements
//========================================================================
template<class ELEMENT>
void PrescribedBoundaryDisplacementProblem<ELEMENT>::set_prescribed_boundary_shape()
{
 // Loop over the elements in the Lagrange multiplier element mesh
 // for elements on the top boundary (boundary 2)
 unsigned n_element=lagrange_multiplier_mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a Lagrange multiplier element
   ImposeDisplacementLagrangeMultiplierElement<ELEMENT> *el_pt = 
    dynamic_cast<ImposeDisplacementLagrangeMultiplierElement<ELEMENT>*>
    (lagrange_multiplier_mesh_pt()->element_pt(i));

   // Set the GeomObject that defines the boundary shape
   el_pt->boundary_shape_geom_object_pt() = 
    &Global_Physical_Variables::Boundary_geom_object;

   // Which bulk boundary are we attached to (needed to extract
   // the boundary coordinate from the bulk nodes)
   unsigned b=2;
   el_pt->set_boundary_number_in_bulk_mesh(b);
  }
 

 // Pin Lagrange multipliers at the end of the to boundary
 for (unsigned b=1;b<4;b=b+2)
  {
   unsigned n_side = solid_mesh_pt()->nboundary_node(b);
   
   //Loop over the nodes
   for(unsigned i=0;i<n_side;i++)
    {
     Node* nod_pt = solid_mesh_pt()->boundary_node_pt(b,i);
     unsigned nval=nod_pt->nvalue();

     // Are there any values? If so they must be Lagrange multipliers
     // and we pin them.
     if (nval>0)
      {
       for (unsigned j=0;j<nval;j++)
        {
         nod_pt->pin(j);
        }
      }
    }  
  }
 

}// end of set boundary shape fct


 
//============start_of_create_lagrange_multiplier_elements===============
/// Create elements that impose the prescribed boundary displacement
//=======================================================================
template<class ELEMENT>
void PrescribedBoundaryDisplacementProblem<ELEMENT>::
create_lagrange_multiplier_elements()
{
 // Lagrange multiplier elements are located on boundary 2:
 unsigned b=2;

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = solid_mesh_pt()->nboundary_element(b);
 
 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    solid_mesh_pt()->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b
   int face_index = solid_mesh_pt()->face_index_at_boundary(b,e);
      
   // Create new element and add to mesh
   Lagrange_multiplier_mesh_pt->add_element_pt(
    new ImposeDisplacementLagrangeMultiplierElement<ELEMENT>(
     bulk_elem_pt,face_index));   
  }  

 // Pass the pointer to the GeomObject that specifies the
 // boundary shape
 set_prescribed_boundary_shape();
 
} // end of create_lagrange_multiplier_elements




//====start_of_delete_lagrange_multiplier_elements=======================
/// Delete elements that impose the prescribed boundary displacement
/// and wipe the associated mesh
//=======================================================================
template<class ELEMENT>
void PrescribedBoundaryDisplacementProblem<ELEMENT>::delete_lagrange_multiplier_elements()
{
 // How many surface elements are in the surface mesh
 unsigned n_element = Lagrange_multiplier_mesh_pt->nelement();
 
 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill surface element
   delete Lagrange_multiplier_mesh_pt->element_pt(e);
  }
 
 // Wipe the mesh
 Lagrange_multiplier_mesh_pt->flush_element_and_node_storage();

} // end of delete_lagrange_multiplier_elements



//==============start_doc===========================================
/// Doc the solution
//==================================================================
template<class ELEMENT>
void PrescribedBoundaryDisplacementProblem<ELEMENT>::doc_solution()
{

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned n_plot = 5; 


 std::cout << "Docing number " <<  Doc_info.number() << std::endl;

 // Output shape of and stress in deformed body
 //--------------------------------------------
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 solid_mesh_pt()->output(some_file,n_plot);
 some_file.close();

 // Output Lagrange multipliers
 //----------------------------
 sprintf(filename,"%s/lagr%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Lagrange_multiplier_mesh_pt->output(some_file);
 some_file.close();

 // Write trace file: Load/displacement characteristics
 Trace_file << Trace_node_pt->x(0) << " " 
            << Trace_node_pt->x(1) << " " 
            << std::endl;

 // Increment label for output files
 Doc_info.number()++;

} //end doc



//=======start_of_main==================================================
/// Driver code
//======================================================================
int main()
{
 
 // Create generalised Hookean constitutive equations
 Global_Physical_Variables::Constitutive_law_pt = 
  new GeneralisedHookean(Global_Physical_Variables::Nu,
                         Global_Physical_Variables::E);
 
 //Set up the problem
 PrescribedBoundaryDisplacementProblem<RefineableQPVDElement<2,3> > problem;
 

 // Max. number of adaptations per solve
 unsigned max_adapt=1;
 
 problem.doc_solution();

 //Parameter incrementation
 unsigned nstep=2; 
 for(unsigned i=0;i<nstep;i++)
  {

   // Increment imposed boundary displacement
   Global_Physical_Variables::Boundary_geom_object.ampl()+=0.1;

   // Solve the problem with Newton's method, allowing
   // up to max_adapt mesh adaptations after every solve.
   problem.newton_solve(max_adapt);
   
   // Doc solution
   problem.doc_solution();

   // For maximum stability: Reset the current nodal positions to be
   // the "stress-free" ones -- this assignment means that the
   // parameter study no longer corresponds to a physical experiment
   // but is what we'd do if we wanted to use the solid solve
   // to update a fluid mesh in an FSI problem, say.
   problem.solid_mesh_pt()->set_lagrangian_nodal_coordinates();
   
  }
 
} //end of main








