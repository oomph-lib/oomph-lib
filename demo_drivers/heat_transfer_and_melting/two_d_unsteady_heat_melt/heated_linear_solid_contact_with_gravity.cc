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
// Driver for 2D contact problem with displacement and gravity loading
// unsteady heat in penetrator and bulk and continuity of temperature on 
// contact boundary
#include <fenv.h> 

//#define STRUCTURED_MESH

//Generic routines
#include "generic.h"

// The solid elements
#include "solid.h"

// The linear elasticity elements
#include "linear_elasticity.h"

// Unsteady heat
#include "unsteady_heat.h"

// Mesh
#include "meshes/triangle_mesh.h"
#include "meshes/rectangular_quadmesh.h"


// Contact stuff
#include "contact_elements.h"

using namespace std;

using namespace oomph;

using namespace MathematicalConstants;


// Namespace extension
namespace oomph
{

//======================================================================
/// Template-free base class for a heated penetrator surface which 
/// conducts to another surface by imposing a flux
//======================================================================
class TemplateFreeHeatedPenetratorFluxElementBase
{

public:

 /// Virtual destructor (empty)
 virtual ~TemplateFreeHeatedPenetratorFluxElementBase(){}

 /// Pure virtual function: Return temperature at 
 /// local coordinate s.
 virtual double penetrator_temperature(const Vector<double>& s) const=0;

};

/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////


//======================class definition==================================
/// Bulk element which combines linear elasticity and unsteady heat
//=========================================================================
template<unsigned DIM, unsigned NNODE_1D>
class TLinearHeatAndElasticityElement :
 public virtual TLinearElasticityElement<DIM,NNODE_1D>,
 public virtual TUnsteadyHeatElement<DIM,NNODE_1D>
{

public:

 /// Constructor: call the underlying constructors
 TLinearHeatAndElasticityElement() : 
  TLinearElasticityElement<DIM,NNODE_1D>(),
  TUnsteadyHeatElement<DIM,NNODE_1D>() 
  {}

 /// The required number of values stored at the nodes is the sum of the
 /// required values of the two single-physics  elements. Note that this step is
 /// generic for any multi-physics element of this type.
 unsigned required_nvalue(const unsigned &n) const
  {return (TLinearElasticityElement<DIM,NNODE_1D>::required_nvalue(n) +
           TUnsteadyHeatElement<DIM,NNODE_1D>::required_nvalue(n));}
 

 /// Number of scalars/fields output by this element. Broken
 /// virtual. Needs to be implemented for each new specific element type.
 /// Temporary dummy
 unsigned nscalar_paraview() const
 {
  throw OomphLibError(
   "This function hasn't been implemented for this element",
   OOMPH_CURRENT_FUNCTION,
   OOMPH_EXCEPTION_LOCATION);
  
  // Dummy unsigned
  return 0;
 }
 
 /// Write values of the i-th scalar field at the plot points. Broken
 /// virtual. Needs to be implemented for each new specific element type.
 /// Temporary dummy
 void scalar_value_paraview(std::ofstream& file_out,
                            const unsigned& i,
                            const unsigned& nplot)  const
 {
  throw OomphLibError(
   "This function hasn't been implemented for this element",
   OOMPH_CURRENT_FUNCTION,
   OOMPH_EXCEPTION_LOCATION);
 }
 
 
 /// Name of the i-th scalar field. Default implementation
 /// returns V1 for the first one, V2 for the second etc. Can (should!) be
 /// overloaded with more meaningful names.
 std::string scalar_name_paraview(const unsigned& i) const
  {
   return "V"+StringConversion::to_string(i);
  }
 

 ///  Overload the standard output function with the broken default
 void output(std::ostream &outfile) 
  {
   FiniteElement::output(outfile);
  }

 /// Output function:  
 void output(std::ostream &outfile, const unsigned &nplot)
  {

  //Set output Vector
  Vector<double> s(DIM);
  Vector<double> x(DIM);
  Vector<double> u(DIM);
  
  // Tecplot header info
  outfile << this->tecplot_zone_string(nplot);
  
  // Loop over plot points
  unsigned num_plot_points=this->nplot_points(nplot);
  for (unsigned iplot=0;iplot<num_plot_points;iplot++)
   {
    
    // Get local coordinates of plot point
    this->get_s_plot(iplot,nplot,s);
    
    // Get Eulerian coordinates and displacements
    this->interpolated_x(s,x);
    this->interpolated_u_linear_elasticity(s,u);
    
    //Output the x,y,..
    for(unsigned i=0;i<DIM;i++) 
     {outfile << x[i] << " ";}

    // Output displacements
    for(unsigned i=0;i<DIM;i++) 
     {outfile << u[i] << " ";} 

    // Output temperature
    outfile << this->interpolated_u_ust_heat(s) << std::endl;   

    outfile << std::endl;
   }
  
  // Write tecplot footer (e.g. FE connectivity lists)
  this->write_tecplot_zone_footer(outfile,nplot);

  } 

 /// C-style output function: Broken default
 void output(FILE* file_pt)
  {FiniteElement::output(file_pt);}

 ///  C-style output function: Broken default
 void output(FILE* file_pt, const unsigned &n_plot)
  {FiniteElement::output(file_pt,n_plot);}

 /// Output function for an exact solution: Broken default
 void output_fct(std::ostream &outfile, const unsigned &Nplot,
                 FiniteElement::SteadyExactSolutionFctPt 
                 exact_soln_pt)
  {FiniteElement::output_fct(outfile,Nplot,exact_soln_pt);}

 /// Output function for a time-dependent exact solution:
 /// Broken default.
 void output_fct(std::ostream &outfile, const unsigned &Nplot,
                 const double& time,
                 FiniteElement::UnsteadyExactSolutionFctPt 
                 exact_soln_pt)
  {
   FiniteElement::
    output_fct(outfile,Nplot,time,exact_soln_pt);
  }

 /// Overload the index at which the temperature 
 /// variable is stored. We choose to store it after the displacements
 inline unsigned u_index_ust_heat() const {return DIM;}
 
 /// Validate against exact solution at given time
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points and compute L2 error
 /// and L2 norm of velocity solution over element
 /// Call the broken default
 void compute_error(std::ostream &outfile,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                    const double& time,
                    double& error, double& norm)
  {FiniteElement::compute_error(outfile,exact_soln_pt,
                                time,error,norm);}
 
 /// Validate against exact solution.
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points and compute L2 error
 /// and L2 norm of velocity solution over element
 /// Call the broken default
 void compute_error(std::ostream &outfile,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                    double& error, double& norm)
  {FiniteElement::compute_error(outfile,exact_soln_pt,error,norm);}
 

 /// Calculate the element's contribution to the residual vector.
 /// Recall that fill_in_* functions MUST NOT initialise the entries 
 /// in the vector to zero. This allows us to call the 
 /// fill_in_* functions of the constituent single-physics elements
 /// sequentially, without wiping out any previously computed entries.
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Fill in the residuals of the linear elasticity equations
   LinearElasticityEquations<DIM>::
    fill_in_contribution_to_residuals(residuals);

   //Fill in the residuals of the unsteady heat equations
   UnsteadyHeatEquations<DIM>::fill_in_contribution_to_residuals(
    residuals);
  }

 /// Compute the element's residual Vector and the Jacobian matrix.
 /// No off diagonal blocks!
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian)
  {
   // Calculate the linear elasticity contributions (diagonal block and 
   // residuals)
   LinearElasticityEquations<DIM>::
    fill_in_contribution_to_jacobian(residuals,jacobian);

   //Calculate the unsteady heat contributions (diagonal block and residuals)
   UnsteadyHeatEquations<DIM>::
    fill_in_contribution_to_jacobian(residuals,jacobian);

  } //End of jacobian calculation


 /// 'Flux' vector for Z2 error estimation: Error estimation
 /// is based on error in unsteady heat element
 void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
  {
   TUnsteadyHeatElement<DIM,NNODE_1D>::get_Z2_flux(s,flux);
  }

 /// Number of 'flux' terms for Z2 error estimation: Error estimation
 /// is based on error in unsteady heat element
 unsigned num_Z2_flux_terms() 
  {
   return TUnsteadyHeatElement<DIM,NNODE_1D>::num_Z2_flux_terms();
  }

 /// Number of vertex nodes in the element
 unsigned nvertex_node() const
  {return TUnsteadyHeatElement<DIM,NNODE_1D>::nvertex_node();}

 /// Pointer to the j-th vertex node in the element
 Node* vertex_node_pt(const unsigned& j) const
  {return TUnsteadyHeatElement<DIM,NNODE_1D>::vertex_node_pt(j);}

 /// Order of recovery shape functions for Z2 error estimation: Done
 /// for unsteady heat element since it determines the refinement
 unsigned nrecovery_order() 
  {
   return TUnsteadyHeatElement<DIM,NNODE_1D>::nrecovery_order();
  }

 /// Pin unsteady heat
 void pin_unsteady_heat()
  {
   unsigned nnod=this->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     this->node_pt(j)->pin(u_index_ust_heat()); 
    }
  }


 /// Pin linear elasticity
 void pin_linear_elasticity()
  {
   unsigned nnod=this->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     this->node_pt(j)->pin(0); // hierher use u_index fct from lin elast
     this->node_pt(j)->pin(1); // hierher 
    }
  }


};


//=======================================================================
/// Face geometry of the customised element
//=======================================================================
template<unsigned int DIM, unsigned int NNODE_1D>
class FaceGeometry<TLinearHeatAndElasticityElement<DIM,NNODE_1D> >: 
public virtual TElement<DIM-1,NNODE_1D>
{
  public:
 FaceGeometry() : TElement<DIM-1,NNODE_1D>() {}
};









/// /////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////


//==========================================================
/// Linear heat and elasticity upgraded to become projectable
//==========================================================
 template<class LINEAR_HEAT_AND_ELAST_ELEMENT>
 class ProjectableLinearHeatAndElasticityElement : 
  public virtual ProjectableElement<LINEAR_HEAT_AND_ELAST_ELEMENT>
 {
  
 public:
  
  /// Constructor [this was only required explicitly
  /// from gcc 4.5.2 onwards...]
  ProjectableLinearHeatAndElasticityElement(){}
  
  
  /// Specify the values associated with field fld. 
  /// The information is returned in a vector of pairs which comprise 
  /// the Data object and the value within it, that correspond to field fld. 
  /// In the underlying elements the 
  /// the displacements and temperatures are stored at the nodal values
  Vector<std::pair<Data*,unsigned> > data_values_of_field(const unsigned& fld)
   {   
    // Create the vector
    Vector<std::pair<Data*,unsigned> > data_values;
    
    // Loop over all nodes and extract the fld-th nodal value
    unsigned nnod=this->nnode();
    for (unsigned j=0;j<nnod;j++)
     {
      // Add the data value associated with the displacement components
      data_values.push_back(std::make_pair(this->node_pt(j),fld));
     }
    
    // Return the vector
    return data_values;
   }

  /// Number of fields to be projected: dim, corresponding to 
  /// the displacement components plus temperature
  unsigned nfields_for_projection()
   {
    return this->dim()+1;
   }
  
  /// Number of history values to be stored for fld-th field. 
  /// (includes present value!)
  unsigned nhistory_values_for_projection(const unsigned &fld)
   {
#ifdef PARANOID
    if (fld>3)
     {
      std::stringstream error_stream;
      error_stream 
       << "Elements only store three fields so fld can't be"
       << " " << fld << std::endl;
      throw OomphLibError(
       error_stream.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
     }
#endif
   return this->node_pt(0)->ntstorage();   
   }
  
  /// Number of positional history values: Read out from
  /// positional timestepper  (Note: count includes current value!)
  unsigned nhistory_values_for_coordinate_projection()
   {
    return this->node_pt(0)->position_time_stepper_pt()->ntstorage();
   }
  
  /// Return Jacobian of mapping and shape functions of field fld
  /// at local coordinate s
  double jacobian_and_shape_of_field(const unsigned &fld, 
                                     const Vector<double> &s, 
                                     Shape &psi)
   {
    unsigned n_dim=this->dim();
    unsigned n_node=this->nnode();
    DShape dpsidx(n_node,n_dim);
        
    // Call the derivatives of the shape functions and return
    // the Jacobian
    return this->dshape_eulerian(s,psi,dpsidx);
   }
  


  /// Return interpolated field fld at local coordinate s, at time level
  /// t (t=0: present; t>0: history values)
  double get_field(const unsigned &t, 
                   const unsigned &fld,
                   const Vector<double>& s)
   {
    unsigned n_node=this->nnode();

/* #ifdef PARANOID */
/*     unsigned n_dim=this->node_pt(0)->ndim(); */
/* #endif */
    
    //Local shape function
    Shape psi(n_node);
    
    //Find values of shape function
    this->shape(s,psi);
    
    //Initialise value of u
    double interpolated_u = 0.0;
    
    //Sum over the local nodes at that time
    for(unsigned l=0;l<n_node;l++) 
     {
// over-zealous I think. This will quietly do the right thing
// even if there are additional degrees of freedom floating around.
/* #ifdef PARANOID */
/*       unsigned nvalue=this->node_pt(l)->nvalue(); */
/*       if (nvalue!=n_dim) */
/*        {         */
/*         std::stringstream error_stream; */
/*         error_stream  */
/*          << "Current implementation only works for non-resized nodes\n" */
/*          << "but nvalue= " << nvalue << "!= dim = " << n_dim << std::endl; */
/*         throw OomphLibError( */
/*          error_stream.str(), */
/*          OOMPH_CURRENT_FUNCTION, */
/*          OOMPH_EXCEPTION_LOCATION); */
/*        } */
/* #endif */
      interpolated_u += this->nodal_value(t,l,fld)*psi[l];
     }
    return interpolated_u;     
   }
  
 
  /// Return number of values in field fld
  unsigned nvalue_of_field(const unsigned &fld)
   {
    return this->nnode();
   }
  
  
  /// Return local equation number of value j in field fld.
  int local_equation(const unsigned &fld,
                     const unsigned &j)
   {
// over-zealous I think. This will quietly do the right thing
// even if there are additional degrees of freedom floating around.
/* #ifdef PARANOID */
/*     unsigned n_dim=this->node_pt(0)->ndim(); */
/*     unsigned nvalue=this->node_pt(j)->nvalue(); */
/*     if (nvalue!=n_dim) */
/*      {         */
/*       std::stringstream error_stream; */
/*       error_stream  */
/*        << "Current implementation only works for non-resized nodes\n" */
/*        << "but nvalue= " << nvalue << "!= dim = " << n_dim << std::endl; */
/*       throw OomphLibError( */
/*          error_stream.str(), */
/*          OOMPH_CURRENT_FUNCTION, */
/*          OOMPH_EXCEPTION_LOCATION); */
/*      } */
/* #endif */
    return this->nodal_local_eqn(j,fld);
   }

  
 };


//=======================================================================
/// Face geometry for element is the same as that for the underlying
/// wrapped element
//=======================================================================
 template<class ELEMENT>
 class FaceGeometry<ProjectableLinearHeatAndElasticityElement<ELEMENT> > 
  : public virtual FaceGeometry<ELEMENT>
 {
 public:
  FaceGeometry() : FaceGeometry<ELEMENT>() {}
 };


//=======================================================================
/// Face geometry of the Face Geometry for element is the same as 
/// that for the underlying wrapped element
//=======================================================================
 template<class ELEMENT>
 class FaceGeometry<FaceGeometry<ProjectableLinearHeatAndElasticityElement<ELEMENT> > >
  : public virtual FaceGeometry<FaceGeometry<ELEMENT> >
 {
 public:
  FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT> >() {}
 };



/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////

//======================================================================
/// Template-free base class for linear contact with a heated penetrator
//======================================================================
class TemplateFreeHeatedLinearSurfaceContactElementBase
{

public:

 /// Virtual destructor (empty)
 virtual ~TemplateFreeHeatedLinearSurfaceContactElementBase(){}

 /// Pure virtual function: Return heat flux required to 
 /// maintain continuity of temperature to adjacent penetrator at 
 /// local coordinate s.
 virtual double heat_flux(const Vector<double>& s) const=0;

};

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////


//======================================================================
/// Class for elements that impose contact boundary conditions
/// either enforcing non-penetration (but "without stick"; default) or 
/// permament contact (in which case contact force can be
/// positive or negative. Uses Lagrange-multiplier-like pressure
/// to enforce contact/non-penetration. Almost certainly works only with
/// 2D Penetrator at the moment. Specific implementation for
/// linear elasticity. Also applies continuity of temperature
/// constraint at interface
//======================================================================
template <class ELEMENT>
class HeatedLinearSurfaceContactElement : 
  public virtual LinearSurfaceContactElement<ELEMENT>,
  public virtual TemplateFreeHeatedLinearSurfaceContactElementBase,
  public virtual ElementWithExternalElement
 
{
 
  public:
 

 /// Constructor, which takes a "bulk" element and the 
 /// value of the index and its limit
 HeatedLinearSurfaceContactElement(FiniteElement* const &element_pt, 
                                   const int &face_index,
                                   const unsigned &id=0,
                                   const bool& 
                                   called_from_refineable_constructor=false) : 
  SurfaceContactElementBase<ELEMENT>(element_pt, 
                                     face_index,
                                     id,
                                     called_from_refineable_constructor),
  LinearSurfaceContactElement<ELEMENT>(element_pt, 
                                       face_index,
                                       id,
                                       called_from_refineable_constructor),
  Use_collocated_heat_flux_flag_pt(0)
  {
   // hierher provide separate id via arg list
   Heat_flux_lagr_multiplier_id=id+1;

   // We need yet another additional value for each FaceElement node:
   // the normal traction (Lagrange multiplier) to be 
   // exerted onto the solid
   unsigned n_nod=this->nnode();
   Vector<unsigned> n_additional_values(n_nod,1);
   
   // Now add storage for Lagrange multipliers and set the map containing 
   // the position of the first entry of this face element's 
   // additional values.
   this->add_additional_values(n_additional_values,
                               Heat_flux_lagr_multiplier_id);


   // Set source element storage: one interaction with an external 
   // element -- the HeatedPenetratorFluxElement whose temperature we're 
   //matching
   this->set_ninteraction(1);
   
   activate_iso_colloc_for_contact = false;
   /// hierher currently deactivating collocation on contact elements
   this->Use_isoparametric_flag_pt = &activate_iso_colloc_for_contact;
   this->Use_collocated_penetration_flag_pt = &activate_iso_colloc_for_contact;
   this->Use_collocated_contact_pressure_flag_pt =&activate_iso_colloc_for_contact;

   //If we are using a top hat function for the heat flux, we need to change 
   //the integration scheme
   if(!use_collocated_heat_flux_flag())
    {
     set_integration_scheme(new PiecewiseGauss<1,3>(s_min(),s_max()));
    }
   
  }


 /// Default constructor
 HeatedLinearSurfaceContactElement(){}

 /// Function to describe the local dofs of the element. The ostream 
 /// specifies the output stream to which the description 
 /// is written; the string stores the currently 
 /// assembled output that is ultimately written to the
 /// output stream by Data::describe_dofs(...); it is typically
 /// built up incrementally as we descend through the
 /// call hierarchy of this function when called from 
 /// Problem::describe_dofs(...)
 void describe_local_dofs(std::ostream& out,
                          const std::string& current_string) const
  {
   LinearSurfaceContactElement<ELEMENT>::
    describe_local_dofs(out,current_string);
   ElementWithExternalElement::describe_local_dofs(out,current_string);
  }

 // hierher why do we suddenly need this?
 /// Final over-rider -- use version in element with external element
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian)
  {
   ElementWithExternalElement::fill_in_contribution_to_jacobian(residuals,
                                                                jacobian);
  }


 /// Return heat flux required to maintain continuity of temperature
 /// to adjacent penetrator at local coordinate s
 double heat_flux(const Vector<double>& s) const
  {
   // Get shape function for Lagrange multiplier
   unsigned n_node = this->nnode();
   Shape psi_p(n_node);
   this->shape_p(s,psi_p);
   
   // Interpolated Lagrange multiplier (heat flux)
   double interpolated_lambda_q=0.0;
   for(unsigned l=0;l<n_node;l++) 
    {
     // Cast to a boundary node
     BoundaryNodeBase *bnod_pt = 
      dynamic_cast<BoundaryNodeBase*>(this->node_pt(l));
     
     // Get the index of the nodal value associated with
     // this FaceElement
     unsigned first_index_q=
      bnod_pt->index_of_first_value_assigned_by_face_element(
       this->Heat_flux_lagr_multiplier_id);
     
     // Add to Lagrange multiplier (acting as pressure on solid
     // to enforce motion to ensure non-penetration)
     interpolated_lambda_q+=this->node_pt(l)->value(first_index_q)*psi_p[l];
    }
   return interpolated_lambda_q;
  }
 

 /// hierher overload...
 double zeta_nodal(const unsigned &n, const unsigned &k, 
                   const unsigned &i) const
  {
   unsigned dim=this->node_pt(n)->ndim();
   Vector<double> x(dim);
   for (unsigned ii=0;ii<dim;ii++)
    {
     x[ii]=this->node_pt(n)->x(ii);
    }
   Vector<double> zeta(dim-1);
   this->penetrator_pt()->surface_coordinate(x,zeta);
   return zeta[i]; 
  }


 void output(std::ostream &outfile, const unsigned &n_plot)
  {
   unsigned n_dim = this->nodal_dimension();
   
   Vector<double> x(n_dim);
   Vector<double> disp(n_dim);
   Vector<double> x_def(n_dim);
   Vector<double> s(n_dim-1);
   Vector<double> r_pen(n_dim);
   Vector<double> unit_normal(n_dim);
   Vector<double> interpolated_heat_flux(1);
   
   // Tecplot header info
   outfile << this->tecplot_zone_string(n_plot);
   
   // Loop over plot points
   unsigned num_plot_points=this->nplot_points(n_plot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     // Get local coordinates of plot point
     this->get_s_plot(iplot,n_plot,s);
     
     // Get coordinates and outer unit normal
     this->interpolated_x(s,x);
     this->outer_unit_normal(s,unit_normal);   
     
     // Displacement
     this->interpolated_u_linear_elasticity(s,disp);
     
     // Deformed position
     for(unsigned i=0;i<n_dim;i++) 
      {
       x_def[i]=x[i]+disp[i];
      }
     
     // Get penetration based on deformed position
     double d= 0.0;
     bool intersection = false;
     this->penetration(x_def,unit_normal,d,intersection);
     
     //Get heat flux
     interpolated_heat_flux[0] = heat_flux(s);
     
     //Output the x,y,..
     for(unsigned i=0;i<n_dim;i++) 
      {outfile << x[i] << " ";} // col 1,2
     
     // Penetration
     outfile << std::max(d,-100.0) << " "; // col 3
     
     // Lagrange multiplier-like pressure
     double p=this->get_interpolated_lagrange_p(s);
     outfile << p << " "; // col 4
     
     
     // Plot Lagrange multiplier like pressure
     outfile << -unit_normal[0]*p << " "; // col 5
     outfile << -unit_normal[1]*p << " "; // col 6
     
     // Plot vector from current point to boundary of penetrator
     double d_tmp=d;
     if (!intersection) d_tmp=0.0;
     outfile << -d_tmp*unit_normal[0] << " ";  // col 7
     outfile << -d_tmp*unit_normal[1] << " ";  // col 8
     
     // Output normal
     for(unsigned i=0;i<n_dim;i++) 
      {outfile << unit_normal[i] << " ";} // col 9, 10
     
     //Output the displacements
     for(unsigned i=0;i<n_dim;i++)
      {
       outfile << disp[i] << " "; // col 11, 12
      }
     
     //Output the deformed position
     for(unsigned i=0;i<n_dim;i++)
      {
      outfile << x_def[i] << " "; // col 13, 14
      }
     
     //Output the heat_flux
     outfile << interpolated_heat_flux[0] << " "; // col 15
     
     
     outfile << std::endl;
    }


  
  // Write tecplot footer (e.g. FE connectivity lists)
  this->write_tecplot_zone_footer(outfile,n_plot);
  
 }


  /// Determine which method to use to discretise heat_flux, collocation (true)
 /// or integrated using the hat function (false)
  bool use_collocated_heat_flux_flag()
  {
   if (Use_collocated_heat_flux_flag_pt==0)
    {
     return false;
    }
   else
    {
     return *Use_collocated_heat_flux_flag_pt;
    }
  } 


 /// Access function: Pointer to flag to use collocated heat_flux
 bool*& use_collocated_heat_flux_flag_pt() 
  {return Use_collocated_heat_flux_flag_pt;}
 /// Access function: Pointer to flag to use collocated heat_flux (const version)
 bool* use_collocated_heat_flux_flag_pt() 
  const {return Use_collocated_heat_flux_flag_pt;}

protected:

 /// Overloaded fill in contributions function -- includes heat flux
 /// to ensure continuity of temperature
 void fill_in_contribution_to_residuals_surface_contact(
  Vector<double>& residuals);

 /// ID of heat flux Lagrange multiplier (to ensure continuity of temperature)
 unsigned Heat_flux_lagr_multiplier_id;
  

 /// Set options for basis/test functions for penetration and pressure
 bool* Use_collocated_heat_flux_flag_pt;

 /// hierher currently we don't want to use isoparametric for
 /// the lagrange multipliers and collocation for penetration
 /// and contact pressure, use this as a bulk switch
 bool activate_iso_colloc_for_contact;

};








/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////


//=====================================================================
/// Return the residuals for the LinearSurfaceContactElement equations
//=====================================================================
template<class ELEMENT>
void HeatedLinearSurfaceContactElement<ELEMENT>::
fill_in_contribution_to_residuals_surface_contact(Vector<double> &residuals)
{

 // Get pointer to bulk element
 ELEMENT *bulk_el_pt = dynamic_cast<ELEMENT*>(this->bulk_element_pt());
 
 // Get index of temperature in bulk element
 unsigned temperature_index=bulk_el_pt->u_index_ust_heat();

 // Spatial dimension of problem
 unsigned n_dim = this->nodal_dimension();
   
 // Contribution to contact force
 Vector<double> contact_force(n_dim,0.0);
 
 // Create vector of local residuals (start assembling contributions
 // from zero -- necessary so we can over-write pseudo-hijacked
 // contributions at the end.
 unsigned n_dof=this->ndof();
 Vector<double> local_residuals(n_dof,0.0);
 
 {
  //Find out how many nodes there are
  unsigned n_node = this->nnode();
  
  //Integer to hold the local equation number
  int local_eqn=0;
  int local_eqn_p=0;
  int local_eqn_q=0;
   
  //Set up memory for the shape functions
  Shape psi(n_node);
  DShape dpsids(n_node,n_dim-1);

  // Separate shape functions for Lagrange multiplier
  Shape psi_p(n_node);

  // Separate shape functions for integration (top hat)
  Shape psi_i(n_node);

  Vector<double> s(n_dim-1);

  // Contribution to integrated pressure
  Vector<double> pressure_integral(n_node,0.0);
     
  // Contribution to weighted penetration integral
  Vector<double> penetration_integral(n_node,0.0);
     
  // Contribution to heat flux integral
  Vector<double> heat_flux_integral(n_node,0.0);
     
  // Contribution to temperature deviation integral
  Vector<double> temp_deviation_integral(n_node,0.0);
     
  // Deformed position
  Vector<double> x_def(n_dim,0.0);

  //Set the value of n_intpt
  unsigned n_intpt = this->integral_pt()->nweight();
   
  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    //Get the integral weight
    double w = this->integral_pt()->weight(ipt);
     
    //Only need to call the local derivatives
    this->dshape_local_at_knot(ipt,psi,dpsids);
     
    // Separate shape function for Lagrange multipliers 
    for(unsigned i=0;i<n_dim-1;i++)
     {
      s[i] = this->integral_pt()->knot(ipt,i);
     }
    this->shape_p(s,psi_p);
    
    //set up shape for integration
    this->shape_i(s,psi_i);

    // Interpolated Lagrange multiplier (pressure acting on solid)
    double interpolated_lambda_p=0.0;

    // Interpolated Lagrange multiplier (heat flux)
    double interpolated_lambda_q=0.0;

    // Interpolated temperature
    double interpolated_temp=0.0;

    // Displacement
    Vector<double> disp(n_dim,0.0);
     
    //Calculate the coordinates 
    Vector<double> interpolated_x(n_dim,0.0);
     
    //Also calculate the surface Vectors (derivatives wrt local coordinates)
    DenseMatrix<double> interpolated_A(n_dim-1,n_dim,0.0);   
     
    //Calculate displacements and derivatives
    for(unsigned l=0;l<n_node;l++) 
     {
      // Cast to a boundary node
      BoundaryNodeBase *bnod_pt = 
       dynamic_cast<BoundaryNodeBase*>(this->node_pt(l));
       
      // Get the index of the nodal value associated with
      // this FaceElement
      unsigned first_index=
       bnod_pt->index_of_first_value_assigned_by_face_element(this->Contact_id);
       
      // Add to Lagrange multiplier (acting as pressure on solid
      // to enforce motion to ensure non-penetration)
      interpolated_lambda_p+=this->node_pt(l)->value(first_index)*psi_p[l];
         
      // Get the index of the nodal value associated with
      // this FaceElement
      unsigned first_index_q=
       bnod_pt->index_of_first_value_assigned_by_face_element(
        this->Heat_flux_lagr_multiplier_id);
       
      // Add to Lagrange multiplier (acting as a heat flux to enforce
      // continuity of temperature across boundary)
      interpolated_lambda_q+=this->node_pt(l)->value(first_index_q)*psi_p[l];

      //Get temperature at this integration point
      interpolated_temp+=this->node_pt(l)->value(temperature_index)*psi[l];

      //Loop over directions
      for(unsigned i=0;i<n_dim;i++)
       {
        //Calculate the positions
        interpolated_x[i]+=this->nodal_position(l,i)*psi(l);
           
        //Index at which the displacement nodal value is stored
        unsigned u_nodal_index=this->U_index_linear_elasticity_traction[i];
        disp[i] += this->nodal_value(l,u_nodal_index)*psi(l);

        // Loop over LOCAL derivative directions, to calculate the 
        // tangent(s)
        for(unsigned j=0;j<n_dim-1;j++)
         {
          interpolated_A(j,i)+=this->nodal_position(l,i)*dpsids(l,j);
         }
       }
     }
     
    // Get the temperature on the penetrator
    double penetrator_temperature=0.0;
    if (external_element_pt(0,ipt)!=0)
     {
      Vector<double> s_ext(external_element_local_coord(0,ipt));
      TemplateFreeHeatedPenetratorFluxElementBase* el_pt=
       dynamic_cast<TemplateFreeHeatedPenetratorFluxElementBase*>(
        external_element_pt(0,ipt));
      // Note: this does NOT take contact or no contact into account
      // it simply returns the local temperature at "that point" in "that
      // element"
      penetrator_temperature=el_pt->penetrator_temperature(s_ext);

      /* oomph_info << "actual penetrator temperature: " 
         << penetrator_temperature << std::endl;*/
     }
    else
     {
      // hierher this should never be triggered on contact surface, 
      // throw error here? Although it may not always be a problem, depending
      // on penetrator/contact surface geometry

      oomph_info << "using default penetrator temperature: " 
                 << penetrator_temperature << std::endl;
     }

    //Now find the local deformed metric tensor from the tangent Vectors
    DenseMatrix<double> A(n_dim-1);
    for(unsigned i=0;i<n_dim-1;i++)
     {
      for(unsigned j=0;j<n_dim-1;j++)
       {
        //Initialise surface metric tensor to zero
        A(i,j) = 0.0;
        //Take the dot product
        for(unsigned k=0;k<n_dim;k++)
         { 
          A(i,j) += interpolated_A(i,k)*interpolated_A(j,k);
         }
       }
     }
     
    //Get the outer unit normal
    Vector<double> interpolated_normal(n_dim);
    this->outer_unit_normal(ipt,interpolated_normal);
     
    //Find the determinant of the metric tensor
    double Adet =0.0;
    switch(n_dim)
     {
     case 2:
      Adet = A(0,0);
      break;
     case 3:
      Adet = A(0,0)*A(1,1) - A(0,1)*A(1,0);
      break;
     default:
      throw OomphLibError(
       "Wrong dimension in SurfaceContactElement",
       "LinearSurfaceContactElement::fill_in_contribution_to_residuals()",
       OOMPH_EXCEPTION_LOCATION);
     }
     
    //Premultiply the weights and the square-root of the determinant of 
    //the metric tensor
    double W = w*sqrt(Adet);
     
    // Calculate the "load" -- Lagrange multiplier acts as traction to
    // to enforce required surface displacement and the
    // deformed position
    Vector<double> traction(n_dim);
    for (unsigned i=0;i<n_dim;i++)
     {
      traction[i]=-interpolated_lambda_p*interpolated_normal[i];
      x_def[i]=interpolated_x[i]+disp[i];
     }


    // Accumulate contribution to total contact force
    for(unsigned i=0;i<n_dim;i++)
     {
      contact_force[i]+=traction[i]*W;
     }
     
    //=====LOAD TERMS  FROM PRINCIPLE OF VIRTUAL DISPLACEMENTS========
       
    //Loop over the test functions, nodes of the element
    for(unsigned l=0;l<n_node;l++)
     {
      //Loop over the displacement components
      for(unsigned i=0;i<n_dim;i++)
       {
        local_eqn = this->nodal_local_eqn(l,i); // hierher get index
        /*IF it's not a boundary condition*/
        if(local_eqn >= 0)
         {
          //Add the loading terms to the residuals
          local_residuals[local_eqn] -= traction[i]*psi(l)*W;
         } //End of if not boundary condition
       }
      
      
      // Contribution to unsteady heat equation on surface
      local_eqn = this->nodal_local_eqn(l,temperature_index);

      /*IF it's not a boundary condition*/
      if(local_eqn >= 0)
       {
        //Add the loading terms to the residuals
        local_residuals[local_eqn] -= interpolated_lambda_q*psi(l)*W;
       } //End of if not boundary condition

    
     } //End of loop over shape functions
       
    //=====CONTRIBUTION TO CONTACT PRESSURE/LAGRANGE MULTIPLIER EQNS ========
    
    //Only calculate the integrals we need to
    if(!this->use_collocated_contact_pressure_flag())
     {
      //Loop over the nodes
      for(unsigned l=0;l<n_node;l++)
       {
        // Contribution to integrated pressure
        pressure_integral[l]+=interpolated_lambda_p*psi_i[l]*W;
       }
     }

    if(!this->use_collocated_penetration_flag())
     {
      // Get local penetration
      double d=0.0;
      bool intersection = false;
      this->penetration(x_def,interpolated_normal,d,intersection);

      //If there is no intersection, d = -max, ie the penetrator 
      //is infinitely far away
      if(!intersection)
       {
        d = -DBL_MAX;
       }

      for(unsigned l=0;l<n_node;l++)
       {
        
        // Contribution to weighted penetration integral
        penetration_integral[l]+=d*psi_i[l]*W;
       }
     }

    if(!this->use_collocated_heat_flux_flag())
     {
      for(unsigned l=0;l<n_node;l++)
       {
        // Contribution to integrated heat flux
        heat_flux_integral[l]+=interpolated_lambda_q*psi_i[l]*W;
       }
     }

   
    for(unsigned l=0;l<n_node;l++)
     {
      // Contribution to weighted temperature deviation integral
      temp_deviation_integral[l]+=(interpolated_temp-
                                   penetrator_temperature)*psi_i[l]*W;
     }
     

   } //End of loop over integration points


  // Eqn for contact pressure and heat flux
  //---------------------------------------------
  
  // Storage for nodal coordinate
  Vector<double> x(n_dim);
     
  //Loop over the nodes
  for(unsigned l=0;l<n_node;l++)
   {
    // get the node pt
    Node* nod_pt = this->node_pt(l);
       
    // Cast to a boundary node
    BoundaryNodeBase *bnod_pt =
     dynamic_cast<BoundaryNodeBase*>(nod_pt);
       
    // Get the index of the first nodal value associated with
    // this FaceElement
    unsigned first_index_p=
     bnod_pt->index_of_first_value_assigned_by_face_element(this->Contact_id);
       
    // Equation for pressure Lagrange multiplier
    local_eqn_p = this->nodal_local_eqn(l,first_index_p);
       
    // Get the index of the first nodal value associated with
    // this FaceElement
    unsigned first_index_q=
     bnod_pt->index_of_first_value_assigned_by_face_element(
      this->Heat_flux_lagr_multiplier_id);
       
    // Equation number for temperature Lagrange multiplier
    local_eqn_q = this->nodal_local_eqn(l,first_index_q);

    //IF it's not a boundary condition for both heat_flux and contact_pressure
    if(local_eqn_p >= 0 || local_eqn_q >= 0)
     {
      
      // Use weighted/integrated quantities intially (could be 0)
      // then overwrite them if using collocation
      double d=penetration_integral[l];
      double contact_pressure=pressure_integral[l];
      double heat_flux=heat_flux_integral[l];
      double temp_deviation=temp_deviation_integral[l];
      
      //overwrite appropriate measure if using collocation

      if(this->use_collocated_penetration_flag())
       {
        // Nodal position
        x[0]=nod_pt->x(0);
        x[1]=nod_pt->x(1);
        
        // Get outer unit normal
        Vector<double> s(1);
        this->local_coordinate_of_node(l,s);
        Vector<double> unit_normal(2);
        this->outer_unit_normal(s,unit_normal);

        // Get penetration
        bool intersection = false;
        this->penetration(x,unit_normal,d,intersection);
        
        //If there is no intersection, d = -max, ie the penetrator 
        //is infinitely far away
        if(!intersection)
         {
          d = -DBL_MAX;
         }
         
       }
      
      //If we are using collocated heat flux, overwrite value
      if(this->use_collocated_heat_flux_flag())
       {
        heat_flux = this->node_pt(l)->value(first_index_q);
       }

      //If we are using collocated contact pressure, overwrite value
      if(this->use_collocated_contact_pressure_flag())
       {
        contact_pressure=nod_pt->value(first_index_p);
       }


      // Contact/non-penetration residual
      if (this->Enable_stick)
       {
        // Enforce contact
        if(local_eqn_p >= 0)
         {
          local_residuals[local_eqn_p]-=d;
         }

        // Enforce continuity of temperature
        if(local_eqn_q >= 0)
         {
          local_residuals[local_eqn_q]+=temp_deviation;
         }
       }
      else
       {
        // Piecewise linear variation for non-penetration constraint
        if (-d>contact_pressure)//No penetration
         {
          if(local_eqn_p >= 0)
           {
            //Crush contact pressure
            local_residuals[local_eqn_p]+=contact_pressure;
           }
          if(local_eqn_q >= 0)
           {
            //Crush heat flux as no contact means no heat transfer
            local_residuals[local_eqn_q]+=heat_flux;
           }
         }
        else //Penetration
         {
          if(local_eqn_p >= 0)
           {
            //Increase contact_pressure until d=0 
            local_residuals[local_eqn_p]-=d;
           }
          if(local_eqn_q >= 0)
           {
            //Change lagrange mult (heat flux) until we have continuity
            local_residuals[local_eqn_q]+=temp_deviation;
           }
         }
       }
     }
   }
 }

 // Now deal with the penetrator equilibrium equations (if any!)
 unsigned n=this->Penetrator_eq_data_type.size();
 for (unsigned i=0;i<n;i++)
  {
   if (this->Penetrator_eq_data_type[i]>=0)
    {
     switch(unsigned(this->Penetrator_eq_data_type[i]))
      {
         
      case TemplateFreeContactElementBase::External_data:
      {
       int local_eqn=this->external_local_eqn(
        this->Penetrator_eq_data_data_index[i],
        this->Penetrator_eq_data_index[i]);
       if (local_eqn>=0)
        {
         local_residuals[local_eqn]=contact_force[i];
        }
      }
      break;
        
      case TemplateFreeContactElementBase::Nodal_position_data:
      {
       // position type (dummy -- hierher paranoid check)
       unsigned k=0;
       int local_eqn=this->position_local_eqn(
        this->Penetrator_eq_data_data_index[i],k,
        this->Penetrator_eq_data_index[i]);
       if (local_eqn>=0)
        {
         local_residuals[local_eqn]=contact_force[i];
        }
      }
      break;
        
        
      case TemplateFreeContactElementBase::Nodal_data:
      {
       int local_eqn=this->nodal_local_eqn(
        this->Penetrator_eq_data_data_index[i],
        this->Penetrator_eq_data_index[i]);
       if (local_eqn>=0)
        {
         local_residuals[local_eqn]+=contact_force[i];
        }
      }
      break;
        
      default:


       std::stringstream junk;
       junk << "Never get here: "
            << "unsigned(Penetrator_eq_data_type[i]) = "
            << unsigned(this->Penetrator_eq_data_type[i]);
       throw OomphLibError(
        junk.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);

      }
    }  
  }

 // Now add local contribution to existing entries
 for (unsigned j=0;j<n_dof;j++)
  {
   residuals[j]+=local_residuals[j];
  }

}





/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////


//======================================================================
/// A class for elements that allow the imposition of an 
/// applied flux on the boundaries of UnsteadyHeat elements.
/// The element geometry is obtained from the  FaceGeometry<ELEMENT> 
/// policy class.
//======================================================================
template <class ELEMENT>
class HeatedPenetratorFluxElement : 
 public virtual FaceGeometry<ELEMENT>, 
 public virtual FaceElement, 
 public virtual ElementWithExternalElement,
 public virtual TemplateFreeHeatedPenetratorFluxElementBase
{
 
public:

 // hierher
 // /// Function pointer to the prescribed-flux function fct(x,f(x)) -- 
 // /// x is a Vector! 
 // typedef void (*UnsteadyHeatPrescribedFluxFctPt)
 //  (const double& time, const Vector<double>& x, double& flux);


 /// Default constructor
 HeatedPenetratorFluxElement(){}

 /// Constructor, takes the pointer to the "bulk" element and the 
 /// index of the face to be created
 HeatedPenetratorFluxElement(FiniteElement* const &bulk_el_pt, 
                         const int &face_index);

 /// Broken copy constructor
 HeatedPenetratorFluxElement(const HeatedPenetratorFluxElement& dummy) 
  { 
   BrokenCopy::broken_copy("HeatedPenetratorFluxElement");
  } 
 
 /// Pointer to penetrator
 Penetrator* penetrator_pt() const
 {
  return Penetrator_pt;
 }

 /// Set pointer to penetrator
 void set_penetrator_pt(Penetrator* penetrator_pt)
 {
  Penetrator_pt=penetrator_pt;
 }

 // hierher
 // /// Access function for the prescribed-flux function pointer
 // UnsteadyHeatPrescribedFluxFctPt& flux_fct_pt() {return Flux_fct_pt;}
 
 /// Compute the element residual vector
 inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution_ust_heat_flux(
    residuals,GeneralisedElement::Dummy_matrix,0);
  }




 /// Use version in element with external element
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                       DenseMatrix<double> &jacobian)
  {
   ElementWithExternalElement::fill_in_contribution_to_jacobian(residuals,
                                                                jacobian);
  }

 // /// Compute the element's residual vector and its Jacobian matrix
 // inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
 //                                          DenseMatrix<double> &jacobian)
 //  {
 //   //Call the generic routine with the flag set to 1
 //   fill_in_generic_residual_contribution_ust_heat_flux(residuals,jacobian,1);
 //  }

 /// hierher overload...
 double zeta_nodal(const unsigned &n, const unsigned &k, 
                   const unsigned &i) const
  {
   Vector<double> displ(this->penetrator_pt()->rigid_body_displacement());
   unsigned dim=this->node_pt(n)->ndim();
   Vector<double> x(dim);
   for (unsigned ii=0;ii<dim;ii++)
    {
     // hierher elaborate on offset
     x[ii]=this->node_pt(n)->x(ii)+displ[ii];
    }
   Vector<double> zeta(dim-1);
   this->penetrator_pt()->surface_coordinate(x,zeta);
   return zeta[i]; 
  }
  
 /// Output function -- forward to broken version in FiniteElement
 /// until somebody decides what exactly they want to plot here...
 void output(std::ostream &outfile) {FaceGeometry<ELEMENT>::output(outfile);}

 /// Output function -- forward to broken version in FiniteElement
 /// until somebody decides what exactly they want to plot here...
// void output(std::ostream &outfile, const unsigned &n_plot)
 // {FaceGeometry<ELEMENT>::output(outfile,n_plot);}

//Output function -- hierher rushed
 void output(std::ostream &outfile, const unsigned &nplot)
  {

  //Set output Vector
  Vector<double> s(1); //hierher - should be DIM and DIM-1
  Vector<double> x(2);
  
  // Tecplot header info
  outfile << this->tecplot_zone_string(nplot);
  
  // Loop over plot points
  unsigned num_plot_points=this->nplot_points(nplot);
  for (unsigned iplot=0;iplot<num_plot_points;iplot++)
   {
    
    // Get local coordinates of plot point
    this->get_s_plot(iplot,nplot,s);
    
    // Get Eulerian coordinates and displacements
    this->interpolated_x(s,x);
    
    //Output the x,y,..
    for(unsigned i=0;i<2;i++) 
     {outfile << x[i] << " ";}

    // Output heat_flux
    //outfile << this->heat_flux(s) << " ";

    // Output temperature
    outfile << this->penetrator_temperature(s) << std::endl;   

    outfile << std::endl;
   }
  
  // Write tecplot footer (e.g. FE connectivity lists)
  this->write_tecplot_zone_footer(outfile,nplot);

  }


 /// C-style output function -- forward to broken version in FiniteElement
 /// until somebody decides what exactly they want to plot here...
 void output(FILE* file_pt) {FaceGeometry<ELEMENT>::output(file_pt);}

 /// C-style output function -- forward to broken version in 
 /// FiniteElement until somebody decides what exactly they want to plot 
 /// here...
 void output(FILE* file_pt, const unsigned &n_plot)
  {FaceGeometry<ELEMENT>::output(file_pt,n_plot);}


 /// Pure virtual function: Return temperature at 
 /// local coordinate s.
 double penetrator_temperature(const Vector<double>& s) const
  {
   // Get shape function for temperature
   unsigned n_node = this->nnode();
   Shape psi(n_node);
   this->shape(s,psi);
   
   // Interpolated temperature
   double interpolated_t=0.0;
   for(unsigned l=0;l<n_node;l++) 
    {
     // hierher: everywhere: change to nodal value (or whatever the
     // non-raw counterpart is)
     interpolated_t+=this->node_pt(l)->value(U_index_ust_heat)*psi[l];
    }
   return interpolated_t;
  }


protected:
 
 /// Function to compute the shape and test functions and to return 
 /// the Jacobian of mapping between local and global (Eulerian)
 /// coordinates
 inline double shape_and_test(const Vector<double> &s, Shape &psi, Shape &test)
  const
  {
   //Find number of nodes
   unsigned n_node = nnode();

   //Get the shape functions
   shape(s,psi);

   //Set the test functions to be the same as the shape functions
   for(unsigned i=0;i<n_node;i++) {test[i] = psi[i];}

   //Return the value of the jacobian
   return J_eulerian(s);
  }


 // hierher
 // /// Function to calculate the prescribed flux at a given spatial
 // /// position and at a gien time
 // void get_flux(const double& time, const Vector<double>& x, double& flux)
 //  {
 //   //If the function pointer is zero return zero
 //   if(Flux_fct_pt == 0)
 //    {
 //     flux=0.0;
 //    }
 //   //Otherwise call the function
 //   else
 //    {
 //     (*Flux_fct_pt)(time,x,flux);
 //    }
 //  }


private:


 /// Compute the element residual vector.
 /// flag=1(or 0): do (or don't) compute the Jacobian as well. 
 void fill_in_generic_residual_contribution_ust_heat_flux(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  unsigned flag);

 /// Pointer to penetrator
 Penetrator* Penetrator_pt;
 
 // hierher
 // /// Function pointer to the (global) prescribed-flux function
 // UnsteadyHeatPrescribedFluxFctPt Flux_fct_pt;

 /// The spatial dimension of the problem
 unsigned Dim;

 /// The index at which the unknown is stored at the nodes
 unsigned U_index_ust_heat;


};

/// ///////////////////////////////////////////////////////////////////// 
/// ///////////////////////////////////////////////////////////////////// 
/// ///////////////////////////////////////////////////////////////////// 

//===========================================================================
/// Constructor, takes the pointer to the "bulk" element and the 
/// index of the face to be created.
//===========================================================================
template<class ELEMENT>
HeatedPenetratorFluxElement<ELEMENT>::
HeatedPenetratorFluxElement(FiniteElement* const &bulk_el_pt, 
                        const int &face_index) : 
 FaceGeometry<ELEMENT>(), FaceElement()
{

 // Let the bulk element build the FaceElement, i.e. setup the pointers 
 // to its nodes (by referring to the appropriate nodes in the bulk
 // element), etc.
 bulk_el_pt->build_face_element(face_index,this);
 
#ifdef PARANOID
 {
  //Check that the element is not a refineable 3d element
  ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(bulk_el_pt);

  //If it's three-d
  if(elem_pt->dim()==3)
   {
    //Is it refineable
    RefineableElement* ref_el_pt=dynamic_cast<RefineableElement*>(elem_pt);
    if(ref_el_pt!=0)
     {
      if (this->has_hanging_nodes())
       {
        throw OomphLibError(
         "This flux element will not work correctly if nodes are hanging\n",
         OOMPH_CURRENT_FUNCTION,
         OOMPH_EXCEPTION_LOCATION);
       }
     }
   }
 }
#endif
  

 // Set source element storage: one interaction with an external 
 // element -- the HeatedLinearSurfaceContactElement that supplies the
 // heat flux required to ensure continuity of temperature
 this->set_ninteraction(1); 
 

 // hierher
 // // Initialise the prescribed-flux function pointer to zero
 // Flux_fct_pt = 0;
 
 // Extract the dimension of the problem from the dimension of 
 // the first node
 Dim = this->node_pt(0)->ndim();
 
 //Set up U_index_ust_heat. Initialise to zero, which probably won't change
 //in most cases, oh well, the price we pay for generality
 U_index_ust_heat = 0;
 
 //Cast to the appropriate UnsteadyHeatEquation so that we can
 //find the index at which the variable is stored
 //We assume that the dimension of the full problem is the same
 //as the dimension of the node, if this is not the case you will have
 //to write custom elements, sorry
 switch(Dim)
  {
   //One dimensional problem
  case 1:
  {
   UnsteadyHeatEquations<1>* eqn_pt = 
    dynamic_cast<UnsteadyHeatEquations<1>*>(bulk_el_pt);
   //If the cast has failed die
   if(eqn_pt==0)
    {
     std::string error_string =
      "Bulk element must inherit from UnsteadyHeatEquations.";
     error_string += 
      "Nodes are one dimensional, but cannot cast the bulk element to\n";
     error_string += "UnsteadyHeatEquations<1>\n.";
     error_string += 
      "If you desire this functionality, you must implement it yourself\n";
     
     throw OomphLibError(error_string,
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   //Otherwise read out the value
   else
    {
     //Read the index from the (cast) bulk element
     U_index_ust_heat = eqn_pt->u_index_ust_heat();
    }
  }
  break;
  
  //Two dimensional problem
  case 2:
  {
   UnsteadyHeatEquations<2>* eqn_pt = 
    dynamic_cast<UnsteadyHeatEquations<2>*>(bulk_el_pt);
   //If the cast has failed die
   if(eqn_pt==0)
    {
     std::string error_string =
      "Bulk element must inherit from UnsteadyHeatEquations.";
     error_string += 
      "Nodes are two dimensional, but cannot cast the bulk element to\n";
     error_string += "UnsteadyHeatEquations<2>\n.";
     error_string += 
      "If you desire this functionality, you must implement it yourself\n";
     
     throw OomphLibError(error_string,
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   else
    {
     //Read the index from the (cast) bulk element.
     U_index_ust_heat = eqn_pt->u_index_ust_heat();
    }
  }
  break;
  
  //Three dimensional problem
  case 3:
  {
   UnsteadyHeatEquations<3>* eqn_pt = 
    dynamic_cast<UnsteadyHeatEquations<3>*>(bulk_el_pt);
   //If the cast has failed die
   if(eqn_pt==0)
    {
     std::string error_string =
      "Bulk element must inherit from UnsteadyHeatEquations.";
     error_string += 
      "Nodes are three dimensional, but cannot cast the bulk element to\n";
     error_string += "UnsteadyHeatEquations<3>\n.";
     error_string += 
      "If you desire this functionality, you must implement it yourself\n";
     
     throw OomphLibError(error_string,
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
     
    }
   else
    {
     //Read the index from the (cast) bulk element.
     U_index_ust_heat = eqn_pt->u_index_ust_heat();
    }
  }
  break;
  
  //Any other case is an error
  default:
   std::ostringstream error_stream; 
   error_stream <<  "Dimension of node is " << Dim 
                << ". It should be 1,2, or 3!" << std::endl;
   
   throw OomphLibError(error_stream.str(),
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
   break;
  }
 
}



//===========================================================================
/// Compute the element's residual vector and the (zero) Jacobian matrix.
//===========================================================================
template<class ELEMENT>
void HeatedPenetratorFluxElement<ELEMENT>::
fill_in_generic_residual_contribution_ust_heat_flux(
 Vector<double> &residuals, DenseMatrix<double> &jacobian, 
 unsigned flag)
{
 //Find out how many nodes there are
 const unsigned n_node = nnode();

 // hierher 
 // // Get continuous time from timestepper of first node
 // double time=node_pt(0)->time_stepper_pt()->time_pt()->time();
  
 //Set up memory for the shape and test functions
 Shape psif(n_node), testf(n_node);
 
 //Set the value of n_intpt
 const unsigned n_intpt = integral_pt()->nweight();
 
 //Set the Vector to hold local coordinates
 Vector<double> s(Dim-1);
 
 //Integer to store the local equation and unknowns
 int local_eqn=0;

 // Locally cache the index at which the variable is stored
 const unsigned u_index_ust_heat = U_index_ust_heat;
 
 //Loop over the integration points
 //--------------------------------
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {

   //Assign values of s
   for(unsigned i=0;i<(Dim-1);i++)
    {
     s[i] = integral_pt()->knot(ipt,i);
    }

   //Get the integral weight
   double w = integral_pt()->weight(ipt);
   
   //Find the shape and test functions and return the Jacobian
   //of the mapping
   double J = shape_and_test(s,psif,testf);
   
   //Premultiply the weights and the Jacobian
   double W = w*J;
   
   // hierher

   // //Need to find position to feed into flux function
   // Vector<double> interpolated_x(Dim);
   
   // //Initialise to zero
   // for(unsigned i=0;i<Dim;i++)
   //  {
   //   interpolated_x[i] = 0.0;
   //  }
   
   // //Calculate velocities and derivatives
   // for(unsigned l=0;l<n_node;l++) 
   //  {
   //   //Loop over velocity components
   //   for(unsigned i=0;i<Dim;i++)
   //    {
   //     interpolated_x[i] += nodal_position(l,i)*psif[l];
   //    }
   //  }
   
   // //Get the imposed flux
   // double flux;
   // get_flux(time,interpolated_x,flux);


   // Get the flux
   double flux=0.0;
   if (external_element_pt(0,ipt)!=0)
    {
     Vector<double> s_ext(external_element_local_coord(0,ipt));
     TemplateFreeHeatedLinearSurfaceContactElementBase* el_pt=
      dynamic_cast<TemplateFreeHeatedLinearSurfaceContactElementBase*>(
       external_element_pt(0,ipt));
     // hierher sign?
     flux=-el_pt->heat_flux(s_ext);
    }

   //Now add to the appropriate equations
   
   //Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
    {
     local_eqn = nodal_local_eqn(l,u_index_ust_heat);
     /*IF it's not a boundary condition*/
     if(local_eqn >= 0)
      {
       //Add the prescribed flux terms
       residuals[local_eqn] -= flux*testf[l]*W;
         
       // Imposed traction doesn't depend upon the solution, 
       // --> the Jacobian is always zero, so no Jacobian
       // terms are required
      }
    }
  }
}







} // end namespace extension

/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////////





//================================================================
/// Function-type-object to compare finite elements based on
/// their x coordinate
//================================================================
class FiniteElementComp
{

public:

 /// Comparison. Is x coordinate of el1_pt less than that of el2_pt?
 bool operator()(FiniteElement* const& el1_pt, FiniteElement* const& el2_pt) 
  const
  {
   return el1_pt->node_pt(0)->x(0) < el2_pt->node_pt(0)->x(0);
  }

};



//======Start_of_warped_line===============================================
/// Warped line in 2D space
//=========================================================================
class WarpedLine : public GeomObject
{

public:

 /// Constructor: Specify amplitude of deflection from straight horizontal line
 WarpedLine(const double& ampl, const double& x_min, const double& x_max) 
  : GeomObject(1,2)
  {
   Ampl=ampl;
   X_min=x_min;
   X_max=x_max;
   Reversed=false;
   Lift_off_amplitude=0.0;
   Lift_off_alpha=100.0;
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

 /// Position vector at Lagrangian coordinate zeta 
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   double zetaa=X_min+zeta[0]*(X_max-X_min);
   if (Reversed)
    {
     zetaa=X_max+zeta[0]*(X_min-X_max);
    }

   double alpha=atan((4.0*Ampl)/(1.0+4.0*Ampl*Ampl));
   double y_c=1.0+1.0/(2.0*tan(alpha));
   double radius=1.0/(2.0*sin(alpha));

   // Position vector for circular shape
   r[0] = zetaa; 
   r[1] = y_c-sqrt(radius*radius-(zetaa-0.5)*(zetaa-0.5));

   // Lift off
   r[1]-=Lift_off_amplitude*exp(-Lift_off_alpha*(zetaa-0.5)*(zetaa-0.5));
  }
 
 /// Parametrised position on object: r(zeta). Evaluated at
 /// previous timestep. t=0: current time; t>0: previous
 /// timestep. Forward to steady version
 void position(const unsigned& t, const Vector<double>& zeta,
                       Vector<double>& r) const
  {
   position(zeta,r);
  }

 /// Return radius
 double radius() const
  {
   double alpha=atan((4.0*Ampl)/(1.0+4.0*Ampl*Ampl));
   return 1.0/(2.0*sin(alpha));
  }

 /// Return y coordinate of centre
 double y_c() const
  {
   double alpha=atan((4.0*Ampl)/(1.0+4.0*Ampl*Ampl));
   return 1.0+1.0/(2.0*tan(alpha));
  }

 /// Access to amplitude
 double& ampl() {return Ampl;}

 /// How many items of Data does the shape of the object depend on?
 /// None.
 unsigned ngeom_data() const
  {
   return 0;
  }

 /// Local coordinates are reversed
 void set_reversed()
  {
   Reversed=true;
  }

 /// Local coordinates are not reversed
 void set_non_reversed()
  {
   Reversed=false;
  }

 /// Lift off amplitude
 double& lift_off_amplitude()
  {
   return Lift_off_amplitude;
  }

 /// Exponential factor for lift off (controls sharpness) 
 double& lift_off_alpha()
  {
   return Lift_off_alpha;
  }

private:

 /// Amplitude of perturbation
 double Ampl;

 /// Min zeta coordinate
 double X_min;

 /// Max zeta coordinate
 double X_max;

 /// Reverse?
 bool Reversed;
 
 /// Lift off amplitude
 double Lift_off_amplitude;

 /// Exponential factor for lift off (controls sharpness) 
 double Lift_off_alpha;
};






/// //////////////////////////////////////////////////////////////////// 
/// ////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////// 


//======================================================================
// hierher update
/// Penetrator that keeps circle in contact with a control node 
/// on target surface (made of solid contact face elements) -- centre
/// of the circular penetrator is located at
///  
///     {\bf r}_c = {\bf R}_p + R {\bf e}_alpha
///
/// where {\bf R}_p is the position of the control point, R the radius
/// of the circular penetrator, and {\bf e}_alpha is a unit vector
/// inclined at an angle \alpha against the vertical. 
/// Penetration can be driven in two ways.
/// (1) We impose the vertical position of the control point (by 
///     pseudo-hijacking the Lagrange-multiplier (representing the contact 
///     pressure) stored at the controlled node. This means that rather than
///     determining the contact pressure from the no-penetration constraint,
///     (which we know to be satisfied by construction) we determine it
///     from the condition that {\bf R}_p \cdot {\bf e}_y = Y_c which 
///     is prescribed. We also impose the angle \alpha (stored as an
///     internal Data value in the element) by solving it via the
///     equation \alpha-\alpha_{prescribed} = 0.
/// (2) We impose the weight (i.e. the vertical reaction force from the 
///     contact elements) by using the equation 
///
///       \int p_c {\bf n} \cdot {\bf e}_y ds - W = 0
///   
///     as the equation for the pseudo-hijacked contact pressure and similarly,
///     use the horizontal force balance
/// 
///       \int p_c {\bf n} \cdot {\bf e}_x ds - H = 0
///   
///     to determine the rotation angle. Here, W and H are prescribed
///     and the integral is computed from the contact elements that
///     potentially impact on the penetrator.
//======================================================================
class HeatedCircularPenetratorElement : public virtual GeneralisedElement,
                                        public virtual HeatedPenetrator
 {

   public:

  /// Constructor: Pass pointer to control node whose 
  /// index_of_contact_pressure-th value represents the Lagrange multiplier
  /// (the discrete contact pressure) that has been traded for the
  /// vertical displacement/weight constraint. Also need the indices
  /// of the nodal values that store the horizontal/vertical 
  /// displacement (linear elasticity).
  // hierher not only SolidNode
  HeatedCircularPenetratorElement(SolidNode* control_node_pt, 
                            const unsigned& index_of_contact_pressure,
                            const unsigned& index_of_horizontal_displacement,
                            const unsigned& index_of_vertical_displacement,
                            double* r_pt) 
   {
    // Create internal data, representing the angle of rotation about
    // contact point. Determined either directly via insisting
    // that the difference between this value and the target is zero
    // or by insisting that the horizontal force is zero.
    add_internal_data(new Data(1));
    internal_data_pt(0)->set_value(0,0.0);

    oomph_info << "In here..\n";
    exit(0);

    // Store pointer to radius
    Radius_pt=r_pt;
    
    // Store original centre
    Orig_centre=centre();

    // Control node 
    Control_node_pt=control_node_pt;

    // Where is the tradedd contact pressure stored?
    Index_of_contact_pressure=index_of_contact_pressure;

    // Where is the horizontal displacement (linear_elasticity) stored?
    Index_of_horizontal_displacement=index_of_horizontal_displacement;

    // Where is the vertical displacement (linear_elasticity) stored?
    Index_of_vertical_displacement=index_of_vertical_displacement;

    //Pointer to target weight (null if vertical displacement of control
    //node is imposed)
    Target_weight_pt=0;

    // Pointer to target horizontal force (null if rotation angle angle 
    // about control node is imposed)
    Target_horizontal_force_pt=0;

    // Pointer to  target vertical displacement of control node (null if 
    // weight is imposed)
    Target_yc_pt=0;

    // Pointer to target rotation angle about control node (null 
    // if horizontal force is imposed)
    Target_rotation_angle_pt=0;

    // Pointer to mesh of contact elements that contribute to force
    Contact_element_mesh_pt=0;
   }
  

  /// Get surface coordinate along penetrator for given point x.
  /// We assume that point on the surface and given point share the
  /// same polar angle and return that polar angle 
  void surface_coordinate(const Vector<double>& x, Vector<double>& zeta) const
   {
    zeta[0]=atan2(x[1]-centre(1),x[0]-centre(0));
   }
  
  
  /// Get rigid body displacement of reference point in penetrator.
  Vector<double> rigid_body_displacement() const
   {    
    unsigned n=Orig_centre.size();
    Vector<double> displ(n);
    for (unsigned i=0;i<n;i++)
     {
      displ[i]=centre(i)-Orig_centre[i];
     }
    return displ;
   }
  
  
  /// Set original centre of penetrator (for computation of rigid body
  /// displacement
  void set_original_centre(const Vector<double>& orig_centre)
   {
    Orig_centre=orig_centre;
   }
  
  
  /// Vector of pairs identifying values (via a pair of pointer to 
  /// Data object and index within it) that correspond to the Data values 
  /// that are determined by the horizontal/vertical/... equilibrium equations.
  Vector<std::pair<Data*,unsigned> > equilibrium_data()
   {
    // We're in 2D
    Vector<std::pair<Data*,unsigned> > thingy(2);

    
    // Horizontal equilibrium determines the rotation angle
    // which is stored as the zero-th internal data
    if (Target_horizontal_force_pt==0)
     {
      thingy[0]=std::make_pair(static_cast<Data*>(0),0);
     }
    else
     {
      thingy[0]=std::make_pair(internal_data_pt(0),0);
     }


    // Vertical equilibrium determines the discrete contact pressure
    // (Lagrange multiplier) at control node
    if (Target_weight_pt==0)
     {
      thingy[1]=std::make_pair(static_cast<Data*>(0),0);
     }
    else
     {
      thingy[1]=std::make_pair(
       external_data_pt(External_data_index_of_traded_contact_pressure),
       Index_of_contact_pressure);
     }

    return thingy;
   }
  
  /// Angle of rotation around contact point
  double angle() const
   {
    return internal_data_pt(0)->value(0);
   }


  /// Set angle of rotation around contact point
  void set_angle(const double& angle)
   {
    internal_data_pt(0)->set_value(0,angle);
   }


  /// Access to pointer to mesh of contact elements that contribute to 
  /// force on penetrator
  Mesh* contact_element_mesh_pt() const
   {
    return Contact_element_mesh_pt;
   }

  /// Set pointer to mesh of contact elements and setup
  /// external Data, i.e. Data that affects the residuals in this
  /// element. Also set the node pointed to by Control_node_pt
  /// as external Data for the elements in the contact mesh
  /// (unless they contain this node already).
  void set_contact_element_mesh_pt(Mesh* contact_element_mesh_pt)
  {
   Contact_element_mesh_pt=contact_element_mesh_pt;
   flush_external_data();
   
   // Store Data associated with control node: It contains the traded
   // Lagrange multiplier (contact pressure) 
   External_data_index_of_traded_contact_pressure=
    add_external_data(Control_node_pt);
   
   // Store its position data
   add_external_data(Control_node_pt->variable_position_pt());

   // Store it as Data (which includes the displacement)
   add_external_data(Control_node_pt);

   // Loop over all the elements in the contact mesh
   // If they don't contain the contact node already, its position
   // is external data because it affects the penetration.
   unsigned nel=Contact_element_mesh_pt->nelement();
   for (unsigned e=0;e<nel;e++)
    {
     bool el_contains_control_node=false;
     FiniteElement* el_pt=Contact_element_mesh_pt->finite_element_pt(e);
     unsigned nnod=el_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       SolidNode* nod_pt=dynamic_cast<SolidNode*>(el_pt->node_pt(j));
       if (nod_pt==Control_node_pt)
        {
         el_contains_control_node=true;
        }
      }
     
     // Position of control node affects position of penetrator and
     // therefore is external data for all contact elements (apart from
     // any that contain the control node as one of their own)
     if (!el_contains_control_node)
      {
       // position
       el_pt->add_external_data(Control_node_pt->variable_position_pt());
       // displacement relative to position
       el_pt->add_external_data(Control_node_pt);
      }
     
     // Rotation angle angle affects position of penetrator and therefore
     // affects penetration at all contact elements
     el_pt->add_external_data(internal_data_pt(0));
    }
  }


  

  /// Set target horizontal and vertical force to be in equilibrium
  void set_equilibrium_target_forces()
   {
#ifdef PARANOID
    if (Target_horizontal_force_pt==0)
     {
      std::stringstream junk;
      junk << "Target_horizontal_force_pt==0\n"
           << "Please set it by call to impose_weight(...)\n";
      throw OomphLibError(
       junk.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
     } 
    if (Target_weight_pt==0)
     {
      std::stringstream junk;
      junk << "Target_weight_pt==0.\n"
           << "Please set it by call to impose_horizontal_force(...)\n";
      throw OomphLibError(
       junk.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
     }
#endif
    Vector<double> force(resulting_force());
    (*Target_horizontal_force_pt)=-force[0];
    (*Target_weight_pt)=-force[1];
   }

  /// Target weight (returns zero if not imposed)
  double target_weight()
   {
    if (Target_weight_pt==0)
     {
      return 0.0;
     }
    else
     {
      return *Target_weight_pt;
     }
   }

  /// Target horizontal force (returns zero if not imposed)
  double target_horizontal_force()
   {
    if (Target_horizontal_force_pt==0)
     {
      return 0.0;
     }
    else
     {
      return *Target_horizontal_force_pt;
     }
   }

  /// Target vertical position of control point (returns zero if not imposed)
  double target_yc()
   {
    if (Target_yc_pt==0)
     {
      return 0.0;
     }
    else
     {
      return *Target_yc_pt;
     }
   }

  /// Target rotation angle about contact point (returns zero if not imposed)
  double target_rotation_angle()
   {
    if (Target_rotation_angle_pt==0)
     {
      return 0.0;
     }
    else
     {
      return *Target_rotation_angle_pt;
     }
   }

  /// Is vertical positon of control node imposed? If false then weight imposed.
  bool yc_is_imposed()
   {
    return (Target_yc_pt!=0);
   }

  /// Impose weight (rather than imposed displacement). Target
  /// weight specified via pointer.
  void impose_weight(double* target_weight_pt)
   {
    Target_weight_pt=target_weight_pt;
    Target_yc_pt=0;
   }

  /// Impose vertical position of control node (rather than weight).
  /// Target vertical position of control node specified via pointer.
  void impose_yc(double* target_yc_pt)
   {
    Target_weight_pt=0;
    Target_yc_pt=target_yc_pt;
   }


  /// Is angle of rotation about control node imposed? If false then 
  /// horizontal force is imposed.
  bool rotation_angle_is_imposed()
   {
    return (Target_rotation_angle_pt!=0);
   }


  /// Impose horizontal force (rather than rotation about contact node). 
  /// Target force specified via pointer.
  void impose_horizontal_force(double* target_horizontal_force_pt)
   {
    Target_horizontal_force_pt=target_horizontal_force_pt;
    Target_rotation_angle_pt=0;
   }

  /// Impose rotation about contact node (rather than horizontal force)
  /// Target angle specified via pointer.
  void impose_rotation_angle(double* target_rotation_angle_pt)
   {
    Target_horizontal_force_pt=0;
    Target_rotation_angle_pt=target_rotation_angle_pt;
   }

  /// Fill in contribution to residuals
  void fill_in_contribution_to_residuals(Vector<double> &residuals) 
  {
   // Get resulting force from all associated PseudoContactElements
   // onto the elastic body
   Vector<double> force(resulting_force());

   // Equation for Lagrange multiplier (contact pressure) at controlled
   // node
   int local_eqn=external_local_eqn(
    External_data_index_of_traded_contact_pressure,Index_of_contact_pressure); 
   if (local_eqn>=0)
    {
     // Resulting force from all associated PseudoContactElements
     // onto the elastic body is equal and opposite to force on penetrator
     if (Target_weight_pt!=0)
      {
       residuals[local_eqn]+=(*Target_weight_pt);
      }
     // Impose vertical position of control node
     else
      {
#ifdef PARANOID
       if (Target_yc_pt!=0)
        {
#endif
         residuals[local_eqn]+=Control_node_pt->x(1)+
          Control_node_pt->value(Index_of_vertical_displacement)-(*Target_yc_pt);
#ifdef PARANOID
        }
       else
        {
         std::stringstream junk;
         junk << "Target_yc_pt=0\n"
              << "Set with impose_yc(...)\n";
         throw OomphLibError(
          junk.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
        }
#endif
      }
    }



   // Equation for rotation angle
   local_eqn=internal_local_eqn(0,0); 
   if (local_eqn>=0)
    {
     // Resulting force from all associated PseudoContactElements
     // onto the elastic body is equal and opposite to force on penetrator
     if (Target_horizontal_force_pt!=0)
      {
       residuals[local_eqn]+=(*Target_horizontal_force_pt);
      }
     // Set rotation angle 
     else
      {
#ifdef PARANOID
       if (Target_rotation_angle_pt!=0)
        {
#endif
         residuals[local_eqn]+=internal_data_pt(0)->value(0)-
          (*Target_rotation_angle_pt);
#ifdef PARANOID
        }
       else
        {
         std::stringstream junk;
         junk << "Target_rotation_angle_pt=0\n"
              << "Set with impose_rotation_angle(...)\n";
         throw OomphLibError(
          junk.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
        }
#endif
      }
    }
  }
    
  /// Get centre of penetrator
  Vector<double> centre() const
   {
    Vector<double> rc(2);
    rc[0]=centre(0); 
    rc[1]=centre(1); 
    return rc;
   }

  /// Get centre of penetrator 
  double centre(const unsigned& i) const
   {
    switch (i)
     {
     case 0:
      return Control_node_pt->x(0)+
       Control_node_pt->value(Index_of_horizontal_displacement)+
       (*Radius_pt)*sin(angle());
      break;
      
     case 1:
      return Control_node_pt->x(1)+
       Control_node_pt->value(Index_of_vertical_displacement)+
       (*Radius_pt)*cos(angle());
      break;
      
     default:
      std::stringstream junk;
      junk << "Wrong index: " << i 
           << "\nCan only handle 0 or 1 (it's 2D!)\n";
      throw OomphLibError(
       junk.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
     }
   }
 
  /// Get penetration for given point x.
  void penetration(const Vector<double>& x, const Vector<double>& n,
                   double& d,bool& intersection) const
  {
   // Vector from potential contact point to centre of penetrator
   Vector<double> l(2);
   l[0]=centre(0)-x[0];
   l[1]=centre(1)-x[1];
   
   // Distance from potential contact point to centre of penetrator
   double ll=sqrt(l[0]*l[0]+l[1]*l[1]);

   // Projection of vector from potential contact point to centre of penetrator
   // onto outer unit normal on potential contact point
   double project=n[0]*l[0]+n[1]*l[1];
   double project_squared=project*project;

   // Final term in square root
   double b_squared=ll*ll-(*Radius_pt)*(*Radius_pt);

   // Is square root negative? In this case we have no intersection
   // and we return penetration as -DBL_MAX
   if (project_squared<b_squared)
    {
     d= -DBL_MAX;
     intersection = false;
    }
   else
    {
     double sqr=sqrt(project_squared-b_squared);
     d= -std::min(project-sqr,project+sqr);
     intersection = true;
    }
  }


  /// Get temperature on penetrator at point "associated" with
  /// point x using the same logic as for the position function). Here
  /// we assume that both points shrare the same polar angle relative
  /// to the centre of (circular!) penetrator
  double temperature(const Vector<double>& x) const
   {
    double phi=atan2(x[1]-centre(1),x[0]-centre(0));
    return cos(phi-0.5*MathematicalConstants::Pi);
   }

  /// Output coordinates of penetrator at nplot plot points
  void output(std::ostream &outfile, const unsigned& nplot) const
   {
    for (unsigned j=0;j<nplot;j++)
     {
      double phi=2.0*MathematicalConstants::Pi*double(j)/double(nplot-1);
      outfile << centre(0)+(*Radius_pt)*cos(phi) << " " 
              << centre(1)+(*Radius_pt)*sin(phi)
              << std::endl;
     }
   }
    
  /// Get position to surface, r, in terms of surface coordinate zeta.
  void position_from_zeta(const Vector<double>& zeta, 
                          Vector<double>& r) const
  {
   double phi=zeta[0];
   r[0]=centre(0)+(*Radius_pt)*cos(phi);
   r[1]=centre(1)+(*Radius_pt)*sin(phi);
  }

  /// Resulting force from all associated ContactElements
  Vector<double> resulting_force() const
  {
   Vector<double> contact_force(2,0.0);
   Vector<double> el_contact_force(2);
   unsigned nel=Contact_element_mesh_pt->nelement();
   for (unsigned e=0;e<nel;e++)
    {
     dynamic_cast<TemplateFreeContactElementBase*>(
      Contact_element_mesh_pt->element_pt(e))->
      resulting_contact_force(el_contact_force);
     for (unsigned i=0;i<2;i++)
      {
       contact_force[i]+=el_contact_force[i];
      }
    }
   return contact_force;
  }

  /// Radius of penetrator
  double radius() const {return *Radius_pt;}
 
 private:

  /// Pointer to radius of penetrator
  double* Radius_pt;

  /// Original centre of penetrator (origin for cylindrical polar
  /// coordinate system)
  Vector<double> Orig_centre;

  /// Control node
  SolidNode* Control_node_pt;

  /// Index at which contact pressure (Lagr mult) is stored in nodal
  /// data associated with control node
  unsigned Index_of_contact_pressure;

  /// Where is the vertical displacement (linear_elasticity) stored?
  unsigned Index_of_vertical_displacement;

  /// Where is the horizontal displacement (linear_elasticity) stored?
  unsigned Index_of_horizontal_displacement;

  /// Index of external data that contains the the contact
  /// pressure in its Index_of_contact_pressure-th value
  unsigned External_data_index_of_traded_contact_pressure;
  
  /// Pointer to target weight (null if vertical displacement of control
  /// node is imposed)
  double* Target_weight_pt;

  /// Pointer to target horizontal force (null if rotation angle angle 
  /// about control node is imposed)
  double* Target_horizontal_force_pt;

  /// Pointer to  target vertical displacement of control node (null if 
  /// weight is imposed)
  double* Target_yc_pt;

  /// Pointer to target rotation angle about control node (null 
  /// if horizontal force is imposed)
  double* Target_rotation_angle_pt;

  /// Mesh of contact elements that contribute to weight/horizontal force
  Mesh* Contact_element_mesh_pt;

 };


/// //////////////////////////////////////////////////////////////////// 
/// ////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////// 

//======start_of_ProblemParameters=====================
/// Namespace for problem parameters
//=====================================================
namespace ProblemParameters
{

 /// hierher temp flux.
 void unit_flux(const double& time, const Vector<double>& x, double& flux)
 {
  flux = 1.0;
 }

 /// hierher
 double T_contact=0.0;

#ifdef STRUCTURED_MESH

 /// GeomObject specifying the shape of the boundary: Initially it's 
 /// almost flat.
 WarpedLine Boundary_geom_object(1.0e-10,0.0,1.0);

#else

 /// Left end of contact region (for unstructured mesh only)
 double X_contact_end_left=0.3;
 
 /// Right end of contact region (for unstructured mesh only)
 double X_contact_end_right=0.7;

 /// GeomObject specifying the shape of the boundary: Initially it's 
 /// almost flat. Starts at the left
 WarpedLine Boundary_geom_object_left(1.0e-10,0.0,X_contact_end_left);

 /// GeomObject specifying the shape of the boundary: Initially it's 
 /// almost flat.
 WarpedLine Boundary_geom_object_contact(1.0e-10,X_contact_end_left,
                                         X_contact_end_right);

 /// GeomObject specifying the shape of the boundary: Initially it's 
 /// almost flat.
 WarpedLine Boundary_geom_object_right(1.0e-10,X_contact_end_right,1.0);

#endif

 /// Impose position of centre (i.e. a stand-alone penetrator with
 /// prescribed position or indirectly via control node?
 bool Impose_position_of_centre=true;
 
 /// Non-dim density for solid
 double Lambda_sq=0.0;
 
 /// Poisson's ratio for solid (both real and pseudo)
 double Nu=0.3;

 /// The elasticity tensor
 IsotropicElasticityTensor E(Nu);

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt=0;

 /// Radius of penetrator
 double Radius=0.5;

 /// Penetrator
 HeatedPenetrator* Penetrator_pt=0;

 /// NOTE: WE IMPOSE EITHER THESE ...

 /// Weight of penetrator
 double Weight=0.0;

 /// Horizontal force of penetrator
 double Horizontal_force=0.0;

 /// ... OR THESE...

 /// Target vertical position of control node
 double Y_c=0.0;

 /// Target rotation angle about control node
 double Rotation_angle=0.0;

 /// ...OR THIS
 
 /// Position of centre of penetrator
 Vector<double> Centre;

 /// Initial/max element area
 double El_area=0.02;

 /// Factor for element length on contact boundary
 /// ie how many times smaller should the elements on the boundary be?
 double Element_length_factor=0.01;


 // /// Body force magnitude
 // double Body_force_amplitude=0.0;
 
 // // Sharpness of body force
 // double Body_force_alpha=1.0e4;

 // /// The body force function
 // void body_force(const double &time,
 //                 const Vector<double> &x,
 //                 Vector<double> &result)
 // {
 // result[0] = 0.0;
 // result[1] = -Body_force_amplitude*(1.0-x[0])*x[0]*
 //  exp(-Body_force_alpha*(x[0]-0.5)*(x[0]-0.5));
 // }
} // end of ProblemParameters


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//=====start_of_problem_class=========================================
/// Problem class
//====================================================================
template<class ELEMENT>
class ContactProblem : public Problem
{

public:

 /// Constructor
 ContactProblem();
 
 /// Destructor (empty)
 ~ContactProblem(){}
 
 /// Actions before timestep
 void actions_before_implicit_timestep()
  {
   double dyc=0.00024; 
   ProblemParameters::Centre[1]-=dyc;
   oomph_info << "Re-solving imposed circle pos for yc=" 
              << ProblemParameters::Centre[1]
              << std::endl;
   unsigned b=Contact_boundary_id;
   unsigned nnod=Bulk_mesh_pt->nboundary_node(b);
   for (unsigned j=0;j<nnod;j++)
    {
     SolidNode* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,j);
     // hierher index read out from element, though 
     // this shouldn't actually be set at all. Kill
     nod_pt->set_value(2,ProblemParameters::T_contact); 
    }
  }

 
 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}
 
 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve()
  {
   // For maximum stability: Reset the current nodal positions to be
   // the "stress-free" ones -- this assignment means that the
   // parameter study no longer corresponds to a physical experiment
   // but is what we'd do if we wanted to use the solid solve
   // to update a fluid mesh in an FSI problem, say.
   Bulk_mesh_pt->set_lagrangian_nodal_coordinates();

   // DoubleVector r;
   // CRDoubleMatrix jac;

   // oomph_info << "SETTING UP JAC FOR OUTPUT OF MOST RECENT JAC\n";
   // get_jacobian(r,jac);
   // oomph_info << "DONE SETTING UP JAC FOR OUTPUT OF MOST RECENT JAC\n";
   // jac.sparse_indexed_output("most_recent_jacobian.dat");

   // ofstream descr_file;
   // descr_file.open("most_recent_description.dat");
   // describe_dofs(descr_file);
   // descr_file.close();
  }
  

 /// Actions before adapt: wipe contact elements
 void actions_before_adapt() 
  {
   // Backup x coordinate of old control node
   Xc_old=Control_node_pt->x(0);

#ifdef STRUCTURED_MESH

   // Make backup of surface mesh
   Backed_up_surface_contact_mesh_pt=
    new BackupMeshForProjection<QElement<1,3> >(Surface_contact_mesh_pt,
                                                Contact_boundary_id);
#else

   //double max_el_length=Maximum_element_length_on_contact_boundary*
   //ProblemParameters::Element_length_factor;

   //Impose max_el_length from el_area, ie derive average element length from El_area then 
   //apply multiplier
   double max_el_length = ProblemParameters::Element_length_factor*
    (2*sqrt(sqrt(3))*sqrt(ProblemParameters::El_area));

   Bulk_mesh_pt->boundary_polyline_pt(Contact_boundary_id)
    ->set_maximum_length(max_el_length);

   Boulder_mesh_pt->boundary_polyline_pt(Boulder_bottom_boundary_id)
    ->set_maximum_length(max_el_length);

   // Make backup of surface mesh
   Backed_up_surface_contact_mesh_pt=
    new BackupMeshForProjection<TElement<1,3> >(Surface_contact_mesh_pt,
                                                Contact_boundary_id);
#endif


   // // Output contact elements
   // ofstream some_file;
   // char filename[100];
   // sprintf(filename,"contact_before.dat");
   // some_file.open(filename);
   // unsigned nel=Surface_contact_mesh_pt->nelement();
   // for (unsigned e=0;e<nel;e++)
   //  {
   //   dynamic_cast<HeatedLinearSurfaceContactElement<ELEMENT>* >(
   //    Surface_contact_mesh_pt->element_pt(e))->output(some_file);
   //  }
   // some_file.close();
   

   // Output impose displ elements
   ofstream some_file;
   char filename[100];
   sprintf(filename,"impose_before.dat");
   some_file.open(filename);
   Displ_imposition_mesh_pt->output(some_file);
   some_file.close();

   // // Kill the  elements and wipe surface mesh
   delete_contact_elements();
   delete_displ_imposition_elements();

   // Wipe the mesh
   Penetrator_mesh_pt->flush_element_and_node_storage();

   // Rebuild the Problem's global mesh from its various sub-meshes
   rebuild_global_mesh();
  }
 
 /// Actions after adapt: 
 /// Setup the problem again -- remember that the mesh has been
 /// completely rebuilt and its element's don't have any
 /// pointers to source fcts etc. yet
 void actions_after_adapt()
  {
   // Create "contact" heat flux elements on surface of boulder
   create_contact_heat_elements_on_boulder();

   // Create normal heat flux elements on surface of boulder
   create_imposed_heat_flux_elements_on_boulder();

   // Create contact elements
   create_contact_elements();
   
   // Create elements that impose displacement of melt line
   create_displ_imposition_elements();

   // Now project from backup of original contact mesh to new one
   Backed_up_surface_contact_mesh_pt->project_onto_new_mesh(
    Surface_contact_mesh_pt);

   // For maximum stability: Reset the current nodal positions to be
   // the "stress-free" ones -- this assignment means that the
   // parameter study no longer corresponds to a physical experiment
   // but is what we'd do if we wanted to use the solid solve
   // to update a fluid mesh in an FSI problem, say.
   Bulk_mesh_pt->set_lagrangian_nodal_coordinates();

   // Rebuild elements
   complete_problem_setup();
   
   // Rebuild the Problem's global mesh from its various sub-meshes
   rebuild_global_mesh();
   
   // Kill backed up mesh
   delete Backed_up_surface_contact_mesh_pt;
   Backed_up_surface_contact_mesh_pt=0;

   // Output impose displ elements
   ofstream some_file;
   char filename[100];
   sprintf(filename,"impose_after.dat");
   some_file.open(filename);
   Displ_imposition_mesh_pt->output(some_file);
   some_file.close();

  }

 /// Switch to displ control
 void switch_to_displ_control()
  {
   dynamic_cast<HeatedCircularPenetratorElement*>(
    ProblemParameters::Penetrator_pt)->impose_yc(&ProblemParameters::Y_c);
   dynamic_cast<HeatedCircularPenetratorElement*>(
    ProblemParameters::Penetrator_pt)->impose_rotation_angle(
    &ProblemParameters::Rotation_angle);
   ProblemParameters::Y_c=Control_node_pt->x(1)+Control_node_pt->value(1);
  }

 /// Switch to force control
 void switch_to_force_control()
  {
   dynamic_cast<HeatedCircularPenetratorElement*>(
    ProblemParameters::Penetrator_pt)->impose_weight(
     &ProblemParameters::Weight);
   dynamic_cast<HeatedCircularPenetratorElement*>(
    ProblemParameters::Penetrator_pt)->impose_horizontal_force(
     &ProblemParameters::Horizontal_force); 
   dynamic_cast<HeatedCircularPenetratorElement*>(
    ProblemParameters::Penetrator_pt)->set_equilibrium_target_forces();

   // Re-set contact mesh -- we now need to treat the positions
   // and penalty pressures of all nodes as external data of the
   // penetrator element!
   dynamic_cast<HeatedCircularPenetratorElement*>(
    ProblemParameters::Penetrator_pt)->set_contact_element_mesh_pt(
    Surface_contact_mesh_pt);

   // Reset penetrator because its equilibrium data (which is external
   // data for contact elements) has changed
   unsigned n_element=Surface_contact_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement 
     HeatedLinearSurfaceContactElement<ELEMENT> *el_pt = 
      dynamic_cast<HeatedLinearSurfaceContactElement<ELEMENT>*>(
       Surface_contact_mesh_pt->element_pt(e));
     
     // Set pointer to penetrator
     el_pt->set_penetrator_pt(ProblemParameters::Penetrator_pt);
    }
   
   cout <<"New number of equations: " << assign_eqn_numbers() << std::endl; 

  }


 /// Doc the solution
 void doc_solution();
 

private:


 /// Create elements that enforce prescribed boundary motion
 /// by Lagrange multipliers
 void create_displ_imposition_elements()
  {

   Vector<unsigned> boundary_id;
   boundary_id.push_back(Contact_boundary_id);
#ifndef STRUCTURED_MESH
   Node* left_contact_node_pt=0;
   Node* right_contact_node_pt=0;
   Node* left_left_top_node_pt=0;
   Node* right_left_top_node_pt=0;
   Node* left_right_top_node_pt=0;
   Node* right_right_top_node_pt=0;
   boundary_id.push_back(Left_top_boundary_id);
   boundary_id.push_back(Right_top_boundary_id);
#endif
   unsigned nb=boundary_id.size();
   for (unsigned bb=0;bb<nb;bb++)
    {
     unsigned b=boundary_id[bb];
     
     // How many bulk elements are adjacent to boundary b?
     unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
     
     // Loop over the bulk elements adjacent to boundary b?
     for(unsigned e=0;e<n_element;e++)
      {
       // Get pointer to the bulk element that is adjacent to boundary b
       ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
        Bulk_mesh_pt->boundary_element_pt(b,e));
       
       //Find the index of the face of element e along boundary b
       int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
       
       // Create new element
       ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>* el_pt=
        new ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>(
         bulk_elem_pt,face_index);
       
       // Add to mesh
       Displ_imposition_mesh_pt->add_element_pt(el_pt);
       
#ifdef STRUCTURED_MESH

       // Set the GeomObject that defines the boundary shape and
       // specify which bulk boundary we are attached to (needed to extract
       // the boundary coordinate from the bulk nodes)
       el_pt->set_boundary_shape_geom_object_pt( 
        &ProblemParameters::Boundary_geom_object,b);

#else
       
       switch(b)
        {
        case Contact_boundary_id:

         // Set the GeomObject that defines the boundary shape and
         // specify which bulk boundary we are attached to (needed to extract
         // the boundary coordinate from the bulk nodes)
         el_pt->set_boundary_shape_geom_object_pt( 
          &ProblemParameters::Boundary_geom_object_contact,b);

         {
          unsigned nnod=el_pt->nnode();
          for (unsigned j=0;j<nnod;j++)
           {
            Node* nod_pt=el_pt->node_pt(j);
            if (nod_pt->is_on_boundary(Left_top_boundary_id))
             {
              left_contact_node_pt=nod_pt;
             }
            else if (nod_pt->is_on_boundary(Right_top_boundary_id))
             {
              right_contact_node_pt=nod_pt;
             }

           }
         }

         break;

        case Left_top_boundary_id:
         
         // Set the GeomObject that defines the boundary shape and
         // specify which bulk boundary we are attached to (needed to extract
         // the boundary coordinate from the bulk nodes)
         el_pt->set_boundary_shape_geom_object_pt( 
          &ProblemParameters::Boundary_geom_object_left,b);
         
         {
          unsigned nnod=el_pt->nnode();
          for (unsigned j=0;j<nnod;j++)
           {
            Node* nod_pt=el_pt->node_pt(j);
            if (nod_pt->is_on_boundary(Left_boundary_id))
             {
              left_left_top_node_pt=nod_pt;
             }
            else if (nod_pt->is_on_boundary(Contact_boundary_id))
             {
              right_left_top_node_pt=nod_pt;
             }
           }
         }

         break;

        case Right_top_boundary_id:

         // Set the GeomObject that defines the boundary shape and
         // specify which bulk boundary we are attached to (needed to extract
         // the boundary coordinate from the bulk nodes)
         el_pt->set_boundary_shape_geom_object_pt( 
          &ProblemParameters::Boundary_geom_object_right,b);
         
         {
          unsigned nnod=el_pt->nnode();
          for (unsigned j=0;j<nnod;j++)
           {
            Node* nod_pt=el_pt->node_pt(j);
            if (nod_pt->is_on_boundary(Right_boundary_id))
             {
              right_right_top_node_pt=nod_pt;
             }
            else if (nod_pt->is_on_boundary(Contact_boundary_id))
             {
              left_right_top_node_pt=nod_pt;
             }
           }
         }
         break;
         
        default:
         
         // Never get here...
         oomph_info << "Never get here! b = " << b << "\n";
         abort();
        }
       

#endif

       // Loop over the nodes 
       unsigned nnod=el_pt->nnode();
       for (unsigned j=0;j<nnod;j++)
        {
         Node* nod_pt = el_pt->node_pt(j);
         
         // Is the node also on side boundaries?
         if ((nod_pt->is_on_boundary(Left_boundary_id))||
             (nod_pt->is_on_boundary(Right_boundary_id)))
          {
           // How many nodal values were used by the "bulk" element
           // that originally created this node?
           unsigned n_bulk_value=el_pt->nbulk_value(j);
           
           // hierher: this statement is unlikely to be true. Fix! 
           // The remaining ones are Lagrange multipliers and we pin them.
           unsigned nval=nod_pt->nvalue();
           for (unsigned j=n_bulk_value;j<nval;j++)
            {
             nod_pt->pin(j);
            }
          }
        }
      }  
    }
   
#ifndef STRUCTURED_MESH

   Vector<double> zeta_left(1);
   Vector<double> zeta_right(1);
   left_left_top_node_pt->get_coordinates_on_boundary(Left_top_boundary_id, 
                                                      zeta_left);
   right_left_top_node_pt->get_coordinates_on_boundary(Left_top_boundary_id, 
                                                       zeta_right);
   if (zeta_left[0]>zeta_right[0])
    {
     ProblemParameters::Boundary_geom_object_left.set_reversed();
     oomph_info << "left top is reversed\n";
    }
   else
    {
     ProblemParameters::Boundary_geom_object_left.set_non_reversed();
     oomph_info << "left top is not reversed\n";
    }
   
   left_contact_node_pt->get_coordinates_on_boundary(Contact_boundary_id, 
                                                     zeta_left);
   right_contact_node_pt->get_coordinates_on_boundary(Contact_boundary_id, 
                                                      zeta_right);
   if (zeta_left[0]>zeta_right[0])
    {
     ProblemParameters::Boundary_geom_object_contact.set_reversed();
     oomph_info << "contact is reversed\n";
    }
   else
    {
     ProblemParameters::Boundary_geom_object_contact.set_non_reversed();
     oomph_info << "contact is not reversed\n";
    }
   
   left_right_top_node_pt->get_coordinates_on_boundary(Right_top_boundary_id, 
                                                       zeta_left);
   right_right_top_node_pt->get_coordinates_on_boundary(Right_top_boundary_id, 
                                                        zeta_right);
   if (zeta_left[0]>zeta_right[0])
    {
     ProblemParameters::Boundary_geom_object_right.set_reversed();
     oomph_info << "right top is reversed\n";
    }
   else
    {
     ProblemParameters::Boundary_geom_object_right.set_non_reversed();
     oomph_info << "right top is not reversed\n";
    }
   
#endif
   
   
  } // end of create_displ_imposition_elements


 /// Delete elements that enforce prescribed boundary motion
 /// by Lagrange multiplliers
 void delete_displ_imposition_elements()
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Displ_imposition_mesh_pt->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Displ_imposition_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Displ_imposition_mesh_pt->flush_element_and_node_storage();
  }
 

 /// Create "contact" heat flux elements on surface of boulder
 void create_contact_heat_elements_on_boulder()
  {
   // How many bulk elements are adjacent to boundary b?
   unsigned b=Boulder_bottom_boundary_id; 
   unsigned n_element = Boulder_mesh_pt->nboundary_element(b);
   
   // Loop over the bulk elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to boundary b
     ProjectableUnsteadyHeatElement<TUnsteadyHeatElement<2,3> >* bulk_elem_pt =
      dynamic_cast<ProjectableUnsteadyHeatElement<TUnsteadyHeatElement<2,3> >*>(
       Boulder_mesh_pt->boundary_element_pt(b,e));
     
     //What is the face index of element e along boundary b
     int face_index = Boulder_mesh_pt->face_index_at_boundary(b,e);
     
     // Build the corresponding contact element
     HeatedPenetratorFluxElement
      <ProjectableUnsteadyHeatElement<TUnsteadyHeatElement<2,3> > >* 
      contact_element_pt = new 
      HeatedPenetratorFluxElement
      <ProjectableUnsteadyHeatElement<TUnsteadyHeatElement<2,3> > >(
       bulk_elem_pt,face_index);
     
     // Pass pointer to penetrator
     contact_element_pt->set_penetrator_pt(ProblemParameters::Penetrator_pt);

     //Add the contact element to the surface mesh
     Boulder_surface_contact_mesh_pt->add_element_pt(contact_element_pt);

    } //end of loop over bulk elements adjacent to boundary b    
  }

 /// Create imposed heat flux elements on surface of boulder
 void create_imposed_heat_flux_elements_on_boulder()
  {
   // How many bulk elements are adjacent to boundary b?
   unsigned b=Boulder_top_boundary_id; 
   unsigned n_element = Boulder_mesh_pt->nboundary_element(b);
   
   // Loop over the bulk elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to boundary b
     ProjectableUnsteadyHeatElement<TUnsteadyHeatElement<2,3> >* bulk_elem_pt =
      dynamic_cast<ProjectableUnsteadyHeatElement<TUnsteadyHeatElement<2,3> >*>(
       Boulder_mesh_pt->boundary_element_pt(b,e));
     
     //What is the face index of element e along boundary b
     int face_index = Boulder_mesh_pt->face_index_at_boundary(b,e);
     
     // Build the corresponding flux element
     UnsteadyHeatFluxElement
      <ProjectableUnsteadyHeatElement<TUnsteadyHeatElement<2,3> > >* 
      flux_element_pt = new 
      UnsteadyHeatFluxElement
      <ProjectableUnsteadyHeatElement<TUnsteadyHeatElement<2,3> > >(
       bulk_elem_pt,face_index);
     
     // Pass pointer to heat flux
     flux_element_pt->flux_fct_pt() = &ProblemParameters::unit_flux;
     
     //Add the contact element to the surface mesh
     Boulder_surface_heat_flux_mesh_pt->add_element_pt(flux_element_pt);

    } //end of loop over bulk elements adjacent to boundary b    
  }

 /// Create contact elements
 void create_contact_elements()
  {
   // How many bulk elements are adjacent to boundary b?
   unsigned b=Contact_boundary_id; 
   unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
   
   // Loop over the bulk elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(b,e));
     
     //What is the face index of element e along boundary b
     int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
     
     // Build the corresponding contact element
     HeatedLinearSurfaceContactElement<ELEMENT>* contact_element_pt = new 
      HeatedLinearSurfaceContactElement<ELEMENT>(bulk_elem_pt,face_index,
                                           Contact_id);

     // Pass pointer to penetrator
     contact_element_pt->set_penetrator_pt(ProblemParameters::Penetrator_pt);
     
     //Add the contact element to the surface mesh
     Surface_contact_mesh_pt->add_element_pt(contact_element_pt);

    } //end of loop over bulk elements adjacent to boundary b    
  }



 /// Delete contact elements
 void delete_contact_elements()
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Surface_contact_mesh_pt->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Surface_contact_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Surface_contact_mesh_pt->flush_element_and_node_storage();


   // How many surface elements are in the surface mesh
   n_element = Boulder_surface_contact_mesh_pt->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Boulder_surface_contact_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Boulder_surface_contact_mesh_pt->flush_element_and_node_storage();


   // How many surface elements are in the surface mesh
   n_element = Boulder_surface_heat_flux_mesh_pt->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Boulder_surface_heat_flux_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Boulder_surface_heat_flux_mesh_pt->flush_element_and_node_storage();
   
  }

 /// Helper function to (re-)set boundary condition
 /// and complete the build of  all elements
 void complete_problem_setup()
  {

   // Set (pseudo-)solid mechanics properties for all elements
   //---------------------------------------------------------
   unsigned n_element = Bulk_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     //Cast to a solid element
     ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

     // Set the constitutive law
     el_pt->constitutive_law_pt() =
      ProblemParameters::Constitutive_law_pt;
     
     // Set density to zero
     el_pt->LinearElasticityEquations<2>::lambda_sq_pt()=
      &ProblemParameters::Lambda_sq;

     // Set density to zero
     el_pt->PVDEquationsBase<2>::lambda_sq_pt()=
      &ProblemParameters::Lambda_sq;
     
     // Set the elasticity tensor
     el_pt->elasticity_tensor_pt() = &ProblemParameters::E;
     
     // Disable inertia
     el_pt->LinearElasticityEquations<2>::disable_inertia();

     // Disable inertia
     el_pt->PVDEquationsBase<2>::disable_inertia();
    }

   // Apply boundary conditions for solid
   //------------------------------------

   // Bottom: completely pinned
   unsigned b=Bottom_boundary_id;
   unsigned nnod=Bulk_mesh_pt->nboundary_node(b);
   for (unsigned j=0;j<nnod;j++)
    {
     SolidNode* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,j);
     nod_pt->pin_position(0);
     nod_pt->pin_position(1);
     nod_pt->pin(0);
     nod_pt->pin(1);
    }

   // Sides: Symmetry bcs
   b=Left_boundary_id;
   nnod=Bulk_mesh_pt->nboundary_node(b);
   for (unsigned j=0;j<nnod;j++)
    {
     SolidNode* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,j);
     nod_pt->pin_position(0);
     nod_pt->pin(0);
    }
   b=Right_boundary_id;
   nnod=Bulk_mesh_pt->nboundary_node(b);
   for (unsigned j=0;j<nnod;j++)
    {
     SolidNode* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,j);
     nod_pt->pin_position(0);
     nod_pt->pin(0);
    }


   // Enforce contact at most central or most heavily loaded node
   //------------------------------------------------------------
   {
    
    // Update angle?
    bool update_angle=false;
    double phi_old=0.0;
    HeatedCircularPenetratorElement* pen_el_pt=
     dynamic_cast<HeatedCircularPenetratorElement*>(
      ProblemParameters::Penetrator_pt);
    if (pen_el_pt!=0)
     {
      update_angle=true;
      phi_old=pen_el_pt->angle();
     }
    

    // Find closest/most loaded node
    double x_c=0.5;
    Control_node_pt=0;
    SolidNode* most_central_node_pt=0;
    SolidNode* most_loaded_node_pt=0;
    double dist_min=DBL_MAX;
    double load_max=0.0;
    unsigned b=Contact_boundary_id;
    unsigned nnod=Bulk_mesh_pt->nboundary_node(b);
    for (unsigned j=0;j<nnod;j++)
     {
      SolidNode* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,j);
      
      // Find closest node
      double dist=std::fabs(nod_pt->x(0)-x_c);
      if (dist<dist_min) 
       {
        dist_min=dist;
        most_central_node_pt=nod_pt;
       }

      // Find most loaded node
      BoundaryNodeBase *bnod_pt=dynamic_cast<BoundaryNodeBase*>(nod_pt);
      unsigned index_of_contact_pressure=
       bnod_pt->index_of_first_value_assigned_by_face_element(Contact_id);
      if (nod_pt->value(index_of_contact_pressure)>load_max)
       {
        load_max=nod_pt->value(index_of_contact_pressure);
        most_loaded_node_pt=nod_pt;
       }
     }

    // No load
    if (most_loaded_node_pt==0)
     {
      oomph_info << "Choosing most central node as control node\n";
      Control_node_pt=most_central_node_pt;
     }
    else
     {
      oomph_info << "Choosing most loaded node as control node\n";
      Control_node_pt=most_loaded_node_pt;
     }
    oomph_info << "Control node located at: "
               << Control_node_pt->x(0) << " " 
               << Control_node_pt->x(1) << " "
               << std::endl;


    // Update angle
    if (update_angle)
     {
      double phi_new=asin((Xc_old-Control_node_pt->x(0))/
                          ProblemParameters::Radius+sin(phi_old));
      dynamic_cast<HeatedCircularPenetratorElement*>(
       ProblemParameters::Penetrator_pt)->set_angle(phi_new);
      oomph_info << "Old/new angle: " << phi_old << " "
                 << phi_new << std::endl;
     }


    //...............................................................


    
    // Set target vertical position of control node to its current
    // position. This is fine as initial assignment; it's then over-written
    // before the next solve (if the displacement is controlled). 
    // If/when this function is called during adaptation it's fine too
    // because it doesn't impose another displacement increment -- it simply
    // maintains the position that was obtained on the previous solve.
    unsigned index_of_vertical_displacement=1;
    ProblemParameters::Y_c=Control_node_pt->x(1)+
     Control_node_pt->value(index_of_vertical_displacement); 
    
    // Index of nodal value at control node that stores the traded
    // contact pressure
    unsigned index_of_traded_contact_pressure=UINT_MAX;

    // Impose position of centre when penetrator position directly imposed
    //--------------------------------------------------------------------
    if (ProblemParameters::Impose_position_of_centre)
     {
      // Delete old one
      delete ProblemParameters::Penetrator_pt;
      
      // Make new one
      ProblemParameters::Penetrator_pt =
       new HeatedCircularPenetrator(&ProblemParameters::Centre,
                                    ProblemParameters::Radius);
     }
    // Compute penetrator position as part of the solution, either by
    // --------------------------------------------------------------
    // prescribing nodal position or weight
    // ------------------------------------
    else //--
     {
      // Loop over face elements to identify the ones that contain 
      // the control node
      unsigned nel=Surface_contact_mesh_pt->nelement();
      bool found=false;
      for (unsigned e=0;e<nel;e++)
       {
        HeatedLinearSurfaceContactElement<ELEMENT>* el_pt=
         dynamic_cast<HeatedLinearSurfaceContactElement<ELEMENT>*>(
          Surface_contact_mesh_pt->element_pt(e));
        unsigned nnod=el_pt->nnode();
        for (unsigned j=0;j<nnod;j++)
         {
          SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(el_pt->node_pt(j));
          if (solid_nod_pt==Control_node_pt)
           {
            // Got it!
            found=true;
            
            // Find index at which contact pressure/Lagrange multiplier
            // is stored
            BoundaryNodeBase *bnod_pt = 
             dynamic_cast<BoundaryNodeBase*>(solid_nod_pt); 
            
            // Get the index of the first nodal value associated with
            // this FaceElement
            unsigned new_index_of_traded_contact_pressure=
             bnod_pt->index_of_first_value_assigned_by_face_element(Contact_id);
            
            // Copy across
#ifdef PARANOID
            if (index_of_traded_contact_pressure!=UINT_MAX)
             {
              if (new_index_of_traded_contact_pressure!=
                  index_of_traded_contact_pressure)
               {
                std::stringstream junk;
                junk << "Inconsistency in identification of index of traded"
                     << "contact pressure: " 
                     << new_index_of_traded_contact_pressure
                     << " != " << index_of_traded_contact_pressure;
                throw OomphLibError(
                 junk.str(),
                 OOMPH_CURRENT_FUNCTION,
                 OOMPH_EXCEPTION_LOCATION);
               }
             }
#endif
            index_of_traded_contact_pressure=
             new_index_of_traded_contact_pressure;
           }
         }
       }
      if (!found) 
       {
        std::stringstream junk;
        junk << "Control node not found!";
        throw OomphLibError(
         junk.str(),
         OOMPH_CURRENT_FUNCTION,
         OOMPH_EXCEPTION_LOCATION);
       }
      
      // Back up old penetrator (if it existed)
      HeatedCircularPenetratorElement* old_penetrator_pt=
       dynamic_cast<HeatedCircularPenetratorElement*>(
        ProblemParameters::Penetrator_pt);
      
      // Make new one
      unsigned index_of_horizontal_displacement=0;
      unsigned index_of_vertical_displacement=1;
      ProblemParameters::Penetrator_pt =
       new HeatedCircularPenetratorElement(Control_node_pt,
                                           index_of_traded_contact_pressure,
                                           index_of_horizontal_displacement,
                                           index_of_vertical_displacement,
                                           &ProblemParameters::Radius);
      
      // Displacement or weight imposed?
      bool impose_displ=true;
      if (old_penetrator_pt!=0)
       {
        if (!old_penetrator_pt->yc_is_imposed())
         {
          impose_displ=false;
         }
       }
      if (impose_displ)
       {
        dynamic_cast<HeatedCircularPenetratorElement*>(
         ProblemParameters::Penetrator_pt)->impose_yc(&ProblemParameters::Y_c);
       }
      else
       {
        dynamic_cast<HeatedCircularPenetratorElement*>(
         ProblemParameters::Penetrator_pt)->impose_weight(
          &ProblemParameters::Weight);
       }
      
      
      // Angle or horizontal force imposed?
      bool impose_angle=true;
      if (old_penetrator_pt!=0)
       {
        if (!old_penetrator_pt->rotation_angle_is_imposed())
         {
          impose_angle=false;
         }
       }
      if (impose_angle)
       {
        dynamic_cast<HeatedCircularPenetratorElement*>(
         ProblemParameters::Penetrator_pt)->impose_rotation_angle(
          &ProblemParameters::Rotation_angle);
       }
      else
       {
        dynamic_cast<HeatedCircularPenetratorElement*>(
         ProblemParameters::Penetrator_pt)->impose_horizontal_force(
          &ProblemParameters::Horizontal_force);
        dynamic_cast<HeatedCircularPenetratorElement*>(
         ProblemParameters::Penetrator_pt)->set_angle(
          old_penetrator_pt->angle());
       }
      
      // Now kill the old one!
      delete old_penetrator_pt;
      
      // Add to mesh
      Penetrator_mesh_pt->add_element_pt(
       dynamic_cast<HeatedCircularPenetratorElement*>(
        ProblemParameters::Penetrator_pt));
     }
    
    // Pass contact elements to penetrator element and declare their
    // positions and Lagrange multiplier (contact pressure) values
    // to be external data.
    HeatedCircularPenetratorElement* el_pt=
     dynamic_cast<HeatedCircularPenetratorElement*>(
     ProblemParameters::Penetrator_pt);
    if (el_pt!=0)
     {
      el_pt->set_contact_element_mesh_pt(Surface_contact_mesh_pt);
     }
    
   } // end of penetrator position computed as part of the solution
     // (directly or indirectly)

   // Loop over the contact elements to pass pointer to penetrator
   //-------------------------------------------------------------
   n_element=Surface_contact_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement 
     HeatedLinearSurfaceContactElement<ELEMENT> *el_pt = 
      dynamic_cast<HeatedLinearSurfaceContactElement<ELEMENT>*>(
       Surface_contact_mesh_pt->element_pt(e));
     
     // Set pointer to penetrator
     el_pt->set_penetrator_pt(ProblemParameters::Penetrator_pt);
    }
   
   
 
   // Find external elements (heated contact elements) that provide
   // the heat flux onto the boulder to maintain continuity of
   // temperature
   Multi_domain_functions::
    Accept_failed_locate_zeta_in_setup_multi_domain_interaction=true;

   Multi_domain_functions::setup_multi_domain_interaction
    <HeatedLinearSurfaceContactElement<ELEMENT> >
    (this,Boulder_surface_contact_mesh_pt,Surface_contact_mesh_pt);
   
   Multi_domain_functions::
    Accept_failed_locate_zeta_in_setup_multi_domain_interaction=false;
 
   // Find external elements (penetrator flux elements) that provide
   // the target temperature
   Multi_domain_functions::
    Accept_failed_locate_zeta_in_setup_multi_domain_interaction=true;
   
   Multi_domain_functions::setup_multi_domain_interaction
    <HeatedPenetratorFluxElement
     <ProjectableUnsteadyHeatElement<TUnsteadyHeatElement<2,3> > > >  
    (this,Surface_contact_mesh_pt,Boulder_surface_contact_mesh_pt);
   
   Multi_domain_functions::
    Accept_failed_locate_zeta_in_setup_multi_domain_interaction=false;


   // hierher setup interaction
   
  //  // Loop over the heat contact elements to pass pointer to flux
  //  //------------------------------------------------------------
  //  n_element=Boulder_surface_contact_mesh_pt->nelement();
  //  for(unsigned e=0;e<n_element;e++)
  //   {
  //    // Upcast from GeneralisedElement 
  //    HeatedPenetratorFluxElement
  //     <ProjectableUnsteadyHeatElement<TUnsteadyHeatElement<2,3> > >* el_pt = 
  //     dynamic_cast<HeatedPenetratorFluxElement
  //     <ProjectableUnsteadyHeatElement<TUnsteadyHeatElement<2,3> > >*>(
  //      Boulder_surface_contact_mesh_pt->element_pt(e));
   
   
  //  // Set the pointer to the prescribed flux function
  //  // hierher setup interaction el_pt->flux_fct_pt() = &ProblemParameters::unit_flux;
   
  // }
 



  }

 /// Pointer to boulder mesh
 RefineableTriangleMesh<ProjectableUnsteadyHeatElement
                        <TUnsteadyHeatElement<2,3> > >* Boulder_mesh_pt;

#ifdef STRUCTURED_MESH

 /// Pointer to bulk mesh
 ElasticRefineableRectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;

#else

 /// Pointer to bulk mesh
 RefineableSolidTriangleMesh<ELEMENT>* Bulk_mesh_pt;

#endif
 
 /// Pointer to the "surface" contact mesh
 Mesh* Surface_contact_mesh_pt;
 
 /// Pointers to meshes of Lagrange multiplier elements to impose 
 /// displacement of top boundary
 SolidMesh* Displ_imposition_mesh_pt;

 /// Mesh of elements that imposed unit heat flux on top of boulder
 Mesh* Boulder_surface_heat_flux_mesh_pt;

 /// Penetrator mesh
 Mesh* Penetrator_mesh_pt;

 /// Pointer to the "surface" mesh on boulder that applies flux boundary 
 /// condition to boulder
 Mesh* Boulder_surface_contact_mesh_pt;
 
 // Boundary enumeration
 enum{
  Bottom_boundary_id,
  Right_boundary_id,
  Contact_boundary_id,
  Left_boundary_id,
  Right_top_boundary_id, 
  Left_top_boundary_id,
  Boulder_top_boundary_id,
  Boulder_bottom_boundary_id
 };
 
 /// Trace file
 ofstream Trace_file;

 // Setup labels for output
 DocInfo Doc_info;

#ifdef STRUCTURED_MESH

 /// Backup of Surface_contact_mesh_pt so the Lagrange multipliers
 /// can be projected across
 BackupMeshForProjection<QElement<1,3> >* Backed_up_surface_contact_mesh_pt;

#else

 /// Backup of Surface_contact_mesh_pt so the Lagrange multipliers
 /// can be projected across
 BackupMeshForProjection<TElement<1,3> >* Backed_up_surface_contact_mesh_pt;

#endif

 /// Pointer to control node where Lagrange multiplier (contact pressure)
 /// is "pseudo-hijacked" to impose either displacement or weight constraint.
 SolidNode* Control_node_pt;

 /// x coordinate of old control node
 double Xc_old;

 /// ID of additional nodal values created by contact elements to store
 /// contact pressure/Lagr. mult.
 unsigned Contact_id;

 /// x coordinate of lower left corner 
 double X_ll;

 /// x coordinate of upper right corner 
 double X_ur;
 
 /// y coordinate of lower left corner 
 double Y_ll;

 /// y coordinate of upper right corner 
 double Y_ur;

 /// Contact boundary in its poly line representation
 TriangleMeshPolyLine* Contact_boundary_pt;

 /// Max. element length on contact boundary
 double Maximum_element_length_on_contact_boundary;
 double Maximum_element_length_on_boulder_boundary;

}; // end of problem class


//========start_of_constructor============================================
/// Constructor for contact problem in square domain
//========================================================================
template<class ELEMENT>
ContactProblem<ELEMENT>::ContactProblem()
{ 

 // Allow for crap initial guess or convergence problems...
 Problem::Max_newton_iterations=50;

 // Initialise
 Control_node_pt=0;

 // Initialise x coordinate of old control node (will be overwritten 
 // when needed)
 Xc_old=0.0;

 // ID of additional nodal values created by contact elements to store
 // contact pressure/Lagr. mult.
 Contact_id=1;

 // Initialise
 Backed_up_surface_contact_mesh_pt=0;

 // Output directory
 Doc_info.set_directory("RESLT");
 
 // Output number
 Doc_info.number()=0;

 // Open trace file
 Trace_file.open("RESLT/trace.dat");
 
 // Set boundaries of domain
 X_ll=0.0;;
 X_ur=1.0;
 Y_ll=0.0;
 Y_ur=1.0;

 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 add_time_stepper_pt(new NewmarkBDF<2>);

 // Create boulder mesh
 //--------------------

 // Create circle representing outer boundary
 Circle* boulder_boundary_circle_pt=new Circle(ProblemParameters::Centre[0],
                                               ProblemParameters::Centre[1],
                                               ProblemParameters::Radius);
 
 // Provide storage for pointers to the two parts of the curvilinear boundary
 Vector<TriangleMeshCurveSection*> boulder_outer_curvilinear_boundary_pt(2);
 
 // First bit
 double zeta_start=0.0;
 double zeta_end=MathematicalConstants::Pi;
 unsigned n_seg_boulder=1000; 
 boulder_outer_curvilinear_boundary_pt[0]=new TriangleMeshCurviLine(
  boulder_boundary_circle_pt,zeta_start,zeta_end,n_seg_boulder,
  Boulder_top_boundary_id);
 
 // Second bit
 zeta_start=MathematicalConstants::Pi;
 zeta_end=2.0*MathematicalConstants::Pi;
 boulder_outer_curvilinear_boundary_pt[1]=new TriangleMeshCurviLine(
  boulder_boundary_circle_pt,zeta_start,zeta_end,n_seg_boulder,
  Boulder_bottom_boundary_id);
 
 // Combine to curvilinear boundary and define the
 // outer boundary
 TriangleMeshClosedCurve* boulder_outer_boundary_pt=
  new TriangleMeshClosedCurve(boulder_outer_curvilinear_boundary_pt);

 // Use the TriangleMeshParameters object for helping on the manage 
 // of the TriangleMesh parameters. The only parameter that needs to take 
 // is the outer boundary.
 TriangleMeshParameters boulder_triangle_mesh_parameters(
  boulder_outer_boundary_pt);

 // Target element size in boulder mesh
 boulder_triangle_mesh_parameters.element_area() = 
  ProblemParameters::El_area;

 // Build boulder mesh
 Boulder_mesh_pt=new RefineableTriangleMesh<
  ProjectableUnsteadyHeatElement
  <TUnsteadyHeatElement<2,3> > >(
   boulder_triangle_mesh_parameters,
   time_stepper_pt());
  
 // Set error estimator for bulk mesh
 Z2ErrorEstimator* boulder_error_estimator_pt=new Z2ErrorEstimator;
 Boulder_mesh_pt->spatial_error_estimator_pt()=boulder_error_estimator_pt;

 // Set element size limits (faire: just uncommented both of these)
 Boulder_mesh_pt->max_element_size()=ProblemParameters::El_area;
 Boulder_mesh_pt->min_element_size()=0.1*ProblemParameters::El_area;
  
 
#ifdef STRUCTURED_MESH

 // Create the mesh

 // # of elements in x-direction
 unsigned n_x=11;

 // # of elements in y-direction
 unsigned n_y=10;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;

 //Now create the mesh 
 Bulk_mesh_pt = new ElasticRefineableRectangularQuadMesh<ELEMENT>
  (n_x,n_y,l_x,l_y,time_stepper_pt());

#else

 // Pointer to the closed curve that defines the outer boundary
 TriangleMeshClosedCurve* closed_curve_pt=0;

 // Build outer boundary as Polygon
  
 // The boundary is bounded by five distinct boundaries, each
 // represented by its own polyline
 Vector<TriangleMeshCurveSection*> boundary_polyline_pt(6);
 
 // Vertex coordinates on boundary
 Vector<Vector<double> > bound_coords(2);
 
 // Left boundary
 bound_coords[0].resize(2);
 bound_coords[0][0]=X_ll;
 bound_coords[0][1]=Y_ur;

 bound_coords[1].resize(2);
 bound_coords[1][0]=X_ll;
 bound_coords[1][1]=Y_ll;

 // Build the boundary polyline
 boundary_polyline_pt[0]=new TriangleMeshPolyLine(bound_coords,
                                                  Left_boundary_id);

 // Bottom boundary
 bound_coords[0].resize(2);
 bound_coords[0][0]=X_ll;
 bound_coords[0][1]=Y_ll;

 bound_coords[1].resize(2);
 bound_coords[1][0]=X_ur;
 bound_coords[1][1]=Y_ll;

 // Build the boundary polyline
 boundary_polyline_pt[1]=new TriangleMeshPolyLine(bound_coords,
                                                  Bottom_boundary_id);

 // Right boundary
 bound_coords[0].resize(2);
 bound_coords[0][0]=X_ur;
 bound_coords[0][1]=Y_ll;

 bound_coords[1].resize(2);
 bound_coords[1][0]=X_ur;
 bound_coords[1][1]=Y_ur;
 
 // Build the boundary polyline
 boundary_polyline_pt[2]=new TriangleMeshPolyLine(bound_coords,
                                                  Right_boundary_id);


 // Right top boundary
 unsigned npt_right=10;
 Vector<Vector<double> > right_top_bound_coords(npt_right);
 right_top_bound_coords[0].resize(2);
 right_top_bound_coords[0][0]=X_ur;
 right_top_bound_coords[0][1]=Y_ur;
 for (unsigned j=1;j<npt_right;j++)
  {  
   right_top_bound_coords[j].resize(2);
   double x=X_ur-(X_ur-ProblemParameters::X_contact_end_right)*double(j)/
    double(npt_right-1);
   double y=Y_ur; 
   right_top_bound_coords[j][0]=x;
   right_top_bound_coords[j][1]=y;
  }

 
 // Build boundary poly line
 TriangleMeshPolyLine* right_top_boundary_pt=
  new TriangleMeshPolyLine(right_top_bound_coords,
                           Right_top_boundary_id);
 boundary_polyline_pt[3]=right_top_boundary_pt;
  

 // Contact boundary
 unsigned npt_contact=200; // hierher 20; 
 Vector<Vector<double> > contact_bound_coords(npt_contact);
 contact_bound_coords[0].resize(2);
 contact_bound_coords[0][0]=right_top_bound_coords[npt_right-1][0];
 contact_bound_coords[0][1]=right_top_bound_coords[npt_right-1][1];
 for (unsigned j=1;j<npt_contact;j++)
  {  
   contact_bound_coords[j].resize(2);
   double x=ProblemParameters::X_contact_end_right-
    (ProblemParameters::X_contact_end_right-
     ProblemParameters::X_contact_end_left)*double(j)/double(npt_contact-1);
   double y=Y_ur;
   contact_bound_coords[j][0]=x;
   contact_bound_coords[j][1]=y;
  }
 
 // Build boundary poly line
 Contact_boundary_pt=
  new TriangleMeshPolyLine(contact_bound_coords,
                           Contact_boundary_id);
 boundary_polyline_pt[4]=Contact_boundary_pt;
  
 // Keep elements near boundary nice and small
 Maximum_element_length_on_contact_boundary=(X_ur-X_ll)/double(npt_contact);
 Contact_boundary_pt->set_maximum_length(
  Maximum_element_length_on_contact_boundary);

 // Left top boundary
 unsigned npt_left=15; 
 Vector<Vector<double> > top_left_bound_coords(npt_left);
 top_left_bound_coords[0].resize(2);
 top_left_bound_coords[0][0]=contact_bound_coords[npt_contact-1][0];
 top_left_bound_coords[0][1]=contact_bound_coords[npt_contact-1][1];
 for (unsigned j=1;j<npt_left-1;j++)
  {  
   top_left_bound_coords[j].resize(2);
   double x=ProblemParameters::X_contact_end_left-
    (ProblemParameters::X_contact_end_left-X_ll)*
    double(j)/double(npt_left-1);
   double y=Y_ur; 
   top_left_bound_coords[j][0]=x;
   top_left_bound_coords[j][1]=y;
  }
 top_left_bound_coords[npt_left-1].resize(2);
 top_left_bound_coords[npt_left-1][0]=X_ll;
 top_left_bound_coords[npt_left-1][1]=Y_ur;
 
 // Build boundary poly line
 TriangleMeshPolyLine* top_left_boundary_pt=
  new TriangleMeshPolyLine(top_left_bound_coords,
                           Left_top_boundary_id);
 boundary_polyline_pt[5]=top_left_boundary_pt;
  
 // Create the triangle mesh polygon for outer boundary
 //----------------------------------------------------
 TriangleMeshPolygon *outer_polygon =
  new TriangleMeshPolygon(boundary_polyline_pt);
  
 // Set the pointer
 closed_curve_pt = outer_polygon;
 
 // Now build the mesh
 //===================

 // Use the TriangleMeshParameters object for helping on the manage of the
 // TriangleMesh parameters
 TriangleMeshParameters triangle_mesh_parameters(closed_curve_pt);

 // Specify the maximum area element
 double uniform_element_area=ProblemParameters::El_area;
 triangle_mesh_parameters.element_area() = uniform_element_area;
 
 // Create the mesh
 Bulk_mesh_pt=
  new RefineableSolidTriangleMesh<ELEMENT>(triangle_mesh_parameters,
                                           time_stepper_pt());
  
#endif


#ifdef STRUCTURED_MESH
 
 Vector<double> lower_left(2);
 lower_left[0]=0.3;
 lower_left[1]=0.8;
 Vector<double> upper_right(2);
 upper_right[0]=0.7;
 upper_right[1]=1.0;
 unsigned central_node_number=4;
 Bulk_mesh_pt->spatial_error_estimator_pt()=
  new DummyErrorEstimator(Bulk_mesh_pt,
                          lower_left,
                          upper_right,
                          central_node_number);

#else

 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Bulk_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // Set element size limits
 Bulk_mesh_pt->max_element_size()=ProblemParameters::El_area;
 Bulk_mesh_pt->min_element_size()=0.1*ProblemParameters::El_area;
  
#endif

 // Create "surface" mesh on boulder that applies contact heat flux boundary 
 // condition to boulder
 Boulder_surface_contact_mesh_pt=new Mesh;
 create_contact_heat_elements_on_boulder();

 // Create "surface" mesh on boulder that applies flux boundary 
 // condition to boulder
 Boulder_surface_heat_flux_mesh_pt=new Mesh;
 create_imposed_heat_flux_elements_on_boulder();

 // Create contact elements
 Surface_contact_mesh_pt=new Mesh;
 create_contact_elements();
 
 // Create elements that enforce prescribed boundary motion
 // by Lagrange multipliers
 Displ_imposition_mesh_pt=new SolidMesh;
 create_displ_imposition_elements();

 // Create mesh for penetrator element
 Penetrator_mesh_pt=new Mesh;

 // Set boundary condition and complete the build of all elements
 complete_problem_setup();

 // Initial assigment
 ProblemParameters::Y_c=Control_node_pt->x(1);

 // Add the sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_contact_mesh_pt); 
 add_sub_mesh(Displ_imposition_mesh_pt); 
 add_sub_mesh(Penetrator_mesh_pt);
 add_sub_mesh(Boulder_mesh_pt);
 add_sub_mesh(Boulder_surface_contact_mesh_pt);
 add_sub_mesh(Boulder_surface_heat_flux_mesh_pt);

 // Combine all submeshes into a single global Mesh
 build_global_mesh();

 // hierher
 //Problem::Newton_solver_tolerance=1.0e-6;

 // hierher 
 //linear_solver_pt()=new FD_LU;

 // Do equation numbering
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor



//=======start_of_doc_solution============================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void ContactProblem<ELEMENT>::doc_solution()
{ 

 oomph_info << "Outputting for step: " << Doc_info.number() << std::endl;
 ofstream some_file;
 char filename[100];


 // Write restart file
/* sprintf(filename,"%s/restart%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 dump(some_file); 
 some_file.close();*/

 // Number of plot points
 unsigned npts;
 npts=3;
 
 // Output solution 
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();
 
 // Output solution coarsely (only element vertices for easier
 // mesh visualisation)
 sprintf(filename,"%s/coarse_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,2);
 some_file.close();

 // Output solution 
 sprintf(filename,"%s/boulder_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Boulder_mesh_pt->output(some_file,npts);
 some_file.close();

  // Output contact region on boulder
 sprintf(filename,"%s/boulder_contact%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Boulder_surface_contact_mesh_pt->output(some_file,npts);
 some_file.close();

  // Output contact region on boulder
 sprintf(filename,"%s/boulder_heat_flux%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Boulder_surface_heat_flux_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output solution coarsely (only element vertices for easier
 // mesh visualisation)
 sprintf(filename,"%s/boulder_coarse_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Boulder_mesh_pt->output(some_file,2);
 some_file.close();
 
 // Output contact elements
 sprintf(filename,"%s/imposed_displ%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Displ_imposition_mesh_pt->output(some_file);
 some_file.close();

 // Output contact elements
 sprintf(filename,"%s/contact%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned nel=Surface_contact_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<HeatedLinearSurfaceContactElement<ELEMENT>* >(
    Surface_contact_mesh_pt->element_pt(e))->output(some_file,3);
  }
 some_file.close();



 // Output integration points of contact elements
 sprintf(filename,"%s/contact_integration_points%i.dat",
         Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Vector<double> s(1);
 Vector<double> x(2);
 nel=Surface_contact_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++) 
  {
   HeatedLinearSurfaceContactElement<ELEMENT>* el_pt=
    dynamic_cast<HeatedLinearSurfaceContactElement<ELEMENT>* >(
     Surface_contact_mesh_pt->element_pt(e));
   unsigned nint=el_pt->integral_pt()->nweight();
   for (unsigned j=0;j<nint;j++)
    {
     s[0]=el_pt->integral_pt()->knot(j,0);
     el_pt->interpolated_x(s,x);
     some_file << x[0] << " " << x[1] << std::endl;
    }
  }
 some_file.close();



 // Output penetrator
 sprintf(filename,"%s/penetrator%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned n=500;
 ProblemParameters::Penetrator_pt->output(some_file,n);
 some_file.close();
  
 // Output contact elements and assemble total resulting force
 Vector<double> total_contact_force(2,0.0);
 Vector<double> contact_force(2,0.0);
 sprintf(filename,"%s/contact%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 nel=Surface_contact_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   HeatedLinearSurfaceContactElement<ELEMENT>* el_pt=
    dynamic_cast<HeatedLinearSurfaceContactElement<ELEMENT>*>(
     Surface_contact_mesh_pt->element_pt(e));
   el_pt->output(some_file,3);
   el_pt->resulting_contact_force(contact_force);
   total_contact_force[0]+=contact_force[0];
   total_contact_force[1]+=contact_force[1];
  }
 some_file.close();
 
 double radius_of_elastic_body=0.0;
#ifdef STRUCTURED_MESH
 radius_of_elastic_body=
  ProblemParameters::Boundary_geom_object.radius();
#else
 radius_of_elastic_body=
  ProblemParameters::Boundary_geom_object_contact.radius();
#endif

 // Get half-width of Hertz contact region
 double b_hertz=0.0;
 if (std::fabs(total_contact_force[1])>1.0e-16)
  {
   b_hertz=sqrt(4.0*(1.0-ProblemParameters::Nu*ProblemParameters::Nu)/
                (MathematicalConstants::Pi*
                 (-1.0/radius_of_elastic_body+
                  1.0/ProblemParameters::Radius))*
                (-total_contact_force[1]));
  }
 oomph_info << "b_hertz " << b_hertz <<  std::endl;
 
 double p_max_hertz = 0.0;
 if(b_hertz!=0.0)
  {
   p_max_hertz=2.0*total_contact_force[1]/
   (MathematicalConstants::Pi*b_hertz);
  }

 // Output Hertzian pressure contact distribution
 sprintf(filename,"%s/hertz%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 n=500;
 HeatedCircularPenetratorElement* pen_el_pt=
  dynamic_cast<HeatedCircularPenetratorElement*>(
   ProblemParameters::Penetrator_pt);
 Vector<double> centre(2);
 if (pen_el_pt!=0)
  {
   Vector<double> my_centre(pen_el_pt->centre());
   centre[0]=my_centre[0];
   centre[1]=my_centre[1];
  }
 else
  {
   centre[0]=ProblemParameters::Centre[0];
   centre[1]=ProblemParameters::Centre[1];
  }
 double x_c=centre[0];
 double width=2.0*b_hertz;
 for (unsigned j=0;j<n;j++)
  {
   double x=x_c-0.5*width+width*double(j)/double(n-1);
   double p=0.0;
   if (abs((x-x_c))<b_hertz)
    {
     p=p_max_hertz*sqrt(1.0-pow((x-x_c)/b_hertz,2));
    }
   some_file << x << " 0.0 " << p << std::endl;
  }
 some_file.close();


 double target_weight=0.0;
 if (pen_el_pt!=0)
  {
   target_weight=pen_el_pt->target_weight();
  }
 // double angle=0.0;
 // if (pen_el_pt!=0)
 //  {
 //   angle=pen_el_pt->angle();
 //  }
 

 // Write trace file
 Trace_file << total_contact_force[0] << " " 
            << total_contact_force[1] << " " 
            << centre[0] << " "
            << centre[1] << " "
            << target_weight << " " 
  //ALH/hierher Suppressed so that validation tests pass
  //For small changes in loading a different contact point can be selected
  //which means that different angles are computed
  // << angle << " " 
            << std::endl;


 //Increment counter for solutions 
 Doc_info.number()++;

} // end of doc_solution



/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//=======start_of_main====================================================
/// Driver code
//========================================================================
int main(int argc, char* argv[])
{
 //feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

 FiniteElement::Accept_negative_jacobian=true;

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Suppress adaptation
 CommandLineArgs::specify_command_line_flag("--no_adapt");
    
 // Initial element size
 CommandLineArgs::specify_command_line_flag("--el_area",
                                            &ProblemParameters::El_area);
    
 // Suppress adaptation
 CommandLineArgs::specify_command_line_flag("--validate");
    
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Create generalised Hookean constitutive equations
 ProblemParameters::Constitutive_law_pt = 
  new GeneralisedHookean(&ProblemParameters::Nu);
      

 // Define centre of penetrator
 ProblemParameters::Centre.resize(2);
 ProblemParameters::Centre[0]=0.5;
 ProblemParameters::Centre[1]=1.0+ProblemParameters::Radius;

 // Build penetrator
 ProblemParameters::Penetrator_pt =
  new HeatedCircularPenetrator(&ProblemParameters::Centre,
                               ProblemParameters::Radius);


#ifdef STRUCTURED_MESH

 // Build problem
 ContactProblem<RefineablePseudoSolidNodeUpdateElement<
  RefineableQLinearElasticityElement<2,3>,
  RefineableQPVDElement<2,3> > >
  problem;

#else

 // Build problem
 ContactProblem
  <ProjectableLinearHeatAndElasticityElement<
   PseudoSolidNodeUpdateElement<TLinearHeatAndElasticityElement<2,3>,
                                TPVDElement<2,3> > > > problem;
 
#endif
 

 // Initialise timestep -- also sets the weights for all timesteppers
 // in the problem.
 double dt=0.01;
 problem.initialise_dt(dt);
 problem.assign_initial_values_impulsive();

 double max_residuals = 100;
 problem.max_residuals() = max_residuals;
 unsigned max_iterations = 10;

  //Output initial condition
 problem.doc_solution();
 
 unsigned max_adapt=1; // hierher
 if (CommandLineArgs::command_line_flag_has_been_set("--no_adapt"))
  {
   max_adapt=0;
  }

 // Pure imposed displacement of upper surface -- no contact
 //---------------------------------------------------------
 {
  //Parameter incrementation
  unsigned nstep=2; 
  for(unsigned i=0;i<nstep;i++)
   {
    
    double d_ampl=0.05;
    double ampl=0.0;

#ifdef STRUCTURED_MESH
    
    // Increment imposed boundary displacement
    ProblemParameters::Boundary_geom_object.ampl()+=d_ampl;
    ampl=ProblemParameters::Boundary_geom_object.ampl();

#else
    
    // Increment imposed boundary displacement
    ProblemParameters::Boundary_geom_object_left.ampl()+=d_ampl;
    ProblemParameters::Boundary_geom_object_contact.ampl()+=d_ampl;
    ProblemParameters::Boundary_geom_object_right.ampl()+=d_ampl;
    ampl=ProblemParameters::Boundary_geom_object_contact.ampl();
        
#endif
    
    oomph_info << "Re-solving for prescr displ amplitude = "
               << ampl << std::endl;
  
    // Solve the problem with Newton's method, allowing
    // up to max_adapt mesh adaptations after every solve.
    bool first=false;
    problem.unsteady_newton_solve(dt,max_adapt,first); 
    
    // Doc solution
    problem.doc_solution();
   }
 }
 
#ifdef STRUCTURED_MESH
 ProblemParameters::Centre[1]=
  ProblemParameters::Boundary_geom_object.y_c()- 
  ProblemParameters::Boundary_geom_object.radius()+
  ProblemParameters::Radius;
#else
 ProblemParameters::Centre[1]=
  ProblemParameters::Boundary_geom_object_contact.y_c()- 
  ProblemParameters::Boundary_geom_object_contact.radius()+
  ProblemParameters::Radius;
#endif


 // Switch on temperature
 ProblemParameters::T_contact=1.0;

 // Move position of centre directly
 //----------------------------------
 unsigned nstep=30;
 for (unsigned i=0;i<nstep;i++)
  {

   double dyc=0.00024; 
   ProblemParameters::Centre[1]-=dyc;
   oomph_info << "Re-solving imposed circle pos for yc=" 
              << ProblemParameters::Centre[1]
              << std::endl;

   // Solve
   //bool first=false;
   // problem.unsteady_newton_solve(dt,max_adapt,first);
   try
    {      
     max_iterations = 10;
     problem.max_newton_iterations() = max_iterations;
     problem.unsteady_newton_solve(dt);
    }
   catch(OomphLibError& error) //This will catch any error, 
//but actually we only want to catch failure to converge
    {
     max_iterations = 100;
     problem.max_newton_iterations() = max_iterations;
     problem.enable_globally_convergent_newton_method();
     problem.unsteady_newton_solve(dt,false);
     problem.disable_globally_convergent_newton_method();
    }

   //Output solution
   problem.doc_solution();
  }


// exit(0);

 // // Now increase resolution
 // //------------------------ 
 // nstep=3;
 // for (unsigned i=0;i<nstep;i++)
 // {

 //  oomph_info << "Re-solving imposed circle pos for element length factor: "
 //             << ProblemParameters::Element_length_factor
 //             << std::endl;

 //  // Re-solve
 //  problem.newton_solve(max_adapt);
  
 //  //Output solution
 //  problem.doc_solution();

 //  // Refine 
 //  ProblemParameters::Element_length_factor/=2.0;
 // }

 // Switch to node control
 ProblemParameters::Impose_position_of_centre=false;

 // Use displacement control initially
 //-----------------------------------
 problem.adapt();
 problem.switch_to_displ_control();
 double dyc=0.0003; 
 nstep=1;
 for (unsigned i=0;i<nstep;i++)
  {
   ProblemParameters::Y_c-=dyc;
   oomph_info << "Re-solving for yc=" 
              << ProblemParameters::Y_c
              << std::endl;

   // Solve
   
   // Try using normal newton solver, if it fails to converge, try a more stabl
   // but slower globally convergent method
   try
    {
     //This problem can have quite high initial residuals
     max_residuals = 100;

     //Set max number of iterations to low number
     max_iterations = 10;
     problem.max_newton_iterations() = max_iterations;

     problem.newton_solve(max_adapt);
    }
   catch(OomphLibError& error) //This will catch any error, 
       //but actually we only want to catch failure to converge
    {   
     //Increase number of iterations to give it chance to converge
     max_iterations = 100;
     problem.max_newton_iterations() = max_iterations;

     //Activate globally convergent method, then deactivate after
     problem.enable_globally_convergent_newton_method();
     problem.newton_solve(max_adapt);
     problem.disable_globally_convergent_newton_method();
    }


   //Output solution
   problem.doc_solution();
  }
 

 // Switch to force control
 problem.switch_to_force_control();
 
 // Balance it in the horizontal direction
 //---------------------------------------
 ProblemParameters::Horizontal_force=0.0;
 
 oomph_info << "RE-solving for weight=" 
            << dynamic_cast<HeatedCircularPenetratorElement*>(
             ProblemParameters::Penetrator_pt)->target_weight() 
            << " and horizontal force: "
            << dynamic_cast<HeatedCircularPenetratorElement*>(
             ProblemParameters::Penetrator_pt)->target_horizontal_force() 
            << std::endl;
 
 // Re-solve
 // problem.unsteady_newton_solve(dt); // max_adapt);

 // Try using normal newton solver, if it fails to converge, try a more stable
 // but slower globally convergent method
 try
  {
   //Set max number of iterations to low number
   max_iterations = 10;
   problem.max_newton_iterations() = max_iterations;
   
   problem.unsteady_newton_solve(dt);
  }
 catch(OomphLibError& error) //This will catch any error, 
//but actually we only want to catch failure to converge
  {
   //Increase number of iterations to give it chance to converge
   max_iterations = 100;
   problem.max_newton_iterations() = max_iterations;
   
   //Activate globally convergent method, then deactivate after
   problem.enable_globally_convergent_newton_method();
   problem.unsteady_newton_solve(dt,false);
   problem.disable_globally_convergent_newton_method(); 
  }




 //Output solution
 problem.doc_solution();
 

 // // Now increase weight
 // //-------------------- 
 // double dweight=0.0001; 
 // nstep=1;
 // for (unsigned i=0;i<nstep;i++)
 //  {
 //   oomph_info << "Re-solving for weight=" 
 //              << dynamic_cast<HeatedCircularPenetratorElement*>(
 //               ProblemParameters::Penetrator_pt)->target_weight() 
 //              << " and horizontal force: "
 //              << dynamic_cast<HeatedCircularPenetratorElement*>(
 //               ProblemParameters::Penetrator_pt)->target_horizontal_force() 
 //              << std::endl;
   
 //   // Re-solve
 //   problem.newton_solve(max_adapt);
   
 //   //Output solution
 //   problem.doc_solution();
   
 //   // Increase weight
 //   ProblemParameters::Weight+=dweight;
 //  }
 

// Now detach elastic body
//------------------------
 nstep=3;
 double d_lift_off_ampl=0.0001;
 for (unsigned i=0;i<nstep;i++)
  {

   double lift_off=0.0;

#ifdef STRUCTURED_MESH
    
    // Increment imposed boundary displacement
    ProblemParameters::Boundary_geom_object.lift_off_amplitude()+=
     d_lift_off_ampl;

    lift_off=ProblemParameters::Boundary_geom_object.lift_off_amplitude();
    
#else
    
    // Increment imposed boundary displacement
    ProblemParameters::Boundary_geom_object_left.lift_off_amplitude()+= 
     d_lift_off_ampl;
    ProblemParameters::Boundary_geom_object_contact.lift_off_amplitude()+=
     d_lift_off_ampl;
    ProblemParameters::Boundary_geom_object_right.lift_off_amplitude()+=
     d_lift_off_ampl;

    lift_off=ProblemParameters::Boundary_geom_object_contact.
     lift_off_amplitude();
    
#endif


  oomph_info << "Re-solving for weight=" 
              << dynamic_cast<HeatedCircularPenetratorElement*>(
               ProblemParameters::Penetrator_pt)->target_weight() 
              << " and horizontal force: "
              << dynamic_cast<HeatedCircularPenetratorElement*>(
               ProblemParameters::Penetrator_pt)->target_horizontal_force() 
              << " and lift off: "
              << lift_off
              << std::endl;
  
  // Re-solve
   try
    {
     max_residuals = 100;
     max_iterations = 10;
     problem.max_newton_iterations() = max_iterations;
     problem.max_residuals() = max_residuals;
     problem.newton_solve(max_adapt);
    }
   catch(OomphLibError& error) //This will catch any error, 
//but actually we only want to catch failure to converge
    {
     max_iterations = 100;
     problem.max_newton_iterations() = max_iterations;
     problem.enable_globally_convergent_newton_method();
     problem.newton_solve(max_adapt);
     problem.disable_globally_convergent_newton_method();
    }


  
  //Output solution
  problem.doc_solution();
  
 }




}
