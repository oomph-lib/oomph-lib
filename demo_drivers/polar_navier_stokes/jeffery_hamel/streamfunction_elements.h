//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
//Header file for Streamfunction elements
#ifndef OOMPH_POLAR_STREAMFUNCTION
#define OOMPH_POLAR_STREAMFUNCTION

namespace oomph
{

//=============================================================
/// A class for solving my poisson like streamfunction problem
//=============================================================
class PolarStreamfunctionEquations : public virtual FiniteElement
{
protected:

 /// Pointer to the angle alpha
 double *Alpha_pt;

public:

 /// Alpha
 const double &alpha() const {return *Alpha_pt;}

 /// Pointer to Alpha
 double* &alpha_pt() {return Alpha_pt;}

 /// Constructor (must initialise the Source_fct_pt to null)
 PolarStreamfunctionEquations(){}
 
 /// Broken copy constructor
 PolarStreamfunctionEquations(const PolarStreamfunctionEquations& dummy) 
  { 
   BrokenCopy::broken_copy("PolarStreamfunctionEquations");
  } 
 

 /// Return the indicies at which the unknown values are stored. 
 virtual inline unsigned u_index_streamfunction() const {return 0;}

 /// Return the indicies at which the (known) velocities are stored.
 virtual inline unsigned u_index_velocity(const unsigned &i) const {return i+1;}

 /// Output with default number of plot points
 void output(std::ostream &outfile) 
  {
   const unsigned n_plot=5;
   output(outfile,n_plot);
  }

 /// Output FE representation of soln: x,y,u or x,y,z,u at 
 /// n_plot^DIM plot points
 void output(std::ostream &outfile, const unsigned &n_plot);
 
 /// C_style output with default number of plot points
 void output(FILE* file_pt)
  {
   const unsigned n_plot=5;
   output(file_pt,n_plot);
  }

 /// C-style output FE representation of soln: x,y,u or x,y,z,u at 
 /// n_plot^DIM plot points
 void output(FILE* file_pt, const unsigned &n_plot);

 /// Get flux: Necessary for Z2 error estimator in refineable version
 /// Doesn't matter what it is though as we're NEVER? going to refine
 void get_flux(const Vector<double>& s, Vector<double>& flux) const
  {
   //Find out how many nodes there are in the element
   const unsigned n_node = nnode();

   //Get the index at which the unknown is stored
   const unsigned u_nodal_index = u_index_streamfunction();

   //Set up memory for the shape and test functions
   Shape psi(n_node);
   DShape dpsidx(n_node,2);
 
   //Call the derivatives of the shape and test functions
   dshape_eulerian(s,psi,dpsidx);
     
   //Initialise to zero
   for(unsigned j=0;j<2;j++)
    {
     flux[j] = 0.0;
    }
   
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     //Loop over derivative directions
     for(unsigned j=0;j<2;j++)
      {                               
       flux[j] += this->nodal_value(l,u_nodal_index)*dpsidx(l,j);
      }
    }
  }

 /// Add the element's contribution to its residual vector (wrapper)
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution(
    residuals,GeneralisedElement::Dummy_matrix,0);
  }

 
 /// Add the element's contribution to its residual vector and 
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution(residuals,jacobian,1);
  }
 

 /// Return FE representation of streamfunction 
 /// at local coordinate s
 inline double interpolated_streamfunction(const Vector<double> &s) const
  {
   //Find number of nodes
   const unsigned n_node = nnode();

   //Get the index at which the poisson unknown is stored
   const unsigned u_nodal_index = u_index_streamfunction();
   
   //Local shape function
   Shape psi(n_node);

   //Find values of shape function
   shape(s,psi);

   //Initialise value of u
   double interpolated_u = 0.0;

   //Loop over the local nodes and sum
   for(unsigned l=0;l<n_node;l++) 
    {
     interpolated_u += this->nodal_value(l,u_nodal_index)*psi[l];
    }

   return(interpolated_u);
  }

 /// Return FE interpolated velocity u[i] at local coordinate s
 double interpolated_velocity(const Vector<double> &s, const unsigned &i) const
  {
   //Find number of nodes
   unsigned n_node = nnode();

   //Get the index at which the poisson unknown is stored
   const unsigned u_nodal_index = u_index_velocity(i);

   //Local shape function
   Shape psi(n_node);
   //Find values of shape function
   shape(s,psi);
   
   //Initialise value of u
   double interpolated_u = 0.0;
   //Loop over the local nodes and sum
   for(unsigned l=0;l<n_node;l++) 
    {
     interpolated_u += this->nodal_value(l,u_nodal_index)*psi[l];
    }
   
   return(interpolated_u);
  }

 /// Return FE interpolated velocity derivative du[i]/dx[j]     
 /// at local coordinate s                                      
 double interpolated_dudx(const Vector<double> &s, const unsigned &i,const unsigned &j) const
  {
   //Find number of nodes
   unsigned n_node = nnode(); 

   //Get the index at which the poisson unknown is stored
   const unsigned u_nodal_index = u_index_velocity(i);

   //Set up memory for the shape and test functions
   Shape psi(n_node);
   DShape dpsidx(n_node,2);

   //double J = 
   dshape_eulerian(s,psi,dpsidx);

   //Initialise to zero
   double interpolated_dudx = 0.0;
   
   //Calculate velocity derivative:

   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {                              
     interpolated_dudx += this->nodal_value(l,u_nodal_index)*dpsidx(l,j);
    }
  
   return(interpolated_dudx);
  }

 /// Construct and Return FE representation of vorticity 
 /// at local coordinate s
 inline double interpolated_vorticity(const Vector<double> &s) const
  {
   //Find number of nodes
   unsigned n_node = nnode(); 

   //Get Alpha
   const double Alpha = alpha();

   //Local shape function
   Shape psi(n_node);
   //Find values of shape function
   shape(s,psi);

   //Allocate and initialise to zero
   double r = 0.0;
 
   // Loop over nodes to calculate r
   for(unsigned l=0;l<n_node;l++) 
    {
       r += nodal_position(l,0)*psi(l);
    }

   //Initialise value of u
   double interpolated_u = 0.0;

   interpolated_u += interpolated_dudx(s,1,0);
   interpolated_u += (1./r)*interpolated_velocity(s,1);
   interpolated_u -= (1./(r*Alpha))*interpolated_dudx(s,0,1);

   return(interpolated_u);
  }

protected:

 /// Shape/test functions and derivs w.r.t. to global coords at 
 /// local coord. s; return  Jacobian of mapping
 virtual double dshape_and_dtest_eulerian_poisson(const Vector<double> &s, 
                                                  Shape &psi, 
                                                  DShape &dpsidx, Shape &test, 
                                                  DShape &dtestdx) const=0;
 

 /// Shape/test functions and derivs w.r.t. to global coords at 
 /// integration point ipt; return  Jacobian of mapping
 virtual double dshape_and_dtest_eulerian_at_knot_poisson(const unsigned &ipt, 
                                                          Shape &psi, 
                                                          DShape &dpsidx,
                                                          Shape &test, 
                                                          DShape &dtestdx) 
  const=0;

 /// Compute element residual Vector only (if flag=and/or element 
 /// Jacobian matrix 
 virtual void fill_in_generic_residual_contribution(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  unsigned flag); 
 
};

/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////

//======================================================================
/// StreamfunctionElements are quadrilateral elements 
/// with isoparametric interpolation for the function.
//======================================================================
 class PolarStreamfunctionElement : public virtual QElement<2,3>,
 public virtual PolarStreamfunctionEquations
{

private:

 /// Static int that holds the number of variables at 
 /// nodes: always the same
 static const unsigned Initial_Nvalue;
 
  public:


 /// Constructor: Call constructors for QElement and 
 /// Streamfunction equations
 PolarStreamfunctionElement() : QElement<2,3>(), PolarStreamfunctionEquations()
  {}
 
 /// Broken copy constructor
 PolarStreamfunctionElement(const PolarStreamfunctionElement& dummy) 
  { 
   BrokenCopy::broken_copy("PolarStreamfunctionElement");
  } 
 
 ///  Required  # of `values' (pinned or dofs) 
 /// at node n
 inline unsigned required_nvalue(const unsigned &n) const 
  {return Initial_Nvalue;}

 /// Output function:  
 ///  x,y,u   or    x,y,z,u
 void output(std::ostream &outfile)
  {PolarStreamfunctionEquations::output(outfile);}

 ///  Output function:  
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(std::ostream &outfile, const unsigned &n_plot)
  {PolarStreamfunctionEquations::output(outfile,n_plot);}

 /// C-style output function:  
 ///  x,y,u   or    x,y,z,u
 void output(FILE* file_pt)
  {PolarStreamfunctionEquations::output(file_pt);}

 ///  C-style output function:  
 ///   x,y,u   or    x,y,z,u at n_plot^DIM plot points
 void output(FILE* file_pt, const unsigned &n_plot)
  {PolarStreamfunctionEquations::output(file_pt,n_plot);}

 /// Output function for an exact solution:
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {PolarStreamfunctionEquations::output_fct(outfile,n_plot,exact_soln_pt);}

 /// Output function for a time-dependent exact solution.
 ///  x,y,u_exact   or    x,y,z,u_exact at n_plot^DIM plot points
 /// (Calls the steady version)
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 const double& time,
                 FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {PolarStreamfunctionEquations::output_fct(outfile,n_plot,time,exact_soln_pt);}


protected:

/// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
 inline double dshape_and_dtest_eulerian_poisson(
  const Vector<double> &s, Shape &psi, DShape &dpsidx, 
  Shape &test, DShape &dtestdx) const;


 /// Shape, test functions & derivs. w.r.t. to global coords. at
 /// integration point ipt. Return Jacobian.
 inline double dshape_and_dtest_eulerian_at_knot_poisson(const unsigned& ipt,
                                                         Shape &psi, 
                                                         DShape &dpsidx, 
                                                         Shape &test,
                                                         DShape &dtestdx) 
  const;

};




//Inline functions:


//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
 double PolarStreamfunctionElement::dshape_and_dtest_eulerian_poisson(
  const Vector<double> &s,
  Shape &psi, 
  DShape &dpsidx,
  Shape &test, 
  DShape &dtestdx) const
{
 //Call the geometrical shape functions and derivatives  
 const double J = this->dshape_eulerian(s,psi,dpsidx);

 //Set the test functions equal to the shape functions
 test = psi;
 dtestdx= dpsidx;
 
 //Return the jacobian
 return J;
}



//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
double PolarStreamfunctionElement::
 dshape_and_dtest_eulerian_at_knot_poisson(
  const unsigned &ipt,
  Shape &psi, 
  DShape &dpsidx,
  Shape &test, 
  DShape &dtestdx) const
{
 //Call the geometrical shape functions and derivatives  
 const double J = this->dshape_eulerian_at_knot(ipt,psi,dpsidx);

 //Set the pointers of the test functions
 test = psi;
 dtestdx = dpsidx;

 //Return the jacobian
 return J;
}

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//=======================================================================
/// Face geometry for the PolarStreamfunctionElement elements: The spatial 
/// dimension of the face elements is one lower than that of the
/// bulk element but they have the same number of points
/// along their 1D edges.
//=======================================================================
template<>
class FaceGeometry<PolarStreamfunctionElement>: public virtual QElement<1,3>
{

  public:
 
 /// Constructor: Call the constructor for the
 /// appropriate lower-dimensional QElement
 FaceGeometry() : QElement<1,3>() {}

};

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//======================================================================
/// Set the data for the number of Variables at each node, always one
/// in every case
//======================================================================
 const unsigned PolarStreamfunctionElement::Initial_Nvalue = 3;


//======================================================================
/// Compute element residual Vector and/or element Jacobian matrix 
/// 
/// flag=1: compute both
/// flag=0: compute only residual Vector
///
/// Pure version without hanging nodes
//======================================================================
void  PolarStreamfunctionEquations::
fill_in_generic_residual_contribution( Vector<double> &residuals, 
                                   DenseMatrix<double> &jacobian, 
                                   unsigned flag) 
{
 //Find out how many nodes there are
 const unsigned n_node = nnode();

 //Get Alpha
 const double Alpha = alpha();

 //Set up memory for the shape and test functions
 Shape psi(n_node), test(n_node);
 DShape dpsidx(n_node,2), dtestdx(n_node,2);

 //Indicies at which the unknowns are stored
 const unsigned s_nodal_index = u_index_streamfunction();

 //Find the indices at which the local velocities are stored
 unsigned u_nodal_index[2];
 for(unsigned i=0;i<2;i++) {u_nodal_index[i] = u_index_velocity(i);}
 
 //Set the value of n_intpt
 const unsigned n_intpt = integral_pt()->nweight();

 //Integers to store the local equation and unknown numbers
 int local_eqn=0, local_unknown=0;

 //Loop over the integration points
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {
   //Get the integral weight
   double w = integral_pt()->weight(ipt);

   //Call the derivatives of the shape and test functions
   double J = dshape_and_dtest_eulerian_at_knot_poisson(ipt,psi,dpsidx,
                                                        test,dtestdx);
       
   //Premultiply the weights and the Jacobian
   double W = w*J;

   //Allocate and initialise to zero
   Vector<double> interpolated_x(2,0.0);
   double interpolated_s=0.0;
   Vector<double> interpolated_dsdx(2,0.0);
   Vector<double> interpolated_u(2);
   DenseMatrix<double> interpolated_dudx(2,2);
   
   //Calculate function value and derivatives:
   //-----------------------------------------
   // Loop over nodes
   for(unsigned l=0;l<n_node;l++) 
    {
     //Get the nodal value of the poisson unknown
     double s_value = raw_nodal_value(l,s_nodal_index);
     interpolated_s += s_value*psi(l);
     // Loop over directions1
     for(unsigned i=0;i<2;i++)
      {
       interpolated_x[i] += raw_nodal_position(l,i)*psi(l);
       interpolated_dsdx[i] += s_value*dpsidx(l,i);
       double u_value = this->nodal_value(l,u_nodal_index[i]);
       interpolated_u[i] += u_value*psi(l);
       // Loop over directions2
       for(unsigned j=0;j<2;j++)
        {                               
         interpolated_dudx(i,j) += u_value*dpsidx(l,j);
        }
      }
    }

   // Assemble residuals and Jacobian
   //--------------------------------
       
   // Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
    {
      //Get the local equation
      local_eqn = nodal_local_eqn(l,s_nodal_index);
      /*IF it's not a boundary condition*/
      if(local_eqn >= 0)
       {
	// Laplacian of the streamfunction (integrated by parts)
	residuals[local_eqn] += interpolated_dsdx[0]*dtestdx(l,0)*interpolated_x[0]*Alpha*W;
	residuals[local_eqn] += (1./(interpolated_x[0]*Alpha))*interpolated_dsdx[1]
	                       *(1./(interpolated_x[0]*Alpha))*dtestdx(l,1)
	  *interpolated_x[0]*Alpha*W;

	// Should equal the vorticity
	residuals[local_eqn] -= (interpolated_dudx(1,0)+(interpolated_u[1]/interpolated_x[0]))*test(l)*interpolated_x[0]*Alpha*W;
	residuals[local_eqn] += (1./(interpolated_x[0]*Alpha))*interpolated_dudx(0,1)*test(l)*interpolated_x[0]*Alpha*W;

	// Calculate the jacobian
	//-----------------------
	if(flag)
	 {
          //Loop over the velocity shape functions again
          for(unsigned l2=0;l2<n_node;l2++)
           {
	     local_unknown = nodal_local_eqn(l2,s_nodal_index);
	     //If at a non-zero degree of freedom add in the entry
	     if(local_unknown >= 0)
              {
               //Add contribution to Elemental Matrix
                 jacobian(local_eqn,local_unknown) 
		   += dpsidx(l2,0)*dtestdx(l,0)*interpolated_x[0]*Alpha*W;
                 jacobian(local_eqn,local_unknown) 
		   += (1./(interpolated_x[0]*Alpha))*dpsidx(l2,1)*(1./(interpolated_x[0]*Alpha))*dtestdx(l,1)
			*interpolated_x[0]*Alpha*W;
	      } 

	   } // End of loop over l2
	 } // End of if Jacobian flag statement
       } // End of Residual if not boundary condition statement
      
    } // End of loop over test functions (l)

  } // End of loop over integration points
}   

//======================================================================
/// Output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
void  PolarStreamfunctionEquations::output(std::ostream &outfile, 
                                    const unsigned &nplot)
{

 //Vector of local coordinates
 Vector<double> s(2);
 
 // Tecplot header info
 outfile << tecplot_zone_string(nplot);
 
 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   
   for(unsigned i=0;i<2;i++) 
    {
     outfile << interpolated_x(s,i) << " ";
    }
   outfile << interpolated_streamfunction(s) << " "
	   << interpolated_vorticity(s) << " " 
	   << interpolated_velocity(s,0) << " "
	   << interpolated_velocity(s,1) << std::endl;   
   
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(outfile,nplot);

}

//======================================================================
/// C-style output function:
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
void  PolarStreamfunctionEquations::output(FILE* file_pt,
                                    const unsigned &nplot)
{
 //Vector of local coordinates
 Vector<double> s(2);
 
 // Tecplot header info
 fprintf(file_pt,"%s",tecplot_zone_string(nplot).c_str());

 // Loop over plot points
 unsigned num_plot_points=nplot_points(nplot);
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   // Get local coordinates of plot point
   get_s_plot(iplot,nplot,s);
   
   for(unsigned i=0;i<2;i++) 
    {
     fprintf(file_pt,"%g ",interpolated_x(s,i));
    }
   fprintf(file_pt,"%g ",interpolated_streamfunction(s));
   fprintf(file_pt,"%g \n",interpolated_vorticity(s));
  }

 // Write tecplot footer (e.g. FE connectivity lists)
 write_tecplot_zone_footer(file_pt,nplot);
}

}

#endif
