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
//Header file for my streamfunction problem class

//Needs to include the refineable streamfunction element class 
#include "refineable_streamfunction_elements.h"

namespace oomph
{

//======================================================================
///  My streamfunction mesh class
//======================================================================
class streamfunction_mesh : public virtual Refineable_r_mesh<RefineablePolarStreamfunctionElement>
{
protected:

 //Storage for pointers to my traction elements
 Vector<PolarStreamfunctionTractionElement<RefineablePolarStreamfunctionElement>*> Inlet_traction_elt_pt;
 Vector<PolarStreamfunctionTractionElement<RefineablePolarStreamfunctionElement>*> Outlet_traction_elt_pt;

public:

 /// Return  pointer to inlet traction element e
 PolarStreamfunctionTractionElement<RefineablePolarStreamfunctionElement>* inlet_traction_elt_pt(unsigned e) 
   {return Inlet_traction_elt_pt[e];}
 /// Return length of inlet traction element vector
 unsigned inlet_traction_elt_length() {return Inlet_traction_elt_pt.size();}
 /// Return  pointer to outlet traction element e
 PolarStreamfunctionTractionElement<RefineablePolarStreamfunctionElement>* outlet_traction_elt_pt(unsigned e) 
   {return Outlet_traction_elt_pt[e];}
 /// Return length of outlet traction element vector
 unsigned outlet_traction_elt_length() {return Outlet_traction_elt_pt.size();}

 /// Constructor, which "builds" the mesh. The arguments are the number
 /// of elements in each direction.
 streamfunction_mesh(const unsigned int &nx,const unsigned int &ny) : 
   Refineable_r_mesh<RefineablePolarStreamfunctionElement>(nx,ny)
  {  
   //Now bolt on traction stuff

   //Assign fluid elements to vector
   this->assign_fluid_element_vector();

   //We need traction elements at both inlet and outlet
   make_traction_elements(false);
   make_traction_elements(true);
  }

 //Function to add the traction boundary elements
 void make_traction_elements(const bool& outlet)
  {
   //Specify inlet/outlet specific quantities
   unsigned ibound; int index;
   if(outlet) {ibound=1;index=1;}
   else {ibound=3;index=-1;}

   unsigned num_elt = this->nboundary_element(ibound);

   //Loop over the number of elements on the boundary
   for(unsigned ielt=0;ielt<num_elt;ielt++)
    {
     PolarStreamfunctionTractionElement<RefineablePolarStreamfunctionElement> *surface_element_pt =
       new PolarStreamfunctionTractionElement<RefineablePolarStreamfunctionElement>
       (this->boundary_element_pt(ibound,ielt),index);
     //Push it back onto the Element_pt Vector
     this->Element_pt.push_back(surface_element_pt);

     if(outlet) { this->Outlet_traction_elt_pt.push_back(surface_element_pt); }
     else { this->Inlet_traction_elt_pt.push_back(surface_element_pt); }

     //Any other information to pass: 
     surface_element_pt->set_boundary(index);
     surface_element_pt->alpha_pt() = &Global_Physical_Variables::Alpha;
  
    }//End of loop over elements

   std::cout << std::endl << "Streamfunction traction elements attached to streamfunction_mesh" << std::endl << std::endl;

  }//End of make traction elements

 //Function to remove the traction boundary elements
 void remove_traction_elements()
  { 
   //Find the number of traction elements
   unsigned Ntraction = this->inlet_traction_elt_length();
   Ntraction += this->outlet_traction_elt_length();

   //The traction elements are ALWAYS? stored at the end
   //So delete and remove them
   for(unsigned e=0;e<Ntraction;e++)
    {
     delete this->Element_pt.back();
     this->Element_pt.pop_back();
    }

   //Now clear the vectors of pointers to traction elements
   Inlet_traction_elt_pt.clear();
   Outlet_traction_elt_pt.clear();

   std::cout << std::endl << "Streamfunction traction elements removed from streamfunction_mesh" << std::endl << std::endl;

  }//End of remove_traction_elements

  //Function to put current fluid elements into a vector of their own
  void assign_fluid_element_vector()
  {
   unsigned check=inlet_traction_elt_length();
   check+=outlet_traction_elt_length();

   this->Fluid_elt_pt.clear();

   if(check!=0)
    {
     std::cout << "Warning, attempting to assemble streamfunction element vector ";
     std::cout << "whilst streamfunction traction elements are attached to the mesh" << std::endl;
    }
   else
    {
     unsigned n_elt = this->Element_pt.size();
     for(unsigned e=0;e<n_elt;e++)
      {
       this->Fluid_elt_pt.push_back(this->Element_pt[e]);
      }
    }
  }//End of assign_fluid_element_vector

}; //End of streamfunction_mesh class

//====================================================================
/// Polar Streamfunction problem class 
//====================================================================

//====== start_of_problem_class=======================================
/// 2D Streamfunction problem on rectangular domain
//====================================================================
class StreamfunctionProblem : public Problem
{
protected:

 /// Storage for whether we're solving for an eigenfunction or not
 bool eigenfunction;

public:

 /// Constructor: Pass pointer to source function
 StreamfunctionProblem(const bool& eigen);

 /// Destructor (empty -- all the cleanup is done in base class)
 ~StreamfunctionProblem(){};

 /// Update the problem specs before solve: Reset boundary conditions
 /// to the values from the exact solution.
 void actions_before_solve();

 /// Update the problem after solve (empty)
 void actions_after_solve(){}

 /// Before adaptation: 
 void actions_before_adapt()
  {
   // Strip traction elements from mesh before refine
   mesh_pt()->remove_traction_elements();
  }

 /// After adaptation: Unpin all pressures and then pin redudant pressure dofs.
 void actions_after_adapt()
  {
   //Reassign fluid elements to storage
   mesh_pt()->assign_fluid_element_vector();

   //Add traction elements back onto the mesh
   mesh_pt()->make_traction_elements(false);
   mesh_pt()->make_traction_elements(true);

   //Repin the velocities stored at all nodes
   pin_velocities();

   //Repin boundaries?
   //pin_boundaries();
  }

 // Access function for the specific mesh
  streamfunction_mesh* mesh_pt() 
  {
   // Upcast from pointer to the Mesh base class to the specific 
   // element type that we're using here.
   return dynamic_cast<streamfunction_mesh* >(Problem::mesh_pt());
  }

 /// Pin boundaries
 void pin_boundaries();

 /// Pin velocities
 void pin_velocities();

 /// Assign velocities to nodes
 void assign_velocities();

 /// Return the value of the streamfunction at given global coordinates
 void get_vels(const Vector<double>& x_to_get, double& stream);

 /// Output the streamfunction on a regular grid
 void my_output(const int radial,
                const int azimuthal,
                bool log_output,string file_name);
 
 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info);

 // Header for doc_solution
 void header( ofstream &some_file );

}; // end of problem class

//=====start_of_constructor===============================================
/// Constructor for Streamfunction problem: Pass pointer to source function.
//========================================================================
StreamfunctionProblem::StreamfunctionProblem(const bool& eigen)
{ 
 // Are we solving for an eigenfunction?
 this->eigenfunction=eigen;

 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=Global_Physical_Variables::xmesh;
 unsigned n_y=Global_Physical_Variables::ymesh;

 // Build and assign mesh 
 //Problem::mesh_pt() = new streamfunction_mesh<RefineablePolarStreamfunctionElement>(n_x,n_y);
 Problem::mesh_pt() = new streamfunction_mesh(n_x,n_y);

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;

 // Pin the relevant boundaries 
 pin_boundaries();

 // Now loop over all nodes and pin the nodal dofs storing velocities
 pin_velocities();

 // Complete the build of all elements so they are fully functional

 //Find number of fluid elements in mesh
 unsigned n_fluid = mesh_pt()->fluid_elt_length();
 for(unsigned i=0;i<n_fluid;i++)
  {
   // Upcast from GeneralsedElement to the present element
   RefineablePolarStreamfunctionElement *el_pt = dynamic_cast<RefineablePolarStreamfunctionElement*>(mesh_pt()->element_pt(i));

   //Set the pointer to the angle Alpha
   el_pt->alpha_pt() = &Global_Physical_Variables::Alpha;
  }

 //Find number of traction elements in mesh

 // Setup equation numbering scheme
 std::cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor

//========================================start_of_actions_before_solve===
/// Update the problem specs before solve: (Re-)set boundary conditions
/// to the values from the exact solution.
//========================================================================
void StreamfunctionProblem::actions_before_solve()
{
 // How many boundaries are there?
 unsigned n_bound = mesh_pt()->nboundary();
 
 //Loop over the side boundaries
 for(unsigned i=0;i<n_bound;i+=2)
  {
   // How many nodes are there on this boundary?
   unsigned n_node = mesh_pt()->nboundary_node(i);

   // Loop over the nodes on boundary
   for (unsigned n=0;n<n_node;n++)
    {
     // Get pointer to node
     Node* nod_pt=mesh_pt()->boundary_node_pt(i,n);

     // Extract nodal coordinates from node:
     Vector<double> x(2);
     x[0]=nod_pt->x(0);
     x[1]=nod_pt->x(1);

     double stream = 0.5*x[1];
     // Zero boundary conditions for eigenproblem
     if(this->eigenfunction) stream=0.0;

     // Assign values to the first nodal value
     nod_pt->set_value(0,stream);

    }
  } 
}  // end of actions before solve

//========================================================================
/// pin boundaries
//========================================================================
void StreamfunctionProblem::pin_boundaries()
{
 unsigned n_bound = mesh_pt()->nboundary();
 //for(unsigned i=0;i<n_bound;i++)
 for(unsigned i=0;i<n_bound;i+=2)
  {
   unsigned n_node = mesh_pt()->nboundary_node(i);
   for (unsigned n=0;n<n_node;n++)
    {
     //Pin both streamfunction and vorcitity at boundaries
     mesh_pt()->boundary_node_pt(i,n)->pin(0); 
    }
  }

} // end of pin_boundaries

//========================================================================
/// pin velocities at all nodes
//========================================================================
void StreamfunctionProblem::pin_velocities()
{
 // Loop over all nodes and asign values for velocites
 unsigned num_nod = mesh_pt()->nnode();
 for(unsigned i=0;i<num_nod;i++)
  {
   mesh_pt()->node_pt(i)->pin(1);
   mesh_pt()->node_pt(i)->pin(2);
  }

} // end of pin_velocities

//========================================================================
/// Assign velocities to nodes
//========================================================================
void StreamfunctionProblem::assign_velocities()
{
 // Loop over all nodes and asign values for velocites
 unsigned num_nod = mesh_pt()->nnode();
 for(unsigned i=0;i<num_nod;i++)
  {
    // Extract nodal coordinates from node:
    Vector<double> x(2);
    x[0]=mesh_pt()->node_pt(i)->x(0);
    x[1]=mesh_pt()->node_pt(i)->x(1);

    //double A = Global_Physical_Variables::Alpha;

    //double u=(x[0]/A)*cos(x[1])-(1./(4.*A));
    //double v=-2.*x[0]*sin(x[1])+0.25*x[1];
    //A=0.0;
    double u=0.0;double v=0.0;

    mesh_pt()->node_pt(i)->set_value(1,u);
    mesh_pt()->node_pt(i)->set_value(2,v);
  }

} // end of assign_velocities

//========================================================================
// Rich's function for locating the element which contains given 
// coordinates and returning u,v there
// Note that this function currently fails at r and theta max
//========================================================================
void StreamfunctionProblem::get_vels(const Vector<double>& x_to_get, double& stream)
{
  // local coord
  Vector<double> s( 2, 0.0 );
  // global coords
  Vector<double> x_lower_left( 2, 0.0 );
  Vector<double> x_centre( 2, 0.0 );
  Vector<double> x_upper_right( 2, 0.0 );
  Vector<double> x( 2, 0.0 );
  double a,b,c,sub;
  //Reset interpolated velocity - check they're being set at some point
  stream=-2.;

  // Loop over all the (fluid) elements 
  unsigned long n_element = mesh_pt() -> nelement();
  for ( unsigned long e = 0; e < n_element; e++ )
	{
	  // Cast to the particular element type, this is necessary because
	  // the base elements don't have the member functions that we're about
	  // to call!
	  RefineablePolarStreamfunctionElement* elt_pt = dynamic_cast<RefineablePolarStreamfunctionElement*>( mesh_pt() -> element_pt(e) );
	  
	  s[0] = -1.0; s[1] = -1.0;
	  elt_pt -> get_x( s, x_lower_left );

	  s[0] = 1.0; s[1] = 1.0;
	  elt_pt -> get_x( s, x_upper_right );
	  
	  // check the elts lower left corner is SW of point to return
	  if ( ( x_lower_left[0] <= x_to_get[0] ) 
		   && ( x_lower_left[1] <= x_to_get[1] ) )
		{
		  // check the elts upper right corner is NE of point to return
		  if ( ( x_upper_right[0] > x_to_get[0] ) 
			   && ( x_upper_right[1] > x_to_get[1] ) )
			{
			  //Quadratic interpolation:
			  a=x_lower_left[0];
			  s[0] = 0.0; s[1] = 0.0;
			  elt_pt -> get_x( s, x_centre );
			  b=x_centre[0];
			  c=x_upper_right[0];
			  sub=pow(c-a,2.0)-8.0*(a+c-2.0*b)*(b-x_to_get[0]);

			  if(abs(a+c-2.*b)<1.e-8)s[0] = -1.0 + 2.*( x_to_get[0] - x_lower_left[0] ) 
						/ ( x_upper_right[0] - x_lower_left[0] );
			  else s[0]=((a-c)+sqrt(sub))/(2.0*(a+c-2.0*b));

			  // our point to get is in the current elt
			  // now work out the local coordinate
			  // This only works for evenly spaced global
			  // coordinates across the element - Not true 
			  // in radial direction. 
			  s[1] = -1.0 + 2.*( x_to_get[1] - x_lower_left[1] ) 
			  	/ ( x_upper_right[1] - x_lower_left[1] );

			  // check on my quadratic interpolation
			  elt_pt -> get_x( s, x );

			  if(abs(x[0]-x_to_get[0])>1.e-12)
			   {
			    if(abs(a+c-2.*b)<1.e-8) std::cout << "Element almost linear: " << abs(a+c-2.*b) << std::endl;
			    std::cout << "error in interpolation: " << (x[0]-x_to_get[0]) << std::endl;
			    std::cout << "r:" << x[0] << "r to get: " << x_to_get[0] << std::endl;
			   }

			  // get the position of this
			  elt_pt -> get_x( s, x );
			  stream = elt_pt -> interpolated_streamfunction( s );
			
			  break;
			}
		}
	  //return vels;
	}
}

//========================================================================
// My function for outputting the streamfunction using 
// Rich's get_vels() function
//========================================================================
void StreamfunctionProblem::my_output(const int radial,const int azimuthal,bool log_output,string file_name)
{
  // Quick adjustment so that I don't ever end up exactly on a boundary
  // (get_vels fails on certain boundaries)
  // (current difference between x_to_get and x_upper_right around 7.e-15 so
  // 1.e-13 should be fine...)
  double inner=Global_Physical_Variables::R_l+1.e-13;
  double outer=Global_Physical_Variables::R_r-1.e-13;
  // And the same for azimuthal direction
  double left=-1.+1.e-13;
  double right=1.-1.e-13;

  // Set up common ratio for log spaced output
  double R = pow((outer/inner),(1.0/(radial - 1.0)));

  // New output file
  ofstream out;
  out.open(file_name.c_str());

  out << "# Streamfunction output at " << radial << " radial and " << azimuthal << " azimuthal points" << std::endl;
  out << "# Re = " << Global_Physical_Variables::Re << " Alpha = " << Global_Physical_Variables::Alpha 
      << " R_l = " << Global_Physical_Variables::R_l << std::endl;
  out << "# Initial xmesh: " << Global_Physical_Variables::xmesh << " Initial ymesh " << Global_Physical_Variables::ymesh << std::endl;
  out << "# Uniform refinements: " << Global_Physical_Variables::uniform << " Adaptive refinements: " << Global_Physical_Variables::adaptive << std::endl;
  out << "# log_mesh = " << Global_Physical_Variables::log_mesh << " new_outlet_region = " << Global_Physical_Variables::new_outlet_region << std::endl;
  out << "# Inlet traction: " << Global_Physical_Variables::inlet_traction << " Outlet traction: " << Global_Physical_Variables::outlet_traction << std::endl;
  out << "# eta_inlet: " << Global_Physical_Variables::eta_inlet << " eta_outlet: " << Global_Physical_Variables::eta_outlet << std::endl;
  out << "# Pin v at inlet: " << Global_Physical_Variables::pinv << std::endl;
  out << std::endl;

  std::cout << endl;
  std::cout << "Outputting streamfunction to file " << file_name << " at " << radial 
	    << " radial plot points and " << azimuthal << " azimuthal points." << std::endl;
  std::cout << std::endl;
  
  //Position of elements
  Vector<double> position(2,0.0);
  double stream=0.0;
  
  //Loops now radial inside azimuthal to agree with my poisson solver
  if(log_output)
    {
      for(int i=0;i<radial;i++)
	{
	  // Update radial position
	   position[0]=inner*pow(R,i);
	   
	   for(int j=0;j<azimuthal;j++)
	     {
	       // Update azimuthal position
	       position[1]=left + ((right-left)*j)/(azimuthal-1.0);
	       
	       // Read off the velocity at this point
	       get_vels(position,stream);
	       
	       // Output to file
	       out << position[0] << " " << position[1] << " " << stream << std::endl;
	       
	     } // End of loop over azimuthal
	   
	   out << std::endl;
	   
	} // End of loop over radial
      
    } // End of log spaced output loop
   else
     {
       for(int i=0;i<radial;i++)
	 {
	   // Update radial position
	   position[0]=inner + ((outer-inner)*i)/(radial-1.0);
	   
	   for(int j=0;j<azimuthal;j++)
	     {
	       // Update azimthual position
	       position[1]=left + ((right-left)*j)/(azimuthal-1.0);
	       
	       // Read off the velocity at this point
	       get_vels(position,stream);
	       
	       // Output to file
	       out << position[0] << " " << position[1] << " " << stream << std::endl;

	     } // End of loop over azimuthal

	   out << std::endl;

	 } // End of loop over radial
        
     } // End of regular spaced output loop

  out.close();
  
}

//===============start_of_doc=============================================
/// Doc the solution: doc_info contains labels/output directory etc.
//========================================================================
void StreamfunctionProblem::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=3;

 // Output solution 
 //-----------------
 sprintf(filename,"%s/streamfunction%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 header(some_file);
 mesh_pt()->output(some_file,npts);
 some_file.close();

} // end of doc

void StreamfunctionProblem::header( ofstream &some_file )
 {
  some_file << "# Refineable streamfunction mesh" << "\n";
  some_file << "# Re = " << Global_Physical_Variables::Re << " Alpha = " << Global_Physical_Variables::Alpha
	    << " R_l = " << Global_Physical_Variables::R_l << "\n";
  some_file << "# Initial xmesh = " << Global_Physical_Variables::xmesh 
	    << " Initial ymesh = " << Global_Physical_Variables::ymesh << "\n";
  some_file << "# Uniform refinements: " << Global_Physical_Variables::uniform 
	    << " Adaptive refinements: " << Global_Physical_Variables::adaptive << "\n";
  some_file << "# inlet_traction = " << Global_Physical_Variables::inlet_traction 
	    << " eta_inlet = " << Global_Physical_Variables::eta_inlet << "\n";
  some_file << "# outlet_traction = " << Global_Physical_Variables::outlet_traction 
	    << " eta_outlet = " << Global_Physical_Variables::eta_outlet << "\n";
  some_file << "# log_mesh = " << Global_Physical_Variables::log_mesh 
	    << " new_outlet_region = " << Global_Physical_Variables::new_outlet_region << "\n";
  some_file << "# pinv = " << Global_Physical_Variables::pinv << "\n";
  some_file << "\n";
 }

//====================================================================
/// End of Polar Streamfunction problem class 
//====================================================================

} //End of namespace oomph
