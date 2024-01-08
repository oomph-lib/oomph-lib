//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
// Driver code for a simple unstructured solid problem using a mesh
// generated from an input file generated by the triangle mesh generator
// Triangle.

//Generic routines
#include "generic.h"
#include "solid.h"
#include "constitutive.h"

// The mesh
#include "meshes/quad_from_triangle_mesh.h"

using namespace std;
using namespace oomph;

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

//================start_mesh===============================================
/// Triangle-based mesh upgraded to become a solid mesh
//=========================================================================
template<class ELEMENT>
class ElasticQuadFromTriangleMesh : public virtual QuadFromTriangleMesh<ELEMENT>, 
				    public virtual SolidMesh 
{
 
public:
 
 /// Constructor: 
 ElasticQuadFromTriangleMesh(const std::string& node_file_name,
			     const std::string& element_file_name,
			     const std::string& poly_file_name,
			     TimeStepper* time_stepper_pt=
			     &Mesh::Default_TimeStepper) : 
  QuadFromTriangleMesh<ELEMENT>(node_file_name,
				element_file_name,
				poly_file_name,
				time_stepper_pt)
  {
   // Assign the Lagrangian coordinates
   set_lagrangian_nodal_coordinates();

   // Identify special boundaries
   set_nboundary(3);

   unsigned n_node=this->nnode();
   for (unsigned j=0;j<n_node;j++)
   {
    Node* nod_pt=this->node_pt(j);

    // Boundary 1 is lower boundary
    if (abs(nod_pt->x(1))<10e-14)
    {
     this->convert_to_boundary_node(nod_pt);
     this->remove_boundary_node(0,nod_pt);
     this->add_boundary_node(1,nod_pt);
    }

    // Boundary 2 is upper boundary
    if (abs(nod_pt->x(1)-3.0)<10e-14)
    {
     this->convert_to_boundary_node(nod_pt);
     this->remove_boundary_node(0,nod_pt);
     this->add_boundary_node(2,nod_pt);
    }
   }

   // Re-setup boundary info, i.e. elements next to boundaries
   QuadFromTriangleMesh<ELEMENT>::setup_boundary_element_info();
  }

 // Calculate the L2 norm of the displacement u=R-r to overload the
 // compute_norm function in the GeneralisedElement base class
 void compute_norm(double& norm)
  {
   // Initialse the norm
   norm=0.0;

   // Per-element norm
   double el_norm=0;

   // Loop over the elements
   unsigned long n_element=this->nelement();
   for (unsigned long e=0;e<n_element;e++)
   {
    // Get a pointer to the e'th element 
    ELEMENT* el_pt=dynamic_cast<ELEMENT*>(this->element_pt(e));

#ifdef OOMPH_HAS_MPI
    // Compute error for each non-halo element
    if (!(el_pt->is_halo()))
#endif
    {
     // Compute the norm of the solution over this element
     el_compute_norm(el_pt,el_norm);
    }
    
    norm+=el_norm;
   } // (unsigned long e=0;e<n_element;e++)
  } // End of compute_norm(...)
 
 // Calculate the L2 norm of the displacement u=R-r to overload the
 // compute_norm function in the GeneralisedElement base class
 void el_compute_norm(ELEMENT* el_pt,double& el_norm)
  {
   // Initialise el_norm to 0.0
   el_norm=0.0;

   // Get the dimension of the element
   unsigned ndim=el_pt->dim();
    
   // Vector of local coordinates
   Vector<double> s(ndim);
    
   // Displacement vector, Lagrangian coordinate vector and Eulerian
   // coordinate vector respectively
   Vector<double> disp(ndim,0.0);
   Vector<double> xi(ndim,0.0);
   Vector<double> x(ndim,0.0);
    
   // Find out how many nodes there are in the element
   unsigned n_node=el_pt->nnode();

   // Construct a shape function
   Shape psi(n_node);

   // Get the number of integration points
   unsigned n_intpt=el_pt->integral_pt()->nweight();
    
   // Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {     
    // Assign values of s
    for(unsigned i=0;i<ndim;i++)
    {
     s[i]=el_pt->integral_pt()->knot(ipt,i);
    }

    // Get the integral weight
    double w=el_pt->integral_pt()->weight(ipt);

    // Get jacobian of mapping
    double J=el_pt->J_eulerian(s);
   
    // Premultiply the weights and the Jacobian
    double W=w*J;
   
    // Get the Lagrangian and Eulerian coordinate at this point, respectively
    el_pt->interpolated_xi(s,xi);
    el_pt->interpolated_x(s,x);

    // Calculate the displacement, u=R-r=x-xi
    for (unsigned idim=0;idim<ndim;idim++)
    {
     disp[idim]=x[idim]-xi[idim];
    }

    // Add to norm
    for (unsigned ii=0;ii<ndim;ii++)
    {
     el_norm+=(disp[ii]*disp[ii])*W;
    }
   }    
  } // End of compute_norm(...)
 

 /// Empty Destructor
 virtual ~ElasticQuadFromTriangleMesh() { }

};

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Physical_Variables
{
 /// Poisson's ratio
 double Nu=0.3;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt=0;

 /// Non-dim gravity
 double Gravity=0.0;

 /// Non-dimensional gravity as body force
 void gravity(const double& time, 
              const Vector<double> &xi, 
              Vector<double> &b)
 {
  b[0]=0.0;
  b[1]=-Gravity;
 }
 
 /// Uniform pressure
 double P = 0.0;

 /// Constant pressure load. The arguments to this function are imposed
 /// on us by the SolidTractionElements which allow the traction to 
 /// depend on the Lagrangian and Eulerian coordinates x and xi, and on the 
 /// outer unit normal to the surface. Here we only need the outer unit
 /// normal.
 void constant_pressure(const Vector<double> &xi, const Vector<double> &x,
                        const Vector<double> &n, Vector<double> &traction)
 {
  unsigned dim = traction.size();
  for(unsigned i=0;i<dim;i++)
  {
   traction[i] = -P*n[i];
  }
 } 

} //end namespace



//==============start_problem=========================================
/// Unstructured solid problem
//====================================================================
template<class ELEMENT> 
class UnstructuredSolidProblem : public Problem
{

public:

 /// Constructor: 
 UnstructuredSolidProblem();

 /// Destructor (empty)
 ~UnstructuredSolidProblem(){}
 
 /// Doc the solution
 void doc_solution();
 
private:
 
 /// Bulk mesh
 ElasticQuadFromTriangleMesh<ELEMENT>* Solid_mesh_pt;
 
 /// Pointer to mesh of traction elements
 SolidMesh* Traction_mesh_pt;

 /// Trace file
 ofstream Trace_file;
 
 /// DocInfo object for output
 DocInfo Doc_info;
};



//===============start_constructor========================================
/// Constructor for unstructured solid problem
//========================================================================
template<class ELEMENT>
UnstructuredSolidProblem<ELEMENT>::UnstructuredSolidProblem()
{  
 // Create solid mesh
 string node_file_name="solid.1.node";
 string element_file_name="solid.1.ele";
 string poly_file_name="solid.1.poly"; 
 Solid_mesh_pt=new ElasticQuadFromTriangleMesh<ELEMENT>(node_file_name,
							element_file_name,
							poly_file_name);
 
 // Traction elements are located on boundary 2:
 unsigned b=2;
  
 // Make traction mesh
 Traction_mesh_pt=new SolidMesh;
 
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element=Solid_mesh_pt->nboundary_element(b);
 
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
 {
  // Get pointer to the bulk element that is adjacent to boundary b
  ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
   Solid_mesh_pt->boundary_element_pt(b,e));
   
  // Find the index of the face of element e along boundary b
  int face_index = Solid_mesh_pt->face_index_at_boundary(b,e);
   
  // Create solid traction element
  SolidTractionElement<ELEMENT> *el_pt = 
   new SolidTractionElement<ELEMENT>(bulk_elem_pt,face_index);   
   
  // Add to mesh
  Traction_mesh_pt->add_element_pt(el_pt);
 
  // Set the traction function
  el_pt->traction_fct_pt() = Global_Physical_Variables::constant_pressure;
 }  
 
 // Add sub meshes
 add_sub_mesh(Solid_mesh_pt);
 add_sub_mesh(Traction_mesh_pt);
 
 // Build global mesh
 build_global_mesh();
 
 // Doc pinned solid nodes
 std::ofstream bc_file("pinned_nodes.dat");

 // Pin both positions at lower boundary (boundary 1)
 unsigned ibound=1;
 unsigned num_nod= mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
 {  
  // Get node
  SolidNode* nod_pt=Solid_mesh_pt->boundary_node_pt(ibound,inod);
   
  // Pin both directions
  for (unsigned i=0;i<2;i++)
  {
   nod_pt->pin_position(i);
     
   // ...and doc it as pinned
   bc_file << nod_pt->x(i) << " ";
  }

  // Make a new line
  bc_file << std::endl;
 }
 // Close file
 bc_file.close();

 // Complete the build of all elements so they are fully functional
 n_element = Solid_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
 {
  // Cast to a solid element
  ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Solid_mesh_pt->element_pt(i));
   
  // Set the constitutive law
  el_pt->constitutive_law_pt() =
   Global_Physical_Variables::Constitutive_law_pt;
   
  // Set the body force
  el_pt->body_force_fct_pt() = Global_Physical_Variables::gravity;
 }
   
 // Setup equation numbering scheme
 cout << "Number of equations: " << assign_eqn_numbers() << std::endl; 

 // Set directory name
 Doc_info.set_directory("RESLT");
 
 // Open trace file
 char filename[100];   
 snprintf(filename, sizeof(filename), "%s/trace.dat",Doc_info.directory().c_str());
 Trace_file.open(filename);
 
} //end constructor


//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredSolidProblem<ELEMENT>::doc_solution()
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;

 // Output solution
 //----------------
 snprintf(filename, sizeof(filename), "%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Solid_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output traction
 //----------------
 snprintf(filename, sizeof(filename), "%s/traction%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Traction_mesh_pt->output(some_file,npts);
 some_file.close();

 // Output boundaries
 //------------------
 snprintf(filename, sizeof(filename), "%s/boundaries.dat",Doc_info.directory().c_str());
 some_file.open(filename);
 Solid_mesh_pt->output_boundaries(some_file);
 some_file.close();

 // Write trace file: Norm values
 //------------------------------
 double norm;
 Solid_mesh_pt->compute_norm(norm);
 Trace_file << Global_Physical_Variables::P << " " << norm << std::endl;
 
 // Increment label for output files
 Doc_info.number()++;
}


//===========start_main===================================================
/// Demonstrate how to solve an unstructured solid problem
//========================================================================
int main(int argc, char **argv)
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define case to be run
 unsigned i_case=0;
 CommandLineArgs::specify_command_line_flag("--case",&i_case);
 
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();
 
 // Create generalised Hookean constitutive equations
 Global_Physical_Variables::Constitutive_law_pt= 
  new GeneralisedHookean(&Global_Physical_Variables::Nu);

 if (i_case==0)
 {
  std::cout << "Running with pure displacement formulation\n";

  // Typedef the element
  typedef QPVDElement<2,3> ELEMENT;
   
  // Set up the problem
  UnstructuredSolidProblem<ELEMENT> problem;
  
  // Output initial configuration
  problem.doc_solution();
  
  // Parameter study
  Global_Physical_Variables::Gravity=2.0e-4;
  Global_Physical_Variables::P=0.0;
  double pressure_increment=-1.0e-4;
  
  unsigned nstep=2; // 10;
  for (unsigned istep=0;istep<nstep;istep++)
  {
   // Solve the problem
   problem.newton_solve();
    
   // Output solution
   problem.doc_solution();
    
   // Bump up suction
   Global_Physical_Variables::P+=pressure_increment;
  }
 } // end_displacement_formulation
 else if (i_case==1)
  // Displacement/pressure formulation 
 {
  std::cout << "Running with pressure/displacement formulation\n";
  
  // Typedef the element
  typedef QPVDElementWithContinuousPressure<2> ELEMENT;
   
  // Set up the problem
  UnstructuredSolidProblem<ELEMENT> problem;
  
  // Output initial configuration
  problem.doc_solution();
  
  // Parameter study
  Global_Physical_Variables::Gravity=2.0e-4;
  Global_Physical_Variables::P=0.0;
  double pressure_increment=-1.0e-4;
  
  unsigned nstep=2;
  for (unsigned istep=0;istep<nstep;istep++)
  {
   // Solve the problem
   problem.newton_solve();
    
   // Output solution
   problem.doc_solution();
    
   // Bump up suction
   Global_Physical_Variables::P+=pressure_increment;
  }
 }
 else if (i_case==2)
  // Displacement/pressure formulation enforcing incompressibility
 {
  std::cout << "Running with pressure/displacement "
	    << "formulation (incompressible) \n";
  
  // Typedef the element
  typedef QPVDElementWithContinuousPressure<2> ELEMENT;
   
  // Set up the problem
  UnstructuredSolidProblem<ELEMENT> problem;
  
  // Output initial configuration
  problem.doc_solution();
  
  // Loop over all elements and set incompressibility flag
  {
   const unsigned n_element = problem.mesh_pt()->nelement();
   for(unsigned e=0;e<n_element;e++)
   {
    //Cast the element to the equation base of our 2D elastiticy elements
    PVDEquationsWithPressure<2> *cast_el_pt =
     dynamic_cast<PVDEquationsWithPressure<2>*>(
      problem.mesh_pt()->element_pt(e));
     
    //If the cast was successful, it's a bulk element, 
    //so set the incompressibilty flag
    if(cast_el_pt)
    {
     cast_el_pt->set_incompressible();
    }
   }
  }

  // Output initial configuration
  problem.doc_solution();
  
  // Parameter study
  Global_Physical_Variables::Gravity=2.0e-4;
  Global_Physical_Variables::P=0.0;
  double pressure_increment=-1.0e-4;
  
  unsigned nstep=2;
  for (unsigned istep=0;istep<nstep;istep++)
  {
   // Solve the problem
   problem.newton_solve();
    
   //Output solution
   problem.doc_solution();
    
   // Bump up suction
   Global_Physical_Variables::P+=pressure_increment;
  }
 }

 

} // end main



