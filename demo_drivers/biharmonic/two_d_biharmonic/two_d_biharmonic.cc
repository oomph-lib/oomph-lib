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
// include required libraries
#include "math.h"
#include "generic.h"
#include "biharmonic.h"

// use oomph namespace
using namespace oomph;


namespace RayParam
{

  // Variables to set up problem.
  bool Unclamp_right_bound = false;
  unsigned Noel = 0;
  int X_min = 0;
  int X_max = 0;
  int Y_min = 0;
  int Y_max = 0;

  bool Print_connectivity_mat = false;
  bool Print_natural_jacobian = false;
  bool Print_subblocks = false;
  bool Print_elemental_jacobian = false;

  std::string Matbase_str = "";

  std::string Connectivity_dir = "";
  std::string Natural_jac_dir = "";
  std::string Blocked_dir = "";
  std::string Elemental_dir = "";
}


//=============================================================================
// two dimensional biharmonic plate problem 1 namespace - contains all the 
// problem functions
//=============================================================================
namespace BiharmonicTestFunctions1
{

 // DIRICHLET BOUNDARY CONDITIONS

 void u_NE(const double& s, double& u)
 {
  double x = (s+1) / 2;
  u = x*x*x;
 }
 void dudn_NE(const double& s, double& dudn)
 {
  double x = (s+1) / 2;
  dudn = 3*x*x*x;
 }
 void u_SW(const double& s, double& u)
 {
  u = 0;
 }
 void dudn_SW(const double& s, double& dudn)
 {
  dudn = 0;
 }  

 // SURFACE LOAD FUNCTION

 void surface_load(const Vector<double>& x, double& f)
 {
   //f = 1.0;
   f = 72*x[0]*x[1];
 }


 // NEUMANN BOUNDARY CONDITIONS

 void flux1_NE(const double& s, double& flux1)
 {
  double x = (s+1) / 2;
  flux1 = 6*x*x*x + 18*x;
 }
 void flux0_NE(const double& s, double& flux0)
 {
  double x = (s+1) / 2;
  flux0 = 6*(x*x*x + x);
 }
 

 // SOLUTION

 void solution(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = x[0]*x[0]*x[0]*x[1]*x[1]*x[1];
 }
}



//=============================================================================
// Two Dimensional Biharmonic Test Problem (square)
// All Edges Clamped
//=============================================================================
class BiharmonicTestProblem1 : public BiharmonicProblem<2>
{

private:

 // Domain pointer
 TopologicallyRectangularDomain* Domain_pt;

public:

 // constructor
 BiharmonicTestProblem1(const unsigned n_element)
  {
   // force linear
   Problem::Problem_is_nonlinear = false;

   // use test functions 1 namespace
   using namespace BiharmonicTestFunctions1;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(1,1);

   // assemble mesh	
   this->build_bulk_mesh(n_element, n_element, Domain_pt);

   // set the bulk element source function
   set_source_function(surface_load);

   // clamp edge on all boundaries
   set_dirichlet_boundary_condition(0,u_SW,dudn_SW);
   set_dirichlet_boundary_condition(1,u_NE,dudn_NE);
   set_dirichlet_boundary_condition(2,u_NE,dudn_NE);
   set_dirichlet_boundary_condition(3,u_SW,dudn_SW);

   // assign equation numbers
   this->build_global_mesh_and_assign_eqn_numbers();
  }

 // constructor - more flexible, for testing purposes.
 BiharmonicTestProblem1(const int x_min, const int x_max,
                        const int y_min, const int y_max,
                        const unsigned n_element)
  {
   // force linear
   Problem::Problem_is_nonlinear = false;

   // use test functions 1 namespace
   using namespace BiharmonicTestFunctions1;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(x_min,x_max,y_min,y_max);

   // assemble mesh	
   this->build_bulk_mesh(n_element, n_element, Domain_pt);

   // set the bulk element source function
   set_source_function(surface_load);

   // clamp edge on all boundaries
   set_dirichlet_boundary_condition(0,u_SW,dudn_SW);
   if(!RayParam::Unclamp_right_bound)
   {
     oomph_info << "Clamping right boundary" << std::endl; 
     
     set_dirichlet_boundary_condition(1,u_NE,dudn_NE);
   }
   set_dirichlet_boundary_condition(2,u_NE,dudn_NE);
   set_dirichlet_boundary_condition(3,u_SW,dudn_SW);

   // assign equation numbers
   this->build_global_mesh_and_assign_eqn_numbers();
  }

 /// Destructor - just deletes domain pt
 virtual ~BiharmonicTestProblem1()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };

 void dump_jacobian(std::string& jac_str)
 {
   DoubleVector res;
   CRDoubleMatrix jac;

   this->get_jacobian(res,jac);

   jac.sparse_indexed_output(jac_str,15,true);
 }
};


//=============================================================================
// TWO DIMENSIONAL BIHARMONIC TEST 2 (PLATE PROBLEM)
//  u = (r-2)^3-r
//   theta in [0,1/4*pi]
//   r in [pi,3pi]
//=============================================================================



// two dimensional biharmonic plate problem 2 namespace - contains all the 
// problem functions
namespace BiharmonicTestFunctions2
{


 // PARAMETERS

 double Pi = MathematicalConstants::Pi;
 double theta = Pi/4;
 double r_min = sqrt(1*Pi);
 double r_max = sqrt(3*Pi);


 // BOUNDARIES
 
 void boundary_N(const double& s, Vector<double>& r)
 {  
  r[0] = (r_min + (0.5*(s+1)*(r_max-r_min)))*cos(theta);  
  r[1] = (r_min + (0.5*(s+1)*(r_max-r_min)))*sin(theta);
 }
 void boundary_E(const double& s, Vector<double>& r)
 {    
  r[0] = r_max * cos((s+1)*theta/2);
  r[1] = r_max * sin((s+1)*theta/2);
 } 
 void boundary_S(const double& s, Vector<double>& r)
 {
  r[0] = r_min + (0.5*(s+1)*(r_max-r_min));  
  r[1] = 0.0;
 }
 void boundary_W(const double& s, Vector<double>& r)
 {
  r[0] = r_min * cos((s+1)*theta/2);
  r[1] = r_min * sin((s+1)*theta/2);
 }


 // NORMALS

 void normal_N(const double& s, Vector<double>& n)
 {
  n[0] = -sin(theta);
  n[1] = cos(theta);
 }
 void normal_E(const double& s, Vector<double>& n)
 {
  double t = (s+1)*theta/2; 
  n[0] = cos(t);
  n[1] = sin(t);
 }
 void normal_S(const double& s, Vector<double>& n)
 {
  n[0] = 0.0;
  n[1] = -1.0;
 }
 void normal_W(const double& s, Vector<double>& n)
 {
  double t = (s+1)*theta/2; 
  n[0] = -cos(t);
  n[1] = -sin(t);
 }


 // DIRICHLET BCs

 void u_N(const double& s, double& u)
 {  
  double r = r_min+0.5*(s+1)*(r_max-r_min);
  u = sin(r*r)*tan(theta);
 }
 void u_E(const double& s, double& u)
 {
  double t = (s+1)*theta/2;
  u = sin(r_max*r_max)*tan(t);
 } 
 void u_S(const double& s, double& u)
 {
  u = 0;
 }
 void u_W(const double& s, double& u)
 {  
  double t = (s+1)*theta/2;
  u = sin(r_min*r_min)*tan(t);
 }
 double dudx_0(const Vector<double> x)
 {
  return(2*cos(x[0]*x[0]+x[1]*x[1])*x[1]
         -sin(x[0]*x[0]+x[1]*x[1])*x[1]/(x[0]*x[0]));
 }
 double dudx_1(const Vector<double> x)
 {
  return(2*cos(x[0]*x[0]+x[1]*x[1])*x[1]*x[1]/x[0]
         +sin(x[0]*x[0]+x[1]*x[1])/x[0]);
 }
 void dudn_N(const double& s, double& dudn)
 {  
  Vector<double> x(2);
  boundary_N(s,x);
  Vector<double> n(2);
  normal_N(s,n);
  dudn = dudx_0(x)*n[0] + dudx_1(x)*n[1];
 }
 void dudn_E(const double& s, double& dudn)
 {  
  Vector<double> x(2);
  boundary_E(s,x);
  Vector<double> n(2);
  normal_E(s,n);
  dudn = dudx_0(x)*n[0] + dudx_1(x)*n[1];
 }
 void dudn_S(const double& s, double& dudn)
 {
  double x = r_min+0.5*(s+1)*(r_max-r_min);
  dudn = -sin(x*x)/x;
 }
 void dudn_W(const double& s, double& dudn)
 {  
  Vector<double> x(2);
  boundary_W(s,x);
  Vector<double> n(2);
  normal_W(s,n);
  dudn = dudx_0(x)*n[0] + dudx_1(x)*n[1];
 } 

 // SURFACE LOAD FUNCTION
 
 void surface_load(const Vector<double>& x, double& f)
 {
  double sinr2 = sin(x[0]*x[0]+x[1]*x[1]);
  double cosr2 = cos(x[0]*x[0]+x[1]*x[1]);
  f = (16*sinr2*x[0]*x[0]*x[0]*x[1]
       -64*cosr2*x[0]*x[1]
       -48*sinr2*x[1]/x[0]
       +24*sinr2*x[1]/(x[0]*x[0]*x[0]*x[0]*x[0])
       +16*sinr2*x[1]*x[1]*x[1]*x[1]*x[1]/x[0]
       -64*cosr2*x[1]*x[1]*x[1]/x[0]
       +32*sinr2*x[1]*x[1]*x[1]*x[0]
       -16*sinr2*x[1]*x[1]*x[1]/(x[0]*x[0]*x[0]));
 }

 // SOLUTION

 void solution(const Vector<double>& x, Vector<double>& u)
 {
  u[0] = sin(x[0]*x[0] + x[1]*x[1])*x[1]/x[0];
 }
}



//=============================================================================
// Two Dimensional Biharmonic Test Problem
// All Edges Clamped
//=============================================================================
class BiharmonicTestProblem2 : public BiharmonicProblem<2>
{

private:

 // Domain pointer
 TopologicallyRectangularDomain* Domain_pt;

public:

 // constructor
 BiharmonicTestProblem2(const unsigned n_element)
  {
   // force linear
   Problem::Problem_is_nonlinear = false;

   // use test functions 1 namespace
   using namespace BiharmonicTestFunctions2;

   // create the domain describing the geometry of the problem
   Domain_pt = new TopologicallyRectangularDomain(boundary_N,boundary_E,
                                                  boundary_S,boundary_W);

   // assemble mesh	
   this->build_bulk_mesh(n_element, n_element, Domain_pt);

   // set the bulk element source function
   set_source_function(surface_load);

   // clamp edge on all boundaries
   set_dirichlet_boundary_condition(0,u_S,dudn_S);
   set_dirichlet_boundary_condition(1,u_E,dudn_E);
   set_dirichlet_boundary_condition(2,u_N,dudn_N);
   set_dirichlet_boundary_condition(3,u_W,dudn_W);

   // assign equation numbers
   this->build_global_mesh_and_assign_eqn_numbers();
  }

 /// Destructor - just deletes domain pt
 virtual ~BiharmonicTestProblem2()
  {
   // delete the domain
   delete Domain_pt;
   Domain_pt = 0;
  };
};

void print_elemental_jacobian(const unsigned& element_number, 
    const Problem* const problem_pt)
{
  AssemblyHandler* const assembly_handler_pt 
    = problem_pt->assembly_handler_pt();

  const Mesh* const mesh_pt = problem_pt->mesh_pt();

  const unsigned n_element = mesh_pt->nelement();
  const unsigned n_ele_1d = sqrt(n_element);

  oomph_info << "Elements: 1D: " << n_ele_1d 
             << ", total: " << n_element << std::endl;

#ifdef PARANOID
  if(element_number >= n_element)
  {
    std::ostringstream err_stream;
    err_stream << "Supplied element number: " << element_number << ",\n"
      << "But number of elements is: " << n_element << std::endl;
    throw OomphLibError(err_stream.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
  }
#endif

  // Get pointer to the element
  GeneralisedElement* elem_pt 
    = mesh_pt->element_pt(element_number);

  // Find number of dofs in the element
  const unsigned n_element_dofs = assembly_handler_pt->ndof(elem_pt);

  // Set up an array
  Vector<double> element_residuals(n_element_dofs);

  // Set up a matrix
  DenseMatrix<double> element_jacobian(n_element_dofs);

  // Fill the array
  assembly_handler_pt->get_jacobian(elem_pt,
      element_residuals,
      element_jacobian);

  std::ostringstream ele_stream;
  ele_stream << "N"<< n_ele_1d << "_ele_" << element_number;
  element_jacobian.sparse_indexed_output(ele_stream.str(),15,true);

  //  // Output the coordinates just to make sure...
  //  FiniteElement* finite_element_pt 
  //    = mesh_pt->finite_element_pt(element_number);
  //
  //  const unsigned n_node = finite_element_pt->nnode();
  //
  //  // Loop over nodes and output coordinates
  //  for (unsigned n = 0; n < n_node; n++) 
  //  {
  //    Node* nod_pt = finite_element_pt->node_pt(n);
  //    const double x = nod_pt->x(0);
  //    const double y = nod_pt->x(1);
  //    oomph_info << "Node: " << n << ", " << x << ", " << y << std::endl;
  //  }
}


void print_connectivity_matrix(const Problem* const problem_pt)
{
  const Mesh* const mesh_pt = problem_pt->mesh_pt();

  const unsigned n_element = mesh_pt->nelement();
//  const unsigned n_ele_1d = sqrt(n_element);

  // Set up out file
  std::ostringstream filenamestream;
  filenamestream << RayParam::Connectivity_dir
    << "/"
    << RayParam::Matbase_str;


  std::ofstream outfile;
  outfile.open(filenamestream.str().c_str());

  // Loop through the number of elements
  for (unsigned ele_i = 0; ele_i < n_element; ele_i++) 
  {
    outfile << "Element number: " << ele_i << std::endl;

    // Get pointer to the element
    FiniteElement* elem_pt 
      = mesh_pt->finite_element_pt(ele_i);

    unsigned nnod=elem_pt->nnode();

    for (unsigned nod_i = 0; nod_i < nnod; nod_i++) 
    {
      outfile << "Node number: " << nod_i << std::endl;

      // Get the node
      Node* nod_pt = elem_pt->node_pt(nod_i);
      const unsigned nval = nod_pt->nvalue();

      outfile << "x: " << nod_pt->x(0) 
              << ", y: " << nod_pt->x(1) 
              << std::endl;

      for (unsigned val_i = 0; val_i < nval; val_i++) 
      {
        outfile << val_i << " ";
        long eqn_num = nod_pt->eqn_number(val_i);
        outfile << eqn_num << " ";
        if(!nod_pt->is_pinned(val_i))
        {
          outfile << elem_pt->local_eqn_number(eqn_num) << " ";
        }

        outfile << "\n";
      } // for values in node
    } // for node
  } // for elements
  outfile.close();
}


//=============================================================================
/// main
//=============================================================================
int main(int argc, char *argv[])
{
  // number of element
  unsigned n_element = 20;
  
  // Set up doc info
  DocInfo doc_info;
  doc_info.set_directory("RESLT");
  doc_info.number()=0;


  // Store commandline arguments
  CommandLineArgs::setup(argc,argv);

  CommandLineArgs::specify_command_line_flag("--unclamp_right_bound");

  // Need to make a string to reflect the type of problem
  CommandLineArgs::specify_command_line_flag("--noel", 
      &RayParam::Noel);
  CommandLineArgs::specify_command_line_flag("--xmin", 
      &RayParam::X_min);
  CommandLineArgs::specify_command_line_flag("--xmax", 
      &RayParam::X_max);
  CommandLineArgs::specify_command_line_flag("--ymin", 
      &RayParam::Y_min);
  CommandLineArgs::specify_command_line_flag("--ymax", 
      &RayParam::Y_max);

  CommandLineArgs::specify_command_line_flag("--connectivity_mat",
      &RayParam::Connectivity_dir);
  CommandLineArgs::specify_command_line_flag("--natural_jacobian",
      &RayParam::Natural_jac_dir);
  CommandLineArgs::specify_command_line_flag("--sub_blocks",
      &RayParam::Blocked_dir);
  CommandLineArgs::specify_command_line_flag("--elemental_jacobian",
      &RayParam::Elemental_dir);

  // Parse the above flags.
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();



  // If there is more than one argument, we Milan's tests.
  if(argc > 1)
  {

  if(CommandLineArgs::command_line_flag_has_been_set(
        "--unclamp_right_bound"))
  {
    RayParam::Unclamp_right_bound = true;
    RayParam::Matbase_str = "twod_bihar_right_bound_neumann";
  }
  else
  {
    RayParam::Unclamp_right_bound = false;
    RayParam::Matbase_str = "twod_bihar_clamp_all_bound";
  }

  if(!CommandLineArgs::command_line_flag_has_been_set("--noel"))
  {
    std::ostringstream err_stream;
    err_stream << "Please set --noel." << std::endl;
    throw OomphLibError(err_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }
  else
  {
    std::ostringstream str_ss;
    str_ss << RayParam::Matbase_str << "_Noel" << RayParam::Noel;
    RayParam::Matbase_str = str_ss.str();
  }

  if(!CommandLineArgs::command_line_flag_has_been_set("--xmin"))
  {
    std::ostringstream err_stream;
    err_stream << "Please set --xmin." << std::endl;
    throw OomphLibError(err_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }
  if(!CommandLineArgs::command_line_flag_has_been_set("--xmax"))
  {
    std::ostringstream err_stream;
    err_stream << "Please set --xmax." << std::endl;
    throw OomphLibError(err_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }
  if(!CommandLineArgs::command_line_flag_has_been_set("--ymin"))
  {
    std::ostringstream err_stream;
    err_stream << "Please set --ymin." << std::endl;
    throw OomphLibError(err_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }
  if(!CommandLineArgs::command_line_flag_has_been_set("--ymax"))
  {
    std::ostringstream err_stream;
    err_stream << "Please set --ymax." << std::endl;
    throw OomphLibError(err_stream.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }

  if(CommandLineArgs::command_line_flag_has_been_set("--connectivity_mat"))
  {
    RayParam::Print_connectivity_mat = true;
  }
  else
  {
    RayParam::Print_connectivity_mat = false;
  }

  if(CommandLineArgs::command_line_flag_has_been_set("--natural_jacobian"))
  {
    RayParam::Print_natural_jacobian = true;
  }
  else
  {
    RayParam::Print_natural_jacobian = false;
  }

  if(CommandLineArgs::command_line_flag_has_been_set("--sub_blocks"))
  {
    RayParam::Print_subblocks = true;
  }
  else
  {
    RayParam::Print_subblocks = false;
  }

  if(CommandLineArgs::command_line_flag_has_been_set("--elemental_jacobian"))
  {
    RayParam::Print_elemental_jacobian = true;
  }
  else
  {
    RayParam::Print_elemental_jacobian = false;
  }



    oomph_info << "n_element: " << RayParam::Noel << std::endl;
    oomph_info << "x_min: " << RayParam::X_min << std::endl;
    oomph_info << "x_max: " << RayParam::X_max << std::endl;
    oomph_info << "y_min: " << RayParam::Y_min << std::endl;
    oomph_info << "y_max: " << RayParam::Y_max << std::endl;

    BiharmonicTestProblem1 problem(RayParam::X_min,
                                   RayParam::X_max,
                                   RayParam::Y_min,
                                   RayParam::Y_max,
                                   RayParam::Noel);

    BiharmonicPreconditioner my_prec;
    my_prec.bulk_element_mesh_pt() = problem.bulk_element_mesh_pt();
    my_prec.preconditioner_type() = 0;

    IterativeLinearSolver* solver_pt = new  CG<CRDoubleMatrix>;
    solver_pt->preconditioner_pt() = &my_prec;

    // Apply the solver
    problem.linear_solver_pt() = solver_pt;

    if(RayParam::Print_natural_jacobian)
    {
      std::ostringstream dump_stream;
      dump_stream << RayParam::Natural_jac_dir 
                  << "/"
                  << RayParam::Matbase_str;
      std::string dump_str = dump_stream.str();
      problem.dump_jacobian(dump_str);
    }

    if(RayParam::Print_connectivity_mat)
    {
      print_connectivity_matrix(&problem);
    }

//    if(RayParam::Print_elemental_jacobian)
//    {
//      //print_elemental_jacobian(&problem);
//    }

    if(RayParam::Print_subblocks)
    {
      std::ostringstream dump_stream;
      dump_stream << RayParam::Blocked_dir 
                  << "/"
                  << RayParam::Matbase_str;
      std::string dump_str = dump_stream.str();

//      my_prec.print_subblocks(dump_str);

      std::ostringstream dump_ss;
      dump_ss << RayParam::Natural_jac_dir
              << "/"
              << RayParam::Matbase_str;
      std::string another_dump_str = dump_ss.str();

//      my_prec.set_fullblock_dir(another_dump_str);

      problem.newton_solve();

    }
  }
  else
  {
    // Biharmonic Problem 1 (square)
    // Exact Biharmonic Preconditioner
    {
      oomph_info 
        << "/////////////////////////////////////////////////////////////////////"
        << std::endl;
      oomph_info << "TESTING: Square 2D Biharmonic Problem w/ "
        << "Exact Preconditioning"
        << std::endl;
      oomph_info 
        << "/////////////////////////////////////////////////////////////////////"
        << std::endl;

      // create the problem
      BiharmonicTestProblem1 problem(n_element);

      // setup the preconditioner
      BiharmonicPreconditioner* prec_pt = new BiharmonicPreconditioner;
      prec_pt->bulk_element_mesh_pt() = problem.bulk_element_mesh_pt();
      prec_pt->preconditioner_type() = 0;

      // setup the solver
      IterativeLinearSolver* solver_pt = new CG<CRDoubleMatrix>;  
      solver_pt->preconditioner_pt() = prec_pt;

      // apply the solver
      problem.linear_solver_pt() = solver_pt;
      problem.newton_solve();

      // ouput the solution
      problem.doc_solution(doc_info,BiharmonicTestFunctions1::solution);
      doc_info.number()++;

      // clean up
      delete solver_pt;
      delete prec_pt;
    }

    // Biharmonic Problem 2 (section of annulus)
    // Inexact Biharmonic Preconditioner w/ SuperLU
    {
      oomph_info 
        << "/////////////////////////////////////////////////////////////////////"
        << std::endl;
      oomph_info << "TESTING: Curved 2D Biharmonic Problem w/ "
        << "Inexact Preconditioning"
        << std::endl;
      oomph_info 
        << "/////////////////////////////////////////////////////////////////////"
        << std::endl;

      // create the problem
      BiharmonicTestProblem2 problem(n_element);

      // setup the preconditioner
      BiharmonicPreconditioner* prec_pt = new BiharmonicPreconditioner;
      prec_pt->bulk_element_mesh_pt() = problem.bulk_element_mesh_pt();
      prec_pt->preconditioner_type() = 1;

      // setup the solver
      IterativeLinearSolver* solver_pt = new CG<CRDoubleMatrix>;  
      solver_pt->preconditioner_pt() = prec_pt;

      // apply the solver
      problem.linear_solver_pt() = solver_pt;
      problem.newton_solve();

      // ouput the solution
      problem.doc_solution(doc_info,BiharmonicTestFunctions2::solution);
      doc_info.number()++;

      // clean up
      delete solver_pt;
      delete prec_pt;
    }
  }
}
