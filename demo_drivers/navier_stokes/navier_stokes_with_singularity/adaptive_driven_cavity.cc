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
// Driver for adaptive 2D rectangular driven cavity. Solved with black
// box adaptation, using Taylor Hood and Crouzeix Raviart elements.

// Generic oomph-lib header
#include "generic.h"

// Navier Stokes headers
#include "navier_stokes.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//==start_of_namespace===================================================
/// Namespace for physical parameters
//=======================================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re=100.0;

 // Max permitted error during adaptation (if negative use default)
 double Max_permitted_error=-1.0;

 // Max permitted error during adaptation (if negative use default)
 double Min_permitted_error=-1.0;

} // end_of_namespace


#include "finite_re_perturbation.h"


//==start_of_problem_class============================================
/// Driven cavity problem in rectangular domain, templated
/// by element type. 
//====================================================================
template<class ELEMENT>
class RefineableDrivenCavityProblem : public Problem
{

public:

 /// Constructor
 RefineableDrivenCavityProblem();

 /// Destructor: Empty
 ~RefineableDrivenCavityProblem() {}

 /// Update the after solve (empty)
 void actions_after_newton_solve() {}

 /// Update the problem specs before solve. 
 /// (Re-)set velocity boundary conditions just to be on the safe side...
 void actions_before_newton_solve()
  { 
  // Setup tangential flow along boundary 0:
  unsigned ibound=0; 
  unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    // Tangential flow
    unsigned i=0;
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,1.0);
    // No penetration
    i=1;
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
   }
  
  // Overwrite with no flow along all other boundaries
  unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned ibound=1;ibound<num_bound;ibound++)
   {
    unsigned num_nod= mesh_pt()->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      for (unsigned i=0;i<2;i++)
       {
        mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
       }
     }
   }
  } // end_of_actions_before_newton_solve


 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableNavierStokesEquations<2>::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());
    
    // Pin redundant pressure dofs
   RefineableNavierStokesEquations<2>::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   // Now set the first pressure dof in the last element to 0.0
   unsigned nel=mesh_pt()->nelement();
   fix_pressure(nel-1,0,0.0);
   
   // Now we have the mesh, let's set up the plot points
   setup_line_plot_points();

  } // end_of_actions_after_adapt
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
 
 /// Build line plot points
 void setup_line_plot_points()
  {
   // Get geom object
   delete Mesh_as_geom_object_pt;
   Mesh_as_geom_object_pt=new MeshAsGeomObject(mesh_pt());

   // Number of radial lines
   unsigned n_phi=5;
   Line_plot_point.resize(n_phi);

   // How many points to you want?
   unsigned npt=100;
   double r_max=0.1;
   Vector<double> zeta(2);   
   Vector<double> s(2);
   GeomObject* geom_object_pt=0;
   for (unsigned i=0;i<n_phi;i++)
    {
     Line_plot_point[i].resize(npt);
     double phi=0.5*MathematicalConstants::Pi*double(i)/double(n_phi-1);
     oomph_info << "setup at phi : " << phi/(0.5*MathematicalConstants::Pi) 
                << " pi/2" << std::endl;
     for (unsigned j=0;j<npt;j++)
      {
       double r=r_max*double(j+1)/double(npt);
       zeta[0]=r*cos(phi);
       zeta[1]=r*sin(phi);
       Mesh_as_geom_object_pt->locate_zeta(zeta,geom_object_pt,s);
        if (geom_object_pt==0)
         {
          oomph_info << "Point : " 
                     << zeta[0] << " " 
                     << zeta[1] << " " 
                     << " not found in setup of line plots" 
                     << std::endl;
         }   
        Line_plot_point[i][j]=std::make_pair(geom_object_pt,s);
      }
    }
  }

private:

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to proper element and fix pressure
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure

 /// Mesh as geom object representation of mesh
 MeshAsGeomObject* Mesh_as_geom_object_pt;

 /// Line_plot_point[i_phi][i_rho]
 Vector<Vector<std::pair<GeomObject*,Vector<double> > > > Line_plot_point;

}; // end_of_problem_class



//==start_of_constructor==================================================
/// Constructor for RefineableDrivenCavity problem 
///
//========================================================================
template<class ELEMENT>
RefineableDrivenCavityProblem<ELEMENT>::RefineableDrivenCavityProblem()
{ 

 // Null out 
 Mesh_as_geom_object_pt=0;
 
 // Setup mesh

 // # of elements in x-direction
 unsigned n_x=10;

 // # of elements in y-direction
 unsigned n_y=10;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;

 // Build and assign mesh
 Problem::mesh_pt() = 
  new RefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y);

 // Now we have the mesh, let's set up the plot points
 setup_line_plot_points();

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(mesh_pt())->
  spatial_error_estimator_pt()=error_estimator_pt;

 // Allow for lots of refinement
 dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(mesh_pt())->
  max_refinement_level()=50;

 // Overwrite default for max. permitted error?
 if (Global_Physical_Variables::Max_permitted_error>0.0)
  {
   dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(mesh_pt())->
    max_permitted_error()=Global_Physical_Variables::Max_permitted_error;
  }

 // Overwrite default for min. permitted error?
 if (Global_Physical_Variables::Min_permitted_error>0.0)
  {
   dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(mesh_pt())->
    min_permitted_error()=Global_Physical_Variables::Min_permitted_error;
  }

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here: All boundaries are Dirichlet boundaries.
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over values (u and v velocities)
     for (unsigned i=0;i<2;i++)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
      }
    }
  } // end loop over boundaries

 //Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to Reynolds
 // number
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
  } // end loop over elements
 
 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 
  // Now set the first pressure dof in the first element to 0.0
 fix_pressure(0,0,0.0);
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} // end_of_constructor



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RefineableDrivenCavityProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file, some_file2;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 


 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();

 // Output solution 
 sprintf(filename,"%s/coarse_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 npts=2;
 mesh_pt()->output(some_file,npts);
 some_file.close();
 
 // Output perturbation solution 
 //-----------------------------
 sprintf(filename,"%s/perturbation_soln_two_term%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 PerturbationSolution::N_terms_for_plot=2;
 some_file.open(filename);
 mesh_pt()->output_fct(
  some_file,npts,
  PerturbationSolution::perturbation_soln_for_plot);
 some_file.close();

 // Output perturbation solution 
 //-----------------------------
 sprintf(filename,"%s/perturbation_soln_one_term%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 PerturbationSolution::N_terms_for_plot=1;
 mesh_pt()->output_fct(
  some_file,npts,
  PerturbationSolution::perturbation_soln_for_plot);
 some_file.close();


 // Output first-order perturbation solution (Stokes)
 //--------------------------------------------------
 sprintf(filename,"%s/first_order_perturbation%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_fct(
  some_file,npts,
  PerturbationSolution::first_order_perturbation_for_plot);
 some_file.close();



 // Output second-order perturbation solution (Stokes)
 //--------------------------------------------------
 sprintf(filename,"%s/second_order_perturbation%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_fct(
  some_file,npts,
  PerturbationSolution::second_order_perturbation_for_plot);
 some_file.close();



 // Do line plots
 sprintf(filename,"%s/line_plot%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 sprintf(filename,"%s/second_order_perturbation_line_plot%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file2.open(filename);
 Vector<double> s(2); 
 Vector<double> x(2); 
 Vector<double> u(2);
 double p=0.0;
 Vector<double> soln(3);
 unsigned nphi=Line_plot_point.size();
 for (unsigned i=0;i<nphi;i++)
  {
   unsigned nr=Line_plot_point[i].size();
   some_file << "ZONE T=\"Re=" << Global_Physical_Variables::Re 
             << "; <GREEK>j</GREEK> = " << i << "/" << nphi-1 
             << " <GREEK>p</GREEK>/2 \"" 
             << std::endl;  
   some_file2 << "ZONE T=\"Re=" << Global_Physical_Variables::Re 
             << "; <GREEK>j</GREEK> = " << i << "/" << nphi-1 
             << " <GREEK>p</GREEK>/2 \"" 
             << std::endl;  
   for (unsigned j=0;j<nr;j++)
    {
     ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Line_plot_point[i][j].first);
     s=Line_plot_point[i][j].second;
     el_pt->interpolated_x(s,x);
     el_pt->interpolated_u_nst(s,u);
     p=el_pt->interpolated_p_nst(s);
     some_file << x[0] << " " 
               << x[1] << " " 
               << sqrt(x[0]*x[0]+x[1]*x[1]) << " " 
               << u[0] << " " 
               << u[1] << " " 
               << p << " " 
               << std::endl;

     PerturbationSolution::second_order_perturbation_for_plot(x,soln);
     some_file2 << x[0] << " " 
                << x[1] << " " 
                << sqrt(x[0]*x[0]+x[1]*x[1]) << " " 
                << soln[0] << " " 
                << soln[1] << " " 
                << soln[2] << " " 
                << std::endl;
    }
  }
 some_file.close();
 some_file2.close();


} // end_of_doc_solution




//==start_of_main======================================================
/// Driver for RefineableDrivenCavity test problem 
//=====================================================================
int main(int argc, char* argv[])
{

 // Allow indefinite refinement without complaints about near singular
 // Jacobian
 FiniteElement::Tolerance_for_singular_jacobian=0.0;

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
 // Reynolds number
 CommandLineArgs::specify_command_line_flag("--re",
                                            &Global_Physical_Variables::Re);

 // Max error
 CommandLineArgs::specify_command_line_flag(
  "--max_permitted_error",
  &Global_Physical_Variables::Max_permitted_error);

 // Min error
 CommandLineArgs::specify_command_line_flag(
  "--min_permitted_error",
  &Global_Physical_Variables::Min_permitted_error);

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Set output directory
 DocInfo doc_info;
 doc_info.set_directory("RESLT");

 // Set max. number of black-box adaptation
 unsigned max_adapt=1;
 
 //Build problem
 RefineableDrivenCavityProblem<RefineableQTaylorHoodElement<2> > problem;
 
 unsigned n_adapt=10;
 for (unsigned i=0;i<n_adapt;i++)
  {
   // Solve the problem with automatic adaptation
   problem.newton_solve(max_adapt);
   
   //Output solution
   problem.doc_solution(doc_info);
   
   // Step number
   doc_info.number()++;
  }
 

} // end_of_main











