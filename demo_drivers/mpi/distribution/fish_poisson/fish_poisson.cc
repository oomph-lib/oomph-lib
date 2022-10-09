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
// Driver for solution of 2D Poisson equation in fish-shaped domain with
// adaptivity

 
// Generic oomph-lib headers
#include "generic.h"

// The Poisson equations
#include "poisson.h"

// The fish mesh 
#include "meshes/fish_mesh.h"

using namespace std;

using namespace oomph;


//============ start_of_namespace=====================================
/// Namespace for const source term in Poisson equation
//====================================================================
namespace ConstSourceForPoisson
{ 
 
 /// Strength of source function: default value -1.0
 double Strength=-1.00;

/// Const source function
 void source_function(const Vector<double>& x, double& source)
 {
  source = Strength;
 }

} // end of namespace




//======start_of_problem_class========================================
/// Refineable Poisson problem in fish-shaped domain.
/// Template parameter identifies the element type.
//====================================================================
template<class ELEMENT>
class RefineableFishPoissonProblem : public Problem
{

public:

 /// Constructor
 RefineableFishPoissonProblem();

 /// Destructor: Empty; all memory gets cleaned up in base destructor
 virtual ~RefineableFishPoissonProblem(){};

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}

 /// Update the problem specs befor solve (also empty)
 void actions_before_newton_solve() {}

 /// Actions before adapt (essentially empty)
 void actions_before_adapt() {oomph_info << "Min/max about to adapt\n";}

 /// Overloaded version of the problem's access function to 
 /// the mesh. Recasts the pointer to the base Mesh object to 
 /// the actual mesh type.
 RefineableFishMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableFishMesh<ELEMENT>*>(Problem::mesh_pt());
  }

 /// Doc the solution. Output directory and labels are specified 
 /// by DocInfo object
 void doc_solution(DocInfo& doc_info);

}; // end of problem class





//===========start_of_constructor=========================================
/// Constructor for adaptive Poisson problem in fish-shaped
/// domain.
//========================================================================
template<class ELEMENT>
RefineableFishPoissonProblem<ELEMENT>::RefineableFishPoissonProblem()
{ 
    
 // Build fish mesh -- this is a coarse base mesh consisting 
 // of four elements. We'll refine/adapt the mesh later.
 Problem::mesh_pt()=new RefineableFishMesh<ELEMENT>;
 
 // Create/set error estimator
 mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;
  
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. Since the boundary values are never changed, we set
 // them here rather than in actions_before_solve(). 
 unsigned n_bound = mesh_pt()->nboundary();
 for(unsigned i=0;i<n_bound;i++)
  {
   unsigned n_node = mesh_pt()->nboundary_node(i);
   for (unsigned n=0;n<n_node;n++)
    {
     // Pin the single scalar value at this node
     mesh_pt()->boundary_node_pt(i,n)->pin(0); 

     // Assign the homogenous boundary condition for the one and only
     // nodal value
     mesh_pt()->boundary_node_pt(i,n)->set_value(0,0.0); 
    }
  }

 // Loop over elements and set pointers to source function
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from FiniteElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //Set the source function pointer
   el_pt->source_fct_pt() = &ConstSourceForPoisson::source_function;
  }

 // Setup the equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor




//=======start_of_doc=====================================================
/// Doc the solution in tecplot format.
//========================================================================
template<class ELEMENT>
void RefineableFishPoissonProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 unsigned min_refinement_level=0;
 unsigned max_refinement_level=0;
 mesh_pt()->get_refinement_levels(min_refinement_level,
                                  max_refinement_level);

 oomph_info << "Min/max actual and max/min allowed refinement level" 
            << min_refinement_level << " " 
            << max_refinement_level << " "
            << mesh_pt()->min_refinement_level() << " "
            << mesh_pt()->max_refinement_level() << "\n";


 ofstream some_file;
 char filename[100];

 // Number of plot points in each coordinate direction.
 unsigned npts;
 npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i_on_proc%i.dat",doc_info.directory().c_str(),
         doc_info.number(),this->communicator_pt()->my_rank());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file << "TEXT X=25,Y=93,F=HELV,HU=POINT,C=BLUE,H=10,T=\"" 
           << doc_info.label() << "; Min/max: "
           << min_refinement_level << " " 
           << max_refinement_level << " "
           << "\"" << std::endl; 
 some_file.close();
 
} // end of doc

 




//=====================start_of_incremental===============================
/// Demonstrate how to solve 2D Poisson problem in 
/// fish-shaped domain with mesh adaptation. First we solve on the original
/// coarse mesh. Next we do a few uniform refinement steps and re-solve.
/// Finally, we enter into an automatic adapation loop.
//========================================================================
void solve_with_incremental_adaptation()
{
 
 oomph_info << "Min/max about to start incremental\n";

 //Set up the problem with nine-node refineable Poisson elements
 RefineableFishPoissonProblem<RefineableQPoissonElement<2,3> > problem;
 
 // Setup labels for output
 //------------------------
 DocInfo doc_info, mesh_doc_info;
 bool report_stats=true;
 
 // Set output directory
 doc_info.set_directory("RESLT_incremental2"); 
 mesh_doc_info.set_directory("RESLT_incr_mesh"); 
 
 // Step number
 doc_info.number()=1;
 mesh_doc_info.number()=10;
  
 // Increase minimum permitted error from default (10e-5) to 10e-4
 problem.mesh_pt()->min_permitted_error()=1.0e-4;
  
 // Doc (default) refinement targets
 //----------------------------------
 problem.mesh_pt()->doc_adaptivity_targets(cout);
 
 // solve
 problem.newton_solve();
  
 //Output solution
 doc_info.label()="initial";
 problem.doc_solution(doc_info);
 
 //Increment counter for solutions 
 doc_info.number()++;

 //Loop for refining uniformly
 for (unsigned mym=0;mym<2;mym++)
  {
   //On the first step refine uniformly twice and then distribute
   if (mym==0) 
    {
     problem.refine_uniformly();
     problem.refine_uniformly();

     oomph_info << "Distributing problem" << std::endl;
     std::ifstream input_file;
     std::ofstream output_file;
     char filename[100];

     // Get the partition to be used from file
     unsigned n_partition=problem.mesh_pt()->nelement();
     Vector<unsigned> element_partition(n_partition,0);
     sprintf(filename,"fish_incremental_partition.dat");
     input_file.open(filename);
     std::string input_string;
     for (unsigned e=0;e<n_partition;e++)
      {
       getline(input_file,input_string,'\n');
       element_partition[e]=atoi(input_string.c_str());
      }

     // Distribute and check halo schemes
     problem.distribute(element_partition,mesh_doc_info,report_stats);
     problem.check_halo_schemes(mesh_doc_info);
    }
  
   // another round of refinement
   oomph_info << "Refine again..." << std::endl;
   problem.refine_uniformly();
 
   // Solve the problem 
   oomph_info << "Solve again..." << std::endl;

   problem.newton_solve();
 
   //Output solution
   doc_info.label()="after first solve after distribution and uniform refinement";
   problem.doc_solution(doc_info);
 
   //Increment counter for solutions 
   doc_info.number()++;
  }

 oomph_info << "Automatic adaptation starting" << std::endl << std::endl;

 // Now do (up to) four rounds of fully automatic adapation in response to 
 //-----------------------------------------------------------------------
 // error estimate
 //---------------
 unsigned max_solve=4;
 for (unsigned isolve=0;isolve<max_solve;isolve++)
  {
   oomph_info << "Adapting problem - isolve=" << isolve << std::endl;

   // Adapt problem/mesh
   problem.adapt(); 

   // Re-solve the problem if the adaptation has changed anything
#ifdef OOMPH_HAS_MPI
   // Make sure all processors know if refinement is taking place
   // Adaptation only converges if ALL the processes have no
   // refinement or unrefinement to perform
   unsigned total_refined;
   unsigned total_unrefined;
   unsigned n_refined=problem.mesh_pt()->nrefined();
   unsigned n_unrefined=problem.mesh_pt()->nunrefined();

   MPI_Allreduce(&n_refined,&total_refined,1,MPI_INT,MPI_SUM,
                 problem.communicator_pt()->mpi_comm());
   n_refined=total_refined;

   MPI_Allreduce(&n_unrefined,&total_unrefined,1,MPI_INT,MPI_SUM,
                 problem.communicator_pt()->mpi_comm());
   n_unrefined=total_unrefined;
#endif

   oomph_info << "---> " << n_refined << " elements to be refined, and " 
              << n_unrefined << " to be unrefined, in total." << std::endl;

   if ((n_refined!=0)||(n_unrefined!=0))
    {
     problem.newton_solve();
    }
   else
    {
     cout << "Mesh wasn't adapted --> we'll stop here" << std::endl;
     break;
    }
   
   //Output solution
   doc_info.label()="solve after adapt";
   problem.doc_solution(doc_info);
   
   //Increment counter for solutions 
   doc_info.number()++;
  }

 // Loop for uniform unrefinement
 for(unsigned myn=0;myn<10;myn++)
  {
   // Test the unrefine_uniformly command and re-solve
   //-------------------------------------------------
   problem.unrefine_uniformly();
 
   // Solve the problem 
   problem.newton_solve();
 
   //Output solution
   doc_info.label()="solve after uniform unrefinement";
   problem.doc_solution(doc_info);
 
   //Increment counter for solutions 
   doc_info.number()++;
  }

 oomph_info << "Min/max end incremental\n";

} // end of incremental



//================================start_black_box=========================
/// Demonstrate how to solve 2D Poisson problem in 
/// fish-shaped domain with fully automatic mesh adaptation
//========================================================================
void solve_with_fully_automatic_adaptation()
{
 //Set up the problem with nine-node refineable Poisson elements
 RefineableFishPoissonProblem<RefineableQPoissonElement<2,3> > problem;
  
 // Setup labels for output
 //------------------------
 DocInfo doc_info, mesh_doc_info;
 bool report_stats=true;
  
 // Set output directory
 doc_info.set_directory("RESLT_fully_automatic"); 
 mesh_doc_info.set_directory("RESLT_adapt_mesh"); 
  
 // Step number
 doc_info.number()=0;
 mesh_doc_info.number()=20;

 // Doc (default) refinement targets
 //----------------------------------
 problem.mesh_pt()->doc_adaptivity_targets(cout);

 // Solve/doc the problem with fully automatic adaptation
 //------------------------------------------------------

 // Refine coarse original mesh first
 //----------------------------------
 problem.refine_uniformly();
 problem.refine_uniformly();

 oomph_info << "-----------------------------------------" << std::endl;
 oomph_info << "--- Distributing problem (fully auto) ---" << std::endl;
 oomph_info << "-----------------------------------------" << std::endl;

 std::ifstream input_file;
 std::ofstream output_file;
 char filename[100];

 // Get the partition to be used from file
 unsigned n_partition=problem.mesh_pt()->nelement();
 Vector<unsigned> element_partition(n_partition,0);
 sprintf(filename,"fish_fully_automatic_partition.dat");
 input_file.open(filename);
 std::string input_string;
 for (unsigned e=0;e<n_partition;e++)
  {
   getline(input_file,input_string,'\n');
   element_partition[e]=atoi(input_string.c_str());
  }

 // Distribute and check halo schemes
 problem.distribute(element_partition,mesh_doc_info,report_stats);
 problem.check_halo_schemes(mesh_doc_info);

 //Maximum number of adaptations:
 unsigned max_adapt=5;

 oomph_info << "Solve with max_adapt=" << max_adapt << std::endl;

 //Solve the problem; perform up to specified number of adaptations.
 problem.newton_solve(max_adapt);

 //Output solution
 oomph_info << "-----------------------" << std::endl;
 oomph_info << "Now output the solution" << std::endl;
 oomph_info << "-----------------------" << std::endl;
 problem.doc_solution(doc_info);   

 //increment doc_info number
 doc_info.number()++;

} // end black box

//========================================================================
/// Solve 2D fish Poisson problem with a selected refinement pattern
//========================================================================
void solve_with_selected_refinement_pattern()
{
 //Set up the problem with nine-node refineable Poisson elements
 RefineableFishPoissonProblem<RefineableQPoissonElement<2,3> > problem;
  
 // Setup labels for output
 //------------------------
 DocInfo doc_info, mesh_doc_info;
 bool report_stats=true;
  
 // Set output directory
 doc_info.set_directory("RESLT_select_refine"); 
 mesh_doc_info.set_directory("RESLT_select_mesh"); 
  
 // Step number
 doc_info.number()=0;
 mesh_doc_info.number()=20;

 // Doc (default) refinement targets
 //----------------------------------
 problem.mesh_pt()->doc_adaptivity_targets(cout);

 // Refine coarse original mesh first
 //----------------------------------
 problem.refine_uniformly();
 problem.refine_uniformly();

 // Distribute the problem and doc mesh info
 oomph_info << "----------------------------------------" << std::endl;
 oomph_info << "--- Distributing problem (selective) ---" << std::endl;
 oomph_info << "----------------------------------------" << std::endl;

 std::ifstream input_file;
 std::ofstream output_file;
 char filename[100];

 // Get the partition to be used from file
 unsigned n_partition=problem.mesh_pt()->nelement();
 Vector<unsigned> element_partition(n_partition,0);
 sprintf(filename,"fish_selective_partition.dat");
 input_file.open(filename);
 std::string input_string;
 for (unsigned e=0;e<n_partition;e++)
  {
   getline(input_file,input_string,'\n');
   element_partition[e]=atoi(input_string.c_str());
  }

 // Distribute and check halo schemes
 problem.distribute(element_partition,mesh_doc_info,report_stats);
 problem.check_halo_schemes(mesh_doc_info);

 // Add some selective refinement on top of this
 unsigned nsoln=3; // number of selective refinements - validate, 3

 for (unsigned isoln=0; isoln<nsoln; isoln++)
  {
   // Refine selected elements based upon y-coordinate
   // of the central node of the element (node 4)
   Vector<unsigned> more_middle_els;
   unsigned nels=problem.mesh_pt()->nelement();
   for (unsigned e=0; e<nels; e++)
    {
     FiniteElement* el_pt=dynamic_cast<FiniteElement*>
      (problem.mesh_pt()->element_pt(e));
     if ((el_pt->node_pt(4)->x(1)>=(-0.35)) &&
         (el_pt->node_pt(4)->x(1)<=0.35))
      {
       more_middle_els.push_back(e);
      }
    }

   oomph_info << "--------------------------------" << std::endl;
   oomph_info << "This time, refine " << more_middle_els.size() 
              << " selected elements" << std::endl;
   oomph_info << "--------------------------------" << std::endl;
   problem.refine_selected_elements(more_middle_els);

   // Now solve the problem
   problem.newton_solve();

   // Doc the solution
   problem.doc_solution(doc_info);   

   // Increment doc_info number
   doc_info.number()++;
  }

}

//=================start_of_main==========================================
/// Demonstrate how to solve 2D Poisson problem in 
/// fish-shaped domain with mesh adaptation.
//========================================================================
int main(int argc, char **argv)
{
 // initialise MPI
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif

 // Solve using a pre-selected refinement pattern
 solve_with_selected_refinement_pattern();

 // Solve with adaptation, docing the intermediate steps
 solve_with_incremental_adaptation();

 // Solve directly, with fully automatic adaptation
 solve_with_fully_automatic_adaptation();

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif

} // end of main
