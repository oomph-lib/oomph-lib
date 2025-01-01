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
//Driver for 3D rectangular driven cavity

#include <fenv.h>

//Generic includes
#include "generic.h"
#include "navier_stokes.h"

#include "meshes/simple_cubic_mesh.h"
#include "meshes/simple_cubic_tet_mesh.h"

using namespace std;

using namespace oomph;



/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

#ifdef OOMPH_HAS_HYPRE

//=============================================================================
/// helper method for the block diagonal F block preconditioner to allow 
/// hypre to be used for as a subsidiary block preconditioner
//=============================================================================
namespace Hypre_Subsidiary_Preconditioner_Helper
{
 Preconditioner* set_hypre_preconditioner()
 {
  return new HyprePreconditioner;
 }
}

#endif

 
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Variables
{

 /// Enumeration for the problem ids
 enum {Driven_cavity, Through_flow};


 /// Reynolds number
 double Re=50.0;

 /// Storage for number of iterations during Newton steps 
 Vector<unsigned> Iterations;

 /// Storage for linear solver times during Newton steps 
 Vector<double> Linear_solver_time;

 /// Traction at the outflow boundary
 void prescribed_traction(const double& t,
                          const Vector<double>& x,
                          const Vector<double> &n,
                          Vector<double>& traction)
 {
  traction.resize(3);
  traction[0]=1.0;
  traction[1]=0.0;
  traction[2]=0.0;
 } 

} // end_of_namespace



/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Test problem for Fp/PCD preconditioner
//====================================================================
class FpTestProblem : public Problem
{

public:


 /// Constructor
 FpTestProblem(const unsigned& n_element,
               const bool& use_tets,
               const bool& use_adaptivity, 
               const bool& use_lsc,
               const bool& use_hypre_for_pressure,
               const bool& use_block_diagonal_for_momentum,
               const bool& use_hypre_for_momentum_diagonals,
               const int& problem_id);

 
 /// Destructor: Cleanup
 ~FpTestProblem()
  {
   delete Solver_pt;
   delete Prec_pt;
   delete P_matrix_preconditioner_pt;
   delete F_matrix_preconditioner_pt;
  }

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to full element type and fix the pressure at that element
   dynamic_cast<NavierStokesEquations<3>*>(Bulk_mesh_pt->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end of fix_pressure


 
 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableNavierStokesEquations<3>::
    unpin_all_pressure_dofs(Bulk_mesh_pt->element_pt());
   
   // Pin redundant pressure dofs
   RefineableNavierStokesEquations<3>::
    pin_redundant_nodal_pressures(Bulk_mesh_pt->element_pt());
   
   // Now set the first pressure dof in the first element to 0.0
   if (Problem_id==Global_Variables::Driven_cavity) fix_pressure(0,0,0.0);
  } // end_of_actions_after_adapt
 

 /// Actions after Newton step record number of iterations
 void actions_after_newton_step() 
  {                               
   Global_Variables::Iterations.push_back(
    dynamic_cast<IterativeLinearSolver*>
    (this->linear_solver_pt())->iterations());
   
   Global_Variables::Linear_solver_time.push_back(
    linear_solver_pt()->linear_solver_solution_time());
  }  

 /// Update the after solve (empty)
 void actions_after_newton_solve(){}

 /// Update the problem specs before solve. 
 void actions_before_newton_solve()
 {
  // Initialise counter for iterations
  Global_Variables::Iterations.clear();
  Global_Variables::Linear_solver_time.clear();

  // Driven cavity bcs
  if (Problem_id==Global_Variables::Driven_cavity)
   {
    // Setup tangential flow along driven boundary
    unsigned ibound=Driven_boundary;
    unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      // Tangential flow
      Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,1.0);
      Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(1,0.0);
      
      // No penetration
      Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(2,0.0);
     }
    
    // Overwrite with no flow along the other boundaries
    unsigned num_bound = Bulk_mesh_pt->nboundary();
    for(unsigned ibound=0;ibound<num_bound;ibound++)
     {
      if (ibound!=Driven_boundary)
       {
        unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
        for (unsigned inod=0;inod<num_nod;inod++)
         {
          for (unsigned i=0;i<3;i++)
           {
            Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(i,0.0);
           }
         }
       }
     }
   }
  // Inflow at left boundary
  else
   {
    // Inflow in upper half of inflow boundary
    unsigned ibound=Inflow_boundary; 
    unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
      double y=nod_pt->x(1);
      double z=nod_pt->x(2);
      if ((y>0.5)&&(z>0.5))
       {
        double u=(y-0.5)*(1.0-y)*(z-0.5)*(1.0-z);
        nod_pt->set_value(0,u);
       }
      else
       {
        nod_pt->set_value(0,0.0);
       }
      nod_pt->set_value(1,0.0);
      nod_pt->set_value(2,0.0);
     }
   }
  
 } // end_of_actions_before_newton_solve

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Create traction elements on outflow boundary
 template<class ELEMENT>
 void create_traction_elements();

 /// Create refineable traction elements on outflow boundary
 template<class ELEMENT>
 void create_refineable_traction_elements();

 /// Validate fp
 template<class FP_ELEMENT>
 void validate_fp()
  {
   DocInfo my_doc_info;
   dynamic_cast<NavierStokesSchurComplementPreconditioner*>(
    Prec_pt)->template validate<FP_ELEMENT>(my_doc_info,this);
   pause("done validation");
  }

 /// Pointer to the "bulk" mesh
 Mesh*& bulk_mesh_pt() {return Bulk_mesh_pt;}

private:

 /// Solver
 IterativeLinearSolver* Solver_pt;

 /// Solver
 Preconditioner* Prec_pt;

 /// Inexact solver for P block
 Preconditioner* P_matrix_preconditioner_pt;

 /// Inexact solver for F block
 Preconditioner* F_matrix_preconditioner_pt;

 /// ID of driven boundary
 unsigned Driven_boundary;

 /// ID of inflow boundary
 unsigned Inflow_boundary;

 /// ID of outflow boundary
 unsigned Outflow_boundary;

 /// Problem id
 unsigned Problem_id;

 /// Pointer to the "bulk" mesh
 Mesh* Bulk_mesh_pt;
 
 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_pt;

}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for DrivenCavity problem 
//========================================================================
FpTestProblem::FpTestProblem(const unsigned& n_el,
                             const bool& use_tets,
                             const bool& use_adaptivity, 
                             const bool& use_lsc,
                             const bool& use_hypre_for_pressure,
                             const bool& use_block_diagonal_for_momentum,
                             const bool& use_hypre_for_momentum_diagonals,
                             const int& problem_id)
{ 
 
 // Store flag
 Problem_id=problem_id;

 // Setup mesh
 
 // # of elements in x-direction
 unsigned n_x=n_el;
 
 // # of elements in y-direction
 unsigned n_y=n_el;

 // # of elements in z-direction
 unsigned n_z=n_el;
 
 // Domain length in x-direction
 double l_x=1.0;
 
 // Domain length in y-direction
 double l_y=1.0;
 
 // Domain length in y-direction
 double l_z=1.0;
 

 // Build and assign mesh
 if (use_tets)
  {
   Bulk_mesh_pt = new SimpleCubicTetMesh<TTaylorHoodElement<3> >
    (n_x,n_y,n_z,l_x,l_y,l_z);
     
   dynamic_cast<SimpleCubicTetMesh<TTaylorHoodElement<3> >*>(Bulk_mesh_pt)
    ->setup_boundary_element_info();

   Driven_boundary=4;
   Inflow_boundary=0;
   Outflow_boundary=1;
  }
 else
  {
   if (use_adaptivity)
    {
     Bulk_mesh_pt = 
      new RefineableSimpleCubicMesh<RefineableQTaylorHoodElement<3> >
      (n_x,n_y,n_z,l_x,l_y,l_z);
     
     Driven_boundary=0;
     Inflow_boundary=4;
     Outflow_boundary=2;     
    }
   else
    {
     Bulk_mesh_pt = 
      new SimpleCubicMesh<QTaylorHoodElement<3> >(n_x,n_y,n_z,l_x,l_y,l_z);
     
     Driven_boundary=0;
     Inflow_boundary=4;
     Outflow_boundary=2;
    }
  } 

 // Create "surface mesh" that will contain only the prescribed-traction 
 // elements.
 Surface_mesh_pt = new Mesh;

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();
 
 // Build preconditoner
 NavierStokesSchurComplementPreconditioner* prec_pt = 
  new NavierStokesSchurComplementPreconditioner(this);
 Prec_pt=prec_pt;
   
 
 // By default, the Schur Complement Preconditioner uses SuperLU as
 // an exact preconditioner (i.e. a solver) for the
 // momentum and Schur complement blocks. 
 // Can overwrite this by passing pointers to 
 // other preconditioners that perform the (approximate)
 // solves of these blocks.


 // Create internal preconditioners used on Schur block
 P_matrix_preconditioner_pt=0;
 if (use_hypre_for_pressure)
  {
#ifdef OOMPH_HAS_HYPRE

   // Create preconditioner
   P_matrix_preconditioner_pt = new HyprePreconditioner;
   
   // Set parameters for use as preconditioner on Poisson-type problem
   Hypre_default_settings::set_defaults_for_2D_poisson_problem(
    static_cast<HyprePreconditioner*>(P_matrix_preconditioner_pt));
   
   // Use Hypre for the Schur complement block
   prec_pt->set_p_preconditioner(P_matrix_preconditioner_pt);
   
   // Shut up!
   static_cast<HyprePreconditioner*>(P_matrix_preconditioner_pt)->
    disable_doc_time();
   
#endif
  }
 
 // Create block-diagonal preconditioner used on momentum block
 F_matrix_preconditioner_pt=0;   
 if (use_block_diagonal_for_momentum)
  {
   
   F_matrix_preconditioner_pt = 
    new BlockDiagonalPreconditioner<CRDoubleMatrix>;
   
   // Use Hypre as block preconditioner
   if (use_hypre_for_pressure)
    {
#ifdef OOMPH_HAS_HYPRE
     dynamic_cast<BlockDiagonalPreconditioner<CRDoubleMatrix>* >
      (F_matrix_preconditioner_pt)->set_subsidiary_preconditioner_function
      (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_preconditioner);
#endif       
    }
   
   // Use Hypre for momentum block 
   prec_pt->set_f_preconditioner(F_matrix_preconditioner_pt);
  }
 
 
   
   // Use LSC?
   if (use_lsc)
    {
     prec_pt->use_lsc();
    }
   else
    {
     prec_pt->use_fp();
    }
   
   // Set Navier Stokes mesh
   prec_pt->set_navier_stokes_mesh(Bulk_mesh_pt);    
   
   
   // Set the boundary conditions for this problem: All nodes are
   // free by default -- just pin the ones that have Dirichlet conditions
   // here. 
   unsigned num_bound = Bulk_mesh_pt->nboundary();
   for(unsigned ibound=0;ibound<num_bound;ibound++)
    {
     unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Loop over values (u, v and w velocities)
       for (unsigned i=0;i<3;i++)
        {
         Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(i); 
        }
      }
    } // end loop over boundaries

 

 // In/outflow bcs
 if (Problem_id==Global_Variables::Through_flow)
  {
   unsigned ibound=Outflow_boundary;
   unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
     // Only free if node is ONLY on a single boundary
     std::set<unsigned>* bnd_pt=0;
     nod_pt->get_boundaries_pt(bnd_pt);
     if (bnd_pt!=0)
      {
       if (bnd_pt->size()<2)
        {
         if (!(nod_pt->is_on_boundary(0)))
          {
           if ((nod_pt->x(1)<0.5)&&nod_pt->x(2)<0.5) nod_pt->unpin(0);
          }
        }
      }
    }
  }

 // Complete the build of all elements so they are fully functional

 //Find number of elements in mesh
 unsigned n_element = Bulk_mesh_pt->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   NavierStokesEquations<3>* el_pt = 
    dynamic_cast<NavierStokesEquations<3>*>(Bulk_mesh_pt->element_pt(e));
   
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Variables::Re;
  } // end loop over elements
 
 
 // Pin redundant pressure dofs
 if (use_adaptivity)
  {
   RefineableNavierStokesEquations<3>::
    pin_redundant_nodal_pressures(Bulk_mesh_pt->element_pt());
  }

 // Now set the first pressure value in element 0 to 0.0
 if (Problem_id==Global_Variables::Driven_cavity) fix_pressure(0,0,0.0);

 // Setup equation numbering scheme
 oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

#ifdef OOMPH_HAS_TRILINOS

 // Build iterative linear solver
 oomph_info << "Using Trilinos GMRES\n"; 
 TrilinosAztecOOSolver* iterative_linear_solver_pt = 
  new TrilinosAztecOOSolver;

 Solver_pt=iterative_linear_solver_pt;

#else

 // Build solve and preconditioner
 Solver_pt = new GMRES<CRDoubleMatrix>;
 dynamic_cast<GMRES<CRDoubleMatrix>*>(Solver_pt)->set_preconditioner_RHS();

#endif

 // Set solver and preconditioner
 Solver_pt->preconditioner_pt() = Prec_pt;
 linear_solver_pt() = Solver_pt;
 
} // end_of_constructor


//============start_of_create_traction_elements==========================
/// Create Navier-Stokes traction elements on outflow boundary
//=======================================================================
template<class ELEMENT>
void FpTestProblem::create_refineable_traction_elements()
{

 unsigned b=Outflow_boundary;

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = Bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->boundary_element_pt(b,e));
   
   //What is the index of the face of element e along boundary b
   int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
   
   // Build the corresponding prescribed-flux element
   RefineableNavierStokesTractionElement<ELEMENT>* flux_element_pt = new 
    RefineableNavierStokesTractionElement<ELEMENT>(bulk_elem_pt,face_index);
   
   //Add the prescribed-flux element to the surface mesh
   Surface_mesh_pt->add_element_pt(flux_element_pt);
   
   // Set the pointer to the prescribed traction function
   flux_element_pt->traction_fct_pt() = &Global_Variables::prescribed_traction;
   
  } //end of loop over bulk elements adjacent to boundary b

 // Now rebuild the global mesh
 rebuild_global_mesh();

 // Reassign equation numbers
 oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of create_traction_elements






//============start_of_create_traction_elements==========================
/// Create Navier-Stokes traction elements on outflow boundary
//=======================================================================
template<class ELEMENT>
void FpTestProblem::create_traction_elements()
{

 unsigned b=Outflow_boundary;

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = Bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->boundary_element_pt(b,e));
   
   //What is the index of the face of element e along boundary b
   int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
   
   // Build the corresponding prescribed-flux element
   NavierStokesTractionElement<ELEMENT>* flux_element_pt = new 
      NavierStokesTractionElement<ELEMENT>(bulk_elem_pt,face_index);

   //Add the prescribed-flux element to the surface mesh
   Surface_mesh_pt->add_element_pt(flux_element_pt);
   
   // Set the pointer to the prescribed traction function
   flux_element_pt->traction_fct_pt() = &Global_Variables::prescribed_traction;
   
  } //end of loop over bulk elements adjacent to boundary b

 // Now rebuild the global mesh
 rebuild_global_mesh();

 // Reassign equation numbers
 oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of create_traction_elements



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
void FpTestProblem::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

} // end_of_doc_solution





/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//==start_of_main======================================================
/// Driver for Fp preconditioner
//=====================================================================
int main(int argc, char **argv)
{


#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif


 //Label for output
 DocInfo doc_info;
 
 //Set output directory
 doc_info.set_directory("RESLT");
 
 //Doc number of gmres iterations
 ofstream out_file;
 
 // Set flags
 bool use_hypre_for_pressure=true;
 bool use_block_diagonal_for_momentum=false;
 bool use_hypre_for_momentum_diagonals=false;
   
 //Loop over problems: Driven cavity and step
 for (unsigned problem_id=0;problem_id<2;problem_id++) 
  {
   
   if (problem_id==0)
    {
     out_file.open("three_d_iter_driven_cavity.dat");
    }
   else 
    {      
     out_file.open("three_d_iter_through_flow.dat");
    }
   
   out_file
    << "VARIABLES=\"nel_1d\","
    << "\"ndof\"," 
    << "\"Re\"," 
    << "\" Newton iteration\","
    << "\"GMRES iterations\","
    << "\"Linear solver time\","
    << "\"doc number\""
    << std::endl;
   
   std::string header1;
   std::string header2;
   
     
   // Loop over preconditioners: iprec=0: LSC
   //                            iprec=1: Fp
   bool use_lsc=true;
   //bool use_robin=true;
   for (unsigned iprec=0;iprec<2;iprec++)
    {
     
     // Loop over three cases (tets, non-refineable/refineable bricks)
     for (unsigned icase=0;icase<=2;icase++) 
      {
       bool use_tets=false;
       bool use_adaptivity=false;
       
       switch(icase)
        {
         
        case 0:
         
         header1=" Tets";
         use_tets=true;
         use_adaptivity=false;
         break;
         
        case 1:
         
         header1=" Bricks";
         use_tets=false;
         use_adaptivity=false;
         break;
         
        case 2:
         
         header1=" Refineable Bricks";
         use_tets=false;
         use_adaptivity=true;
         break;
         
        default:
         break;
        }
       
       oomph_info << "Doing it with " << header1 << std::endl;
       
       // Set preconditioner
       if (iprec==0)
        {
         use_lsc=true;
         header2=", LSC";
        }
       else if (iprec==1)
        {
         use_lsc=false;
         //use_robin=true;
         header2=", Fp";
        }   
       
        oomph_info << "Doing it with " << header2 << " preconditioner\n";

        // Write tecplot header
        string header="ZONE T=\""+header1+header2+"\"\n";
        out_file << header;
                 
        // Number of elements in x/y directions (reduced for validaton)
        unsigned max_nel_1d=4; 
        if (argc>1) max_nel_1d=2;
        for (unsigned nel_1d = 2; nel_1d <= max_nel_1d; nel_1d*=2) 
         {
           
          // Build the problem 
          FpTestProblem problem(
           nel_1d,use_tets,use_adaptivity,use_lsc,
           use_hypre_for_pressure,use_block_diagonal_for_momentum,
           use_hypre_for_momentum_diagonals,problem_id);
           
          
          
          // Refine a few times
          if (use_adaptivity)
           {             
            // Note: This manually chosen refinement pattern introduces
            // doubly hanging nodes on boundary and in the interior
            Vector<unsigned> elements_to_be_refined;
            unsigned nel=problem.bulk_mesh_pt()->nelement();
            unsigned e=0;
            while (e<nel)
             {
              elements_to_be_refined.push_back(e);
              e+=7;
             }
            problem.refine_selected_elements(0,elements_to_be_refined); 
            elements_to_be_refined.clear();
            elements_to_be_refined.push_back(1);
            problem.refine_selected_elements(0,elements_to_be_refined); 
            elements_to_be_refined.clear();
            elements_to_be_refined.push_back(4);
            problem.refine_selected_elements(0,elements_to_be_refined); 
           }
          // Attach traction elements now with the problem in its
          // most refined state
          if (use_tets)
           {
            problem.create_traction_elements<TTaylorHoodElement<3> >();
           }
          else 
           {
            if (use_adaptivity)
             {
              problem.create_refineable_traction_elements
               <RefineableQTaylorHoodElement<3> >();
             }
            else
             {
              problem.create_traction_elements<QTaylorHoodElement<3> >();
             }
           }
          
          
          
          // Loop over Reynolds numbers (limited during validation)
          double start_re = 0.0; 
          if (argc>1) start_re=50;
          double end_re = 50.0; 
          for (double re = start_re; re <= end_re; re+=50.0)
           {
            
            // Set Reynolds
            Global_Variables::Re=re;
            
            // Solve the problem 
            problem.newton_solve();
            
            // Doc solution
            problem.doc_solution(doc_info);
            doc_info.number()++;
            
            // Doc iteration counts for each Newton iteration
            unsigned ndof = problem.ndof();
            unsigned iter = Global_Variables::Iterations.size();

            // Doc for all Newton iterations or only the last one?
            unsigned j_lo=0;
            //j_lo=iter-1;
            for (unsigned j = j_lo; j < iter; j++)
             {
              out_file
               << nel_1d << " "
               << ndof << " "
               << re << " " 
               << j << " "
               << Global_Variables::Iterations[j] << " "
               << Global_Variables::Linear_solver_time[j] << " "
               << doc_info.number()-1 << " "
               << std::endl;
             }
            
           }
         }     
      }
     }
   
   out_file.close();
   
  }
 

#ifdef OOMPH_HAS_MPI
   MPI_Helpers::finalize();
#endif
   
   
} // end_of_main


