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
// Driver for adaptive 2D quarter circle driven cavity. Solved with
// hp-adaptation, using Crouzeix-Raviart elements.

// Generic oomph-lib header
#include "generic.h"

// Navier Stokes headers
#include "navier_stokes.h"

// The mesh
#include "meshes/quarter_circle_sector_mesh.h"

using namespace std;

using namespace oomph;

//==start_of_namespace===================================================
/// Namespace for physical parameters
//=======================================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re=100;

 /// Reynolds/Froude number
 double Re_invFr=100;

 /// Gravity vector
 Vector<double> Gravity(2);

 /// Functional body force
 void body_force(const double& time, const Vector<double>& x, 
                 Vector<double>& result)
 {
  result[0]=0.0;
  result[1]=-Re_invFr;
 }

 /// Zero functional body force
 void zero_body_force(const double& time, const Vector<double>& x, 
                      Vector<double>& result)
 {
  result[0]=0.0;
  result[1]=0.0;
 }

} // end_of_namespace


//==start_of_problem_class============================================
/// Driven cavity problem in quarter circle domain, templated
/// by element type. 
//====================================================================
template<class ELEMENT>
class QuarterCircleDrivenCavityProblem : public Problem
{

public:

 /// Constructor
 QuarterCircleDrivenCavityProblem(
  NavierStokesEquations<2>::NavierStokesBodyForceFctPt body_force_fct_pt,
  const unsigned &driven_bdry);

 /// Destructor: Empty
 ~QuarterCircleDrivenCavityProblem() {}

 /// Update the after solve (empty)
 void actions_after_newton_solve() {}

 /// Update the problem specs before solve. 
 /// (Re-)set velocity boundary conditions just to be on the safe side...
 void actions_before_newton_solve()
  {
   if(Driven_boundary == 0)
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
    }
   else
    {
     // Setup tangential flow along boundary 1:
     unsigned ibound=1; 
     unsigned num_nod= mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // get coordinates
       double x=mesh_pt()->boundary_node_pt(ibound,inod)->x(0);
       double y=mesh_pt()->boundary_node_pt(ibound,inod)->x(1);
       // find Lagrangian coordinate (the angle)
       double zeta=0.0;
       if (x!=0.0)
        {
         zeta=atan(y/x);
        }
       // Tangential flow u0
       unsigned i=0;
       mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,-sin(zeta));
       // Tangential flow u1
       i=1;
       mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,cos(zeta));
      }

     // Overwrite with no flow along all boundaries
     unsigned num_bound = mesh_pt()->nboundary();
     for(unsigned ibound=0;ibound<num_bound;ibound++)
      {
       if (ibound!=1)
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
   
   // Now pin the first pressure dof in the first element and set it to 0.0
   fix_pressure(0,0,0.0);
  } // end_of_actions_after_adapt
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
private:

 /// Pointer to body force function
 NavierStokesEquations<2>::NavierStokesBodyForceFctPt Body_force_fct_pt;

 /// Which boundary is driven?
 unsigned Driven_boundary;

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to proper element and fix pressure
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure

}; // end_of_problem_class



//==start_of_constructor==================================================
/// Constructor for driven cavity problem in quarter circle domain
//========================================================================
template<class ELEMENT>
QuarterCircleDrivenCavityProblem<ELEMENT>::QuarterCircleDrivenCavityProblem(
 NavierStokesEquations<2>::NavierStokesBodyForceFctPt body_force_fct_pt,
 const unsigned &driven_bdry) :
 Body_force_fct_pt(body_force_fct_pt), Driven_boundary(driven_bdry)
{ 

 // Build geometric object that parametrises the curved boundary
 // of the domain

 // Half axes for ellipse
 double a_ellipse=1.0;
 double b_ellipse=1.0;

 // Setup elliptical ring 
 GeomObject* Wall_pt=new Ellipse(a_ellipse,b_ellipse);

 // End points for wall
 double xi_lo=0.0;
 double xi_hi=2.0*atan(1.0);

 //Now create the mesh
 double fract_mid=0.5;
 Problem::mesh_pt() = new 
  RefineableQuarterCircleSectorMesh<ELEMENT>(
   Wall_pt,xi_lo,fract_mid,xi_hi);

 // Raise the max refinement level
 dynamic_cast<TreeBasedRefineableMeshBase*>(mesh_pt())->max_refinement_level()=7;
 
 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 dynamic_cast<RefineableQuarterCircleSectorMesh<ELEMENT>*>(
  mesh_pt())->spatial_error_estimator_pt()=error_estimator_pt;
 
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
   //Set the Re/Fr
   el_pt->re_invfr_pt() = &Global_Physical_Variables::Re_invFr;
   //Set Gravity vector
   el_pt->g_pt() = &Global_Physical_Variables::Gravity;
   //set body force function
   el_pt->body_force_fct_pt() = Body_force_fct_pt;

  } // end loop over elements

 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 
 // Now pin the first pressure dof in the first element and set it to 0.0
 fix_pressure(0,0,0.0);
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} // end_of_constructor



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void QuarterCircleDrivenCavityProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 


 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
 
} // end_of_doc_solution




//==start_of_main======================================================
/// Driver for QuarterCircleDrivenCavityProblem test problem 
//=====================================================================
int main()
{
 //Drive straight and curved boundaries in turn
 for(unsigned driven_bdry=0; driven_bdry<2; driven_bdry++)
  {
   // Set output directory
   DocInfo doc_info;
   if(driven_bdry==0)
    doc_info.set_directory("RESLT_circ_bdry0");
   else
    doc_info.set_directory("RESLT_circ_bdry1");

   // Solve problem 1 with hp-refineable Crouzeix-Raviart elements
   //-------------------------------------------------------------
   {
    // Set up downwards-Gravity vector
    Global_Physical_Variables::Gravity[0] = 0.0;
    Global_Physical_Variables::Gravity[1] = -1.0;
  
    // Set up Gamma vector for stress-divergence form
    NavierStokesEquations<2>::Gamma[0]=1;
    NavierStokesEquations<2>::Gamma[1]=1;

    // Build problem with Gravity vector in stress divergence form, 
    // using zero body force function
    QuarterCircleDrivenCavityProblem<PRefineableQCrouzeixRaviartElement<2> > 
     problem(&Global_Physical_Variables::zero_body_force,driven_bdry);

    //Set initial doc_info number
    doc_info.number()=10;
 
    // Initial refinement level
    problem.refine_uniformly();
    problem.refine_uniformly();
    problem.p_refine_uniformly();

    problem.newton_solve();
    cout << "output file: " << doc_info.number() << endl;
    problem.doc_solution(doc_info);
    doc_info.number()++;
  
    problem.adapt();
    problem.newton_solve();
    cout << "output file: " << doc_info.number() << endl;
    problem.doc_solution(doc_info);
    doc_info.number()++;
  
    problem.p_adapt();
    problem.newton_solve();
    cout << "output file: " << doc_info.number() << endl;
    problem.doc_solution(doc_info);
    doc_info.number()++;
  
    problem.adapt();
    problem.newton_solve();
    cout << "output file: " << doc_info.number() << endl;
    problem.doc_solution(doc_info);
    doc_info.number()++;
  
    //Doc again for check
    doc_info.number()=0;
    problem.doc_solution(doc_info);

   } // end of problem 1


   // Solve problem 2 with hp-refineable Crouzeix-Raviart elements
   //-------------------------------------------------------------
   {
    // Set up zero-Gravity vector
    Global_Physical_Variables::Gravity[0] = 0.0;
    Global_Physical_Variables::Gravity[1] = 0.0;

    // Set up Gamma vector for simplified form
    NavierStokesEquations<2>::Gamma[0]=0;
    NavierStokesEquations<2>::Gamma[1]=0;

    // Build problem with body force function and simplified form,
    // using body force function
    QuarterCircleDrivenCavityProblem<PRefineableQCrouzeixRaviartElement<2> >
     problem(&Global_Physical_Variables::body_force,driven_bdry);

    //Set initial doc_info number
    doc_info.number()=20;

    // Initial refinement level
    problem.refine_uniformly();
    problem.refine_uniformly();
    problem.p_refine_uniformly();

    problem.newton_solve();
    cout << "output file: " << doc_info.number() << endl;
    problem.doc_solution(doc_info);
    doc_info.number()++;
  
    problem.adapt();
    problem.newton_solve();
    cout << "output file: " << doc_info.number() << endl;
    problem.doc_solution(doc_info);
    doc_info.number()++;
  
    problem.p_adapt();
    problem.newton_solve();
    cout << "output file: " << doc_info.number() << endl;
    problem.doc_solution(doc_info);
    doc_info.number()++;
  
    problem.adapt();
    problem.newton_solve();
    cout << "output file: " << doc_info.number() << endl;
    problem.doc_solution(doc_info);
    doc_info.number()++;

    //Doc again for check
    doc_info.number()=1;
    problem.doc_solution(doc_info);
  
   } // end of problem 2

  }//end of both boundaries

} // end_of_main


