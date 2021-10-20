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

// Generic oomph-lib stuff
#include "generic.h"

//Specific mesh
#include "half_rectangle_with_hole_mesh.h"

// Navier Stokes
#include "multi_domain_axisym_boussinesq_elements.h"

using namespace std;

using namespace oomph;



//===============================================
/// Global parameters
//===============================================
namespace Global_Parameters
{

 /// Reynolds number
 double Re = 1.0;

 /// Prandtl number
 double Pr = 0.73;

 /// Peclet number
 double Pe = 1.0;

 /// Rayleigh number
 double Ra = 10.0;
 
 /// Gravity
 Vector<double> G;

 /// Location of the centre of the sphere on the axis
 double Sphere_centre_z= 0.0;

}


///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////




//===================================================================
/// Flow around a cylinder in rectangular domain
//===================================================================
template<class NST_ELEMENT,class AD_ELEMENT>
class FlowAroundHalfCylinderProblem : public Problem
{

public:

 /// Constructor: Pass geometric object that represents
 /// central cylinder, and length and height of domain.
 FlowAroundHalfCylinderProblem(GeomObject* cylinder_pt, 
                               const double &radius, 
                               const double &length);
 
 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}

 ///  Update the problem specs before solve (empty; all prescribed
 /// velocities are constant along their respective boundares, therefore
 /// their FE interpolation onto the newly created nodes is sufficiently
 /// accurate)
 void actions_before_newton_solve() {}

 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableAxisymmetricNavierStokesEquations::
    unpin_all_pressure_dofs(nst_mesh_pt()->element_pt());
    
   // Pin redundant pressure dofs
   RefineableAxisymmetricNavierStokesEquations::
    pin_redundant_nodal_pressures(nst_mesh_pt()->element_pt());

   //Pin all swirl velocities to zero
   unsigned n_node = nst_mesh_pt()->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     nst_mesh_pt()->node_pt(n)->pin(2);
     nst_mesh_pt()->node_pt(n)->set_value(2,0.0);
    }
   
   //Plug flow everywhere
   for(unsigned ibound=0;ibound<5;++ibound)
    {
     unsigned num_nod= nst_mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       //No radial flow on symmetry boundary
       if(ibound==3) 
        {
         nst_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
        }
       
       //No radial flow on outlet
       if(ibound==2) 
        {
         nst_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
        }


       //No flow on side boundary
       if(ibound==1)
        {
         nst_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
         nst_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
        }

       //No flow on inlet
       if(ibound==0)
        {
         nst_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
         nst_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
        }

       //No slip on sphere
       if(ibound==4)
        {
         nst_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
         nst_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
       }
      }
    }

   //Fixed concentration on sphere
   for(unsigned ibound=0;ibound<5;++ibound)
    {
     unsigned num_nod= adv_diff_mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       //No concentration on inlet
       if(ibound==0)
        {
         adv_diff_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
        }
       //No concentration on side boundary
       if(ibound==1)
        {
         adv_diff_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
        }
       //Fixed concentration on sphere
       if(ibound==4)
        {
         adv_diff_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,1.0);
       }
      }
    }
   // Set external elements for the multi-domain solution.
   Multi_domain_functions::
    setup_multi_domain_interactions<NST_ELEMENT,AD_ELEMENT>
    (this,nst_mesh_pt(),adv_diff_mesh_pt());

  } // end_of_actions_after_adapt


 /// Access function for the specific mesh
 RefineableHalfRectangleWithHoleMesh<NST_ELEMENT>* nst_mesh_pt() 
  {
   return dynamic_cast<RefineableHalfRectangleWithHoleMesh<NST_ELEMENT>*>
    (Nst_mesh_pt);
  }

 /// Access function for the specific mesh
 RefineableHalfRectangleWithHoleMesh<AD_ELEMENT>* adv_diff_mesh_pt() 
  {
   return dynamic_cast<RefineableHalfRectangleWithHoleMesh<AD_ELEMENT>*>
    (Adv_diff_mesh_pt);
  }


 private:
 
  /// Height of the domain
  double Domain_radius;

 /// Length of the domain
 double Domain_length;
  
 /// Navier Stokes mesh
 RefineableQuadMesh<NST_ELEMENT>* Nst_mesh_pt;

 /// Advection diffusion mesh
 RefineableQuadMesh<AD_ELEMENT>* Adv_diff_mesh_pt;

};




//========================================================================
/// Constructor 
//========================================================================
template<class NST_ELEMENT,class AD_ELEMENT>
FlowAroundHalfCylinderProblem<NST_ELEMENT,AD_ELEMENT>::FlowAroundHalfCylinderProblem(
 GeomObject* cylinder_pt, const double &radius, const double &length) 
{ 
 
 Domain_radius = radius;
 Domain_length = length;
 
 // Build meshes
 Nst_mesh_pt =
  new RefineableHalfRectangleWithHoleMesh<NST_ELEMENT>(cylinder_pt,
                                                       radius,length,
                                                       4.0,2,4.0,2,1.0,2);
 Adv_diff_mesh_pt =
  new RefineableHalfRectangleWithHoleMesh<AD_ELEMENT>(cylinder_pt,
                                                      radius,length,
                                                      4.0,2,4.0,2,1.0,2);
 
 
 // Set error estimatora
 Nst_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 Adv_diff_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

  
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 

 //Pin all swirl velocities
 unsigned n_node = nst_mesh_pt()->nnode();
 for(unsigned n=0;n<n_node;n++)
  {
   nst_mesh_pt()->node_pt(n)->pin(2);
  }

 // Navier--Stokes boundary conditions
 //Pin both velocities at all boundaries
 unsigned num_bound = nst_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= nst_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //Pin boundary conditions on cylinder side boundary and inlet
     if((ibound==0) || (ibound==1) || (ibound==4))
      {
       nst_mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
       nst_mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
      }
     //Otherwise pin radial flow
     else
      {
       nst_mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
      }
    }
  }

 // Concentration boundary conditions
 // Pin the concentration on the sphere, inlet and side wall
 num_bound = adv_diff_mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= adv_diff_mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //Pin boundary conditions on cylinder, side boundary  and inlet
     if((ibound==0) || (ibound==1) || (ibound==4))
      {
       adv_diff_mesh_pt()->boundary_node_pt(ibound,inod)->pin(0);
      }
    }
  }
    

  // Pin redudant pressure dofs
  RefineableAxisymmetricNavierStokesEquations::
   pin_redundant_nodal_pressures(nst_mesh_pt()->element_pt());
  
  // Pass pointer to Reynolds number to elements
  unsigned n_nst_elem=nst_mesh_pt()->nelement();
  for (unsigned e=0;e<n_nst_elem;e++)
   {
    NST_ELEMENT* el_pt = 
     dynamic_cast<NST_ELEMENT*>(nst_mesh_pt()->element_pt(e));
    
    
    // Set the Reynolds number (1/Pr in our non-dimensionalisation)
    el_pt->re_pt() = &Global_Parameters::Re;
    
    // Set ReSt (also 1/Pr in our non-dimensionalisation)
    el_pt->re_st_pt() = &Global_Parameters::Re;
    
    // Set the Rayleigh number
    el_pt->ra_pt() = &Global_Parameters::Ra;
    
    //Set Gravity vector
    el_pt->g_pt() = &Global_Parameters::G;
    
    // We can ignore the external geometric data in the "external"
    // advection diffusion element when computing the Jacobian matrix
    // because the interaction does not involve spatial gradients of 
    // the temperature (and also because the mesh isn't moving!)
    el_pt->ignore_external_geometric_data();
   }
  
  unsigned n_ad_elem=adv_diff_mesh_pt()->nelement();
  for (unsigned e=0;e<n_ad_elem;e++)
   {
    AD_ELEMENT* el_pt = 
     dynamic_cast<AD_ELEMENT*>(adv_diff_mesh_pt()->element_pt(e));
    
    // Set the Peclet number
    el_pt->pe_pt() = &Global_Parameters::Pe;
    
    // Set the Peclet number multiplied by the Strouhal number
    el_pt->pe_st_pt() =&Global_Parameters::Pe;
    
    // We can ignore the external geometric data in the "external"
    // Navier Stokes element when computing the Jacobian matrix
    // because the interaction does not involve spatial gradients of 
    // the velocities (and also because the mesh isn't moving!)
    el_pt->ignore_external_geometric_data();
   }

  //Combine the sub meshes
  add_sub_mesh(Nst_mesh_pt);
  add_sub_mesh(Adv_diff_mesh_pt);
  build_global_mesh();
  
  // Set external elements for the multi-domain solution.
  Multi_domain_functions::
   setup_multi_domain_interactions<NST_ELEMENT,AD_ELEMENT>
   (this,nst_mesh_pt(),adv_diff_mesh_pt());
  
  //Attach the boundary conditions to the mesh
  cout <<"Number of equations: " << assign_eqn_numbers() << endl; 
  
  //Set swirl velocity to zero
  for(unsigned n=0;n<n_node;n++)
   {
    nst_mesh_pt()->node_pt(n)->set_value(2,0.0);
   }
  
   //Plug flow everywhere
  for(unsigned ibound=0;ibound<5;++ibound)
    {
     unsigned num_nod= nst_mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      { 
       //No radial flow on symmetry boundary
       if(ibound==3) 
        {
         nst_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
        }
      
       //No radial flow on outlet
       if(ibound==2) 
        {
         nst_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
        }

       //No flow on side boundary
       if(ibound==1)
        {
         nst_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
         nst_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
        }
       
       //No flow on inlet
       if(ibound==0) 
        {
         nst_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
         nst_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
        }

       //No slip on sphere
       if(ibound==4)
        {
         nst_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
         nst_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
       }
      }
    }

   //Fixed concentration on sphere
   for(unsigned ibound=0;ibound<5;++ibound)
    {
     unsigned num_nod= adv_diff_mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       //No concentration on inlet
       if(ibound==0) 
        {
         adv_diff_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,0.0);
        }

       //Fixed concentration on sphere
       if(ibound==4)
        {
         adv_diff_mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,1.0);
       }
      }
    }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


//=====================================================================
/// Driver
//=====================================================================
int main()
{
 //Multi_domain_functions::Doc_stats=true;
 //Multi_domain_functions::Doc_full_stats=true;

 Global_Parameters::G.resize(3);
 Global_Parameters::G[0] = 0.0;
 Global_Parameters::G[1] = -1.0;
 Global_Parameters::G[2] = 0.0;
 
 // radius and Length of domain
 double radius=5.0;
 double length=10.0;
 Global_Parameters::Sphere_centre_z = 5.0;


 //Create a new ellipse object as the central cylinder
 HalfEllipse* cylinder_pt = 
  new HalfEllipse(Global_Parameters::Sphere_centre_z,0.5,0.5);

 // Create Problem
 FlowAroundHalfCylinderProblem 
  <RefineableQAxisymCrouzeixRaviartBoussinesqElement,
   RefineableQAxisymAdvectionDiffusionBoussinesqElement> 
  problem(cylinder_pt,radius,length);

 //Refine the problem a couple of times
 //problem.refine_uniformly();
 //problem.refine_uniformly();
 

 // Solve adaptively with up to max_adapt rounds of refinement
 unsigned max_adapt=1;

 // Output filename
 char filename[100];
 // Trace file
 //std::ofstream trace("trace.dat");

 //Step up in Reynolds number (and therefore Peclet number)
 for(unsigned i=0;i<2;i++)
  {
   problem.newton_solve(max_adapt);
   
   //Open an output file
   sprintf(filename,"vel_soln_Re%g_Ra%g.dat", Global_Parameters::Re,
           Global_Parameters::Ra);
   //Doc result
   ofstream outfile(filename);
   problem.nst_mesh_pt()->output(outfile,5);
   outfile.close();
   
   sprintf(filename,"conc_soln_Re%g_Ra%g.dat", Global_Parameters::Re,
           Global_Parameters::Ra);
   outfile.open(filename);
   problem.adv_diff_mesh_pt()->output(outfile,5);
   outfile.close();

   Global_Parameters::Ra += 40.0;
  }

 //trace.close();
}

