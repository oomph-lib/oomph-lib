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
// Driver for Solid deformation -- driven by boundary motion which
// is imposed directly (without Lagrange multipliers)

//Oomph-lib includes
#include "generic.h"
#include "solid.h"
#include "constitutive.h"

//The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;


//======Start_of_warped_line===============================================
/// Warped line in 2D space
//=========================================================================
class WarpedLine : public GeomObject
{

public:

 /// Constructor: Specify amplitude of deflection from straight horizontal line
 WarpedLine(const double& ampl) : GeomObject(1,2)
  {
   Ampl=ampl;
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

 ///  Position vector at Lagrangian coordinate zeta 
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   // Position vector
   r[0] = zeta[0]+5.0*Ampl*zeta[0]*(zeta[0]-1.0)*(zeta[0]-0.7);
   r[1] = 1.0+Ampl*0.5*(1.0-cos(2.0*MathematicalConstants::Pi*zeta[0]));
  }
 
 ///  Parametrised position on object: r(zeta). Evaluated at
 /// previous timestep. t=0: current time; t>0: previous
 /// timestep. Forward to steady version
 void position(const unsigned& t, const Vector<double>& zeta,
                       Vector<double>& r) const
  {
   position(zeta,r);
  }

 /// Access to amplitude
 double& ampl() {return Ampl;}
 
private:

 /// Amplitude of perturbation
 double Ampl;

};



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//=======start_namespace==========================================
/// Global parameters
//================================================================
namespace Global_Physical_Variables
{

 /// GeomObject specifying the shape of the boundary: Initially it's flat.
 WarpedLine Boundary_geom_object(0.0);

 /// Poisson's ratio
 double Nu=0.3;

 // Generalised Hookean constitutive equations
 GeneralisedHookean Constitutive_law(&Global_Physical_Variables::Nu);
 
} //end namespace



//=============begin_problem============================================ 
/// Problem class for deformation of elastic block by prescribed
/// boundary motion.
//====================================================================== 
template<class ELEMENT>
class PrescribedBoundaryDisplacementProblem : public Problem
{

public:

 /// Constructor:
 PrescribedBoundaryDisplacementProblem();
 
 /// Update function (empty)
 void actions_after_newton_solve() {}

 /// Update boundary position directly
 void actions_before_newton_solve()
  {

   // Loop over all nodes on top boundary (boundary 2)
   unsigned b=2;
   unsigned n_nod = solid_mesh_pt()->nboundary_node(b);
   for(unsigned i=0;i<n_nod;i++)
    {
     Node* nod_pt= solid_mesh_pt()->boundary_node_pt(b,i);

     // Get boundary coordinate associated with boundary 2
     Vector<double> zeta(1);
     nod_pt->get_coordinates_on_boundary(b,zeta);

     // Get prescribed position from GeomObject
     Vector<double> r(2);
     Global_Physical_Variables::Boundary_geom_object.position(zeta,r);

     // Update position
     nod_pt->x(0)=r[0];
     nod_pt->x(1)=r[1];
    }

  } // end actions_before_newton_solve
 
 /// Access function for the solid mesh
 ElasticRefineableRectangularQuadMesh<ELEMENT>*& solid_mesh_pt() 
  {return Solid_mesh_pt;} 

 /// Actions after adapt: Pin the redundant solid pressures (if any)
 void actions_after_adapt()
  {
   PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
    solid_mesh_pt()->element_pt());
  }
 
 /// Doc the solution
 void doc_solution();

private:

 /// Pointer to solid mesh
 ElasticRefineableRectangularQuadMesh<ELEMENT>* Solid_mesh_pt;

 /// DocInfo object for output
 DocInfo Doc_info;
 
};


//===========start_of_constructor======================================= 
/// Constructor: 
//====================================================================== 
template<class ELEMENT>
PrescribedBoundaryDisplacementProblem<ELEMENT>::PrescribedBoundaryDisplacementProblem() 
{

 // Create the mesh

 // # of elements in x-direction
 unsigned n_x=5;

 // # of elements in y-direction
 unsigned n_y=5;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;

 //Now create the mesh 
 solid_mesh_pt() = new ElasticRefineableRectangularQuadMesh<ELEMENT>(
  n_x,n_y,l_x,l_y);

 // Copy across
 Problem::mesh_pt()=solid_mesh_pt();

 // Set error estimator
 solid_mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 //Assign the physical properties to the elements before any refinement
 //Loop over the elements in the main mesh
 unsigned n_element =solid_mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(i));
   
   // Set the constitutive law
   el_pt->constitutive_law_pt()=&Global_Physical_Variables::Constitutive_law;
  }

 // Refine the mesh uniformly
 solid_mesh_pt()->refine_uniformly();


 
 // Pin nodal positions on all boundaries including the top one 
 for (unsigned b=0;b<4;b++)
  {
   unsigned n_side = solid_mesh_pt()->nboundary_node(b);
   
   //Loop over the nodes
   for(unsigned i=0;i<n_side;i++)
    {
     solid_mesh_pt()->boundary_node_pt(b,i)->pin_position(0);
     solid_mesh_pt()->boundary_node_pt(b,i)->pin_position(1);
    }
  }
 

// Pin the redundant solid pressures (if any)
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
  solid_mesh_pt()->element_pt());
 
 // Setup equation numbering scheme
 cout << "Number of dofs: " << assign_eqn_numbers() << std::endl; 
 
 // Set output directory
 Doc_info.set_directory("RESLT");

} //end of constructor



//==============start_doc===========================================
/// Doc the solution
//==================================================================
template<class ELEMENT>
void PrescribedBoundaryDisplacementProblem<ELEMENT>::doc_solution()
{

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned n_plot = 5; 

 // Output shape of deformed body
 //------------------------------
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 solid_mesh_pt()->output(some_file,n_plot);
 some_file.close();

 // Increment label for output files
 Doc_info.number()++;

} //end doc



//=======start_of_main==================================================
/// Driver code
//======================================================================
int main()
{
 
 //Set up the problem
 PrescribedBoundaryDisplacementProblem<RefineableQPVDElement<2,3> > problem;
 
 // Doc initial domain shape
 problem.doc_solution();

 // Max. number of adaptations per solve
 unsigned max_adapt=1;

 //Parameter incrementation
 unsigned nstep=2; // 16;  
 for(unsigned i=0;i<nstep;i++)
  {
   // Increment imposed boundary displacement
   Global_Physical_Variables::Boundary_geom_object.ampl()+=0.0125;

   // Solve the problem with Newton's method, allowing
   // up to max_adapt mesh adaptations after every solve.
   problem.newton_solve(max_adapt);
   
   // Doc solution
   problem.doc_solution();

   // For maximum stability: Reset the current nodal positions to be
   // the "stress-free" ones -- this assignment means that the
   // parameter study no longer corresponds to a physical experiment
   // but is what we'd do if we wanted to use the solid solve
   // to update a fluid mesh in an FSI problem, say.
   problem.solid_mesh_pt()->set_lagrangian_nodal_coordinates();
   
  }
 
} //end of main








