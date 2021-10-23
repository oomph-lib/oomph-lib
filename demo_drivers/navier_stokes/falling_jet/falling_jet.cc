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
/// Driver for a jet falling from a circular orifice under the action
///of gravity

//Generic routines
#include "generic.h"
#include "navier_stokes.h"
#include "fluid_interface.h"
#include "solid.h"
#include "constitutive.h"

// The mesh
#include "meshes/quarter_tube_mesh.h"

using namespace std;

using namespace oomph;


//=start of mesh class===============================================
template<class ELEMENT,class INTERFACE_ELEMENT>
class ElasticQuarterTubeMesh : public RefineableQuarterTubeMesh<ELEMENT>,
public SolidMesh
{
 /// Vector of pointers to element in the fluid layer
 Vector <GeneralisedElement *> Bulk_element_pt;

 /// Vector of pointers to interface elements
 Vector<FiniteElement *> Interface_element_pt;

  /// Vector of pointers to interface elements
 Vector<FiniteElement *> Interface_edge_element_pt;


public:

  /// Access functions for pointers to interface elements
 FiniteElement* &interface_element_pt(const unsigned long &i) 
  {return Interface_element_pt[i];}

 /// Number of elements on interface
 unsigned long ninterface_element() const {return Interface_element_pt.size();}
 
 /// Access functions for pointers to interface elements
 FiniteElement* &interface_edge_element_pt(const unsigned long &i) 
  {return Interface_edge_element_pt[i];}

 /// Number of elements on interface
 unsigned long ninterface_edge_element() const 
  {return Interface_edge_element_pt.size();}
 
 Vector<GeneralisedElement*> &bulk_element_pt()
  {return Bulk_element_pt;}

 ///Access functions for pointers to elements in bulk
 GeneralisedElement* &bulk_element_pt(const unsigned long &i) 
  {return Bulk_element_pt[i];}

 ///Number of elements in bulk 
 unsigned long nbulk() const {return Bulk_element_pt.size();}

 ///Constructor
 ElasticQuarterTubeMesh(GeomObject* wall_pt,
                           const Vector<double>& xi_lo,
                           const double& fract_mid,
                           const Vector<double>& xi_hi,
                           const unsigned& nlayer,
                           TimeStepper* time_stepper_pt=
                           &Mesh::Default_TimeStepper):
  QuarterTubeMesh<ELEMENT>(wall_pt,xi_lo,fract_mid,xi_hi,
                           nlayer,time_stepper_pt),
  
  RefineableQuarterTubeMesh<ELEMENT>(wall_pt,xi_lo,fract_mid,xi_hi,
                                      nlayer,time_stepper_pt),
  SolidMesh()
  {
   //Perform a singleuniform refinement
   this->refine_uniformly();
   //this->refine_uniformly();

   //Update the boundary information
   this->setup_boundary_element_info();

   //Add the elements to the bulk element
   Bulk_element_pt = this->Element_pt;

   //Let's create the interface elements
   unsigned b=3;
   unsigned n_element = this->nboundary_element(b);
   unsigned np = this->boundary_element_pt(b,0)->nnode_1d();
   for(unsigned e=0;e<n_element;e++)
    {
     INTERFACE_ELEMENT *interface_element_pt =
      new INTERFACE_ELEMENT(this->boundary_element_pt(b,e),
                            this->face_index_at_boundary(b,e));

     //Push it back onto the stack
     this->Element_pt.push_back(interface_element_pt);

     //Add it to the interface
     this->Interface_element_pt.push_back(interface_element_pt);

     //Now check whether the corners are on boundary 0 (the bottom)
     //or boundary 4 (the top)
     for(unsigned b=0;b<5;b+=4)
      {
       std::vector<bool> corner(4);
       corner[0] = interface_element_pt->node_pt(0)->is_on_boundary(b);
       corner[1] = interface_element_pt->node_pt(np-1)->is_on_boundary(b);
       corner[2] = interface_element_pt->node_pt(np*(np-1))->is_on_boundary(b);
       corner[3] = interface_element_pt->node_pt(np*np-1)->is_on_boundary(b);

       //Work out which edge is on the boundary
       bool bottom = corner[0] && corner[1];
       bool top = corner[2] && corner[3];
       bool left = corner[0] && corner[2];
       bool right = corner[1] && corner[3];
       
       //Create the edge elements, don't check for repetition
       //two edges of an element *could* be on the line
       if(bottom)
        {
         Interface_edge_element_pt.push_back(interface_element_pt
                                             ->make_bounding_element(-2));
        }
       
       if(top)
        {
         Interface_edge_element_pt.push_back(interface_element_pt
                                             ->make_bounding_element(2));
        }
       
       if(left)
        {
         Interface_edge_element_pt.push_back(interface_element_pt
                                             ->make_bounding_element(-1));
        }
       
       if(right)
        {
         Interface_edge_element_pt.push_back(interface_element_pt
                                             ->make_bounding_element(1));
        }
      } //End of loop over the top and the bottom
    }

   //Add all the edge elements to the standard element list
   for(Vector<FiniteElement*>::iterator it
        = Interface_edge_element_pt.begin();
       it != Interface_edge_element_pt.end(); ++it)
    {
     this->Element_pt.push_back(*it);
    }
   
   //Now let's print the things shall we
   /*{
    ofstream surface("surface.dat");
    for(Vector<FiniteElement*>::iterator it = Interface_element_pt.begin();
        it!=Interface_element_pt.end();++it)
     {
      (*it)->output(surface,5);
     }
    surface.close();
   }

   {
    ofstream line("line.dat");
    for(Vector<FiniteElement*>::iterator it = 
         Interface_edge_element_pt.begin();
        it!=Interface_edge_element_pt.end();++it)
     {
      (*it)->output(line,5);
     }
    line.close();
    }*/


   //Make the current configuration the undeformed one
   this->set_lagrangian_nodal_coordinates();
  }

 
};




//=start_of_namespace================================================
/// Namespace for physical parameters
//===================================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re=0.0;

 /// Capillary number
 double Ca = 1.0;


 /// Stokes number
 double St = 0.0;

 /// Gravity direction
 Vector<double> G(3);

 /// Pseudo-solid Poisson ratio
 double Nu=0.1;

} // end_of_namespace


//=start_of_problem_class=============================================
/// Entry flow problem in quarter tube domain
//====================================================================
template<class ELEMENT>
class EntryFlowProblem : public Problem
{
 //Constitutive law used to determine the mesh deformation
 ConstitutiveLaw *Constitutive_law_pt;

 //Pointer to external pressure data
 Data *Pext_pt;

public:

 /// Constructor: Pass DocInfo object and target errors
 EntryFlowProblem(DocInfo& doc_info, const double& min_error_target,
                  const double& max_error_target);

 /// Destructor (empty)
 ~EntryFlowProblem() {}

 /// Doc the solution after solve
 void actions_after_newton_solve() 
  {
   // Increment label for output files
   Doc_info.number()++;
   // Doc solution after solving
   doc_solution();
  }

 /// Update the problem specs before solve 
 void actions_before_newton_solve() { } 

 /// After adaptation: Pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Pin redudant pressure dofs
   RefineableNavierStokesEquations<3>::
    pin_redundant_nodal_pressures(mesh_pt()->bulk_element_pt());
  } 

 /// Doc the solution
 void doc_solution();

 /// Overload generic access function by one that returns
 /// a pointer to the specific  mesh
 ElasticQuarterTubeMesh<ELEMENT,ElasticSurfaceFluidInterfaceElement<ELEMENT> >* 
 mesh_pt() 
  {
   return dynamic_cast<ElasticQuarterTubeMesh<ELEMENT,
    ElasticSurfaceFluidInterfaceElement<ELEMENT> >*>(Problem::mesh_pt());
  }

private:

 /// Exponent for bluntness of velocity profile
 int Alpha;
 
 /// Doc info object
 DocInfo Doc_info;

}; // end_of_problem_class




//=start_of_constructor===================================================
/// Constructor: Pass DocInfo object and error targets
//========================================================================
template<class ELEMENT>
EntryFlowProblem<ELEMENT>::EntryFlowProblem(DocInfo& doc_info,
                                            const double& min_error_target,
                                            const double& max_error_target) 
 : Doc_info(doc_info)
{
 //Set the external pressure
 Pext_pt = new Data(1);
 Pext_pt->set_value(0,-1.0);
 add_global_data(Pext_pt);
 Pext_pt->pin(0);

 // Setup mesh:
 //------------

 // Create geometric objects: Elliptical tube with half axes = radius = 1.0
 double radius=1.0;
 GeomObject* Wall_pt=new EllipticalTube(radius,radius);

 // Boundaries on object
 Vector<double> xi_lo(2);
 // height of inflow
 xi_lo[0]=0.0;
 // start of Wall_pt
 xi_lo[1]=0.0;

 Vector<double> xi_hi(2);
 // height of outflow
 xi_hi[0]=7.0;
 // end of Wall_pt
 xi_hi[1]=2.0*atan(1.0);

 // # of layers
 unsigned nlayer=3;

 //Radial divider is located half-way along the circumference
 double frac_mid=0.5;

 // Build and assign mesh
 Problem::mesh_pt()=
  new ElasticQuarterTubeMesh<ELEMENT,
  ElasticSurfaceFluidInterfaceElement<ELEMENT> >(
   Wall_pt,xi_lo,frac_mid,xi_hi,nlayer);
 

 // Set error estimator 
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
 
 // Error targets for adaptive refinement
 mesh_pt()->max_permitted_error()=max_error_target; 
 mesh_pt()->min_permitted_error()=min_error_target; 



 //Doc the boundaries
 /*ofstream some_file;
 char filename[100];
 sprintf(filename,"boundaries.dat");
 some_file.open(filename);
 mesh_pt()->output_boundaries(some_file);
 some_file.close();*/
 
 
 // Set the boundary conditions for this problem: All nodal values are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. 

 //We need to find whether the pressure is stored at the nodes or not
 //Wrap this up?
 //If it is not the Lagrange index is 3
 unsigned lagrange_index=3;
 //If the pressure is stored at a node, the lagrangian index is 4
 if(dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(0))
    ->p_nodal_index_nst() >= 0) {lagrange_index=4;}
 
 //Boundary 4 is the inflow, pin x, y and z
 {
  unsigned b=4;
  unsigned n_node = mesh_pt()->nboundary_node(b);
  for(unsigned n=0;n<n_node;n++)
   {
    for(unsigned i=0;i<3;i++) {mesh_pt()->boundary_node_pt(b,n)->pin(i);}
    mesh_pt()->boundary_node_pt(b,n)->pin_position(0);
    mesh_pt()->boundary_node_pt(b,n)->pin_position(1);
    mesh_pt()->boundary_node_pt(b,n)->pin_position(2);
    if(mesh_pt()->boundary_node_pt(b,n)->nvalue() > lagrange_index)
     {
      mesh_pt()->boundary_node_pt(b,n)->pin(lagrange_index);
     }
   }
 }

 //Boundary three is the wall (free)

 //Boundary zero is the outflow (free)/
 //Need to pin the upward positions
 {
  unsigned b=0;
  unsigned n_node = mesh_pt()->nboundary_node(b);
  for(unsigned n=0;n<n_node;n++)
   {
    mesh_pt()->boundary_node_pt(b,n)->pin_position(2);
   }
 }


 //Boundary 1 is the vertical symmetry boundary
 //only allow flow in the y-direction.
 {
  unsigned b=1;
  unsigned n_node = mesh_pt()->nboundary_node(b);
  for(unsigned n=0;n<n_node;n++)
   {
    mesh_pt()->boundary_node_pt(b,n)->pin(0);
    mesh_pt()->boundary_node_pt(b,n)->pin(2);
    mesh_pt()->boundary_node_pt(b,n)->pin_position(0);
   }
 }

 //Boundary 2 is the horizontal symmetry boundary
 //only allow flow in the x-direction.
 {
  unsigned b=2;
  unsigned n_node = mesh_pt()->nboundary_node(b);
  for(unsigned n=0;n<n_node;n++)
   {
    mesh_pt()->boundary_node_pt(b,n)->pin(1);
    mesh_pt()->boundary_node_pt(b,n)->pin(2);
    mesh_pt()->boundary_node_pt(b,n)->pin_position(1);
   }
 }



 //Set the constituive law
 Constitutive_law_pt = new GeneralisedHookean(&Global_Physical_Variables::Nu);


 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 unsigned n_bulk = mesh_pt()->nbulk();
 for(unsigned e=0;e<n_bulk;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //Set the constitutive law
   el_pt->constitutive_law_pt() = Constitutive_law_pt;

   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
   el_pt->re_invfr_pt() = &Global_Physical_Variables::St;
   el_pt->g_pt() = &Global_Physical_Variables::G;
  }

 //Loop over the interface elements and set the capillary number
 unsigned n_interface = mesh_pt()->ninterface_element();
 for(unsigned e=0;e<n_interface;e++)
  {
   ElasticSurfaceFluidInterfaceElement<ELEMENT>* el_pt
    = dynamic_cast<ElasticSurfaceFluidInterfaceElement<ELEMENT>*>(
     mesh_pt()->interface_element_pt(e));
   
   //set the capillary number
   el_pt->ca_pt() =  &Global_Physical_Variables::Ca;
   //Set the external pressure
   el_pt->set_external_pressure_data(Pext_pt);
  }

 //Loop over the interface edge element
 unsigned n_interface_edge = mesh_pt()->ninterface_edge_element();
 for(unsigned e=0;e<n_interface_edge;e++)
  {
   ElasticLineFluidInterfaceBoundingElement<ELEMENT>* el_pt
    = dynamic_cast<ElasticLineFluidInterfaceBoundingElement<ELEMENT>*>(
     mesh_pt()->interface_edge_element_pt(e));

   //Set the capillary number
   el_pt->ca_pt() = &Global_Physical_Variables::Ca;

  }


 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<3>::
  pin_redundant_nodal_pressures(mesh_pt()->bulk_element_pt());

 // Set Plug flow as initial for solution
 // Inflow will be overwritten in actions_before_solve()
 unsigned n_nod=mesh_pt()->nnode();
 for (unsigned j=0;j<n_nod;j++)
  {
   SolidNode* node_pt=mesh_pt()->node_pt(j);
   // Recover coordinates
   //double x=node_pt->x(0);
   //double y=node_pt->x(1);
   //double r=sqrt(x*x+y*y );  
   
   // Poiseuille flow
   node_pt->set_value(0,0.0);
   node_pt->set_value(1,0.0);
   //node_pt->set_value(2,(1.0-r*r));
   node_pt->set_value(2,-1.0);
  }

 // Set the exponent for bluntness: Alpha=2 --> Poisseuille; anything
 // larger makes the inflow blunter
 Alpha=20;

 //Attach the boundary conditions to the mesh
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end_of_constructor


//=start_of_actions_before_solve==========================================
/// Set the inflow boundary conditions
//========================================================================
/*template<class ELEMENT>
void EntryFlowProblem<ELEMENT>::actions_before_solve()
{

 // (Re-)assign velocity profile at inflow values
 //--------------------------------------------

 // Setup bluntish parallel inflow on boundary 0:
 unsigned ibound=0; 
 unsigned num_nod= mesh_pt()->nboundary_node(ibound); 
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   // Recover coordinates
   double x=mesh_pt()->boundary_node_pt(ibound,inod)->x(0);
   double y=mesh_pt()->boundary_node_pt(ibound,inod)->x(1);
   double r=sqrt(x*x+y*y);  
   
   // Bluntish profile for axial velocity (component 2)
   mesh_pt()->boundary_node_pt(ibound,inod)->
    set_value(2,(1.0-pow(r,Alpha)));
  }

} // end_of_actions_before_solve
*/

//=start_of_doc_solution==================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void EntryFlowProblem<ELEMENT>::doc_solution()
{ 
 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 //Only output the bulk elements, because the surface elements are
 //assigned in a different order on different machines
 unsigned n_bulk = mesh_pt()->nbulk(); 
 for(unsigned e=0;e<n_bulk;e++)
  {
   dynamic_cast<FiniteElement*>(mesh_pt()->bulk_element_pt(e))
    ->output(some_file,npts);
  }
 some_file.close();

} // end_of_doc_solution

 


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//=start_of_main=======================================================
/// Driver for 3D entry flow into a quarter tube. If there are
/// any command line arguments, we regard this as a validation run
/// and perform only a single adaptation
//=====================================================================
int main(int argc, char* argv[]) 
{
 //Set the gravity direction
 Global_Physical_Variables::G[0] = 0.0;
 Global_Physical_Variables::G[1] = 0.0;
 Global_Physical_Variables::G[2] = -1.0;


 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Allow (up to) five rounds of fully automatic adapation in response to 
 //-----------------------------------------------------------------------
 // error estimate
 //---------------
 //unsigned max_adapt;
 double max_error_target,min_error_target;

 // Set max number of adaptations in black-box Newton solver and
 // error targets for adaptation
 if (CommandLineArgs::Argc==1)
  {
   // Up to five adaptations
   //max_adapt=5;

   // Error targets for adaptive refinement
   max_error_target=0.005;
   min_error_target=0.0005;
  } 
 // Validation run: Only one adaptation. Relax error targets
 // to ensure that not all elements are refined so we're getting
 // some hanging nodes.
 else
  {
   // Validation run: Just one round of adaptation
   //max_adapt=1;
   
   // Error targets for adaptive refinement
   max_error_target=0.02;
   min_error_target=0.002;
  }
 // end max_adapt setup
 

 // Set up doc info
 DocInfo doc_info;
 
 // Do Taylor-Hood elements
 //------------------------
 {
  // Set output directory
  doc_info.set_directory("RESLT_TH");
  
  // Step number
  doc_info.number()=0;
  
  // Build problem
  EntryFlowProblem<
   RefineablePseudoSolidNodeUpdateElement<
   RefineableQTaylorHoodElement<3>,
   RefineableQPVDElement<3,3> > >
   problem(doc_info,min_error_target,max_error_target);
  
  cout << " Doing Taylor-Hood elements " << std::endl;
  
  // Doc solution before solving
  problem.doc_solution();
  
  // Solve the problem 
  for(unsigned i=0;i<2;i++)
   {
    problem.newton_solve();//max_adapt);
    Global_Physical_Variables::St += 0.5;
   }
}


 // Do Crouzeix-Raviart elements
 //-----------------------------
 {
  // Set output directory
  doc_info.set_directory("RESLT_CR");
  
  // Step number
  doc_info.number()=0;

  //Reset the Stokes number to zero
  Global_Physical_Variables::St = 0.0;
  
  // Build problem
  EntryFlowProblem<
   RefineablePseudoSolidNodeUpdateElement<
   RefineableQCrouzeixRaviartElement<3>,
   RefineableQPVDElement<3,3> > >
   problem(doc_info,min_error_target,max_error_target);
  
  cout << " Doing Crouzeix-Raviart elements " << std::endl;
  
  // Doc solution before solving
  problem.doc_solution();
  
  // Solve the problem 
  for(unsigned i=0;i<2;i++)
   {
    problem.newton_solve();//max_adapt);
    Global_Physical_Variables::St += 0.5;
   }

  }

} // end_of_main


