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
///Driver for a 3D navier stokes entry flow problem in quarter tube domain

//Generic routines
#include "generic.h"
#include "navier_stokes.h"
#include "fluid_interface.h"

// The mesh
#include "meshes/quarter_tube_mesh.h"

using namespace std;

using namespace oomph;

//=start_of_namespace================================================
/// Namespace for physical parameters
//===================================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re= 0.0;

 ///Capillary number
 double Ca = 1.0;

 ///Rotation rate
 double Omega = 0.0;

 //Inverse Froude number
 double ReInvFr = 0.0;

 Vector<double> G(3); 

 double Angle;

} // end_of_namespace

namespace WallFunction
{
 void normal(const Vector<double> &x, Vector<double> &normal)
 {
  //Calculate theta
  double theta = std::atan2(x[1],x[0]);
  //Return the normal
  normal[0] = cos(theta);
  normal[1] = sin(theta);
  normal[2] = 0.0;
 }
}


//=start of mesh class===============================================
template<class ELEMENT,class INTERFACE_ELEMENT>
class AxialSpineQuarterTubeMesh : public RefineableQuarterTubeMesh<ELEMENT>,
public SpineMesh
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
 AxialSpineQuarterTubeMesh(GeomObject* wall_pt,
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
  SpineMesh()
  {
   //Perform a singleuniform refinement
   this->refine_uniformly();
   this->refine_uniformly();

   unsigned n_layer=4;

   //Now worry about the spines the trick is to loop over each
   //layer in turn. We determine the layer based on the z coordinate
   //of the first node in each element.
   Bulk_element_pt = Element_pt;

   //Loop over all the elements
   unsigned n_element = this->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     //If the element is on the lowest layer
     if(std::abs(finite_element_pt(e)->node_pt(0)->x(2)) < 1.0e-15)
      {
       //Find the number of 1d nodes
       unsigned n_p = finite_element_pt(e)->nnode_1d();
       //Loop over the lowest layer
       for(unsigned n=0;n<(n_p*n_p);n++)
        {
         //Get pointer to the node
         SpineNode* nod_pt = element_node_pt(e,n);
         //If the node has no spine, create one
         if(nod_pt->spine_pt() == 0)
          {
           //Create a new spine of length 2
           Spine* new_spine_pt = new Spine(2.0);
           //Add it to the mesh
           Spine_pt.push_back(new_spine_pt);
           nod_pt->spine_pt() = new_spine_pt;
           //The fraction is zeo
           nod_pt->fraction() = 0.0;
           //Pointer to the mesh the implements the update
           nod_pt->spine_mesh_pt() = this;
          }
         
         //Loop up the spine and set the other nodes in the element
         for(unsigned m=1;m<n_p;m++)
          {
           SpineNode* nod2_pt = element_node_pt(e,(n_p*n_p)*m + n);
           nod2_pt->spine_pt() = nod_pt->spine_pt();
           nod2_pt->fraction() = nod2_pt->x(2)/2.0;
           nod2_pt->spine_mesh_pt() = this;
          }
        }
      }
    }

   double layer_width = 2.0/(double)n_layer;
   FiniteElement *interface_element_pt=0;

   //Loop over the remaining layers
   for(unsigned l=1;l<n_layer;l++)
    {
     double boundary = layer_width*l;
     //Treat the second layer
     for(unsigned e=0;e<n_element;e++)
      {
       //If the element is on the second layer
       if(std::abs(finite_element_pt(e)->node_pt(0)->x(2)-boundary) < 1.0e-15)
        {
         //If on the last layer
         if(l==n_layer-1)
          {
           //We can build the interface elements at the max
           //of the element (face 3)
           interface_element_pt
            = new INTERFACE_ELEMENT(finite_element_pt(e),3);
           
           //push it back onto the stack
           Element_pt.push_back(interface_element_pt);
           
           //Add it to the stack of interface elements
           Interface_element_pt.push_back(interface_element_pt);
          }
         
         //Find the number of 1d nodes
         unsigned n_p = finite_element_pt(e)->nnode_1d();
         //Loop over the lowest layer
         for(unsigned n=0;n<(n_p*n_p);n++)
          {
           //Get pointer to the spine of the nodes
           Spine* spine_pt = element_node_pt(e,n)->spine_pt();
           
           //Loop up the spine and set the other nodes in the element
           for(unsigned m=1;m<n_p;m++)
            {
             SpineNode* nod_pt = element_node_pt(e,(n_p*n_p)*m + n);
             nod_pt->spine_pt() = spine_pt;
             nod_pt->fraction() = nod_pt->x(2)/2.0;
             nod_pt->spine_mesh_pt() = this;
            }
          }

         if(l==n_layer-1)
          {
           //Is the interface adjacent to the wall
           unsigned n_p = interface_element_pt->nnode_1d();
           //Which boundary
           std::vector<int> adjacent(4,true);
           for(unsigned n=0;n<n_p;n++)
            {
             //Bottom
             adjacent[0] &= 
              interface_element_pt->node_pt(n)->is_on_boundary(3);
             //Right
             adjacent[1] &= 
              interface_element_pt->node_pt(n*n_p + n_p-1)->is_on_boundary(3);
             //Top
             adjacent[2] &= 
              interface_element_pt->node_pt(n_p*(n_p-1) + n)->is_on_boundary(3);
             //Left
             adjacent[3] &=
              interface_element_pt->node_pt(n*n_p)->is_on_boundary(3);
            }

           
           //Now if we have a match make the element
           FluidInterfaceBoundingElement* interface_edge_element_pt=0;
           if(adjacent[0])
            {
             interface_edge_element_pt
              = dynamic_cast<INTERFACE_ELEMENT*>(interface_element_pt)
              ->make_bounding_element(-2);
            }
           else if(adjacent[1])
            {
             interface_edge_element_pt
              = dynamic_cast<INTERFACE_ELEMENT*>(interface_element_pt)
              ->make_bounding_element(1);
            }
           else if(adjacent[2])
            {
             interface_edge_element_pt
              = dynamic_cast<INTERFACE_ELEMENT*>(interface_element_pt)
              ->make_bounding_element(2);
           }
          else if(adjacent[3])
           {
            interface_edge_element_pt
             = dynamic_cast<INTERFACE_ELEMENT*>(interface_element_pt)
             ->make_bounding_element(-1);
            
           }
          
          if(interface_edge_element_pt!=0)
           {
            Element_pt.push_back(interface_edge_element_pt);
            Interface_edge_element_pt.push_back(interface_edge_element_pt);
            dynamic_cast<FluidInterfaceBoundingElement*>
             (interface_edge_element_pt)->wall_unit_normal_fct_pt() 
             = WallFunction::normal;
            //Set the contact angle
            dynamic_cast<FluidInterfaceBoundingElement*>
             (interface_edge_element_pt)->
             set_contact_angle(&Global_Physical_Variables::Angle);
            //and pass the capillary number
            dynamic_cast<FluidInterfaceBoundingElement*>
             (interface_edge_element_pt)->ca_pt() = 
             &Global_Physical_Variables::Ca;
           }
          } //End of final layer stuff
        }
      }
    }
   
   /*  //Treat the third layer
   for(unsigned e=0;e<n_element;e++)
    {
     //If the element is on the third layer
     if(std::abs(finite_element_pt(e)->node_pt(0)->x(2)-1.0) < 1.0e-15)
      {
       //Find the number of 1d nodes
       unsigned n_p = finite_element_pt(e)->nnode_1d();
       //Loop over the lowest layer
       for(unsigned n=0;n<(n_p*n_p);n++)
        {
         //Get pointer to the spine of the nodes
         Spine* spine_pt = element_node_pt(e,n)->spine_pt();
         
         //Loop up the spine and set the other nodes in the element
         for(unsigned m=1;m<n_p;m++)
          {
           SpineNode* nod_pt = element_node_pt(e,(n_p*n_p)*m + n);
           nod_pt->spine_pt() = spine_pt;
           nod_pt->fraction() = nod_pt->x(2)/2.0;
           nod_pt->spine_mesh_pt() = this;
          }
        }
      }
      }*/

   /*  //Treat the final layer
   for(unsigned e=0;e<n_element;e++)
    {
     //If the element is on the fourth layer
     if(std::abs(finite_element_pt(e)->node_pt(0)->x(2)-1.5) < 1.0e-15)
      {

       //We can build the interface elements at the max
       //of the element
       FiniteElement *interface_element_pt
        = new INTERFACE_ELEMENT(finite_element_pt(e),2,1);
       
       //push it back onto the stack
       Element_pt.push_back(interface_element_pt);
       
       //Add it to the stack of interface elements
       Interface_element_pt.push_back(interface_element_pt);
       
       //Find the number of 1d nodes
       unsigned n_p = finite_element_pt(e)->nnode_1d();
       //Loop over the lowest layer
       for(unsigned n=0;n<(n_p*n_p);n++)
        {
         //Get pointer to the spine of the nodes
         Spine* spine_pt = element_node_pt(e,n)->spine_pt();
         
         //Loop up the spine and set the other nodes in the element
         for(unsigned m=1;m<n_p;m++)
          {
           SpineNode* nod_pt = element_node_pt(e,(n_p*n_p)*m + n);
           nod_pt->spine_pt() = spine_pt;
           nod_pt->fraction() = nod_pt->x(2)/2.0;
           nod_pt->spine_mesh_pt() = this;
          }
        }
      }
      } */

  }


 ///Update nodal positions in response to spine changes
 virtual void spine_node_update(SpineNode* spine_node_pt)
  {
   //Get fraction along the spine
   double W = spine_node_pt->fraction();
   //Get spine height
   double H = spine_node_pt->h();
   //Set the value of z
   spine_node_pt->x(2) = this->Xi_lo[0] + W*H;
  }

 
};

//=start_of_problem_class=============================================
/// Entry flow problem in quarter tube domain
//====================================================================
template<class ELEMENT>
class FreeSurfaceRotationProblem : public Problem
{

public:

 /// Constructor: Pass DocInfo object and target errors
 FreeSurfaceRotationProblem(DocInfo& doc_info, const double& min_error_target,
                  const double& max_error_target);

 /// Destructor (empty)
 ~FreeSurfaceRotationProblem() {}

 /// Doc the solution after solve
 void actions_after_newton_solve() 
  {
   // Doc solution after solving
   doc_solution();

   // Increment label for output files
   Doc_info.number()++;
  }

 /// Update the problem specs before solve 
 void actions_before_newton_solve() {}

 /// Spine heights/lengths are unknowns in the problem so their
 /// values get corrected during each Newton step. However,
 /// changing their value does not automatically change the
 /// nodal positions, so we need to update all of them
 void actions_before_newton_convergence_check()
  {
   mesh_pt()->node_update();
  }


 /// Doc the solution
 void doc_solution();

 /// Overload generic access function by one that returns
 /// a pointer to the specific  mesh
 AxialSpineQuarterTubeMesh<ELEMENT,SpineSurfaceFluidInterfaceElement<ELEMENT> >* mesh_pt() 
  {
   return 
    dynamic_cast<AxialSpineQuarterTubeMesh<ELEMENT,
    SpineSurfaceFluidInterfaceElement<ELEMENT> >*>
    (Problem::mesh_pt());
  }

private:
 
 /// Doc info object
 DocInfo Doc_info;

}; // end_of_problem_class




//=start_of_constructor===================================================
/// Constructor: Pass DocInfo object and error targets
//========================================================================
template<class ELEMENT>
FreeSurfaceRotationProblem<ELEMENT>::FreeSurfaceRotationProblem(DocInfo& doc_info,
                                            const double& min_error_target,
                                            const double& max_error_target) 
 : Doc_info(doc_info)
{ 
 Data *Pext_pt = new Data(1);
 Pext_pt->pin(0);
 Pext_pt->set_value(0,0.0);

 //Set the contact angle
 const double pi = MathematicalConstants::Pi;
 Global_Physical_Variables::Angle = 0.5*pi;

 Newton_solver_tolerance = 2.0e-5;
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
 xi_hi[0]=2.0;
 // end of Wall_pt
 xi_hi[1]=0.5*MathematicalConstants::Pi;

 // # of layers
 unsigned nlayer=1;

 //Radial divider is located half-way along the circumference
 double frac_mid=0.5;

 // Build and assign mesh
 Problem::mesh_pt()=
  new AxialSpineQuarterTubeMesh<ELEMENT,SpineSurfaceFluidInterfaceElement<ELEMENT> >
  (Wall_pt,xi_lo,frac_mid,xi_hi,nlayer);
 

 // Set error estimator 
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
 
 // Error targets for adaptive refinement
 mesh_pt()->max_permitted_error()=max_error_target; 
 mesh_pt()->min_permitted_error()=min_error_target; 

 //Doc the boundaries
 ofstream some_file;
 char filename[100];
 sprintf(filename,"boundaries.dat");
 some_file.open(filename);
 mesh_pt()->output_boundaries(some_file);
 some_file.close();
 
 //Set the boundary conditions

 //Boundary 0 is the bottom of the domain (pinned)
 {
  unsigned b = 0;
  unsigned n_node = mesh_pt()->nboundary_node(b);
  for(unsigned n=0;n<n_node;n++)
   {
    for(unsigned i=0;i<3;i++) {mesh_pt()->boundary_node_pt(b,n)->pin(i);}
   }
 }
 
 //Boundary 1 is the boundary x is zero, so pin the x-velocity
 {
  unsigned b = 1;
  unsigned n_node = mesh_pt()->nboundary_node(b);
  for(unsigned n=0;n<n_node;n++)
   {
    mesh_pt()->boundary_node_pt(b,n)->pin(0);
   }
 }


 //Boundary 2 is the boundary y=0, so pin the y-velocity
 //otherwise free
 {
  unsigned b = 2;
  unsigned n_node = mesh_pt()->nboundary_node(b);
  for(unsigned n=0;n<n_node;n++)
   {
    mesh_pt()->boundary_node_pt(b,n)->pin(1);
   }
 }

 //Boundary 3 is the wall, so pinned in all coordinates
 {
  unsigned b = 3;
  unsigned n_node = mesh_pt()->nboundary_node(b);
  for(unsigned n=0;n<n_node;n++)
   {
    for(unsigned i=0;i<3;i++)
     {
      mesh_pt()->boundary_node_pt(b,n)->pin(i);
     }
    if(n==0)
     {
      //Pin one spine heights on the wall
      static_cast<SpineNode*>(mesh_pt()->boundary_node_pt(b,n))->
       spine_pt()->spine_height_pt()->pin(0);
     }
   }
 }
  

 //Boundary 4 must be the top, so it's traction free


 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 unsigned n_element = mesh_pt()->nbulk();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;

   //Set the gravitational force
   el_pt->re_invfr_pt() = &Global_Physical_Variables::ReInvFr; 
   el_pt->g_pt() = &Global_Physical_Variables::G;
  }

 //Loop over the interface elements and set the capillary number
 unsigned n_interface = mesh_pt()->ninterface_element();
 for(unsigned e=0;e<n_interface;e++)
  {
   SpineSurfaceFluidInterfaceElement<ELEMENT>* el_pt
    = dynamic_cast<SpineSurfaceFluidInterfaceElement<ELEMENT>*>(
     mesh_pt()->interface_element_pt(e));
   
   //set the capillary number
   el_pt->ca_pt() =  &Global_Physical_Variables::Ca;
   //Set the 
   el_pt->set_external_pressure_data(Pext_pt);
  }

 
 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<3>::
  pin_redundant_nodal_pressures(mesh_pt()->bulk_element_pt());

 // Set solid body rotation as initial solution
 unsigned n_nod=mesh_pt()->nnode();
 for (unsigned j=0;j<n_nod;j++)
  {
   using namespace Global_Physical_Variables;
   Node* node_pt=mesh_pt()->node_pt(j);
   // Recover coordinates
   double x=node_pt->x(0);
   double y=node_pt->x(1);
   double r=sqrt(x*x+y*y);  
   double theta = atan2(y,x);
   
   // Solid body rotation
   node_pt->set_value(0,-Omega*r*sin(theta));
   node_pt->set_value(1,Omega*r*cos(theta));
   node_pt->set_value(2,0.0);
  }

 //Attach the boundary conditions to the mesh
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end_of_constructor



//=start_of_doc_solution==================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void FreeSurfaceRotationProblem<ELEMENT>::doc_solution()
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
 mesh_pt()->output(some_file,npts);
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

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Allow (up to) five rounds of fully automatic adapation in response to 
 //-----------------------------------------------------------------------
 // error estimate
 //---------------
 //unsigned max_adapt;
 double max_error_target,min_error_target;

 using namespace Global_Physical_Variables;

 G[0] = 0.0; G[1] = 0.0; G[2] = -1.0;


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

 const double pi = MathematicalConstants::Pi;
 // Do Taylor-Hood elements
 //------------------------
 {
  // Set output directory
  doc_info.set_directory("RESLT_TH");
  
  // Step number
  doc_info.number()=0;
  
  // Build problem
  FreeSurfaceRotationProblem<SpineElement<RefineableQTaylorHoodElement<3> > >
   problem(doc_info,min_error_target,max_error_target);
  
  cout << " Doing Taylor-Hood elements " << std::endl;


  // Doc solution after solving
  problem.doc_solution();
  
  // Solve the problem  (DO NOT ADAPT)
  problem.newton_solve();

  
  for(unsigned i=0;i<5;i++)
   {
    Angle -= 0.05*pi;
    problem.newton_solve();
   }
 }


 // Do Crouzeix-Raviart elements
 //-----------------------------
 {
  // Set output directory
  doc_info.set_directory("RESLT_CR");
  
  // Step number
  doc_info.number()=0;
  
  // Build problem
  FreeSurfaceRotationProblem<
   SpineElement<RefineableQCrouzeixRaviartElement<3> > >
   problem(doc_info,min_error_target,max_error_target);
  
  cout << " Doing Crouzeix-Raviart elements " << std::endl;
  
  // Solve the problem  DO NOT ADAPT
  problem.newton_solve();
  }

} // end_of_main


