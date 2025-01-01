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
// Driver for Solid deformation -- driven by boundary motion which
// is imposed via Lagrange multipliers


//Oomph-lib includes
#include "generic.h"
#include "solid.h"
#include "constitutive.h"

//The mesh
#include "meshes/simple_cubic_mesh.h"


using namespace std;

using namespace oomph;


// #define NOREFINE

//================================================================
/// Function-type-object to compare finite elements based on
/// the x and y coordinates of their first node
//================================================================
class FiniteElementComp
{

public:

 /// Comparison. Is x/y coordinate of el1_pt less than that of el2_pt?
 bool operator()(FiniteElement* const& el1_pt, FiniteElement* const& el2_pt) 
  const
  {

   if (el1_pt->node_pt(0)->x(0) < el2_pt->node_pt(0)->x(0))
    {
     return true;
    }
   else if (el1_pt->node_pt(0)->x(0) == el2_pt->node_pt(0)->x(0))
    {
     if (el1_pt->node_pt(0)->x(1) < el2_pt->node_pt(0)->x(1))
      {
       return true;
      }
     else
      {
       return false;
      }
    }
   else
    {
     return false;
    }
  }

};



//======Start_of_warped_plane==============================================
/// Warped plane in 3D space
//=========================================================================
class WarpedPlane : public GeomObject
{

public:

 /// Constructor: Specify amplitude of deflection from flat horizontal plane
 WarpedPlane(const double& ampl) : GeomObject(2,3)
  {
   Ampl=ampl;
  }

 /// Broken copy constructor
 WarpedPlane(const WarpedPlane& dummy) 
  { 
   BrokenCopy::broken_copy("WarpedPlane");
  } 
 
 /// Broken assignment operator
 void operator=(const WarpedPlane&) 
  {
   BrokenCopy::broken_assign("WarpedPlane");
  }


 /// Empty Destructor
 ~WarpedPlane(){}

 /// Position vector at Lagrangian coordinate zeta 
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   // Position vector
   r[0] = zeta[0]+5.0*Ampl*zeta[0]*(zeta[0]-1.0)*(zeta[0]-0.7)*
                           zeta[1]*(zeta[1]-1.0);
   r[1] = zeta[1]+5.0*Ampl*zeta[1]*(zeta[1]-1.0)*(zeta[1]-0.7)*
                           zeta[0]*(zeta[0]-1.0);
   r[2] = 1.0+Ampl*0.25*(1.0-cos(2.0*MathematicalConstants::Pi*zeta[0]))*
    (1.0-cos(2.0*MathematicalConstants::Pi*zeta[1]));
  }
 
 /// Parametrised position on object: r(zeta). Evaluated at
 /// previous timestep. t=0: current time; t>0: previous
 /// timestep. Forward to steady version
 void position(const unsigned& t, const Vector<double>& zeta,
                       Vector<double>& r) const
  {
   position(zeta,r);
  }

 /// Access to amplitude
 double& ampl() {return Ampl;}

 /// How many items of Data does the shape of the object depend on?
 /// None.
 unsigned ngeom_data() const
  {
   return 0;
  }
 
private:

 /// Amplitude of perturbation
 double Ampl;

};



/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////


//=======start_namespace==========================================
/// Global parameters
//================================================================
namespace Global_Physical_Variables
{

 /// GeomObject specifying the shape of the boundary: Initially it's flat.
 WarpedPlane Boundary_geom_object(0.0);

 /// Poisson's ratio
 double Nu=0.3;

 // Generalised Hookean constitutive equations
 GeneralisedHookean Constitutive_law(&Global_Physical_Variables::Nu);
 
} //end namespace


/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////

#ifdef NOREFINE

//=========================================================================
/// Simple cubic mesh upgraded to become a solid mesh
//=========================================================================
template<class ELEMENT>
class ElasticCubicMesh : public virtual SimpleCubicMesh<ELEMENT>, 
                         public virtual SolidMesh
{

public:

 /// Constructor: 
 ElasticCubicMesh(const unsigned &nx, const unsigned &ny, 
                  const unsigned &nz,
                  const double &a, const double &b, 
                  const double &c,
                  TimeStepper* time_stepper_pt = 
                  &Mesh::Default_TimeStepper) :
  SimpleCubicMesh<ELEMENT>(nx,ny,nz,0.0,a,0.0,b,0.0,c,time_stepper_pt),
  SolidMesh()
  {
   
   //Assign the initial lagrangian coordinates
   set_lagrangian_nodal_coordinates();
   
   // Loop over boundary nodes on boundary 5 and setup boundary coordinates
   {
    unsigned b=5;
    unsigned n_node = this->nboundary_node(b);
    Vector<double> zeta(2);
    for (unsigned j=0;j<n_node;j++)
     {
      Node* nod_pt=this->boundary_node_pt(b,j);
      zeta[0]=nod_pt->x(0);
      zeta[1]=nod_pt->x(1);
      nod_pt->set_coordinates_on_boundary(b,zeta);
     }
    this->Boundary_coordinate_exists[b]=true;
   }

  }

 /// Empty Destructor
 virtual ~ElasticCubicMesh() { }

};



#else

//=========================================================================
/// Refineable cubic mesh upgraded to become a solid mesh
//=========================================================================
template<class ELEMENT>
class RefineableElasticCubicMesh : public virtual SimpleCubicMesh<ELEMENT>, 
                                   public virtual RefineableBrickMesh<ELEMENT>,
                                   public virtual SolidMesh
{

public:

 /// Constructor: 
 RefineableElasticCubicMesh(const unsigned &nx, const unsigned &ny, 
                            const unsigned &nz,
                            const double &a, const double &b, 
                            const double &c,
                            TimeStepper* time_stepper_pt = 
                            &Mesh::Default_TimeStepper) :
  SimpleCubicMesh<ELEMENT>(nx,ny,nz,0.0,a,0.0,b,0.0,c,time_stepper_pt),
  RefineableBrickMesh<ELEMENT>(),
  SolidMesh()
  {
   
   this->setup_octree_forest();

   //Assign the initial lagrangian coordinates
   set_lagrangian_nodal_coordinates();
   
   // Loop over boundary nodes on boundary 5 and setup boundary coordinates
   {
    unsigned b=5;
    unsigned n_node = this->nboundary_node(b);
    Vector<double> zeta(2);
    for (unsigned j=0;j<n_node;j++)
     {
      Node* nod_pt=this->boundary_node_pt(b,j);
      zeta[0]=nod_pt->x(0);
      zeta[1]=nod_pt->x(1);
      nod_pt->set_coordinates_on_boundary(b,zeta);
     }
    this->Boundary_coordinate_exists[b]=true;
   }

  }

 /// Empty Destructor
 virtual ~RefineableElasticCubicMesh() { }

};


#endif



/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



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

 /// Update function (empty)
 void actions_before_newton_solve() {}

#ifdef NOREFINE

 /// Access function for the solid mesh
 ElasticCubicMesh<ELEMENT>*& solid_mesh_pt() 
  {return Solid_mesh_pt;} 

#else

 /// Access function for the solid mesh
 RefineableElasticCubicMesh<ELEMENT>*& solid_mesh_pt() 
  {return Solid_mesh_pt;} 

#endif


 /// Actions before adapt: Wipe the mesh of Lagrange multiplier elements
 void actions_before_adapt();

 /// Actions after adapt: Rebuild the mesh of Lagrange multiplier elements
 void actions_after_adapt();

 /// Doc the solution
 void doc_solution();

private:

 /// Create elements that enforce prescribed boundary motion
 /// by Lagrange multiplilers
 void create_lagrange_multiplier_elements();

 /// Delete elements that enforce prescribed boundary motion
 /// by Lagrange multiplilers
 void delete_lagrange_multiplier_elements();

#ifdef NOREFINE

 /// Pointer to solid mesh
 ElasticCubicMesh<ELEMENT>* Solid_mesh_pt;

#else

 /// Pointer to solid mesh
 RefineableElasticCubicMesh<ELEMENT>* Solid_mesh_pt;

#endif


 /// Pointers to meshes of Lagrange multiplier elements
 SolidMesh* Lagrange_multiplier_mesh_pt;

 /// DocInfo object for output
 DocInfo Doc_info;
 
};


//===========start_of_constructor======================================= 
/// Constructor: 
//====================================================================== 
template<class ELEMENT>
PrescribedBoundaryDisplacementProblem<ELEMENT>::
PrescribedBoundaryDisplacementProblem() 
{

 // Create the mesh

 unsigned nel=4;
 if (CommandLineArgs::Argc!=1) nel=2;

 // # of elements in x-direction
 unsigned n_x=nel;

 // # of elements in y-direction
 unsigned n_y=nel;

 // # of elements in y-direction
 unsigned n_z=nel;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;

 // Domain length in z-direction
 double l_z=1.0;


#ifdef NOREFINE

 //Now create the mesh 
 solid_mesh_pt() = new ElasticCubicMesh<ELEMENT>(
  n_x,n_y,n_z,l_x,l_y,l_z);

#else

 //Now create the mesh 
 solid_mesh_pt() = new RefineableElasticCubicMesh<ELEMENT>(
  n_x,n_y,n_z,l_x,l_y,l_z);

 // Set error estimator
 solid_mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 solid_mesh_pt()->max_permitted_error()=0.02;
 solid_mesh_pt()->min_permitted_error()=0.001;

#endif

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

 
 // Construct the mesh of elements that enforce prescribed boundary motion
 // by Lagrange multipliers
 Lagrange_multiplier_mesh_pt=new SolidMesh;
 create_lagrange_multiplier_elements();
 
 // Solid mesh is first sub-mesh
 add_sub_mesh(solid_mesh_pt());

 // Add Lagrange multiplier sub-mesh
 add_sub_mesh(Lagrange_multiplier_mesh_pt);

 // Build combined "global" mesh
 build_global_mesh();
 

 // Pin nodal positions on all boundaries apart from the top one (5) 
 for (unsigned b=0;b<6;b++)
  {
   if (b!=5)
    {
     unsigned n_side = solid_mesh_pt()->nboundary_node(b);
     
     //Loop over the nodes
     for(unsigned i=0;i<n_side;i++)
      {
       solid_mesh_pt()->boundary_node_pt(b,i)->pin_position(0);
       solid_mesh_pt()->boundary_node_pt(b,i)->pin_position(1);
       solid_mesh_pt()->boundary_node_pt(b,i)->pin_position(2);
      }
    }
  }

 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(
  solid_mesh_pt()->element_pt());

 // Setup equation numbering scheme
 cout << "Number of dofs: " << assign_eqn_numbers() << std::endl; 

 // Set output directory
 Doc_info.set_directory("RESLT");

} //end of constructor


//=====================start_of_actions_before_adapt======================
/// Actions before adapt: Wipe the mesh of elements that impose
/// the prescribed boundary displacements
//========================================================================
template<class ELEMENT>
void PrescribedBoundaryDisplacementProblem<ELEMENT>::actions_before_adapt()
{
 // Kill the  elements and wipe surface mesh
 delete_lagrange_multiplier_elements();
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();

}// end of actions_before_adapt



//=====================start_of_actions_after_adapt=======================
///  Actions after adapt: Rebuild the mesh of elements that impose
/// the prescribed boundary displacements
//========================================================================
template<class ELEMENT>
void PrescribedBoundaryDisplacementProblem<ELEMENT>::actions_after_adapt()
{
 // Create the elements that impose the displacement constraint 
 // and attach them to the bulk elements that are
 // adjacent to boundary 5 
 create_lagrange_multiplier_elements();
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();
 
 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(
  solid_mesh_pt()->element_pt());
 
}// end of actions_after_adapt


 
//============start_of_create_lagrange_multiplier_elements===============
/// Create elements that impose the prescribed boundary displacement
//=======================================================================
template<class ELEMENT>
void PrescribedBoundaryDisplacementProblem<ELEMENT>::
create_lagrange_multiplier_elements()
{
 // Lagrange multiplier elements are located on boundary 5:
 unsigned b=5;

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = solid_mesh_pt()->nboundary_element(b);
 
 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    solid_mesh_pt()->boundary_element_pt(b,e));
   
   //Find the index of the face of element e along boundary b
   int face_index = solid_mesh_pt()->face_index_at_boundary(b,e);
      
#ifdef NOREFINE

   // Create new element and add to mesh
   Lagrange_multiplier_mesh_pt->add_element_pt(
    new ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>(
     bulk_elem_pt,face_index));   
   
#else

   // Create new element and add to mesh
   Lagrange_multiplier_mesh_pt->add_element_pt(
    new RefineableImposeDisplacementByLagrangeMultiplierElement<ELEMENT>(
     bulk_elem_pt,face_index));   

#endif


  }  

 
 // Loop over the elements in the Lagrange multiplier element mesh
 // for elements on the top boundary (boundary 5)
 n_element=Lagrange_multiplier_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {

#ifdef NOREFINE

   //Cast to a Lagrange multiplier element
   ImposeDisplacementByLagrangeMultiplierElement<ELEMENT> *el_pt = 
    dynamic_cast<ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>*>
    (Lagrange_multiplier_mesh_pt->element_pt(i));

#else

   //Cast to a Lagrange multiplier element
   RefineableImposeDisplacementByLagrangeMultiplierElement<ELEMENT> *el_pt = 
    dynamic_cast<
    RefineableImposeDisplacementByLagrangeMultiplierElement<ELEMENT>*>
    (Lagrange_multiplier_mesh_pt->element_pt(i));

#endif


   // Set the GeomObject that defines the boundary shape and
   // specify which bulk boundary we are attached to (needed to extract
   // the boundary coordinate from the bulk nodes)
   el_pt->set_boundary_shape_geom_object_pt( 
    &Global_Physical_Variables::Boundary_geom_object,b);
   
   // Loop over the nodes 
   unsigned nnod=el_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = el_pt->node_pt(j);
     
     // Is the node also on boundaries 1, 2, 3 or 4?
     if ((nod_pt->is_on_boundary(1))||
         (nod_pt->is_on_boundary(2))||
         (nod_pt->is_on_boundary(3))||
         (nod_pt->is_on_boundary(4)))
      {
       // How many nodal values were used by the "bulk" element
       // that originally created this node?
       unsigned n_bulk_value=el_pt->nbulk_value(j);
       
       // The remaining ones are Lagrange multipliers and we pin them.
       unsigned nval=nod_pt->nvalue();
       for (unsigned j=n_bulk_value;j<nval;j++)
        {
         nod_pt->pin(j);
        }
      }     
    }
  }
  
} // end of create_lagrange_multiplier_elements




//====start_of_delete_lagrange_multiplier_elements=======================
/// Delete elements that impose the prescribed boundary displacement
/// and wipe the associated mesh
//=======================================================================
template<class ELEMENT>
void PrescribedBoundaryDisplacementProblem<ELEMENT>::delete_lagrange_multiplier_elements()
{
 // How many surface elements are in the surface mesh
 unsigned n_element = Lagrange_multiplier_mesh_pt->nelement();
 
 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill surface element
   delete Lagrange_multiplier_mesh_pt->element_pt(e);
  }
 
 // Wipe the mesh
 Lagrange_multiplier_mesh_pt->flush_element_and_node_storage();

} // end of delete_lagrange_multiplier_elements



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

 // Output Lagrange multipliers
 //----------------------------
 sprintf(filename,"%s/lagr%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);

 // This makes sure the elements are ordered in same way every time
 // the code is run -- necessary for validation tests.
 std::vector<FiniteElement*> el_pt;
 unsigned nelem=Lagrange_multiplier_mesh_pt->nelement();
 for (unsigned e=0;e<nelem;e++)
  {
   el_pt.push_back(Lagrange_multiplier_mesh_pt->finite_element_pt(e));
  }
 std::sort(el_pt.begin(),el_pt.end(),FiniteElementComp());
 for (unsigned e=0;e<nelem;e++)
  {
   el_pt[e]->output(some_file);
  }
 some_file.close();

 // Increment label for output files
 Doc_info.number()++;

} //end doc



//=======start_of_main==================================================
/// Driver code
//======================================================================
int main(int argc, char **argv)
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
 
#ifdef NOREFINE

 //Set up the problem
 PrescribedBoundaryDisplacementProblem<QPVDElement<3,3> > problem;

#else

 //Set up the problem
 PrescribedBoundaryDisplacementProblem<RefineableQPVDElement<3,3> > problem;

 //Refine according to a pattern
 Vector<unsigned> refine_pattern(2);
 refine_pattern[0] = problem.solid_mesh_pt()->nelement()-1; 
 refine_pattern[1] = problem.solid_mesh_pt()->nelement()-4; 
 problem.refine_selected_elements(0,refine_pattern);
 
#endif

 
 // Doc initial domain shape
 problem.doc_solution();

 //Parameter incrementation
 unsigned nstep=10; 
 if (CommandLineArgs::Argc!=1)
  {
   std::cout << "Validation -- only doing one step" << std::endl;
   nstep=1;
  }
 for(unsigned i=0;i<nstep;i++)
  {
   // Increment imposed boundary displacement
   Global_Physical_Variables::Boundary_geom_object.ampl()+=0.02;

#ifdef NOREFINE

   // Solve the problem with Newton's method
   problem.newton_solve();

#else

   // Solve the problem with Newton's method, allowing
   // up to max_adapt mesh adaptations after every solve.
   unsigned max_adapt=1;
   problem.newton_solve(max_adapt);

#endif

   
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








