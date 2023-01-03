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
// Driver for Solid deformation -- driven by boundary motion which
// is imposed via Lagrange multipliers. Slightly adjusted version
// to test procedure for resizing hanging nodes.


//Oomph-lib includes
#include "generic.h"
#include "solid.h"
#include "constitutive.h"

//The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;


//================================================================
/// Function-type-object to compare finite elements based on
/// their x coordinate
//================================================================
class FiniteElementComp
{

public:

 /// Comparison. Is x coordinate of el1_pt less than that of el2_pt?
 bool operator()(FiniteElement* const& el1_pt, FiniteElement* const& el2_pt) 
  const
  {
   return el1_pt->node_pt(0)->x(0) < el2_pt->node_pt(0)->x(0);
  }

};



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

 /// Position vector at Lagrangian coordinate zeta 
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
   // Position vector
   r[0] = zeta[0]+5.0*Ampl*zeta[0]*(zeta[0]-0.4)*(zeta[0]-0.7);
   r[1] = 1.0+Ampl*0.5*(1.0-cos(2.0*MathematicalConstants::Pi*zeta[0]/0.4));
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

 /// Actually attach elements?
 bool Actually_attach_face_elements=true;

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

 /// Update function (empty)
 void actions_before_newton_solve() {}

 /// Access function for the solid mesh
 ElasticRefineableRectangularQuadMesh<ELEMENT>*& solid_mesh_pt() 
  {return Solid_mesh_pt;} 

 /// Actions before adapt: Wipe the mesh of Lagrange multiplier elements
 void actions_before_adapt();

 /// Actions after adapt: Rebuild the mesh of Lagrange multiplier elements
 void actions_after_adapt();

 /// Actions before distribute: Wipe the mesh of Lagrange multiplier elements
 void actions_before_distribute()
  {
   actions_before_adapt();
  }

 /// Actions after distribute: Rebuild the mesh of Lagrange multiplier elements
 void actions_after_distribute()
  {
   actions_after_adapt();
  }

 /// Helper function to build mesh (needed during load balancing and
 /// also used by constructor)
 void build_mesh();

 /// Doc the solution
 void doc_solution();

private:

 /// Create elements that enforce prescribed boundary motion
 /// by Lagrange multiplilers
 void create_lagrange_multiplier_elements();

 /// Delete elements that enforce prescribed boundary motion
 /// by Lagrange multiplilers
 void delete_lagrange_multiplier_elements();

 /// Pointer to solid mesh
 ElasticRefineableRectangularQuadMesh<ELEMENT>* Solid_mesh_pt;

 /// Pointers to meshes of Lagrange multiplier elements
 SolidMesh* Lagrange_multiplier_mesh_pt;

 /// DocInfo object for output
 DocInfo Doc_info;
 
};


//===========start_of_constructor======================================= 
/// Constructor: 
//====================================================================== 
template<class ELEMENT>
PrescribedBoundaryDisplacementProblem<ELEMENT>::PrescribedBoundaryDisplacementProblem() 
{
 // Set output directory
 Doc_info.set_directory("RESLT_resize_test");

 // Create the mesh
 build_mesh();

 // Setup equation numbering scheme
 cout << "Number of dofs: " << assign_eqn_numbers() << std::endl; 

} //end of constructor



//===========start_of_build_mesh======================================== 
/// Helper function to build mesh. Required during load balancing and
/// also called from problem constructor.
//====================================================================== 
template<class ELEMENT>
void PrescribedBoundaryDisplacementProblem<ELEMENT>::build_mesh()
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


 // Demonstrate that resizing halo nodes matters by uncommenting this
 // "optimisation"
 //solid_mesh_pt()->disable_resizing_of_halo_nodes();

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

 // Pin nodal positions on all boundaries apart from the top one (2) 
 for (unsigned b=0;b<4;b++)
  {
   if (b!=2)
    {
     unsigned n_side = solid_mesh_pt()->nboundary_node(b);
     
     //Loop over the nodes
     for(unsigned i=0;i<n_side;i++)
      {
       solid_mesh_pt()->boundary_node_pt(b,i)->pin_position(0);
       solid_mesh_pt()->boundary_node_pt(b,i)->pin_position(1);
      }
    }
  }

 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
  solid_mesh_pt()->element_pt());
 
} //end of build mesh function


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
 // adjacent to boundary 2 
 create_lagrange_multiplier_elements();
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();
 
 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
  solid_mesh_pt()->element_pt());
 
}// end of actions_after_adapt


 
//============start_of_create_lagrange_multiplier_elements===============
/// Create elements that impose the prescribed boundary displacement
//=======================================================================
template<class ELEMENT>
void PrescribedBoundaryDisplacementProblem<ELEMENT>::
create_lagrange_multiplier_elements()
{

 // Actually attach elements?
 if (!Global_Physical_Variables::Actually_attach_face_elements) return;

 // Lagrange multiplier elements are located on boundary 2:
 unsigned b=2;

 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = solid_mesh_pt()->nboundary_element(b);
 
 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    solid_mesh_pt()->boundary_element_pt(b,e));
 
   // Only attach FaceElements to elements that are not haloed
   // (the way this is implemented is a bit hacky and relies on
   // the fact that we're using 10x10 elements in the original mesh,
   // are distributing it over 2 procs along the vertical line of
   // symmetry and are not re-setting the Lagrangian coordinates.
   if (dynamic_cast<SolidNode*>(bulk_elem_pt->node_pt(7))->xi(0)<=0.4)
    {
     //Find the index of the face of element e along boundary b
     int face_index = solid_mesh_pt()->face_index_at_boundary(b,e);
     
     // Create new element and add to mesh
     Lagrange_multiplier_mesh_pt->add_element_pt(
      new ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>(
       bulk_elem_pt,face_index));   
    }
  }  

 // This is crucial: Lagrange multiplier elements have created additional
 // unknowns at their nodes. Corresponding halo nodes may have to
 // resized for consistency (they're not necessarily touched by
 // a Lagrange multiplier element that would do the resizing for us).
 solid_mesh_pt()->resize_halo_nodes();

 // Loop over the elements in the Lagrange multiplier element mesh
 // for elements on the top boundary (boundary 2)
 n_element=Lagrange_multiplier_mesh_pt->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a Lagrange multiplier element
   ImposeDisplacementByLagrangeMultiplierElement<ELEMENT> *el_pt = 
    dynamic_cast<ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>*>
    (Lagrange_multiplier_mesh_pt->element_pt(i));

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
     
     oomph_info << "Lagr mult node at: " 
                << dynamic_cast<SolidNode*>(nod_pt)->xi(0) << " " 
                << dynamic_cast<SolidNode*>(nod_pt)->xi(1) << "\n"; 

     // Is the node also on boundary 1 or 3?
     if ((nod_pt->is_on_boundary(1))||(nod_pt->is_on_boundary(3)))
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

 oomph_info << "Outputting for step: " << Doc_info.number() << "\n";


 // Output nodes and number of values
 //----------------------------------
 sprintf(filename,"%s/nodes%i_on_proc%i.dat",Doc_info.directory().c_str(),
         Doc_info.number(),this->communicator_pt()->my_rank());
 some_file.open(filename);
 unsigned nnod=solid_mesh_pt()->nnode();
 for (unsigned j=0;j<nnod;j++)
  {
   Node* nod_pt = solid_mesh_pt()->node_pt(j);
   some_file << nod_pt->x(0) << " " 
             << nod_pt->x(1) << " "
             << nod_pt->nvalue() << "\n";
  }
 some_file.close();

 // Output shape of deformed body
 //------------------------------
 sprintf(filename,"%s/soln%i_on_proc%i.dat",Doc_info.directory().c_str(),
         Doc_info.number(),this->communicator_pt()->my_rank());
 some_file.open(filename);
 solid_mesh_pt()->output(some_file,n_plot);
 some_file.close();

 // Output Lagrange multipliers
 //----------------------------
 sprintf(filename,"%s/lagr%i_on_proc%i.dat",Doc_info.directory().c_str(),
         Doc_info.number(),this->communicator_pt()->my_rank());
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
int main(int argc, char* argv[])
{
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif

 // Create mayhem by initially building the mesh without attached
 // FaceElements
 Global_Physical_Variables::Actually_attach_face_elements=false;

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Validation run?
 CommandLineArgs::specify_command_line_flag("--validation");
 
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 //Set up the problem
 PrescribedBoundaryDisplacementProblem<RefineableQPVDElement<2,3> > problem;

 // Validation run
 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   // Fake but repeatable load balancing for self-test
   problem.set_default_partition_in_load_balance();
  }

#ifdef OOMPH_HAS_MPI

 // Distribute the problem
 DocInfo mesh_doc_info;
 mesh_doc_info.set_directory("RESLT");
 mesh_doc_info.number()=38;

 /// Setup partitioning along vertical symmetry line
 unsigned n_partition=problem.mesh_pt()->nelement();
 Vector<unsigned> element_partition(n_partition,0);
 std::string input_string;
 for (unsigned e=0;e<n_partition;e++)
  {
   // Base decision on where element is on central node on bottom edge
   if (dynamic_cast<SolidNode*>(problem.mesh_pt()->finite_element_pt(e)->
                                node_pt(1))->xi(0)<0.5)
    {
     element_partition[e]=1;
    }
  }
 
 // Distribute and check halo schemes
 bool report_stats=false;
 problem.distribute(element_partition,mesh_doc_info,report_stats);
 problem.check_halo_schemes(mesh_doc_info);

#endif
 
 // Now start the trouble: The mesh has been distributed and
 // we now attach the face elements only to the non-haloed elements
 // in the left half of the domain. This means that one node gets
 // resized there while its halo counterpart on the other processor
 // doesn't. This will be reconciled by Mesh::resize_halo_nodes()
 Global_Physical_Variables::Actually_attach_face_elements=true;

 // Quick hacky way of attaching the face elements
 problem.actions_after_adapt(); 

 // Must now re-assign eqn numbers because additional unknowns have been
 // created
 problem.assign_eqn_numbers(); 

 // Doc initial domain shape
 problem.doc_solution();

 // Max. number of adaptations per solve
 unsigned max_adapt=1;

 //Parameter incrementation
 unsigned nstep=2; 
 for(unsigned i=0;i<nstep;i++)
  {
   // Increment imposed boundary displacement (more gently than
   // in other code because we're not resetting Lagrangian coordinates
   // to faciliate identification of elements in left and right half
   // of domain
   Global_Physical_Variables::Boundary_geom_object.ampl()+=0.02;

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

   // NOT DOING THIS HERE TO FACILITATE IDENTIFICATION OF ELEMENTS
   // AS BEING IN LEFT OR RIGHT HALF OF DOMAIN
   // problem.solid_mesh_pt()->set_lagrangian_nodal_coordinates();
   
  }

 oomph_info << "Load balancing \n";

 // Load balance and re-solve without further adaptation
 problem.load_balance();
 problem.newton_solve();
 
 // Doc solution
 problem.doc_solution();


 // Increment imposed boundary displacement yet again...
 Global_Physical_Variables::Boundary_geom_object.ampl()+=0.02;

 // ...and re-solve for a final time
 problem.newton_solve();
 
 // Doc solution
 problem.doc_solution();
 
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif

} //end of main








