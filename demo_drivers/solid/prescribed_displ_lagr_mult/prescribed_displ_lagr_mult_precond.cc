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
// is imposed via Lagrange multipliers


//Oomph-lib includes
#include "generic.h"
#include "solid.h"
#include "constitutive.h"

//The mesh
#include "meshes/rectangular_quadmesh.h"

// The preconditioner
#include "multi_physics/pseudo_elastic_preconditioner.h"


using namespace std;

using namespace oomph;

namespace oomph {

//=============start_wrapper==================================================
/// Pseudo-Elastic Solid element class to overload the block preconditioner
/// methods ndof_types() and get_dof_numbers_for_unknowns() to differentiate
/// between DOFs subject to Lagrange multiplier and those that are not.
//============================================================================
template <class ELEMENT>
class PseudoElasticBulkElement : 
 public virtual ELEMENT
{

public:

 /// Constructor
 PseudoElasticBulkElement() : ELEMENT() {}
 
 ///  Returns the number of DOF types associated with this element: Twice
 /// the number of spatial dimensions (for the constrained and 
 /// unconstrained nodal positions).
 unsigned ndof_types() const
  {
   return 2*ELEMENT::dim();
  }
 
 ///  Create a list of pairs for all unknowns in this element,
 /// so that the first entry in each pair contains the global equation
 /// number of the unknown, while the second one contains the number
 /// of the "DOF" that this unknown is associated with.
 /// (Function can obviously only be called if the equation numbering
 /// scheme has been set up.)\n
 /// E.g. in a 3D problem there are 6 types of DOF:\n
 /// 0 - x displacement (without lagr mult traction)\n
 /// 1 - y displacement (without lagr mult traction)\n
 /// 2 - z displacement (without lagr mult traction)\n
 /// 4 - x displacement (with lagr mult traction)\n
 /// 5 - y displacement (with lagr mult traction)\n
 /// 6 - z displacement (with lagr mult traction)\n
 void get_dof_numbers_for_unknowns(
    std::list<std::pair<unsigned long,unsigned> >& dof_lookup_list) const
  {
   // temporary pair (used to store dof lookup prior to being added to list
   std::pair<unsigned,unsigned> dof_lookup;
   
   // number of nodes
   const unsigned n_node = this->nnode();
   
   //Get the number of position dofs and dimensions at the node
   const unsigned n_position_type = ELEMENT::nnodal_position_type();
   const unsigned nodal_dim = ELEMENT::nodal_dimension();
   
   //Integer storage for local unknown
   int local_unknown=0;
   
   //Loop over the nodes
   for(unsigned n=0;n<n_node;n++)
    {
     unsigned offset = 0;
     if (this->node_pt(n)->nvalue() != this->required_nvalue(n))
      {
       offset = ELEMENT::dim();
      }
     
     //Loop over position dofs
     for(unsigned k=0;k<n_position_type;k++)
      {
       //Loop over dimension
       for(unsigned i=0;i<nodal_dim;i++)
        {
         //If the variable is unpinned
         local_unknown = ELEMENT::position_local_eqn(n,k,i);
         if (local_unknown >= 0)
          {
           // store dof lookup in temporary pair: First entry in pair
           // is global equation number; second entry is dof type
           dof_lookup.first = this->eqn_number(local_unknown);
           dof_lookup.second = offset+i;
           
           // add to list
           dof_lookup_list.push_front(dof_lookup);
          }
        }
      }
    }
  }
};





//===========start_face_geometry==============================================
/// FaceGeometry of wrapped element is the same as the underlying element
//============================================================================
template<class ELEMENT>
class FaceGeometry<PseudoElasticBulkElement<ELEMENT> > :
 public virtual FaceGeometry<ELEMENT>
{

public: 

 /// Constructor -- required for more recent versions of gcc
 FaceGeometry() {}

};



}




//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



//================================================================
/// Namespace to "hide" global function that creates an instance
/// of oomph-lib's diagonal preconditioner.
//================================================================
namespace DiagonalPreconditionerHelper
{

///  Create a matrix-based diagonal preconditioner for
/// subsidiary linear systems
 Preconditioner* get_diagonal_preconditioner()
 {
  MatrixBasedDiagPreconditioner* prec_pt=
   new MatrixBasedDiagPreconditioner;
  return prec_pt;
 }

}


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



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

 ///  How many items of Data does the shape of the object depend on?
 /// None.
 unsigned ngeom_data() const
  {
   return 0;
  }
 
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
/// Problem class for deformation of elastic DOF type by prescribed
/// boundary motion.
//====================================================================== 
template<class ELEMENT>
class PrescribedBoundaryDisplacementProblem : public Problem
{

public:

 /// Constructor: Pass in number of elements along axes
 PrescribedBoundaryDisplacementProblem(const unsigned& nel_1d);
 
 /// Update function (empty)
 void actions_after_newton_solve() {}

 ///  Update function (empty)
 void actions_before_newton_solve() {}

 /// Access function for the solid mesh
 ElasticRefineableRectangularQuadMesh<ELEMENT>*& solid_mesh_pt() 
  {return Solid_mesh_pt;} 

 /// Actions before adapt: Wipe the mesh of Lagrange multiplier elements
 void actions_before_adapt();

 /// Actions after adapt: Rebuild the mesh of Lagrange multiplier elements
 void actions_after_adapt();

 /// Doc the solution
 void doc_solution();

private:

 ///  Create elements that enforce prescribed boundary motion
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
/// Constructor:  Pass in number of elements along axes
//====================================================================== 
template<class ELEMENT>
PrescribedBoundaryDisplacementProblem<ELEMENT>::
PrescribedBoundaryDisplacementProblem(const unsigned& nel_1d) 
{

 // Create the mesh

 // # of elements in x-direction
 unsigned n_x=nel_1d;

 // # of elements in y-direction
 unsigned n_y=nel_1d;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;

 //Now create the mesh 
 solid_mesh_pt() = new ElasticRefineableRectangularQuadMesh<ELEMENT>(
  n_x,n_y,l_x,l_y);

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

 // Setup equation numbering scheme
 cout << "Number of dofs: " << assign_eqn_numbers() << std::endl; 

 // Set output directory
 Doc_info.set_directory("RESLT");

 // Create the linear solver
 IterativeLinearSolver* solver_pt=0;

 // If we have trilinos, use it
#ifdef OOMPH_HAS_TRILINOS

 // Create solver
 solver_pt = new TrilinosAztecOOSolver;

 // Use GMRES
 dynamic_cast<TrilinosAztecOOSolver*>(solver_pt)->solver_type() 
  = TrilinosAztecOOSolver::GMRES;
 
#else

 // Use oomph-lib's own GMRES
 solver_pt = new GMRES<CRDoubleMatrix>;

#endif

 // Set solver
 linear_solver_pt() = solver_pt;

 // Create the preconditioner
 PseudoElasticPreconditioner * prec_pt = new PseudoElasticPreconditioner;

 // Set solid and lagrange multiplier meshes
 prec_pt->set_elastic_mesh(Solid_mesh_pt);
 prec_pt->set_lagrange_multiplier_mesh(Lagrange_multiplier_mesh_pt);
 
 // Set the preconditioner
 solver_pt->preconditioner_pt() = prec_pt;

 // Use upper block triangular preconditioner for elastic block?
 if (CommandLineArgs::command_line_flag_has_been_set
     ("--block_upper_for_elastic_block"))
  {
   prec_pt->elastic_preconditioner_type()=
    PseudoElasticPreconditioner::Block_upper_triangular_preconditioner;
  }

#ifdef OOMPH_HAS_HYPRE

 // Use Hypre as subsidiary preconditioner (inexact solver) for
 // linear (sub-)systems to be solved in the elastic block?
 if (CommandLineArgs::command_line_flag_has_been_set
     ("--hypre_for_elastic_blocks"))
  {
   prec_pt->set_elastic_subsidiary_preconditioner
    (Pseudo_Elastic_Preconditioner_Subsidiary_Operator_Helper::
     get_elastic_preconditioner_hypre);
  }

#endif


#ifdef OOMPH_HAS_TRILINOS

 // Use Trilinos CG as subsidiary preconditioner (inexact solver) for
 // linear (sub-)systems to be solved in the Lagrange multiplier block?
 if (CommandLineArgs::command_line_flag_has_been_set
     ("--trilinos_cg_for_lagrange_multiplier_blocks"))
  {
   prec_pt->set_lagrange_multiplier_subsidiary_preconditioner
    (Pseudo_Elastic_Preconditioner_Subsidiary_Operator_Helper::
     get_lagrange_multiplier_preconditioner);
  }

#endif

 // Use diagonal scaling as subsidiary preconditioner (inexact solver) for
 // linear (sub-)systems to be solved in the elastic block?
 if (CommandLineArgs::command_line_flag_has_been_set
     ("--diagonal_scaling_for_elastic_blocks"))
  {
   prec_pt->set_elastic_subsidiary_preconditioner(
    DiagonalPreconditionerHelper::get_diagonal_preconditioner);
  }
 
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
   
   //Find the index of the face of element e along boundary b
   int face_index = solid_mesh_pt()->face_index_at_boundary(b,e);
      
   // Create new element and add to mesh
   Lagrange_multiplier_mesh_pt->add_element_pt(
    new ImposeDisplacementByLagrangeMultiplierElement<ELEMENT>(
     bulk_elem_pt,face_index));   
  }  

 
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
int main(int argc, char* argv[])
{


#ifdef OOMPH_HAS_MPI

 // Start up mpi if oomph-lib has been compiled with parallel support
 // because parallel versions of trilinos and hypre which are then called
 // by default) need it.
 MPI_Helpers::init(argc,argv,false);

#endif

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified

 // Number of elements along axes
 unsigned nel_1d=5;
 CommandLineArgs::specify_command_line_flag("--nel_1d",&nel_1d);

 // Suppress adaptation (for study of iteration count vs uniform mesh 
 // refinement, say)
 CommandLineArgs::specify_command_line_flag("--no_adapt");

 // Use block upper triangular preconditioner for elasticity block
 CommandLineArgs::specify_command_line_flag("--block_upper_for_elastic_block");

 // Use Hypre as subsidiary preconditioner (inexact solver) for
 // linear (sub-)systems to be solved in the elastic block?
 CommandLineArgs::specify_command_line_flag("--hypre_for_elastic_blocks");

 // Use Trilinos CG as subsidiary preconditioner (inexact solver) for
 // linear (sub-)systems to be solved in the Lagrange multiplier block?
 CommandLineArgs::specify_command_line_flag
  ("--trilinos_cg_for_lagrange_multiplier_blocks");

 // Use diagonal scaling as subsidiary preconditioner (inexact solver) for
 // linear (sub-)systems to be solved in the elastic block?
 // [only used for exercise]
 CommandLineArgs::specify_command_line_flag
  ("--diagonal_scaling_for_elastic_blocks");
  
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 //Set up the problem with specified number of elements along axes.
 PrescribedBoundaryDisplacementProblem<
  PseudoElasticBulkElement<RefineableQPVDElement<2,3> > > problem(nel_1d);
 
 // Doc initial domain shape
 problem.doc_solution();

 // Max. number of adaptations per solve
 unsigned max_adapt=1;
 if (CommandLineArgs::command_line_flag_has_been_set("--no_adapt"))
  {
   max_adapt=0;
  }

 //Parameter incrementation
 unsigned nstep=2; 
 for(unsigned i=0;i<nstep;i++)
  {
   // Increment imposed boundary displacement
   Global_Physical_Variables::Boundary_geom_object.ampl()+=0.1;

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








