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

#include<fenv.h>

//Generic routines
#include "generic.h"

// Poisson
#include "poisson.h"

// The mesh
#include "meshes/triangle_mesh.h"
 
// Get the mesh
#include "meshes/tetgen_mesh.h" 
#include "meshes/refineable_tetgen_mesh.h"
#include "meshes/gmsh_tet_mesh.h"

// Get the faceted surfaces
#include "tetmesh_faceted_surfaces.h"

// Tetgen or Gmsh
// #define DO_TETGEN


using namespace std;

using namespace oomph;



/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////



//=============================================================
/// Namespace for problem parameters
//=============================================================
namespace Global_Parameters
{
 /// (Half-)width of the box
 double Box_half_width = 1.5;

 /// (Half)height of the box
 double Box_half_length = 1.0;

 /// Specify how to call gmsh from the command line
 std::string Gmsh_command_line_invocation="/home/mheil/gmesh/bin/bin/gmsh";

}

/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////




//====================================================================
/// Demo class solves Poisson problem using Gmsh mesh
//====================================================================
template<class ELEMENT> 
class TetmeshPoissonProblem : public Problem
{

public:

 /// Constructor
 TetmeshPoissonProblem();
  
 /// Destructor (empty)
 ~TetmeshPoissonProblem()
  {
   //Delete the objects
   unsigned nh = Inner_boundary_pt.size();
   for(unsigned h=0;h<nh;++h)
    {
     delete Inner_boundary_pt[h];
    }
   delete Outer_boundary_pt;
  }
 
 /// Actions before adapt (empty)
 void actions_before_adapt()
  {}

 /// Totally new mesh; build elements and apply boundary conditions
 void actions_after_adapt()
  {
   // Complete problem setup
   complete_problem_setup();
  }
 
 /// Update the problem specs before solve: (empty)
 void actions_before_newton_solve(){}

 /// Update the problem specs before solve (empty)
 void actions_after_newton_solve(){}
 
 /// Doc the solution
 void doc_solution(const unsigned& nplot, DocInfo& doc_info);

private:
 
 /// Apply BCs and make elements functional
 void complete_problem_setup();

 /// Helper function to apply boundary conditions
 void apply_boundary_conditions();
 
#ifdef DO_TETGEN

 /// Bulk mesh
 RefineableTetgenMesh<ELEMENT>* Bulk_mesh_pt;

#else

 /// Bulk mesh
 RefineableGmshTetMesh<ELEMENT>* Bulk_mesh_pt;

#endif


 /// Storage for the outer boundary object
 TetMeshFacetedClosedSurface* Outer_boundary_pt;

 /// Inner boundary
 Vector<TetMeshFacetedSurface*> Inner_boundary_pt;

 /// First boundary ID for outer boundary
 unsigned First_boundary_id_for_outer_boundary;
 

 // Disk by itself
 //---------------

 /// First boundary ID for warped disk
 unsigned First_disk_boundary_id;

 /// Last boundary ID for warped disk
 unsigned Last_disk_boundary_id;


 // Disk with torus round the edges
 //--------------------------------

 /// Region ID for torus around edge of warped disk
 unsigned Torus_region_id;

 /// First boundary ID for disk that is surrounded by torus
 unsigned First_disk_with_torus_boundary_id;
 
 /// Last boundary ID for disk that is surrounded by torus
 unsigned Last_disk_with_torus_boundary_id;

 /// First boundary ID for torus surrounding edge of disk
 unsigned First_torus_boundary_id;

 /// Last boundary ID for torus surrounding edge of disk
 unsigned Last_torus_boundary_id;
 
 /// Storage for one-based boundary IDs for boundaries on disk within
 ///  the torus region
 Vector<unsigned> One_based_boundary_id_for_disk_within_torus;

 /// Storage for one-based boundary IDs for boundaries on disk 
 /// outside the torus region
 Vector<unsigned> One_based_boundary_id_for_disk_outside_torus;


 // Disk with layers
 //-----------------

 /// Region ID for layer above disk
 unsigned Above_disk_region_id;

 /// Region ID for layer below disk
 unsigned Below_disk_region_id;

 /// First boundary ID for disk with layers
 unsigned First_disk_with_layers_boundary_id;
 
 /// Last boundary ID for disk with layers
 unsigned Last_disk_with_layers_boundary_id;

 // Volumes
 //--------

 /// Sanity check: Exact bounded volume
 double Exact_bounded_volume;

};



//========================================================================
/// Constructor
//========================================================================
template<class ELEMENT>
TetmeshPoissonProblem<ELEMENT>::TetmeshPoissonProblem()
{ 

 // OUTER BOUNDARY
 //===============

 // Start boundary IDs for outer boundary from some crazy offset
 // (just for testing). By default the one-based boundary IDs go from
 // 1 to 6; let's start from 1001.
 unsigned outer_boundary_id_offset=1000;

 //Make the outer boundary object
 Outer_boundary_pt = new CubicTetMeshFacetedSurface(
  Global_Parameters::Box_half_width,
  Global_Parameters::Box_half_length,
  outer_boundary_id_offset);

 // Look, we can visualise the faceted surface!
 Outer_boundary_pt->output("outer_faceted_surface.dat");

 // First oomph-lib (zero-based!) boundary ID for outer boundary
 First_boundary_id_for_outer_boundary=outer_boundary_id_offset;
 
 // For sanity check:
 Exact_bounded_volume=
  2.0*Global_Parameters::Box_half_width*
  2.0*Global_Parameters::Box_half_width*
  2.0*Global_Parameters::Box_half_length;

 
 // INTERNAL BOUNDARIES
 //====================



 // A warped disk
 //--------------

 // Warped disk with specified amplitude and wavenumber for warping
 double epsilon=0.1; 
 unsigned n=4;
 double z_offset=0.5;
 WarpedCircularDisk* disk_pt=new WarpedCircularDisk(epsilon,n,z_offset);

 // (Half) number of segments used to represent the disk perimeter
 unsigned half_nsegment=30; 
 
 // Start enumeration from here
 unsigned first_one_based_boundary_id=9001; 
 
 // Provide storage for last boundary ID on disk (don't know
 // in advance how many facets (with distinct IDs) the disk
 // is going to be represented by
 unsigned last_one_based_boundary_for_disk_id=0;
 
 // Create faceted representation of warped disk
 DiskTetMeshFacetedSurface* disk_faceted_surface_pt=
  new DiskTetMeshFacetedSurface(disk_pt,
                                half_nsegment,
                                first_one_based_boundary_id,
                                last_one_based_boundary_for_disk_id);

 // Keep track of the (zero-based) boundary IDs
 First_disk_boundary_id=first_one_based_boundary_id-1;
 Last_disk_boundary_id=last_one_based_boundary_for_disk_id-1;

 // Look, we can visualise the faceted surface!
 disk_faceted_surface_pt->output("warped_disk_faceted_surface.dat");
 
 // Add as inner boundary for mesh
 Inner_boundary_pt.push_back(disk_faceted_surface_pt);



 // A warped disk surrounded by a torus
 //------------------------------------
 
 // Radius of torus region
 double r_torus=0.1;

 // Warped disk with specified amplitude and wavenumber for warping
 epsilon=0.1;
 n=5;

 // Thickness of annular region on disk = radius of torus surrounding the
 // edge
 double h_annulus=r_torus;
 WarpedCircularDiskWithAnnularInternalBoundary* disk2_pt =
  new WarpedCircularDiskWithAnnularInternalBoundary(h_annulus,epsilon,n);

 // Number of vertices around perimeter of torus
 unsigned nvertex_torus=20;

 // Enumerate the boundaries making up the disk starting with this
 // one-based ID
 unsigned first_one_based_disk_with_torus_boundary_id=
  last_one_based_boundary_for_disk_id+200;

 
 // These get returned
 unsigned last_one_based_disk_with_torus_boundary_id=0;
 unsigned first_one_based_torus_boundary_id=0;
 unsigned last_one_based_torus_boundary_id=0;

 // One-based region ID for torus
 unsigned one_based_torus_region_id=4;


 // Build disk with torus around the edge
 DiskWithTorusAroundEdgeTetMeshFacetedSurface* disk_with_torus_pt=
  new DiskWithTorusAroundEdgeTetMeshFacetedSurface(
   disk2_pt,
   half_nsegment,
   r_torus,
   nvertex_torus,
   first_one_based_disk_with_torus_boundary_id,
   one_based_torus_region_id, 
   last_one_based_disk_with_torus_boundary_id,
   first_one_based_torus_boundary_id,
   last_one_based_torus_boundary_id,
   One_based_boundary_id_for_disk_within_torus,
   One_based_boundary_id_for_disk_outside_torus);


 /// Keep track of (zero-based) IDs
 Torus_region_id=one_based_torus_region_id-1; 
 First_disk_with_torus_boundary_id=
  first_one_based_disk_with_torus_boundary_id-1;
 Last_disk_with_torus_boundary_id=
  last_one_based_disk_with_torus_boundary_id-1;
 First_torus_boundary_id=first_one_based_torus_boundary_id-1;
 Last_torus_boundary_id=last_one_based_torus_boundary_id-1;

 // Look, we can visualise the faceted surface!
 disk_with_torus_pt->output("warped_disk_with_torus_faceted_surface.dat");
 
 // Add as inner boundary for mesh
 Inner_boundary_pt.push_back(disk_with_torus_pt);
 


 // A warped disk bounded by a layer above and below
 //-------------------------------------------------

 // Warped disk with specified amplitude and wavenumber for warping
 epsilon=0.1;
 n=6;
 z_offset=-0.5;
 WarpedCircularDisk* disk3_pt=new WarpedCircularDisk(epsilon,n,z_offset);

 // Start enumeration from here
 first_one_based_boundary_id=last_one_based_torus_boundary_id+200;
 
 // Provide storage for last boundary ID on disk (don't know
 // in advance how many facets (with distinct IDs) the disk
 // is going to be represented by
 unsigned last_one_based_boundary_for_disk_with_layers_id=0;
 unsigned last_one_based_boundary_id=0;
  
 unsigned one_based_region_id_above_disk=2408;
 unsigned one_based_region_id_below_disk=2409;


 // Create faceted representation of warped disk
 DiskWithTwoLayersTetMeshFacetedSurface* disk_with_layers_faceted_surface_pt=
  new DiskWithTwoLayersTetMeshFacetedSurface
  (disk3_pt,
   half_nsegment,
   first_one_based_boundary_id,
   one_based_region_id_above_disk,
   one_based_region_id_below_disk,
   last_one_based_boundary_for_disk_with_layers_id,
   last_one_based_boundary_id);

 // // Keep track of the (zero-based) boundary IDs
 First_disk_with_layers_boundary_id=first_one_based_boundary_id-1;
 Last_disk_with_layers_boundary_id=
  last_one_based_boundary_for_disk_with_layers_id-1;

 // Zero-based region IDs
 Above_disk_region_id=one_based_region_id_above_disk-1; 
 Below_disk_region_id=one_based_region_id_below_disk-1; 


 // Look, we can visualise the faceted surface!
 disk_with_layers_faceted_surface_pt->output
  ("warped_disk_with_layers_faceted_surface.dat");
 
 // Add as inner boundary for mesh
 Inner_boundary_pt.push_back(disk_with_layers_faceted_surface_pt);


 // Build the mesh
 //--------------- 

 // Initial element volume
 double initial_element_volume=1.0;

 // Setup parameters for gmsh
 GmshParameters* gmsh_parameters_pt=
  new GmshParameters(Outer_boundary_pt,
                     Global_Parameters::Gmsh_command_line_invocation);

 // Element volume
 gmsh_parameters_pt->element_volume()=initial_element_volume;


 // Specify inner boundaries
 gmsh_parameters_pt->internal_surface_pt()=Inner_boundary_pt;

 // Filename for file in which target element size is stored
 // (for disk-based operation of gmsh)
 gmsh_parameters_pt->stem_for_filename_gmsh_size_transfer()=
  "target_size_on_grid";
 gmsh_parameters_pt->counter_for_filename_gmsh_size_transfer()=0;

 // Problem is linear so we don't need to transfer the solution to the
 // new mesh; we keep it on for self-test purposes...
 // gmsh_parameters_pt->disable_projection();

 // Redirect gmsh on-screen output
 gmsh_parameters_pt->gmsh_onscreen_output_file_name()=
  "RESLT/gmsh_on_screen_output.dat";

 // Not needed, of course, but here to test out the handling
 // of timesteppers
 add_time_stepper_pt(new Steady<1>);

#ifdef DO_TETGEN

 // And now build it... 
 Bulk_mesh_pt =
  new RefineableTetgenMesh<ELEMENT>(Outer_boundary_pt,
                                    Inner_boundary_pt,
                                    initial_element_volume,
                                    this->time_stepper_pt());

 // Problem is linear so we don't need to transfer the solution to the
 // new mesh; we keep it on for self-test purposes...
 Bulk_mesh_pt->disable_projection();

#else

 // And now build it...
 Bulk_mesh_pt = new RefineableGmshTetMesh<ELEMENT>(gmsh_parameters_pt,
                                                   this->time_stepper_pt());

#endif

 // Add sub-mesh
 Problem::mesh_pt()=Bulk_mesh_pt;

 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Bulk_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // Set targets for spatial adaptivity
 Bulk_mesh_pt->max_permitted_error()=0.0005; 
 Bulk_mesh_pt->min_permitted_error()=0.00001;
 
 // Complete problem setup
 complete_problem_setup();
 
#ifdef OOMPH_HAS_HYPRE

 // Create a new Hypre linear solver
 HypreSolver* hypre_linear_solver_pt = new HypreSolver;
 
 // Set the linear solver for problem
 linear_solver_pt() = hypre_linear_solver_pt;
 
 // Set some solver parameters
 hypre_linear_solver_pt->max_iter() = 100;
 hypre_linear_solver_pt->tolerance() = 1e-10;
 hypre_linear_solver_pt->amg_simple_smoother() = 1;
 hypre_linear_solver_pt->disable_doc_time();
 hypre_linear_solver_pt->enable_hypre_error_messages();
 hypre_linear_solver_pt->amg_print_level() = 0;
 hypre_linear_solver_pt->krylov_print_level() = 0;
 hypre_linear_solver_pt->hypre_method() = HypreSolver::BoomerAMG;
   
#endif

 // Setup equation numbering scheme
 oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

}


//========================================================================
/// Complete problem setup
//========================================================================
template<class ELEMENT>
void TetmeshPoissonProblem<ELEMENT>::complete_problem_setup()
{
 // Apply bcs
 apply_boundary_conditions();
}

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void TetmeshPoissonProblem<ELEMENT>::apply_boundary_conditions()
{
 

 ofstream pin_file;
 pin_file.open("pinned_nodes.dat");

 // Identify boundary ids of pinned nodes


 Vector<unsigned> pinned_boundary_id;
 for (unsigned ibound=First_disk_boundary_id;
      ibound<=Last_disk_boundary_id;ibound++)
  {
   pinned_boundary_id.push_back(ibound);
  }
 for (unsigned ibound=First_disk_with_torus_boundary_id;
      ibound<=Last_disk_with_torus_boundary_id;ibound++)
  {
   pinned_boundary_id.push_back(ibound);
  }
 for (unsigned ibound=First_disk_with_layers_boundary_id;
      ibound<=Last_disk_with_layers_boundary_id;ibound++)
  {
   pinned_boundary_id.push_back(ibound);
  }
 for (unsigned ibound=First_boundary_id_for_outer_boundary;
      ibound<First_boundary_id_for_outer_boundary+6;ibound++)
  {
   pinned_boundary_id.push_back(ibound);
  }

 // Loop over pinned boundaries
 unsigned num_pin_bnd=pinned_boundary_id.size();
 for (unsigned bnd=0;bnd<num_pin_bnd;bnd++)
  {
   unsigned ibound=pinned_boundary_id[bnd];
   unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   if (num_nod==0)
    {
     std::ostringstream error_message;
     error_message << "No boundary nodes on boundary " 
                   << ibound << "! Something's gone wrong!\n";
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
     Vector<double> x(3);
     x[0]=Bulk_mesh_pt->boundary_node_pt(ibound,inod)->x(0);
     x[1]=Bulk_mesh_pt->boundary_node_pt(ibound,inod)->x(1);
     x[2]=Bulk_mesh_pt->boundary_node_pt(ibound,inod)->x(2);

     // Zero BC on outer boundary; unit BC elsewhere
     double value=0.0;
     if ((ibound>=First_disk_with_torus_boundary_id)&&
         (ibound<=Last_disk_with_torus_boundary_id))
      {
       value=1.0;
      }
     if ((ibound>=First_disk_boundary_id)&&
         (ibound<=Last_disk_boundary_id))
      {
       value=1.0;
      }
     if ((ibound>=First_disk_with_layers_boundary_id)&&
         (ibound<=Last_disk_with_layers_boundary_id))
      {
       value=1.0;
      }

     Bulk_mesh_pt->boundary_node_pt(ibound,inod)->set_value(0,value);
     pin_file << x[0] << " " 
              << x[1] << " " 
              << x[2] << " " 
              << std::endl;
    }
  }
 
 pin_file.close();



} // end set bc



//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void TetmeshPoissonProblem<ELEMENT>::doc_solution(const unsigned& nplot,
                                                DocInfo& doc_info)
{ 


 bool do_bulk_output=true;
 if (CommandLineArgs::command_line_flag_has_been_set("--suppress_bulk_output"))
  {
   do_bulk_output=false;
  }

 ofstream some_file;
 ofstream some_file2;
 ofstream face_some_file;
 ofstream coarse_some_file;
 char filename[100];


 // Doc mesh quality (Ratio of max. edge length to min. height,
 /// so if it's very large it's BAAAAAD)
 sprintf(filename,"%s/mesh_quality%i.dat",
         doc_info.directory().c_str(),
         doc_info.number()); 
 ofstream quality_file;
 quality_file.open(filename);
 if (do_bulk_output) Bulk_mesh_pt->assess_mesh_quality(quality_file);
 quality_file.close();

 
 // Output elements adjacent to outer boundary
 //-------------------------------------------
 sprintf(filename,"%s/elements_next_to_outer_boundary%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 for (unsigned ibound=First_boundary_id_for_outer_boundary;
      ibound<First_boundary_id_for_outer_boundary+6;ibound++)
  {
   unsigned n_el=Bulk_mesh_pt->nboundary_element(ibound);
   for (unsigned e=0;e<n_el;e++)
    {
     if (do_bulk_output) 
      {
       Bulk_mesh_pt->boundary_element_pt(ibound,e)->
        output(some_file,nplot);
      }
    }
  }
 some_file.close();


 // Output boundary coordinates on outer boundary
 //-----------------------------------------------
 sprintf(filename,"%s/boundary_coordinates_outer_boundary%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 for (unsigned ibound=First_boundary_id_for_outer_boundary;
      ibound<First_boundary_id_for_outer_boundary+6;ibound++)
  {
   Bulk_mesh_pt->Mesh::template 
    doc_boundary_coordinates<ELEMENT>(ibound,some_file);
  }
 some_file.close();

 // Output boundary coordinates on outer boundary
 //-----------------------------------------------
 unsigned n_b=Bulk_mesh_pt->nboundary();
 oomph_info << "number of boundaries in bulk mesh: " << n_b << std::endl;
 sprintf(filename,"%s/boundary_coordinates%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 for (unsigned ibound=0;ibound<n_b;ibound++)
  {
   if (Bulk_mesh_pt->boundary_coordinate_exists(ibound))
    {
     Bulk_mesh_pt->Mesh::template 
      doc_boundary_coordinates<ELEMENT>(ibound,some_file);
    }
  }
 some_file.close();


 // Output boundaries
 //------------------
 sprintf(filename,"%s/boundaries%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output_boundaries(some_file);
 some_file.close();


 // Output volumes and areas
 //-------------------------
 std::ofstream volumes_and_areas_file;
 sprintf(filename,"%s/volumes%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 volumes_and_areas_file.open(filename);


 // Output bulk elements in torus region
 //-------------------------------------
 double volume_in_torus_region=0.0;
 sprintf(filename,"%s/soln_in_torus_region%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 unsigned region_id=Torus_region_id;
 unsigned n_el=Bulk_mesh_pt->nregion_element(region_id);
 for (unsigned e=0;e<n_el;e++)
  {
   if (do_bulk_output) 
    {
     Bulk_mesh_pt->region_element_pt(region_id,e)->output(some_file,nplot);
    }
   volume_in_torus_region+=Bulk_mesh_pt->
    region_element_pt(region_id,e)->size();
  }
 some_file.close();

 // Output bulk elements above disk
 //--------------------------------
 double volume_in_above_disk_region=0.0;
 sprintf(filename,"%s/soln_in_above_disk_region%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 region_id=Above_disk_region_id;
 n_el=Bulk_mesh_pt->nregion_element(region_id);
 for (unsigned e=0;e<n_el;e++)
  {
   if (do_bulk_output) 
    {
     Bulk_mesh_pt->region_element_pt(region_id,e)->output(some_file,nplot);
    }
   volume_in_above_disk_region+=Bulk_mesh_pt->
    region_element_pt(region_id,e)->size();
  }
 some_file.close();


 // Output bulk elements below disk
 //--------------------------------
 double volume_in_below_disk_region=0.0;
 sprintf(filename,"%s/soln_in_below_disk_region%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 region_id=Above_disk_region_id;
 n_el=Bulk_mesh_pt->nregion_element(region_id);
 for (unsigned e=0;e<n_el;e++)
  {
   if (do_bulk_output) 
    {
     Bulk_mesh_pt->region_element_pt(region_id,e)->output(some_file,nplot);
    }
   volume_in_below_disk_region+=Bulk_mesh_pt->
    region_element_pt(region_id,e)->size();
  }
 some_file.close();

 // Output bulk elements in region 0
 //--------------------------------- 
 double volume_in_region0=0.0;
 sprintf(filename,"%s/soln_in_zero_region%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 region_id=0;
 n_el=Bulk_mesh_pt->nregion_element(region_id);
 for (unsigned e=0;e<n_el;e++)
  {
   if (do_bulk_output) 
    {
     Bulk_mesh_pt->region_element_pt(region_id,e)->output(some_file,nplot);
    }
   volume_in_region0+= Bulk_mesh_pt->region_element_pt(region_id,e)->size();
   }
 some_file.close();


 // Get total mesh volume
 double total_mesh_volume=0.0;
 n_el=Bulk_mesh_pt->nelement();
 for (unsigned e=0;e<n_el;e++)
  {
   total_mesh_volume+=Bulk_mesh_pt->finite_element_pt(e)->size();
  }


 // Check volumes:
 oomph_info << "Error in total region volume balance: " <<
  abs(total_mesh_volume-(volume_in_above_disk_region+
                         volume_in_below_disk_region+
                         volume_in_torus_region+
                         volume_in_region0))/total_mesh_volume*100.0 
            << " % " << std::endl;

 oomph_info << "Error in above/below disk region volume balance: " <<
  abs(volume_in_above_disk_region-volume_in_below_disk_region)/
  volume_in_above_disk_region*100.0 << " % " << std::endl;

 
 // Attach face elements to boundary of torus
 //------------------------------------------
 sprintf(filename,"%s/face_elements_on_boundary_of_torus%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 double torus_surface_area=0.0;
 region_id=Torus_region_id;
 for (unsigned b=First_torus_boundary_id;
      b<=Last_torus_boundary_id;b++)
  {
   unsigned nel=Bulk_mesh_pt->nboundary_element_in_region(b,region_id);
   for (unsigned e=0;e<nel;e++)
    {
     FiniteElement* el_pt=
      Bulk_mesh_pt->boundary_element_in_region_pt(b,region_id,e);
     
     // What is the index of the face of the bulk element at the boundary
     int face_index = Bulk_mesh_pt->
      face_index_at_boundary_in_region(b,region_id,e);
     
     // Build the corresponding flux jump element
     PoissonFluxElement<ELEMENT>* flux_element_pt 
      = new PoissonFluxElement<ELEMENT>(el_pt,face_index);
     
     // Get surface area
     torus_surface_area+=flux_element_pt->size();
     
     // Output
     flux_element_pt->output(some_file);
     
     // ...and we're done!
     delete flux_element_pt;
    }
  }
 some_file.close();
 oomph_info << "Torus surface area: " <<  torus_surface_area << std::endl;

 
 // Attach face elements to part of disk inside torus
 //--------------------------------------------------
 sprintf(filename,"%s/face_elements_on_disk_in_torus%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 double disk_in_torus_surface_area=0.0;
 region_id=Torus_region_id;
 unsigned nb=One_based_boundary_id_for_disk_within_torus.size();
 for (unsigned i=0;i<nb;i++)
  {
   unsigned b=One_based_boundary_id_for_disk_within_torus[i]-1;
   unsigned nel=Bulk_mesh_pt->nboundary_element_in_region(b,region_id);
   for (unsigned e=0;e<nel;e++)
    {
     FiniteElement* el_pt=
      Bulk_mesh_pt->boundary_element_in_region_pt(b,region_id,e);
     
     // What is the index of the face of the bulk element at the boundary
     int face_index = Bulk_mesh_pt->
      face_index_at_boundary_in_region(b,region_id,e);
     
     // Build the corresponding flux jump element
     PoissonFluxElement<ELEMENT>* flux_element_pt 
      = new PoissonFluxElement<ELEMENT>(el_pt,face_index);
     
     // Get surface area
     disk_in_torus_surface_area+=flux_element_pt->size();
     
     // Output
     flux_element_pt->output(some_file);
     
     // ...and we're done!
     delete flux_element_pt;
    }
  }
 some_file.close();
 oomph_info << "Disk in torus surface area: "
            <<  disk_in_torus_surface_area << std::endl;

 
 // Attach face elements to part of disk outside torus
 //--------------------------------------------------
 sprintf(filename,"%s/face_elements_on_disk_outside_torus%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 double disk_outside_torus_surface_area=0.0;
 region_id=0; 
 nb=One_based_boundary_id_for_disk_outside_torus.size();
 for (unsigned i=0;i<nb;i++)
  {
   unsigned b=One_based_boundary_id_for_disk_outside_torus[i]-1;
   unsigned nel=Bulk_mesh_pt->nboundary_element_in_region(b,region_id);
   for (unsigned e=0;e<nel;e++)
    {
     FiniteElement* el_pt=
      Bulk_mesh_pt->boundary_element_in_region_pt(b,region_id,e);
     
     // What is the index of the face of the bulk element at the boundary
     int face_index = Bulk_mesh_pt->
      face_index_at_boundary_in_region(b,region_id,e);
     
     // Build the corresponding flux jump element
     PoissonFluxElement<ELEMENT>* flux_element_pt 
      = new PoissonFluxElement<ELEMENT>(el_pt,face_index);
     
     // Get surface area
     disk_outside_torus_surface_area+=flux_element_pt->size();
     
     // Output
     flux_element_pt->output(some_file);
     
     // ...and we're done!
     delete flux_element_pt;
    }
  }
 some_file.close();
 oomph_info << "Disk outside torus surface area: "
            <<  disk_outside_torus_surface_area << std::endl;

 oomph_info << "Total surface area of disk with torus: "
            <<  disk_in_torus_surface_area+disk_outside_torus_surface_area 
            << std::endl;


 // Attach face elements to free-standing disk
 //-------------------------------------------
 sprintf(filename,"%s/face_elements_on_free_standing_disk%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 double free_standing_disk_surface_area=0.0;
 region_id=0; 
 for (unsigned b=First_disk_boundary_id;
      b<Last_disk_boundary_id;b++)
  {
   unsigned nel=Bulk_mesh_pt->nboundary_element_in_region(b,region_id);
   for (unsigned e=0;e<nel;e++)
    {
     FiniteElement* el_pt=
      Bulk_mesh_pt->boundary_element_in_region_pt(b,region_id,e);
     
     // What is the index of the face of the bulk element at the boundary
     int face_index = Bulk_mesh_pt->
      face_index_at_boundary_in_region(b,region_id,e);
     
     // Build the corresponding flux jump element
     PoissonFluxElement<ELEMENT>* flux_element_pt 
      = new PoissonFluxElement<ELEMENT>(el_pt,face_index);
     
     // Get surface area
     free_standing_disk_surface_area+=flux_element_pt->size();
     
     // Output
     flux_element_pt->output(some_file);
     
     // ...and we're done!
     delete flux_element_pt;
    }
  }
 some_file.close();
 oomph_info << "Free-standing disk surface area: "
            << free_standing_disk_surface_area << std::endl;

 
 // Attach face elements to disk surrounded by layers: above
 //---------------------------------------------------------
 sprintf(filename,"%s/face_elements_on_disk_upper_layer%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 double disk_upper_layer_surface_area=0.0;
 region_id=Above_disk_region_id;
 for (unsigned b=First_disk_with_layers_boundary_id;
      b<=Last_disk_with_layers_boundary_id;b++)
  {
   unsigned nel=Bulk_mesh_pt->nboundary_element_in_region(b,region_id);
   for (unsigned e=0;e<nel;e++)
    {
     FiniteElement* el_pt=
      Bulk_mesh_pt->boundary_element_in_region_pt(b,region_id,e);
     
     // What is the index of the face of the bulk element at the boundary
     int face_index = Bulk_mesh_pt->
      face_index_at_boundary_in_region(b,region_id,e);
     
     // Build the corresponding flux jump element
     PoissonFluxElement<ELEMENT>* flux_element_pt 
      = new PoissonFluxElement<ELEMENT>(el_pt,face_index);
     
     // Get surface area
     disk_upper_layer_surface_area+=flux_element_pt->size();
     
     // Output
     flux_element_pt->output(some_file);
     
     // ...and we're done!
     delete flux_element_pt;
    }
  }
 some_file.close();
 oomph_info << "Disk upper layer surface area: "
            <<  disk_upper_layer_surface_area << std::endl;


 // Attach face elements to disk surrounded by layers: below
 //---------------------------------------------------------
 sprintf(filename,"%s/face_elements_on_disk_lower_layer%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 double disk_lower_layer_surface_area=0.0;
 region_id=Below_disk_region_id;
 for (unsigned b=First_disk_with_layers_boundary_id;
      b<=Last_disk_with_layers_boundary_id;b++)
  {
   unsigned nel=Bulk_mesh_pt->nboundary_element_in_region(b,region_id);
   for (unsigned e=0;e<nel;e++)
    {
     FiniteElement* el_pt=
      Bulk_mesh_pt->boundary_element_in_region_pt(b,region_id,e);
     
     // What is the index of the face of the bulk element at the boundary
     int face_index = Bulk_mesh_pt->
      face_index_at_boundary_in_region(b,region_id,e);
     
     // Build the corresponding flux jump element
     PoissonFluxElement<ELEMENT>* flux_element_pt 
      = new PoissonFluxElement<ELEMENT>(el_pt,face_index);
     
     // Get surface area
     disk_lower_layer_surface_area+=flux_element_pt->size();
     
     // Output
     flux_element_pt->output(some_file);
     
     // ...and we're done!
     delete flux_element_pt;
    }
  }
 some_file.close();
 oomph_info << "Disk lower layer surface area: "
            <<  disk_lower_layer_surface_area << std::endl;


 // Doc volumes and areas
 volumes_and_areas_file << volume_in_above_disk_region << " " 
                        << volume_in_below_disk_region << " " 
                        << volume_in_torus_region << " " 
                        << total_mesh_volume << " " 
                        << volume_in_region0 << " "
                        << torus_surface_area << " " 
                        << disk_in_torus_surface_area << " " 
                        << disk_outside_torus_surface_area << " " 
                        << disk_upper_layer_surface_area << " " 
                        << disk_lower_layer_surface_area << " " 
                        << free_standing_disk_surface_area << " " 
                        << std::endl;
 volumes_and_areas_file.close();


 // Output solution
 //----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 if (do_bulk_output) Bulk_mesh_pt->output(some_file,nplot);
 some_file.close();

 // Output solution showing element outlines
 //-----------------------------------------
 sprintf(filename,"%s/coarse_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 if (do_bulk_output) Bulk_mesh_pt->output(some_file,2);
 some_file.close();


 // Get norm of solution
 //---------------------
 sprintf(filename,"%s/norm%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 double norm_soln=0.0;
 Bulk_mesh_pt->compute_norm(norm_soln);  
 some_file << sqrt(norm_soln) << std::endl;
 oomph_info << "Norm of computed solution: "   << sqrt(norm_soln)  << endl;
 some_file.close();

} // end of doc




//========================================================================
/// Driver
//========================================================================
int main(int argc, char* argv[])
{

 MPI_Helpers::init(argc,argv);

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
  
 // length of downstream region occupied by impedance elements
 CommandLineArgs::specify_command_line_flag("--suppress_bulk_output");

#ifndef DO_TETGEN

 // Gmsh command line invocation
 CommandLineArgs::specify_command_line_flag
  ("--gmsh_command_line",
   &Global_Parameters::Gmsh_command_line_invocation);

#endif

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

#ifndef DO_TETGEN

 // Are you suicidal?
 if (!CommandLineArgs::command_line_flag_has_been_set("--gmsh_command_line"))
  {
   std::string error_msg
    ("You haven't specified how gmsh is invoked on the command line\n");
   error_msg += "on your computer, so I'll use the default\n\n" + 
    Global_Parameters::Gmsh_command_line_invocation
    + "\n\nwhich, unless you're mheil, is unlikely to work and I will "
    + "now die...\n";
   throw OomphLibError(error_msg, 
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
  }

#endif

 // Note that this can make tetgen die!
 //feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

 // Shut up prefix
 oomph_info.output_modifier_pt()=&default_output_modifier;

 // Label for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT");
  
 // Number of output points per edge
 unsigned nplot=5;

 // Build problem
 TetmeshPoissonProblem<ProjectablePoissonElement<
  TPoissonElement<3,3> > > problem;


 //Output initial guess
 problem.doc_solution(nplot,doc_info);
 doc_info.number()++;

 unsigned max_adapt=1; 
 for (unsigned i=0;i<=max_adapt;i++)
  {
   // Solve the bastard!
   problem.newton_solve();

   //Output solution
   problem.doc_solution(nplot,doc_info);

   // oomph_info << "tmp stop\n";
   // exit(0);

   //Increment counter for solutions 
   doc_info.number()++;

   if (i!=max_adapt)
    {
     problem.adapt();
    }
  }

}



