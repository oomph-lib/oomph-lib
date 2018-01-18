//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision$
//LIC//
//LIC// $LastChangedDate$
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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

//Generic routines
#include "generic.h"

// Navier--Stokes equations
#include "navier_stokes.h"

// Get the mesh
#include "meshes/tetgen_mesh.h"

//#include "../../../external_src/oomph_tetgen/tetgen.h"
#include "oomph_tetgen.h"

using namespace std;

using namespace oomph;

namespace Global_Parameters
{
 double Re = 0.0;
 double Visc_Ratio = 2.0;
 double Box_width = 3.0;
 double Box_length = 10.0;
}

namespace oomph
{


//=========================================================================
// Unstructured refineable Triangle Mesh 
//=========================================================================
template<class ELEMENT>
 class RefineableTetgenMesh : public virtual TetgenMesh<ELEMENT>,
  public virtual RefineableMeshBase
  {
   
    public:

   /// \short Build mesh, based on a TetgenMeshClosedSurface that specifies
   /// the outer boundary of the domain and any number of internal
   /// closed curves, also specified by TriangleMeshClosedSurfaces.
   /// Also specify target area for uniform element size.
   RefineableTetgenMesh(
    TetgenMeshFacetedSurface* const &outer_boundary_pt,
    Vector<TetgenMeshFacetedSurface*>& internal_closed_surface_pt,
    const double &element_volume,
    TimeStepper* time_stepper_pt=&Mesh::Default_TimeStepper,
    const bool &use_attributes=false) :
    TetgenMesh<ELEMENT>(outer_boundary_pt,
                        internal_closed_surface_pt,
                        element_volume,
                        time_stepper_pt,
                        use_attributes)
    {
     // Initialise the data associated with adaptation
      initialise_adaptation_data();
    }
   
   /// \short Build mesh from specified triangulation and
   /// associated target volumes for elements in it
   RefineableTetgenMesh(const Vector<double> &target_volume,
                        tetgenio* const &tetgen_io_pt,
                        TimeStepper* time_stepper_pt=
                        &Mesh::Default_TimeStepper,
                        const bool &use_attributes=false)  
    {
     // Initialise the data associated with adaptation
     initialise_adaptation_data();
     
     // Store Timestepper used to build elements
     this->Time_stepper_pt=time_stepper_pt;
     
     // Triangulation has been created -- remember to wipe it!
     this->Tetgenio_exists =true;
     this->Tetgenio_pt = new tetgenio;

     // Add the volume constraints to the tetgenio data object
     // which may be bad because it's actually modifying things in the base
     //mesh
     //Create a local copy
     tetgenio *tetgen_input_pt = new tetgenio;;
     this->deep_copy_of_tetgenio(tetgen_io_pt,tetgen_input_pt);
     //Add volume constraints
     tetgen_input_pt->tetrahedronvolumelist = 
      new double[tetgen_input_pt->numberoftetrahedra];
     for(int e=0;e<tetgen_input_pt->numberoftetrahedra;++e)
      {
       tetgen_input_pt->tetrahedronvolumelist[e] = target_volume[e];
      }
     
     // Input string for triangle
     std::stringstream input_string_stream;
     input_string_stream<< "Vqra"; 
     
     // Convert to a *char required by the triangulate function
     char tetswitches[100];
     sprintf(tetswitches,"%s",input_string_stream.str().c_str());
     
     // Build triangulateio refined object
     tetrahedralize(tetswitches, tetgen_input_pt, this->Tetgenio_pt);       
     // Build scaffold
     this->Tmp_mesh_pt=new TetgenScaffoldMesh(*this->Tetgenio_pt);

     // Convert mesh from scaffold to actual mesh
     this->build_from_scaffold(time_stepper_pt,use_attributes);
     
     // Kill the scaffold
     delete this->Tmp_mesh_pt;
     this->Tmp_mesh_pt=0;

     //delete the input
     delete tetgen_input_pt;

     // Setup boundary coordinates for boundaries
     //unsigned nb=nboundary();
     //for (unsigned b=0;b<nb;b++)
     // {
     //  this->setup_boundary_coordinates(b);
     // }       
    }   
   
   /// Empty Destructor
   virtual ~RefineableTetgenMesh() {}
   
   /// \short Problem pointer (needed for multi-domain machinery during
   /// adaptation)
   Problem*& problem_pt(){return Problem_pt;}
   
   /// Max element size allowed during adaptation
   double& max_element_size(){return Max_element_size;}
   
   /// Min element size allowed during adaptation
   double& min_element_size(){return Min_element_size;}
   
   /// Min angle before remesh gets triggered
   double& max_permitted_edge_ratio(){return Max_permitted_edge_ratio;}
   
   /// Doc the targets for mesh adaptation
   void doc_adaptivity_targets(std::ostream &outfile)
   {
    outfile << std::endl;
    outfile << "Targets for mesh adaptation: " << std::endl;
    outfile << "---------------------------- " << std::endl;
    outfile << "Target for max. error: " << Max_permitted_error << std::endl;
    outfile << "Target for min. error: " << Min_permitted_error << std::endl;
    outfile << "Target max edge ratio: " 
            << Max_permitted_edge_ratio << std::endl;
    outfile << "Min. allowed element size: " << Min_element_size << std::endl;
    outfile << "Max. allowed element size: " << Max_element_size << std::endl;
    outfile << "Don't unrefine if less than " << Max_keep_unrefined 
            << " elements need unrefinement." << std::endl;
    outfile << std::endl;
   }
   
   
   /// Refine mesh uniformly and doc process
   void refine_uniformly(DocInfo& doc_info)
   {
    throw OomphLibError("refine_uniformly() not implemented yet",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION); 
   }    
   
   
   /// \short Unrefine mesh uniformly: Return 0 for success,
   /// 1 for failure (if unrefinement has reached the coarsest permitted
   /// level)
   unsigned unrefine_uniformly()
   {
    throw OomphLibError("unrefine_uniformly() not implemented yet",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION); 
    // dummy return
    return 0;
   }
   
   /// Adapt mesh, based on elemental error provided
   void adapt(const Vector<double>& elem_error);
   
   
   //protected:
   
   /// \short Helper function that updates the input polygon's PSLG
   /// by using the end-points of elements from FaceMesh(es) that are
   /// constructed for the boundaries associated with the segments of the
   /// polygon.
   //void update_polygon_using_face_mesh(TriangleMeshPolygon* polygon_pt);
   
   /// \short Generate a new PSLG representation of the inner hole
   /// boundaries
   //virtual void surface_remesh_for_inner_hole_boundaries(
   // Vector<Vector<double> > &internal_point_coord);
   
   
   /// \short Generate a new PSLG representation of the outer boundary
   //virtual void surface_remesh_for_outer_boundary();
   
   
  /// \short Snap the boundary nodes onto any curvilinear boundaries
  //void snap_nodes_onto_boundary(RefineableTriangleMesh<ELEMENT>* &new_mesh_pt,
  //                              const unsigned &b);

   
   /// Helper function to initialise data associated with adaptation
   void initialise_adaptation_data()
   {
    // Set max and min targets for adaptation
    this->Max_element_size=1.0;
    this->Min_element_size=0.001;
    this->Max_permitted_edge_ratio=2.0;
    
    // Initialise problem pointer
    this->Problem_pt=0;
   }
   
   /// \short Build a new tetgenio object from previous TriangulateIO
   /// based on target area for each element
   //void refine_triangulateio(tetgenio& tetgen_io, 
   //                          const Vector<double> &target_volume,
   //                          tetgenio &tetgen_refine);
   

   /// \short Compute target volume based on the element's error and the
   /// error target; return max edge ratio
   double compute_volume_target(const Vector<double> &elem_error,
                                Vector<double> &target_volume)
    {
     double max_edge_ratio=0.0;
     unsigned count_unrefined=0;
     unsigned count_refined=0;
     this->Nrefinement_overruled=0;
     
     unsigned nel=this->nelement();
     for (unsigned e=0;e<nel;e++)
      {
       // Get element
       FiniteElement* el_pt=this->finite_element_pt(e);
       
       // Calculate the volume of the element
       double volume=el_pt->size();
       
       //Find the vertex coordinates
       // (vertices are enumerated first)
       double vertex[4][3];
       for(unsigned n=0;n<4;++n)
        {
         for(unsigned i=0;i<3;++i)
          {
           vertex[n][i] = el_pt->node_pt(n)->x(i);
          }
        }
       
       //Compute the radius of the circumsphere of the tetrahedron
       //Algorithm stolen from tetgen for consistency
       DenseDoubleMatrix A(3);
       for(unsigned i=0;i<3;++i)
        {
         A(0,i) = vertex[1][i] - vertex[0][i];
         A(1,i) = vertex[2][i] - vertex[0][i];
         A(2,i) = vertex[3][i] - vertex[0][i];
        }
       
       Vector<double> rhs(3);
       // Compute the right hand side vector b (3x1).
       for(unsigned i=0;i<3;++i)
        {
         rhs[i] = 0.0;
         for(unsigned k=0;k<3;++k)
          {
           rhs[i] += A(i,k)*A(i,k);
          }
         rhs[i] *= 0.5;
        }
       
       //Solve the linear system, in which the rhs is over-written with 
       //the solution
       A.solve(rhs);
       //Calculate the circum-radius
       double circum_radius = 
        sqrt(rhs[0]*rhs[0] + rhs[1]*rhs[1] + rhs[2]*rhs[2]);
       
       //Now find the shortest edge length
       Vector<double> edge(3);
       double min_length = DBL_MAX;
       for(unsigned start=0;start<4;++start)
        {
         for(unsigned end=start+1;end<4;++end)
          {
           for(unsigned i=0;i<3;++i)
            {
             edge[i] = vertex[start][i] - vertex[end][i];
            }
           double length = 
            sqrt(edge[0]*edge[0] + edge[1]*edge[1] + edge[2]*edge[2]);
           if(length < min_length) {min_length = length;}
          }
        }
       
       //Now calculate the minimum edge ratio for this element
       double local_max_edge_ratio = circum_radius/min_length;
       if(local_max_edge_ratio > max_edge_ratio)
        {
         max_edge_ratio = local_max_edge_ratio;
        }
       
       // Mimick refinement in tree-based procedure: Target volumes
       // for elements that exceed permitted error is 1/4 of their
       // current volume, corresponding to a uniform sub-division.
       if (elem_error[e] > this->max_permitted_error())
        {
         target_volume[e]=std::max(volume/4.0,Min_element_size);
         if (target_volume[e]!=Min_element_size)
          {
           count_refined++;
          }
         else
          {
           this->Nrefinement_overruled++;
          }
        }
       else if (elem_error[e] < this->min_permitted_error())
        {
         target_volume[e]=std::min(4.0*volume,Max_element_size);
         if (target_volume[e]!=Max_element_size)
          {
           count_unrefined++;
          }
        }
       else
        {
         // Leave it alone
         target_volume[e] = std::max(volume,Min_element_size); 
        }
       
      } //End of loop over elements
     
     // Tell everybody
     this->Nrefined=count_refined;
     this->Nunrefined=count_unrefined;
     
     return max_edge_ratio;
   }

   
   /// \short Problem pointer (needed for multi-domain machinery during
   /// adaptation
   Problem* Problem_pt;
   
   /// Max permitted element size
   double Max_element_size;
   
   /// Min permitted element size
   double Min_element_size;
   
   /// Max edge ratio before remesh gets triggered
   double Max_permitted_edge_ratio;
   
  }; 

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

//======================================================================
/// Adapt problem based on specified elemental error estimates
//======================================================================
template <class ELEMENT>
void RefineableTetgenMesh<ELEMENT>::adapt(const Vector<double>& elem_error)
 {    
  // Get refinement targets
  Vector<double> target_size(elem_error.size());
  double max_edge_ratio=compute_volume_target(elem_error,
                                              target_size);
  // Get maximum target volume
  unsigned n=target_size.size();
  double max_size=0.0;
  double min_size=DBL_MAX;
  for (unsigned e=0;e<n;e++)
   {
    if (target_size[e]>max_size) max_size=target_size[e];
    if (target_size[e]<min_size) min_size=target_size[e];
   }
  
  oomph_info << "Maximum target size: " << max_size << std::endl;
  oomph_info << "Minimum target size: " << min_size << std::endl;
  oomph_info << "Number of elements to be refined " 
             << this->Nrefined << std::endl;
  oomph_info << "Number of elements to be unrefined "
             << this->Nunrefined << std::endl;
  oomph_info << "Max edge ratio "<< max_edge_ratio << std::endl;

  double orig_max_size, orig_min_size;
  this->max_and_min_element_size(orig_max_size, orig_min_size);
  oomph_info << "Max/min element size in original mesh: " 
             << orig_max_size  << " "
             << orig_min_size << std::endl;    

  // Should we bother to adapt?
  if ( (Nrefined > 0) || (Nunrefined > this->max_keep_unrefined()) ||
       (max_edge_ratio > this->max_permitted_edge_ratio()) )
   {

    if (! ( (Nrefined > 0) || (Nunrefined > max_keep_unrefined()) ) )
     {
      oomph_info 
       << "Mesh regeneration triggered by edge ratio criterion\n";
     }

    //Generate a new 1D mesh representation of the inner hole boundaries
    //unsigned nhole=this->Internal_polygon_pt.size();
    //Vector<Vector<double> > internal_point_coord(nhole);
    //this->surface_remesh_for_inner_hole_boundaries(internal_point_coord);

    //Update the representation of the outer boundary
    //this->surface_remesh_for_outer_boundary();

    //If there is not a geometric object associated with the boundary
    //the reset the boundary coordinates so that the lengths are consistent
    //in the new mesh and the old mesh.
    //const  unsigned n_boundary = this->nboundary();
    //for(unsigned b=0;b<n_boundary;++b)
    // {
    //  if(this->boundary_geom_object_pt(b)==0)
    //   {
    //    this->setup_boundary_coordinates(b);
    //   }
    // }

    // Are we dealing with a solid mesh?
    SolidMesh* solid_mesh_pt=dynamic_cast<SolidMesh*>(this);

    // Build temporary uniform background mesh
    //----------------------------------------
    // with volume set by maximum required volume
    //---------------------------------------
    RefineableTetgenMesh<ELEMENT>* tmp_new_mesh_pt=0;
    /*  if (solid_mesh_pt!=0)
     {
      tmp_new_mesh_pt=new RefineableSolidTriangleMesh<ELEMENT>
       (closed_curve_pt,
        hole_pt,
        max_size,
        this->Time_stepper_pt,
        this->Use_attributes);
     }
     else*/
    {
     tmp_new_mesh_pt=new RefineableTetgenMesh<ELEMENT>
      (this->Outer_boundary_pt,
       this->Internal_surface_pt,
       max_size,
       this->Time_stepper_pt,
       this->Use_attributes);
    }



    // Snap to curvilinear boundaries (some code duplication as this
    // is repeated below but helper function would take so many
    // arguments that it's nearly as messy...
    
    //Pass the boundary  geometric objects to the new mesh
    //tmp_new_mesh_pt->boundary_geom_object_pt() = 
    // this->boundary_geom_object_pt();
    
    //Reset the boundary coordinates if there is
    //a geometric object associated with the boundary
    //tmp_new_mesh_pt->boundary_coordinate_limits() = 
    // this->boundary_coordinate_limits();
    //for (unsigned b=0;b<n_boundary;b++)
    // {
    //  if(tmp_new_mesh_pt->boundary_geom_object_pt(b)!=0)
    //   {
    //    tmp_new_mesh_pt->setup_boundary_coordinates(b);
    //   }
    // }
    
    //Output the mesh before any snapping takes place
    //tmp_new_mesh_pt->output("pre_mesh_nodes_snapped_0.dat");
    
    //Move the nodes on the new boundary onto the 
    //old curvilinear boundary
    //If the boundary is straight this will do precisely nothing
    //but will be somewhat inefficient
    //for(unsigned b=0;b<n_boundary;b++)
    // {
    //  this->snap_nodes_onto_boundary(tmp_new_mesh_pt,b);
    // }
    
    //Output the mesh after the snapping has taken place
    //tmp_new_mesh_pt->output("mesh_nodes_snapped_0.dat"); 
    
    // Get the tetgenio object associated with that mesh
    tetgenio *tmp_new_tetgenio_pt = tmp_new_mesh_pt->tetgenio_pt();

    
#ifdef PARANOID
    if (this->Problem_pt==0) 
     {
      throw OomphLibError("Problem pointer must be set with problem_pt()",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif

    RefineableTetgenMesh<ELEMENT>* new_mesh_pt=0;

    // Map storing target sizes for elements in temporary 
    // tetgenio mesh
    std::map<GeneralisedElement*,double> target_size_map;


    //////////////////////////////////////////////////////////////
    // NOTE: Repeated setup of multidomain interaction could
    // be avoided by setting up a sufficiently fine bin
    // for the original mesh and reading out the target
    // area information from there
    //////////////////////////////////////////////////////////////

    // Now start iterating to refine mesh recursively
    //-----------------------------------------------
    bool done=false;
    unsigned iter=0;
    while (!done)
     {
      
      // "Project" target volumes from current mesh onto uniform
      //------------------------------------------------------
      // background mesh
      //----------------
      
      // Temporarily switch on projection capabilities to allow
      // storage of pointer to external element.
      // Need to do this for both meshes to ensure that 
      // matching is done based on Eulerian coordinates for both
      // (in case we're dealing with solid meshes where the
      // locate_zeta would otherwise use the Lagrangian coordintes).
      unsigned nelem=this->nelement();
      for (unsigned e=0;e<nelem;e++)
       {
        dynamic_cast<ELEMENT*>(this->element_pt(e))->enable_projection();
       }
      unsigned nelem2=tmp_new_mesh_pt->nelement();
      for (unsigned e=0;e<nelem2;e++)
       {
        dynamic_cast<ELEMENT*>(tmp_new_mesh_pt->element_pt(e))->
         enable_projection();
       }

      // Set up multi domain interactions so we can figure out
      // which element in the intermediate uniform mesh is co-located
      // with given element in current mesh (which is to be refined)
      Multi_domain_functions::setup_multi_domain_interaction
       <ELEMENT>(this->Problem_pt,this,tmp_new_mesh_pt);
      
      target_size_map.clear();
      for (unsigned e=0;e<nelem;e++)
       {
        ELEMENT* el_pt=dynamic_cast<ELEMENT*>(this->element_pt(e));
        unsigned nint=el_pt->integral_pt()->nweight();
        for (unsigned ipt=0;ipt<nint;ipt++)
         {
          GeneralisedElement* ext_el_pt=el_pt->external_element_pt(0,ipt);

          // Use max. rather than min area of any element overlapping the
          // the current element, otherwise we get a rapid outward diffusion
          // of small elements
          target_size_map[ext_el_pt]=std::max(target_size_map[ext_el_pt],
                                              target_size[e]);
         }

        // Switch off projection capability          
        dynamic_cast<ELEMENT*>(this->element_pt(e))->disable_projection();
       }
      for (unsigned e=0;e<nelem2;e++)
       {
        dynamic_cast<ELEMENT*>(tmp_new_mesh_pt->element_pt(e))->
         disable_projection();
       }      

      // Now copy into target area for temporary mesh but limit to
      // the equivalent of one sub-division per iteration
      done=true;
      unsigned nel_new=tmp_new_mesh_pt->nelement();
      Vector<double> new_target_size(nel_new);
      for (unsigned e=0;e<nel_new;e++)
       {
        // No target area found for this element -- keep its size
        // by setting target area to -1 for tetrahedron
        double new_size=target_size_map[tmp_new_mesh_pt->element_pt(e)];
        if (new_size<=0.0) 
         {
          new_target_size[e]=-1.0; 
         }
        else 
         {
          // Limit target area to the equivalent of uniform
          // refinement during this stage of the iteration
          new_target_size[e]=new_size;
          if (new_target_size[e]<
              tmp_new_mesh_pt->finite_element_pt(e)->size()/4.0)
           {
            //ALH: It seems that tetgen "enlarges" the volume constraint
            //so this criterion can never be meet unless dividing by 1.2 
            //as well.
            new_target_size[e]=
             tmp_new_mesh_pt->finite_element_pt(e)->size()/4.0;
            //This is the tetgen adjustment
            new_target_size[e] /= 1.2;
         
            // We'll need to give it another go later
            done=false;
           }
         }
       }
      
      

      // Now create the new mesh from TriangulateIO structure
      //-----------------------------------------------------
      // associated with uniform background mesh and the
      //------------------------------------------------
      // associated target element sizes.
      //---------------------------------
      
      // Solid mesh?
      /*if (solid_mesh_pt!=0)
       {
        new_mesh_pt=new RefineableSolidTriangleMesh<ELEMENT>
         (new_target_area,
          tmp_new_triangulateio,
          this->Time_stepper_pt,
          this->Use_attributes);
       }      
      // No solid mesh
      else */
       { 
        new_mesh_pt=new RefineableTetgenMesh<ELEMENT>
         (new_target_size,
          tmp_new_tetgenio_pt,
          this->Time_stepper_pt,
          this->Use_attributes);
       }    
      
      
      // Snap to curvilinear boundaries (some code duplication as this
      // is repeated below but helper function would take so many
      // arguments that it's nearly as messy...
/*      
      //Pass the boundary  geometric objects to the new mesh 
      new_mesh_pt->boundary_geom_object_pt() = 
       this->boundary_geom_object_pt();
      
      
      // Reset the boundary coordinates if there is
      // a geometric object associated with the boundary
      new_mesh_pt->boundary_coordinate_limits() = 
       this->boundary_coordinate_limits();
      for (unsigned b=0;b<n_boundary;b++)
       {
        if(new_mesh_pt->boundary_geom_object_pt(b)!=0)
         {
          new_mesh_pt->setup_boundary_coordinates(b);
         }
       }
      
      //Output the mesh before any snapping takes place
      //new_mesh_pt->output("pre_mesh_nodes_snapped_1.dat"); 
      
      //Move the nodes on the new boundary onto the 
      //old curvilinear boundary
      //If the boundary is straight this will do precisely nothing
      //but will be somewhat inefficient
      for(unsigned b=0;b<n_boundary;b++)
       {
        this->snap_nodes_onto_boundary(new_mesh_pt,b);
       }
      
      //Output the mesh after the snapping has taken place
      //new_mesh_pt->output("mesh_nodes_snapped_1.dat"); 
      */
      
      // Not done: get ready for another iteration
      iter++;
      //Delete the temporary mesh
      delete tmp_new_mesh_pt;
      //Now transfer over the pointers
      if (!done)
       {
        tmp_new_mesh_pt=new_mesh_pt;
        tmp_new_tetgenio_pt=new_mesh_pt->tetgenio_pt();
       }
      
     } // end of iteration
    

    // Project current solution onto new mesh
    //---------------------------------------
    ProjectionProblem<ELEMENT>* project_problem_pt=
     new ProjectionProblem<ELEMENT>;
    project_problem_pt->mesh_pt()=new_mesh_pt;
    project_problem_pt->project(this);
    
    //this->output("pre_proj",5);
    //new_mesh_pt->output("post_proj.dat",5);
    
    //Flush the old mesh 
    unsigned nnod=nnode();
    for(unsigned j=nnod;j>0;j--)  
     { 
      delete Node_pt[j-1];  
      Node_pt[j-1] = 0; 
     } 
    unsigned nel=nelement(); 
    for(unsigned e=nel;e>0;e--)  
     { 
      delete Element_pt[e-1];  
      Element_pt[e-1] = 0; 
     } 
    
    // Now copy back to current mesh
    //------------------------------
    nnod=new_mesh_pt->nnode();
    Node_pt.resize(nnod);
    nel=new_mesh_pt->nelement();
    Element_pt.resize(nel);  
    for(unsigned j=0;j<nnod;j++)
     { 
      Node_pt[j] = new_mesh_pt->node_pt(j);
     } 
    for(unsigned e=0;e<nel;e++)
     { 
      Element_pt[e] = new_mesh_pt->element_pt(e);
     } 
    
    //Copy the boundary schemes
    unsigned nbound=new_mesh_pt->nboundary();
    Boundary_element_pt.resize(nbound);
    Face_index_at_boundary.resize(nbound);
    Boundary_node_pt.resize(nbound);
    for (unsigned b=0;b<nbound;b++)
     {
      unsigned nel=new_mesh_pt->nboundary_element(b);
      Boundary_element_pt[b].resize(nel);
      Face_index_at_boundary[b].resize(nel);
      for (unsigned e=0;e<nel;e++)
       {
        Boundary_element_pt[b][e]=new_mesh_pt->boundary_element_pt(b,e);
        Face_index_at_boundary[b][e]=new_mesh_pt->face_index_at_boundary(b,e);
       }
      unsigned nnod=new_mesh_pt->nboundary_node(b);
      Boundary_node_pt[b].resize(nnod);
      for (unsigned j=0;j<nnod;j++)
       {
        Boundary_node_pt[b][j]=new_mesh_pt->boundary_node_pt(b,j);
       }
     }

    //Also copy over the new boundary and region information
    unsigned n_region = new_mesh_pt->nregion();
    //Only bother if we have regions
    if(n_region > 1)
     {
      //Deal with the region information first
      this->Region_element_pt.resize(n_region);
      this->Region_attribute.resize(n_region);
      for(unsigned r=0;r<n_region;r++)
       {
        this->Region_attribute[r] = new_mesh_pt->region_attribute(r);
        //Find the number of elements in the region
        unsigned n_region_element = new_mesh_pt->nregion_element(r);
        this->Region_element_pt[r].resize(n_region_element);
        for(unsigned e=0;e<n_region_element;e++)
         {
          this->Region_element_pt[r][e] = new_mesh_pt->region_element_pt(r,e);
         }
       }

      //Now the boundary region information
      this->Boundary_region_element_pt.resize(nbound);
      this->Face_index_region_at_boundary.resize(nbound);
      
      //Now loop over the boundaries
      for(unsigned b=0;b<nbound;++b)
       {
        //Loop over the regions
        for(unsigned r=0;r<n_region;++r)
         {
          unsigned n_boundary_el_in_region = 
           new_mesh_pt->nboundary_element_in_region(b,r);
          
          if(n_boundary_el_in_region > 0)
           {
            //Allocate storage in the map
            this->Boundary_region_element_pt[b][r].
             resize(n_boundary_el_in_region);
            this->Face_index_region_at_boundary[b][r].
             resize(n_boundary_el_in_region);

            //Copy over the information
            for(unsigned e=0;e<n_boundary_el_in_region;++e)
             {
              this->Boundary_region_element_pt[b][r][e]
               = new_mesh_pt->boundary_element_in_region_pt(b,r,e);
              this->Face_index_region_at_boundary[b][r][e] 
               = new_mesh_pt->face_index_at_boundary_in_region(b,r,e);
             }
           }
         }
       } //End of loop over boundaries

     } //End of case when more than one region

    //Snap the newly created nodes onto any geometric objects
    //this->snap_nodes_onto_geometric_objects();

    // Copy the IDs of the vertex nodes
    //this->Oomph_vertex_nodes_id=new_mesh_pt->oomph_vertex_nodes_id();
    
    // Copy TriangulateIO representation
    //TriangleHelper::clear_triangulateio(this->Triangulateio);
    //bool quiet=true;
    //this->Triangulateio=
    // TriangleHelper::deep_copy_of_triangulateio_representation(
    //  new_mesh_pt->triangulateio_representation(),quiet);
    
    this->set_deep_copy_tetgenio_pt(new_mesh_pt->tetgenio_pt());
     //this->Tetgenio_pt = new_mesh_pt->tetgenio_pt();

    // Flush the mesh
    new_mesh_pt->flush_element_and_node_storage();
    
    // Delete the mesh and the problem
    delete new_mesh_pt;
    delete project_problem_pt;

    // Solid mesh?
    if (solid_mesh_pt!=0)
     {
      // Warning
      std::stringstream error_message;
      error_message 
       << "Lagrangian coordinates are currently not projected but are\n"
       << "are re-set during adaptation. This is not appropriate for\n"
       << "real solid mechanics problems!\n";
      OomphLibWarning(error_message.str(),
                      "RefineableTriangleMesh::adapt()",
                      OOMPH_EXCEPTION_LOCATION);
      
      // Reset Lagrangian coordinates
      dynamic_cast<SolidMesh*>(this)->set_lagrangian_nodal_coordinates();
     }
    
    double max_area;
    double min_area;
    this->max_and_min_element_size(max_area, min_area);
    oomph_info << "Max/min element size in adapted mesh: " 
               << max_area  << " "
               << min_area << std::endl;    
   }
  else
   {
    oomph_info << "Not enough benefit in adaptation.\n";
    Nrefined=0;
    Nunrefined=0;
   }
 }

}

//====================================================================
/// Micky mouse  problem.
//====================================================================
template<class ELEMENT> 
class FallingBlockProblem : public Problem
{

public:


 /// Constructor
 FallingBlockProblem();
  
 /// Destructor (empty)
 ~FallingBlockProblem()
  {
   //Delete the objects
   unsigned nh = Inner_boundary_pt.size();
   for(unsigned h=0;h<nh;++h)
    {
     delete Inner_boundary_pt[h];
    }
   delete Outer_boundary_pt;
  }

 /// Snap the boundary nodes onto the sphere
 void snap_onto_sphere()
  {
   this->mesh_pt()->output("pre_moved.dat",5);
   std::ofstream nodes("moved_nodes.dat");

   unsigned n_bound = mesh_pt()->nboundary_node(0);
   for(unsigned n=0;n<n_bound;++n)
    {
     Node* nod_pt = mesh_pt()->boundary_node_pt(0,n);
     double x = nod_pt->x(0);
     double y = nod_pt->x(1);
     double z = nod_pt->x(2);

     nodes << x << " " << y << " " << z << "  ";

     //Now let's snap by calculating the angle
     double r = sqrt(x*x + y*y + z*z);
     double theta = acos(z/r);
     double phi = atan2(y,x);
     
     nodes << r << " " << theta << " " << phi << " ";

     //Do the snapping
     double R_new = sqrt(1.0 + 0.5*0.5*(1.0 + sqrt(5.0))*(1.0 + sqrt(5.0)));
     nod_pt->x(0) = R_new*sin(theta)*cos(phi);
     nod_pt->x(1) = R_new*sin(theta)*sin(phi);
     nod_pt->x(2) = R_new*cos(theta);

     nodes << nod_pt->x(0) << " " << nod_pt->x(1) << " " << nod_pt->x(2) << "\n";
    }
   nodes.close();
   this->mesh_pt()->output("post_moved.dat",5);
  }
      
     
 /// Totally new mesh, need to fix it
 void actions_after_adapt()
  {
   // Set the boundary conditions for this problem 
   // Only on the "outer boundaries"
   unsigned num_bound = mesh_pt()->nboundary();
   for(unsigned ibound=0;ibound<num_bound;ibound++)
    {
     unsigned final_index = 3;
     //Do no pin the outlet z-velocity
     if(ibound==3) {final_index = 2;}
     unsigned num_nod= mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       for(unsigned i=0;i<final_index;++i)
        {
         mesh_pt()->boundary_node_pt(ibound,inod)->pin(i);
        }
      }
    }

   // Complete the build of all elements so they are fully functional
   
   //Find number of elements in mesh
   unsigned n_element = mesh_pt()->nelement();
   
   // Loop over the elements to set up element-specific 
   // things that cannot be handled by constructor
   for(unsigned i=0;i<n_element;i++)
    {
     // Upcast from GeneralElement to the present element
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
     
     //Set the source function pointer
     el_pt->re_pt() = &Global_Parameters::Re;
    }
  }
 


 /// Update the problem specs before solve: (Re)set boundary conditions
 void actions_before_newton_solve()
  {
   //Loop over the boundaries 
   unsigned num_bound = mesh_pt()->nboundary();
   for(unsigned ibound=0;ibound<num_bound;ibound++)
    {
     //Don't boundary 3's z-coordinates
     unsigned final_index = 3;
     if(ibound==3) {final_index = 2;}

     // Loop over the nodes on boundary
     unsigned num_nod=mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
       for(unsigned i=0;i<final_index;++i) {nod_pt->set_value(i,0.0);}
      }
    }

   //Now set a Poiseuille like inflow
   {
    using namespace Global_Parameters;
    
    // Loop over the nodes on boundary
    unsigned num_nod=mesh_pt()->nboundary_node(6);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
       Node* nod_pt=mesh_pt()->boundary_node_pt(6,inod);
       double x = nod_pt->x(0);
       double y = nod_pt->x(1);
       double u = (Box_width-x)*(Box_width+x)*(Box_width-y)*(Box_width+y);

       nod_pt->set_value(2,u);
     }
   }
  }

 /// Update the problem specs before solve (empty)
 void actions_after_newton_solve(){}


 //Access function for the specific mesh
RefineableTetgenMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableTetgenMesh<ELEMENT>*>(Problem::mesh_pt());
  }

  /// Doc the solution
 void doc_solution(const unsigned& nplot, DocInfo& doc_info);

 
 /// Calculate the fluid dissipation
 double get_dissipation()
  {
   double dissipation=0.0;
   const unsigned n_element = this->mesh_pt()->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     //Cast to a fluid element
     ELEMENT *el_pt = 
      dynamic_cast<ELEMENT*>(this->mesh_pt()->element_pt(e));
     //Add to the dissipation
     dissipation += el_pt->dissipation(); 
    }
   return dissipation;
  }

 /// Storage for the outer boundary object
 TetgenMeshFacetedSurface* Outer_boundary_pt;

 Vector<TetgenMeshFacetedSurface*> Inner_boundary_pt;

};



//========================================================================
/// Constructor for FallingBlock problem
//========================================================================
template<class ELEMENT>
FallingBlockProblem<ELEMENT>::FallingBlockProblem()
{ 
 //Let's have stupidly high tolerance
 //Newton_solver_tolerance = 1000;

 //Add a steady time stepper
 this->add_time_stepper_pt(new Steady<0>);

 //Make the external box
 const double box_width = Global_Parameters::Box_width;
 const double box_length = Global_Parameters::Box_length;
 Vector<Vector<double> > box_point(8);
 box_point[0].resize(3);
 box_point[0][0] = -box_width;
 box_point[0][1] = -box_width;
 box_point[0][2] = -box_length;

 box_point[1].resize(3);
 box_point[1][0] = -box_width;
 box_point[1][1] =  box_width;
 box_point[1][2] = -box_length;

 box_point[2].resize(3);
 box_point[2][0] = -box_width;
 box_point[2][1] =  box_width;
 box_point[2][2] =  box_length;

 box_point[3].resize(3);
 box_point[3][0] = -box_width;
 box_point[3][1] = -box_width;
 box_point[3][2] =  box_length;

 box_point[4].resize(3);
 box_point[4][0] =  box_width;
 box_point[4][1] = -box_width;
 box_point[4][2] = -box_length;

 box_point[5].resize(3);
 box_point[5][0] =  box_width;
 box_point[5][1] =  box_width;
 box_point[5][2] = -box_length;

 box_point[6].resize(3);
 box_point[6][0] =  box_width;
 box_point[6][1] =  box_width;
 box_point[6][2] =  box_length;

 box_point[7].resize(3);
 box_point[7][0] =  box_width;
 box_point[7][1] =  -box_width;
 box_point[7][2] =  box_length;

 Vector<Vector<unsigned> > box_facet(6);
 box_facet[0].resize(4);
 box_facet[0][0] = 0;
 box_facet[0][1] = 4;
 box_facet[0][2] = 7;
 box_facet[0][3] = 3;

 box_facet[1].resize(4);
 box_facet[1][0] = 4;
 box_facet[1][1] = 5;
 box_facet[1][2] = 6;
 box_facet[1][3] = 7;

 box_facet[2].resize(4);
 box_facet[2][0] = 3;
 box_facet[2][1] = 7;
 box_facet[2][2] = 6;
 box_facet[2][3] = 2;

 box_facet[3].resize(4);
 box_facet[3][0] = 6;
 box_facet[3][1] = 5;
 box_facet[3][2] = 1;
 box_facet[3][3] = 2;

 box_facet[4].resize(4);
 box_facet[4][0] = 3;
 box_facet[4][1] = 2;
 box_facet[4][2] = 1;
 box_facet[4][3] = 0;

 box_facet[5].resize(4);
 box_facet[5][0] = 5;
 box_facet[5][1] = 1;
 box_facet[5][2] = 0;
 box_facet[5][3] = 4;

 Vector<unsigned> box_facet_boundary_id(6);
 box_facet_boundary_id[0] = 2;
 box_facet_boundary_id[1] = 3;
 box_facet_boundary_id[2] = 4;
 box_facet_boundary_id[3] = 5;
 box_facet_boundary_id[4] = 6;
 box_facet_boundary_id[5] = 7;

 //Make the outer boundary object
 Outer_boundary_pt = 
  new TetgenMeshFacetedSurface(box_point,box_facet,box_facet_boundary_id);
 

  //Set basic icosahedron points
 Vector<Vector<double> > icosa_point(12);
 
 //Golden ratio
 const double phi = 0.5*(1.0 + sqrt(5.0));
 icosa_point[0].resize(3);
 icosa_point[0][0] = 0.0;
 icosa_point[0][1] = 1.0;
 icosa_point[0][2] = phi;
 icosa_point[1].resize(3);
 icosa_point[1][0] = 0.0;
 icosa_point[1][1] = -1.0;
 icosa_point[1][2] = phi;
 icosa_point[2].resize(3);
 icosa_point[2][0] = 0.0;
 icosa_point[2][1] = 1.0;
 icosa_point[2][2] = -phi;
 icosa_point[3].resize(3);
 icosa_point[3][0] = 0.0;
 icosa_point[3][1] = -1.0;
 icosa_point[3][2] = -phi;
 icosa_point[4].resize(3);
 icosa_point[4][0] = 1.0;
 icosa_point[4][1] = phi;
 icosa_point[4][2] = 0.0;
 icosa_point[5].resize(3);
 icosa_point[5][0] = -1.0;
 icosa_point[5][1] = phi;
 icosa_point[5][2] = 0.0;
 icosa_point[6].resize(3);
 icosa_point[6][0] = 1.0;
 icosa_point[6][1] = -phi;
 icosa_point[6][2] = 0.0;
 icosa_point[7].resize(3);
 icosa_point[7][0] = -1.0;
 icosa_point[7][1] = -phi;
 icosa_point[7][2] = 0.0;
 icosa_point[8].resize(3);
 icosa_point[8][0] = phi;
 icosa_point[8][1] = 0.0;
 icosa_point[8][2] = 1.0;
 icosa_point[9].resize(3);
 icosa_point[9][0] = phi;
 icosa_point[9][1] = 0.0;
 icosa_point[9][2] = -1.0;
 icosa_point[10].resize(3);
 icosa_point[10][0] = -phi;
 icosa_point[10][1] = 0.0;
 icosa_point[10][2] = 1.0;
 icosa_point[11].resize(3);
 icosa_point[11][0] = -phi;
 icosa_point[11][1] = 0.0;
 icosa_point[11][2] = -1.0;

 //Set up the connectivity
 Vector<Vector<unsigned> > icosa_facet(20);
 icosa_facet[0].resize(3);
 icosa_facet[0][0] = 0;
 icosa_facet[0][1] = 1;
 icosa_facet[0][2] = 8;
 
 icosa_facet[1].resize(3);
 icosa_facet[1][0] = 0;
 icosa_facet[1][1] = 10;
 icosa_facet[1][2] = 1;

 icosa_facet[2].resize(3);
 icosa_facet[2][0] = 0;
 icosa_facet[2][1] = 5;
 icosa_facet[2][2] = 10;

 icosa_facet[3].resize(3);
 icosa_facet[3][0] = 0;
 icosa_facet[3][1] = 4;
 icosa_facet[3][2] = 5;

 icosa_facet[4].resize(3);
 icosa_facet[4][0] = 0;
 icosa_facet[4][1] = 8;
 icosa_facet[4][2] = 4;

 icosa_facet[5].resize(3);
 icosa_facet[5][0] = 5;
 icosa_facet[5][1] = 11;
 icosa_facet[5][2] = 10;

 icosa_facet[6].resize(3);
 icosa_facet[6][0] = 5;
 icosa_facet[6][1] = 2;
 icosa_facet[6][2] = 11;

 icosa_facet[7].resize(3);
 icosa_facet[7][0] = 4;
 icosa_facet[7][1] = 2;
 icosa_facet[7][2] = 5;

 icosa_facet[8].resize(3);
 icosa_facet[8][0] = 4;
 icosa_facet[8][1] = 9;
 icosa_facet[8][2] = 2;

 icosa_facet[9].resize(3);
 icosa_facet[9][0] = 8;
 icosa_facet[9][1] = 9;
 icosa_facet[9][2] = 4;

 icosa_facet[10].resize(3);
 icosa_facet[10][0] = 6;
 icosa_facet[10][1] = 9;
 icosa_facet[10][2] = 8;

 icosa_facet[11].resize(3);
 icosa_facet[11][0] = 1;
 icosa_facet[11][1] = 6;
 icosa_facet[11][2] = 8;

 icosa_facet[12].resize(3);
 icosa_facet[12][0] = 1;
 icosa_facet[12][1] = 7;
 icosa_facet[12][2] = 6;

 icosa_facet[13].resize(3);
 icosa_facet[13][0] = 10;
 icosa_facet[13][1] = 7;
 icosa_facet[13][2] = 1;

 icosa_facet[14].resize(3);
 icosa_facet[14][0] = 10;
 icosa_facet[14][1] = 11;
 icosa_facet[14][2] = 7;

 icosa_facet[15].resize(3);
 icosa_facet[15][0] = 11;
 icosa_facet[15][1] = 3;
 icosa_facet[15][2] = 7;

 icosa_facet[16].resize(3);
 icosa_facet[16][0] = 7;
 icosa_facet[16][1] = 3;
 icosa_facet[16][2] = 6;

 icosa_facet[17].resize(3);
 icosa_facet[17][0] = 6;
 icosa_facet[17][1] = 3;
 icosa_facet[17][2] = 9;

 icosa_facet[18].resize(3);
 icosa_facet[18][0] = 9;
 icosa_facet[18][1] = 3;
 icosa_facet[18][2] = 2;

 icosa_facet[19].resize(3);
 icosa_facet[19][0] = 2;
 icosa_facet[19][1] = 3;
 icosa_facet[19][2] = 11;

 Vector<unsigned> icosa_facet_boundary_id(20,1);
 
 //Create the inner boundary object
 Inner_boundary_pt.resize(1);
 Inner_boundary_pt[0] = new TetgenMeshFacetedSurface(icosa_point,icosa_facet,
                                                     icosa_facet_boundary_id);

 Vector<double> inner_point(3,0.0);
 Inner_boundary_pt[0]->set_hole(inner_point);

 /*in.numberofregions = 1;
 in.regionlist = new double[5*in.numberofregions];
 in.regionlist[0] = 0.0;
 in.regionlist[1] = 0.0;
 in.regionlist[2] = 0.0;
 in.regionlist[3] = 1;
 in.regionlist[4] = 1;*/

 Problem::mesh_pt() = 
  new RefineableTetgenMesh<ELEMENT>(Outer_boundary_pt,
                                    Inner_boundary_pt,2.0,
                                    this->time_stepper_pt(),
                                    true);

 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;

 // Set targets for spatial adaptivity
 mesh_pt()->max_permitted_error()=0.005;
 mesh_pt()->min_permitted_error()=0.001; 
 mesh_pt()->max_element_size()=1.0;
 mesh_pt()->min_element_size()=0.001; 

 // Use coarser mesh during validation
 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   mesh_pt()->max_element_size()=2.0;
   mesh_pt()->min_element_size()=0.1; 
  }

 // Set the problem pointer
 mesh_pt()->problem_pt()=this;


 //Problem::mesh_pt() = new TetgenMesh<ELEMENT>(out,
 //                                             this->time_stepper_pt(),
 //                                             true);


 // Set the boundary conditions for this problem 
 // Only on the "outer boundaries"
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
 {
  unsigned final_index = 3;
  //Do no pin the outlet z-velocity
  if(ibound==3) {final_index = 2;}
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
   {
    for(unsigned i=0;i<final_index;++i)
     {
      mesh_pt()->boundary_node_pt(ibound,inod)->pin(i);
     }
   }
 }

 // Complete the build of all elements so they are fully functional

 //Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the source function pointer
   el_pt->re_pt() = &Global_Parameters::Re;
  }

 //Loop over the elements in region 1
 /*unsigned n_inner = mesh_pt()->nregion_element(1);
 for(unsigned e=0;e<n_inner;++e)
  {
   // Upcast from GeneralElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->region_element_pt(1,e));

   //Set the source function pointer
   el_pt->viscosity_ratio_pt() = &Global_Parameters::Visc_Ratio;
   }*/

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

}



//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void FallingBlockProblem<ELEMENT>::doc_solution(const unsigned& nplot,
                                           DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Doc local node numbering
 //-------------------------
 sprintf(filename,"%s/node_numbering%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 FiniteElement* el_pt=mesh_pt()->finite_element_pt(0);
 unsigned nnode=el_pt->nnode();
 unsigned ndim=el_pt->node_pt(0)->ndim();
 for (unsigned j=0;j<nnode;j++)
  {
   for (unsigned i=0;i<ndim;i++)
    {
     some_file << el_pt->node_pt(j)->x(i) << " " ;
    }
   some_file << j << std::endl;
  }
 some_file.close();

 // Output boundaries
 //------------------
 sprintf(filename,"%s/boundaries%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_boundaries(some_file);
 some_file.close();


 // Output solution
 //----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,nplot);
 some_file.close();



} // end of doc




//========================================================================
/// Demonstrate how to solve Poisson problem
//========================================================================
int main(int argc, char* argv[])
{
 // Label for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT");
  
 // Number of output points per edge
 unsigned nplot=2;

 //Output trace file
 std::ofstream dissipation("RESLT/diss.dat");

 // Do the problem with quadratic elements
 //---------------------------------------
 {
  FallingBlockProblem<ProjectableCrouzeixRaviartElement<
  TCrouzeixRaviartElement<3> > > problem;
  //problem.snap_onto_sphere();

  for(unsigned n=0;n<1;++n)
   {
    // Solve the problem
    problem.steady_newton_solve(1);
    
    //Output solution with 5 points per edge
    nplot=5;
    problem.doc_solution(nplot,doc_info);
    
    //Increment counter for solutions 
    doc_info.number()++;
    //Output the dissipation
    dissipation << " " << Global_Parameters::Re << "  "
                << problem.get_dissipation() << std::endl;

    //Increase the Reynolds number
    Global_Parameters::Re += 0.1;
   }

  dissipation.close();
 }


}



