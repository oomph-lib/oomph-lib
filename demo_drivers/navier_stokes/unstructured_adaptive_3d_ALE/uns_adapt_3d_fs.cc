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

//Generic routines
#include "generic.h"

//The required governing equations
#include "navier_stokes.h"
#include "solid.h"
#include "constitutive.h"
#include "fluid_interface.h"

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
 double Box_length = 3.0;

 /// Set the Capillary number
 double Ca = 1.0;
 
 /// Pseudo-solid Poisson ratio
 double Nu=0.3;

 /// Set bubble radius
 double Radius = 1.0;

 /// Volume is negative becuase the volume is enclosed by the bulk
 /// fluid
 double Volume = -16.0*atan(1.0)*Radius*Radius*Radius/3.0;

 /// Gravity
 Vector<double> G(3);

 /// The Reynolds INverse Froude number
 double ReInvFr = 0.0;
 
 /// Constitutive law used to determine the mesh deformation
 ConstitutiveLaw* Constitutive_law_pt= new GeneralisedHookean(&Nu);
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

   ///  Build mesh, based on a TetgenMeshClosedSurface that specifies
   /// the outer boundary of the domain and any number of internal
   /// closed curves, also specified by TriangleMeshClosedSurfaces.
   /// Also specify target area for uniform element size.
   RefineableTetgenMesh(
    TetMeshFacetedClosedSurface* const &outer_boundary_pt,
    Vector<TetMeshFacetedSurface*>& internal_closed_surface_pt,
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
   
   ///  Build mesh from specified triangulation and
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
   
   ///  Problem pointer (needed for multi-domain machinery during
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
   
   
   ///  Unrefine mesh uniformly: Return 0 for success,
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

   ///Overload the standard setup of boundary coordinates to use area
   ///coordinates for triangles
   void setup_boundary_coordinates_generic(const unsigned& b,
                                           const bool& switch_normal,
                                           std::ofstream& outfile)
    {
     std::cout << "Calling this local version\n";
    }



   
   ///  Helper function that updates the input faceted surface
   /// by using the flattened elements from FaceMesh(es) that are
   /// constructed for the boundaries associated with the segments of the
   /// polygon.
   void update_faceted_surface_using_face_mesh(
    TetMeshFacetedClosedSurface* faceted_surface_pt)
    {
     //The easiest thing to do is simply to update the
     //positions of the key control nodes, leaving the connectivity alone,
     //but that doesn't allow for any surface remeshing
     
     ///List of vertex nodes
     std::list<Node*> new_vertex_nodes;
     //List of facets and boundary ids
     std::list<std::pair<Vector<unsigned>, unsigned>  > new_facet_list;

     //Loop over the number of old facets
     unsigned n_old_facet = faceted_surface_pt->nfacet();
     for(unsigned f=0;f<n_old_facet;f++)
      {
       //Get the boundary id of the facet. Need to subtract one, 
       //which is confusing now I think about it.
       //ALH: Should fix this.
       unsigned bound = faceted_surface_pt->one_based_facet_boundary_id(f) -1;
       
       // Create a face mesh adjacent to the fluid mesh's bound-th boundary. 
       // The face mesh consists of FaceElements that may also be 
       // interpreted as GeomObjects
       Mesh* face_mesh_pt = new Mesh;
       this->template build_face_mesh
        <ELEMENT,FaceElementAsGeomObject>(bound, face_mesh_pt);
       
       // Loop over these new face elements and tell them the boundary number
       // from the bulk fluid mesh -- this is required to they can
       // get access to the boundary coordinates!
       unsigned n_face_element = face_mesh_pt->nelement();
       for(unsigned e=0;e<n_face_element;e++)
        {
         //Cast the element pointer to the correct thing!
         FaceElementAsGeomObject<ELEMENT>* el_pt=
          dynamic_cast<FaceElementAsGeomObject<ELEMENT>*>
          (face_mesh_pt->element_pt(e));
         
         // Set bulk boundary number
         el_pt->set_boundary_number_in_bulk_mesh(bound);
        }

       //In order to set up the new faceted representation
       //Need to know the positions of the corner nodes
       //and the connectivity

       //Storage for the connectivity information
       Vector<unsigned> new_local_facet(3);

       //Now we have the face mesh loop over the face elements and 
       //print out the end points
       for(unsigned e=0;e<n_face_element;++e)
        {
         //Cache pointer to the element
         FiniteElement const *elem_pt = face_mesh_pt->finite_element_pt(e);

         //Just use the three primary (corner) nodes to define the facet
         unsigned n_vertex_node = 3;
         for(unsigned n=0;n<n_vertex_node;n++)
          {
           //Cache the pointer to the node
           Node* const nod_pt = elem_pt->node_pt(n);
           //If the vertex node is in the list return the number
           unsigned counter=0; bool found_it=false;
           for(std::list<Node*>::iterator it = new_vertex_nodes.begin();
               it!=new_vertex_nodes.end();++it,++counter)
            {
             //If we have found the node then store it
             if(*it==nod_pt) 
              {
               new_local_facet[n] = counter; 
               found_it=true; 
               break;
              }
            }

           //If we haven't found it
           //then add the node to the list and fill in the counter
           if(!found_it)
            {
             new_vertex_nodes.push_back(nod_pt);
             //We know where it is
             new_local_facet[n] = counter;
            }
          }

         //Add the new facet connectivity to the list
         new_facet_list.push_back(std::make_pair(new_local_facet,bound+1));
        }
      } //End of loop over old facets

     
     //Probably want the surface mesh in a better data structure so
     //that we can perform refinement or coarsening on it
     //i.e. want neighbours, edge flips all that fun stuff
     //That will go here!
     
     //Now sort out the facet nodes
     unsigned n_facet_vertex = new_vertex_nodes.size();
     Vector<Vector<double> > facet_point(n_facet_vertex);
     unsigned count=0;
     for(std::list<Node*>::iterator it=new_vertex_nodes.begin();
         it!=new_vertex_nodes.end();++it)
      {
       facet_point[count].resize(3);
       for(unsigned i=0;i<3;i++)
        {
         facet_point[count][i] = (*it)->x(i);
        }
       ++count;
      }
     
     //And also the facets
     unsigned n_facet = new_facet_list.size();
     Vector<Vector<unsigned> > new_facet(n_facet);
     Vector<unsigned> new_facet_boundary_id(n_facet);
     count=0;
     for(std::list<std::pair<Vector<unsigned>, unsigned> >::iterator 
          it=new_facet_list.begin();it!=new_facet_list.end();++it)
      {
       new_facet[count] = (*it).first;
       new_facet_boundary_id[count] = (*it).second;
       ++count;
      }

     for(unsigned f=0;f<n_facet;f++)
      {
       for(unsigned i=0;i<3;i++)
        {
         oomph_info << new_facet[f][i] << " " ;
        }
       oomph_info << " : ";
       oomph_info << new_facet_boundary_id[f] << "\n";
      }

     //Now just create the new boundary
     delete faceted_surface_pt;
     faceted_surface_pt = new TetMeshFacetedClosedSurface(
      facet_point,new_facet,new_facet_boundary_id);


     //Take average to get a new hole position (Won't always work)
     Vector<double> inner_point(3,0.0);
     for(unsigned n=0;n<n_facet_vertex;n++)
      {
       for(unsigned i=0;i<3;i++)
        {
         inner_point[i] += facet_point[n][i];
        }
      }

     for(unsigned i=0;i<3;i++) {inner_point[i] /= n_facet_vertex;}

     //Now set the hole if required
     faceted_surface_pt->set_hole(inner_point);

    }


   ///  Generate a new faceted representation of the inner hole
   /// boundaries
   virtual void surface_remesh_for_inner_hole_boundaries()
    {
     //Loop over the number of internal boundarys
     unsigned n_hole = this->Internal_surface_pt.size();
     for(unsigned ihole=0;ihole<n_hole;ihole++)
      {
       //Cache the pointer to the representation
       TetgenMeshClosedFacetedSurface* const faceted_surface_pt =
        dynamic_cast<TetgenMeshClosedFacetedSurface*>(this->Internal_surface_pt[ihole]);

       //Now can the surface update its own representation goes in here
       
       //If not we have to generate it from the new mesh
       {
        //Update the faceted surface associated with the ihole-th hole
        this->update_faceted_surface_using_face_mesh(faceted_surface_pt);
        
        //Then sort out the hole coordinates
       }
      }
    }
   
   
   ///  Generate a new PSLG representation of the outer boundary
   //virtual void surface_remesh_for_outer_boundary();
   
   
   ///  Snap the boundary nodes onto any curvilinear boundaries
   void snap_nodes_onto_boundary(RefineableTetgenMesh<ELEMENT>* &new_mesh_pt,
                                 const unsigned &b)
    {

     // Quick return
     if (!Boundary_coordinate_exists[b])
      {
       oomph_info << "Not snapping nodes on boundary " << b 
                  << " because no boundary coordinate has been set up.\n";
       return;
      }
     
     //Firstly we set the boundary coordinates of the new nodes
     //In case the mapping between the geometric object's intrinsic coordiante
     //and the arc-length coordinate is nonlinear. 
     //This is only an approximation, 
     //but it will ensure that the nodes that were input to triangle will
     //retain exactly the same boundary coordinates and 
     //then linear interpolation
     //is used between those values for any newly created nodes.

     
     // Create a face mesh adjacent to the fluid mesh's b-th boundary. 
     // The face mesh consists of FaceElements that may also be 
     // interpreted as GeomObjects
     Mesh* face_mesh_pt = new Mesh;
     this->template build_face_mesh
      <ELEMENT,FaceElementAsGeomObject>(b, face_mesh_pt);
     
     // Loop over these new face elements and tell them the boundary number
     // from the bulk fluid mesh -- this is required to they can
     // get access to the boundary coordinates!
     unsigned n_face_element = face_mesh_pt->nelement();
     for(unsigned e=0;e<n_face_element;e++)
      {
       //Cast the element pointer to the correct thing!
       FaceElementAsGeomObject<ELEMENT>* el_pt=
        dynamic_cast<FaceElementAsGeomObject<ELEMENT>*>
        (face_mesh_pt->element_pt(e));
       
       // Set bulk boundary number
       el_pt->set_boundary_number_in_bulk_mesh(b);
      }   
     
     //Now make the mesh as geometric object
     MeshAsGeomObject* mesh_geom_obj_pt = new MeshAsGeomObject(face_mesh_pt);
     
     //Now assign the new nodes positions based on the old meshes
     //potentially curvilinear boundary (its geom object incarnation)
     Vector<double> new_x(3);
     Vector<double> b_coord(2);
     unsigned n_new_boundary_node = new_mesh_pt->nboundary_node(b);
     for(unsigned n=0;n<n_new_boundary_node;n++)
      {
       //Get the boundary coordinate of all new nodes
       Node* const nod_pt = new_mesh_pt->boundary_node_pt(b,n);
       nod_pt->get_coordinates_on_boundary(b,b_coord);
       //Let's find boundary coordinates of the new node
       mesh_geom_obj_pt->position(b_coord,new_x);
       //Now snap to the boundary
       for(unsigned i=0;i<3;i++)
        {
         nod_pt->x(i) = new_x[i];
        }
      }

     //Delete the allocated memory for the geometric object and face mesh
     delete mesh_geom_obj_pt;
     face_mesh_pt->flush_element_and_node_storage();
     delete face_mesh_pt;
     
     //Loop over the elements adjacent to the boundary
     unsigned n_element = new_mesh_pt->nboundary_element(b);
     if(n_element > 0)
      {
       //Make a dummy simplex element
       TElement<3,2> dummy_four_node_element;
       //Make a dummy quadratic element
       TElement<3,3> dummy_ten_node_element;
       Vector<double> s(3);
       Vector<double> x_new(3);
       for(unsigned n=0;n<4;n++) {dummy_four_node_element.construct_node(n);}
       for(unsigned n=0;n<10;n++) {dummy_ten_node_element.construct_node(n);}
       for(unsigned e=0;e<n_element;e++)
        {
         //Cache the element pointer
         ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(
          new_mesh_pt->boundary_element_pt(b,e));
         
         //Find the number of nodes
         const unsigned n_node = elem_pt->nnode();
         //Only do something if not simplex element
         if(n_node > 4)
          {
           //Copy the nodes into the dummy
           for(unsigned n=0;n<4;n++)
            {
             for(unsigned i=0;i<3;i++)
              {
               dummy_four_node_element.node_pt(n)->x(i) = 
                elem_pt->node_pt(n)->x(i);
              }
            }
           
           //Now sort out the mid-side nodes
           for(unsigned n=4;n<10;n++)
            {
             //If it's not on a boundary then reset to be interpolated
             //from the simplex
             if(!elem_pt->node_pt(n)->is_on_boundary())
              {
               elem_pt->local_coordinate_of_node(n,s);
               dummy_four_node_element.interpolated_x(s,x_new);
               for(unsigned i=0;i<3;i++) {elem_pt->node_pt(n)->x(i) = x_new[i];}
              }
            }
          }
         
         //If we have more than 10 nodes interpolate from the quadratic shape
         if(n_node > 10)
          {
           //Copy the nodes into the dummy
           for(unsigned n=0;n<10;n++)
            {
             for(unsigned i=0;i<3;i++)
              {
               dummy_ten_node_element.node_pt(n)->x(i) = 
                elem_pt->node_pt(n)->x(i);
              }
            }
           
           //Now sort out the mid-face and central nodes
           for(unsigned n=10;n<n_node;n++)
            {
             //If it's not on a boundary then reset to be interpolated
             //from the simplex
             if(!elem_pt->node_pt(n)->is_on_boundary())
              {
               elem_pt->local_coordinate_of_node(n,s);
               dummy_ten_node_element.interpolated_x(s,x_new);
               for(unsigned i=0;i<3;i++) {elem_pt->node_pt(n)->x(i) = x_new[i];}
              }
            }
          }
        }
      } //End of fix up of elements
    }
   
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
   
   ///  Build a new tetgenio object from previous TriangulateIO
   /// based on target area for each element
   //void refine_triangulateio(tetgenio& tetgen_io, 
   //                          const Vector<double> &target_volume,
   //                          tetgenio &tetgen_refine);
   

   ///  Compute target volume based on the element's error and the
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

   
   ///  Problem pointer (needed for multi-domain machinery during
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
    this->surface_remesh_for_inner_hole_boundaries();

    //Update the representation of the outer boundary
    //this->surface_remesh_for_outer_boundary();

    //If there is not a geometric object associated with the boundary
    //the reset the boundary coordinates so that the lengths are consistent
    //in the new mesh and the old mesh.
    const  unsigned n_boundary = this->nboundary();
    for(unsigned b=0;b<n_boundary;++b)
     {
      //if(this->boundary_geom_object_pt(b)==0)
       {
        this->setup_boundary_coordinates(b);
       }
     }

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
    tmp_new_mesh_pt->output("pre_mesh_nodes_snapped_0.dat");
    
    //Move the nodes on the new boundary onto the 
    //old curvilinear boundary
    //If the boundary is straight this will do precisely nothing
    //but will be somewhat inefficient
    for(unsigned b=0;b<n_boundary;b++)
     {
      this->snap_nodes_onto_boundary(tmp_new_mesh_pt,b);
     }
    
    //Output the mesh after the snapping has taken place
    tmp_new_mesh_pt->output("mesh_nodes_snapped_0.dat"); 
    
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
      
      //Pass the boundary  geometric objects to the new mesh 
      //new_mesh_pt->boundary_geom_object_pt() = 
      // this->boundary_geom_object_pt();
      
      
      // Reset the boundary coordinates if there is
      // a geometric object associated with the boundary
      //new_mesh_pt->boundary_coordinate_limits() = 
      // this->boundary_coordinate_limits();
      //for (unsigned b=0;b<n_boundary;b++)
      // {
      //  if(new_mesh_pt->boundary_geom_object_pt(b)!=0)
      //   {
      //    new_mesh_pt->setup_boundary_coordinates(b);
      //   }
      // }
      
      //Output the mesh before any snapping takes place
      new_mesh_pt->output("pre_mesh_nodes_snapped_1.dat"); 
      
      //Move the nodes on the new boundary onto the 
      //old curvilinear boundary
      //If the boundary is straight this will do precisely nothing
      //but will be somewhat inefficient
      for(unsigned b=0;b<n_boundary;b++)
       {
        this->snap_nodes_onto_boundary(new_mesh_pt,b);
       }
      
      //Output the mesh after the snapping has taken place
      new_mesh_pt->output("mesh_nodes_snapped_1.dat"); 
      
      
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

//=========================================================================
// Unstructured refineable Tetgen Mesh upgraded to solid mesh
//=========================================================================
 template<class ELEMENT>
  class RefineableSolidTetgenMesh : 
  public virtual RefineableTetgenMesh<ELEMENT>,
public virtual SolidMesh
 {
  
   public:
  
  ///  Build mesh, based on closed curve that specifies
  /// the outer boundary of the domain and any number of internal
  /// closed curves. Specify target area for uniform element size.
  RefineableSolidTetgenMesh(
   TetgenMeshOLDFacetedSurface* const &outer_boundary_pt,
   Vector<TetgenMeshOLDFacetedSurface*>& internal_closed_surface_pt,
   const double &element_volume,
   TimeStepper* time_stepper_pt=&Mesh::Default_TimeStepper,
   const bool &use_attributes=false) :
   TetgenMesh<ELEMENT>(outer_boundary_pt,
                       internal_closed_surface_pt,
                       element_volume,
                       time_stepper_pt,
                       use_attributes),
   RefineableTetgenMesh<ELEMENT>(outer_boundary_pt,
                                 internal_closed_surface_pt,
                                 element_volume,
                                 time_stepper_pt,
                                 use_attributes)
   
   {
    //Assign the Lagrangian coordinates
    set_lagrangian_nodal_coordinates();
   }
  
     
   ///  Build mesh from specified triangulation and
   /// associated target areas for elements in it.
   RefineableSolidTetgenMesh(const Vector<double> &target_volume,
                             tetgenio* const &tetgen_io,
                             TimeStepper* time_stepper_pt=
                             &Mesh::Default_TimeStepper,
                             const bool &use_attributes=false)  :
    RefineableTetgenMesh<ELEMENT>(target_volume,
                                  tetgen_io,
                                  time_stepper_pt,
                                  use_attributes)
   {
    //Assign the Lagrangian coordinates
    set_lagrangian_nodal_coordinates();
   }
  
  /// Empty Destructor
  virtual ~RefineableSolidTetgenMesh() {}
  
 };

}

//====================================================================
/// Micky mouse  problem.
//====================================================================
template<class ELEMENT> 
class RisingBubbleProblem : public Problem
{

public:


 /// Constructor
 RisingBubbleProblem();
  
 /// Destructor (empty)
 ~RisingBubbleProblem()
  {
   //Delete the objects
   unsigned nh = Inner_boundary_pt.size();
   for(unsigned h=0;h<nh;++h)
    {
     delete Inner_boundary_pt[h];
    }
   delete Outer_boundary_pt;
  }

  ///  Create free surface elements
 void create_free_surface_elements()
  {  
   ///  Volume constraint element stores the Data item that stores
   /// the bubble pressure that is adjusted/traded to allow for
   /// volume conservation. Which value is the pressure stored in?
   unsigned p_traded_index=Vol_constraint_el_pt->index_of_traded_pressure();
   
   //Loop over the free surface boundaries
   //which is actually quite a large number
   unsigned nb=Fluid_mesh_pt->nboundary();
   for(unsigned b=6;b<nb;b++)
    {
     // How many bulk fluid elements are adjacent to boundary b?
     unsigned n_element = Fluid_mesh_pt->nboundary_element(b);
   
     // Loop over the bulk fluid elements adjacent to boundary b?
     for(unsigned e=0;e<n_element;e++)
      {
       // Get pointer to the bulk fluid element that is 
     // adjacent to boundary b
       ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
        Fluid_mesh_pt->boundary_element_pt(b,e));
       
       //Find the index of the face of element e along boundary b
       int face_index = Fluid_mesh_pt->face_index_at_boundary(b,e);
       
       // Create new element
       ElasticSurfaceFluidInterfaceElement<ELEMENT>* el_pt =
        new ElasticSurfaceFluidInterfaceElement<ELEMENT>(
         bulk_elem_pt,face_index);   
       
       // Add it to the mesh
       Free_surface_mesh_pt->add_element_pt(el_pt);
       
       //Add the appropriate boundary number
       el_pt->set_boundary_number_in_bulk_mesh(b);
       
       //Specify the capillary number
       el_pt->ca_pt() = &Global_Parameters::Ca;
       
       // Specify the bubble pressure (pointer to Data object and 
       // index of value within that Data object that corresponds
       // to the traded pressure
       el_pt->set_external_pressure_data(
        Vol_constraint_el_pt->p_traded_data_pt(),p_traded_index); 
      } 
    }
  }

 ///  Delete free surface elements 
 void delete_free_surface_elements()
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Free_surface_mesh_pt->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Free_surface_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Free_surface_mesh_pt->flush_element_and_node_storage();
   
  } // end of delete_free_surface_elements

 
 /// Totally new mesh, need to fix it
 void actions_after_adapt()
  {
   // Set the boundary conditions for this problem 
   // Only on the "outer boundaries"
   for(unsigned ibound=0;ibound<6;ibound++)
    {
     unsigned final_index = 3;
     //Do no pin the outlet z-velocity
     if(ibound==2) {final_index = 2;}
     unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       for(unsigned i=0;i<final_index;++i)
        {
         Fluid_mesh_pt->boundary_node_pt(ibound,inod)->pin(i);
        }
      }
    }

   // Complete the build of all elements so they are fully functional
   
   //Find number of elements in mesh
   unsigned n_element = Fluid_mesh_pt->nelement();
   
   // Loop over the elements to set up element-specific 
   // things that cannot be handled by constructor
   for(unsigned i=0;i<n_element;i++)
    {
     // Upcast from GeneralElement to the present element
     ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(i));
     
     //Set the source function pointer
     el_pt->re_pt() = &Global_Parameters::Re;
     
     //Set the gravity term
     el_pt->re_invfr_pt() = &Global_Parameters::ReInvFr;

     el_pt->g_pt() = &Global_Parameters::G;

     //Set the contitutive law
     el_pt->constitutive_law_pt() = Global_Parameters::Constitutive_law_pt;
    }
  }
 


 /// Update the problem specs before solve: (Re)set boundary conditions
 void actions_before_newton_solve()
  {
   //Loop over the boundaries 
   for(unsigned ibound=0;ibound<6;ibound++)
    {
     //Don't boundary 3's z-coordinates
     unsigned final_index = 3;
     if(ibound==2) {final_index = 2;}

     // Loop over the nodes on boundary
     unsigned num_nod=Fluid_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       Node* nod_pt=Fluid_mesh_pt->boundary_node_pt(ibound,inod);
       for(unsigned i=0;i<final_index;++i) {nod_pt->set_value(i,0.0);}
      }
    }

   //Now set a Poiseuille like inflow
   {
    using namespace Global_Parameters;
    
    // Loop over the nodes on boundary
    unsigned num_nod=Fluid_mesh_pt->nboundary_node(5);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      //Node* nod_pt=Fluid_mesh_pt->boundary_node_pt(5,inod);
       //double x = nod_pt->x(0);
       //double y = nod_pt->x(1);
       //double u = (Box_width-x)*(Box_width+x)*(Box_width-y)*(Box_width+y);
       // nod_pt->set_value(2,u);
     }
   }
  }

 /// Update the problem specs before solve (empty)
 void actions_after_newton_solve(){}

  /// Doc the solution
 void doc_solution(const unsigned& nplot, DocInfo& doc_info);


/// Create elements that impose volume constraint on the bubble
 void create_volume_constraint_elements()
  {
   // Add volume constraint element to the mesh
   Volume_constraint_mesh_pt->add_element_pt(Vol_constraint_el_pt);
   
   //Loop over the free surface boundaries
   unsigned nb=Fluid_mesh_pt->nboundary();
   for(unsigned b=6;b<nb;b++)
    {
     // How many bulk fluid elements are adjacent to boundary b?
     unsigned n_element = Fluid_mesh_pt->nboundary_element(b);
     
     // Loop over the bulk fluid elements adjacent to boundary b?
     for(unsigned e=0;e<n_element;e++)
      {
       // Get pointer to the bulk fluid element that is 
       // adjacent to boundary b
       ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
        Fluid_mesh_pt->boundary_element_pt(b,e));
       
       //Find the index of the face of element e along boundary b
       int face_index = Fluid_mesh_pt->face_index_at_boundary(b,e);
       
       // Create new element
       ElasticSurfaceVolumeConstraintBoundingElement<ELEMENT>* el_pt =
        new ElasticSurfaceVolumeConstraintBoundingElement<ELEMENT>(
         bulk_elem_pt,face_index);

       //Set the "master" volume control element
       el_pt->set_volume_constraint_element(Vol_constraint_el_pt);   
       
       // Add it to the mesh
       Volume_constraint_mesh_pt->add_element_pt(el_pt);     
      } 
    }
  }
 
 ///  Delete volume constraint elements
 void delete_volume_constraint_elements()
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Volume_constraint_mesh_pt->nelement();
   
   // Loop over the surface elements (but don't kill the volume constraint
   // element (element 0))
   unsigned first_el_to_be_killed=1;
   for(unsigned e=first_el_to_be_killed;e<n_element;e++) 
    {
     delete Volume_constraint_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Volume_constraint_mesh_pt->flush_element_and_node_storage();
   
  } // end of delete_volume_constraint_elements
 

 /// Snap the boundary nodes onto the sphere
 void snap_onto_sphere()
  {
   const unsigned n_boundary = Fluid_mesh_pt->nboundary();
   for(unsigned b=6;b<n_boundary;b++)
    {
     unsigned n_node = Fluid_mesh_pt->nboundary_node(b);
     for(unsigned n=0;n<n_node;++n)
      {
       Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(b,n);
       double x = nod_pt->x(0);
       double y = nod_pt->x(1);
       double z = nod_pt->x(2);
       
       //Now let's snap by calculating the angle
       double r = sqrt(x*x + y*y + z*z);
       double theta = acos(z/r);
       double phi = atan2(y,x);
       
       //Do the snapping
       double R_new = Global_Parameters::Radius;
       nod_pt->x(0) = R_new*sin(theta)*cos(phi);
       nod_pt->x(1) = R_new*sin(theta)*sin(phi);
       nod_pt->x(2) = R_new*cos(theta);
      }
     
     //Loop over the elements adjacent to the boundary
     unsigned n_element = this->Fluid_mesh_pt->nboundary_element(b);
     if(n_element > 0)
      {
       //Make a dummy simplex element
       TElement<3,2> dummy_four_node_element;
       //Make a dummy quadratic element
       TElement<3,3> dummy_ten_node_element;
       Vector<double> s(3);
       Vector<double> x_new(3);
       for(unsigned n=0;n<4;n++) {dummy_four_node_element.construct_node(n);}
       for(unsigned n=0;n<10;n++) {dummy_ten_node_element.construct_node(n);}
       for(unsigned e=0;e<n_element;e++)
        {
         //Cache the element pointer
         ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(
          this->Fluid_mesh_pt->boundary_element_pt(b,e));
         
         //Find the number of nodes
         const unsigned n_node = elem_pt->nnode();
         //Only do something if not simplex element
         if(n_node > 4)
          {
           //Copy the nodes into the dummy
           for(unsigned n=0;n<4;n++)
            {
             for(unsigned i=0;i<3;i++)
              {
               dummy_four_node_element.node_pt(n)->x(i) = 
                elem_pt->node_pt(n)->x(i);
              }
            }
           
           //Now sort out the mid-side nodes
           for(unsigned n=4;n<10;n++)
            {
             //If it's not on a boundary then reset to be interpolated
             //from the simplex
             if(!elem_pt->node_pt(n)->is_on_boundary())
              {
               elem_pt->local_coordinate_of_node(n,s);
               dummy_four_node_element.interpolated_x(s,x_new);
               for(unsigned i=0;i<3;i++) {elem_pt->node_pt(n)->x(i) = x_new[i];}
              }
            }
          }
         
         //If we have more than 10 nodes interpolate from the quadratic shape
         if(n_node > 10)
          {
           //Copy the nodes into the dummy
           for(unsigned n=0;n<10;n++)
            {
             for(unsigned i=0;i<3;i++)
              {
               dummy_ten_node_element.node_pt(n)->x(i) = 
                elem_pt->node_pt(n)->x(i);
              }
            }
           
           //Now sort out the mid-face and central nodes
           for(unsigned n=10;n<n_node;n++)
            {
             //If it's not on a boundary then reset to be interpolated
             //from the simplex
             if(!elem_pt->node_pt(n)->is_on_boundary())
              {
               elem_pt->local_coordinate_of_node(n,s);
               dummy_ten_node_element.interpolated_x(s,x_new);
               for(unsigned i=0;i<3;i++) {elem_pt->node_pt(n)->x(i) = x_new[i];}
              }
            }
          }
        }
      } //End of fix up of elements
    }

  }

 /// Pointer to data that will store the bubble pressure
 Data *Bubble_pressure_data_pt;

public:

 /// Pointer to the fluid mesh
 RefineableSolidTetgenMesh<ELEMENT>* Fluid_mesh_pt;
 
 /// Pointers to mesh of free surface elements
 Mesh* Free_surface_mesh_pt;
 
 /// Pointer to mesh containing elements that impose volume constraint
 Mesh* Volume_constraint_mesh_pt;

 /// Pointer to element that imposes volume constraint for bubble
 VolumeConstraintElement* Vol_constraint_el_pt;

 /// Storage for the outer boundary object
 TetgenMeshOLDFacetedSurface* Outer_boundary_pt;

 Vector<TetgenMeshOLDFacetedSurface*> Inner_boundary_pt;

};



//========================================================================
/// Constructor for RisingBubble problem
//========================================================================
template<class ELEMENT>
RisingBubbleProblem<ELEMENT>::RisingBubbleProblem()
{ 

 //Add a time stepper
 this->add_time_stepper_pt(new BDF<2>);

 // Create bubble pressure as global Data
 Bubble_pressure_data_pt = new Data(1);
 unsigned index_of_traded_pressure=0;
 this->add_global_data(Bubble_pressure_data_pt);

 Vol_constraint_el_pt = 
  new VolumeConstraintElement(&Global_Parameters::Volume,
                              Bubble_pressure_data_pt,
                              index_of_traded_pressure);
 
 //Provide a reasonable initial guess for bubble pressure (hydrostatics):
 Bubble_pressure_data_pt->set_value(index_of_traded_pressure,
                                    2.0*Global_Parameters::Ca/
                                    Global_Parameters::Radius);

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
 box_facet_boundary_id[0] = 1;
 box_facet_boundary_id[1] = 2;
 box_facet_boundary_id[2] = 3;
 box_facet_boundary_id[3] = 4;
 box_facet_boundary_id[4] = 5;
 box_facet_boundary_id[5] = 6;

 //Make the outer boundary object
 Outer_boundary_pt = 
  new TetgenMeshOLDFacetedSurface(box_point,box_facet,box_facet_boundary_id);
 

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

 //Scale the points by the radius
 for(unsigned p=0;p<12;p++)
  {
   double r2 = 0.0;
   for(unsigned i=0;i<3;i++) {r2 += icosa_point[p][i]*icosa_point[p][i];}
   
   //Now scale all of the points
   for(unsigned i=0;i<3;i++) 
    {
     icosa_point[p][i] *= Global_Parameters::Radius/sqrt(r2);
    }
  }

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

 //Set each icosahedral facet to be a separate boundary
 //This is required so that we can set up a unique surface 
 //coordinate system
 Vector<unsigned> icosa_facet_boundary_id(20);
 for(unsigned i=0;i<20;i++)
  {
   icosa_facet_boundary_id[i] = 7+i;
  }


 //Create the inner boundary object
 Inner_boundary_pt.resize(1);
 Inner_boundary_pt[0] = 
  new TetgenMeshClosedFacetedSurface(icosa_point,icosa_facet,
                                     icosa_facet_boundary_id);

 Vector<double> inner_point(3,0.0);
 dynamic_cast<TetgenMeshClosedFacetedSurface*>(Inner_boundary_pt[0])->
  set_hole(inner_point);

 /*in.numberofregions = 1;
 in.regionlist = new double[5*in.numberofregions];
 in.regionlist[0] = 0.0;
 in.regionlist[1] = 0.0;
 in.regionlist[2] = 0.0;
 in.regionlist[3] = 1;
 in.regionlist[4] = 1;*/

 Fluid_mesh_pt = 
  new RefineableSolidTetgenMesh<ELEMENT>(Outer_boundary_pt,
                                         Inner_boundary_pt,2.0,
                                         this->time_stepper_pt(),
                                         true);

 //Must split the elements in the corners
 Fluid_mesh_pt->split_elements_in_corners(this->time_stepper_pt());
 
 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Fluid_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // Set targets for spatial adaptivity
 Fluid_mesh_pt->max_permitted_error()=0.005;
 Fluid_mesh_pt->min_permitted_error()=0.001; 
 Fluid_mesh_pt->max_element_size()=1.0;
 Fluid_mesh_pt->min_element_size()=0.001; 

 // Use coarser mesh during validation
 if (CommandLineArgs::command_line_flag_has_been_set("--validation"))
  {
   Fluid_mesh_pt->max_element_size()=2.0;
   Fluid_mesh_pt->min_element_size()=0.1; 
  }

 // Set the problem pointer
 Fluid_mesh_pt->problem_pt()=this;

 // Set the boundary conditions for this problem 
 // Only on the "outer boundaries"
 for(unsigned ibound=0;ibound<6;ibound++)
 {
  unsigned final_index = 3;
  //Do no pin the outlet z-velocity
  if(ibound==2) {final_index = 2;}
   unsigned num_nod= Fluid_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
   {
    Node* nod_pt = Fluid_mesh_pt->boundary_node_pt(ibound,inod);
    //Pin the fluid velocities
    for(unsigned i=0;i<final_index;++i) {nod_pt->pin(i);}
    
    //Now pin all positions not on bubble
    for(unsigned i=0;i<3;i++) 
     {dynamic_cast<SolidNode*>(nod_pt)->pin_position(i);}
   }
 }

 // Complete the build of all elements so they are fully functional

 //Find number of elements in mesh
 unsigned n_element = Fluid_mesh_pt->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(i));

   //Set the source function pointer
   el_pt->re_pt() = &Global_Parameters::Re;

   //Set the gravity term
   el_pt->re_invfr_pt() = &Global_Parameters::ReInvFr;

   //Set gravity
   el_pt->g_pt() = &Global_Parameters::G;
   
   //Set the contitutive parameter
   el_pt->constitutive_law_pt() = Global_Parameters::Constitutive_law_pt;
  }

 //Loop over the elements in region 1
 /*unsigned n_inner = Fluid_mesh_pt->nregion_element(1);
 for(unsigned e=0;e<n_inner;++e)
  {
   // Upcast from GeneralElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->region_element_pt(1,e));

   //Set the source function pointer
   el_pt->viscosity_ratio_pt() = &Global_Parameters::Visc_Ratio;
   }*/

 //Now create the other meshes
 Free_surface_mesh_pt = new Mesh;
 create_free_surface_elements();

 Volume_constraint_mesh_pt = new Mesh;
 create_volume_constraint_elements();

 //Combine the meshes
 //Bubble_pressure_data_pt->pin(0);

 // Add volume constraint sub mesh
 this->add_sub_mesh(this->Volume_constraint_mesh_pt);

 // Add Fluid_mesh_pt sub meshes
 this->add_sub_mesh(Fluid_mesh_pt);

 // Add Free_surface sub meshes
 this->add_sub_mesh(this->Free_surface_mesh_pt);
 
 // Build global mesh
 this->build_global_mesh();


 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

}



//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RisingBubbleProblem<ELEMENT>::doc_solution(const unsigned& nplot,
                                           DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Doc local node numbering
 //-------------------------
 sprintf(filename,"%s/node_numbering%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 FiniteElement* el_pt=Fluid_mesh_pt->finite_element_pt(0);
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

 //Tell me how many elements there are
 std::cout << "Surface " << Free_surface_mesh_pt->nelement() << "\n";

 // Output boundaries
 //------------------
 sprintf(filename,"%s/surface%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Free_surface_mesh_pt->output(some_file,nplot);

 //Tell me what we have
 some_file.close();


 // Output solution
 //----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Fluid_mesh_pt->output(some_file,nplot);
 some_file.close();



} // end of doc




//========================================================================
/// Demonstrate how to solve Poisson problem
//========================================================================
int main(int argc, char* argv[])
{

 //Set the direction of gravity
 Global_Parameters::G[0] = 0.0;
 Global_Parameters::G[1] = 0.0;
 Global_Parameters::G[2] = -1.0;

 // Label for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT");
  
 // Number of output points per edge
 unsigned nplot=2;

 // Do the problem with quadratic elements
 //---------------------------------------
 {
  RisingBubbleProblem<ProjectableTaylorHoodElement<
  PseudoSolidNodeUpdateElement<TTaylorHoodElement<3>,
   TPVDElement<3,3> > > > problem;
  problem.snap_onto_sphere();
  problem.Fluid_mesh_pt->set_lagrangian_nodal_coordinates();
  problem.doc_solution(5,doc_info);
  doc_info.number()++;

  Global_Parameters::Volume -= 0.1255;

  problem.steady_newton_solve();
  problem.doc_solution(5,doc_info);
  doc_info.number()++;

  //Assume it's always been this way
  double dt = 0.1;
  problem.Fluid_mesh_pt->set_lagrangian_nodal_coordinates();
  problem.assign_initial_values_impulsive(dt);
  //Turn on gravity
  Global_Parameters::ReInvFr = 1.0;

  for(unsigned t=0;t<10;t++)
   {
    problem.unsteady_newton_solve(dt);
    problem.doc_solution(5,doc_info);
    doc_info.number()++;
    //Reset the lagrangian
    problem.Fluid_mesh_pt->set_lagrangian_nodal_coordinates();
   }

  //Now can I adapt
  problem.adapt();
  
  problem.unsteady_newton_solve(dt);
  problem.doc_solution(5,doc_info);
  doc_info.number()++;
  //Reset the lagrangian
  problem.Fluid_mesh_pt->set_lagrangian_nodal_coordinates();

  exit(1);

  for(unsigned n=0;n<1;++n)
   {
    // Solve the problem
    problem.steady_newton_solve(1);
    
    //Output solution with 5 points per edge
    nplot=5;
    problem.doc_solution(nplot,doc_info);
    
    //Increment counter for solutions 
    doc_info.number()++;
    
    Global_Parameters::Re += 0.1;
   }

 }


}



