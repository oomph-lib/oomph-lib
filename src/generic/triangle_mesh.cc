//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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

#include <algorithm>
#include "map_matrix.h"
#include "triangle_mesh.h"


namespace oomph
{



/// Namespace that allows the specification of a tolerance 
/// between vertices at the ends of polylines that are supposed
/// to be at the same position.
namespace ToleranceForVertexMismatchInPolygons
{
 
 ///  Acceptable discrepancy for mismatch in vertex coordinates.
 /// In paranoid mode, the code will die if the beginning/end of
 /// two adjacent polylines differ by more than that. If the
 /// discrepancy is smaller (but nonzero) one of the vertex coordinates
 /// get adjusted to match perfectly; without paranoia the vertex
 /// coordinates are taken as they come...
 double Tolerable_error=1.0e-14;

}


//=========================================================================
/// Constructor: Specify vector of pointers to TriangleMeshPolyLines
/// that define the boundary of the segments of the polygon.
/// Each TriangleMeshPolyLine has its own boundary ID and can contain
/// multiple (straight-line) segments. If there is just a single
/// polyline, the first and last vertices should not coincide -- we 
/// will close the polygon for you! However, if there multiple
/// polylines their joint vertices must be specified in both
/// polylines (since the polylines may be used in isolation). 
//=========================================================================
TriangleMeshPolygon::TriangleMeshPolygon(const Vector<TriangleMeshPolyLine*>& 
                                         boundary_polyline_pt) :
 Boundary_polyline_pt(boundary_polyline_pt)
  {
   
#ifdef PARANOID
   
   // Check that the polylines are contiguous
   bool contiguous=true;
   unsigned i_offensive=0;
   unsigned nbound=Boundary_polyline_pt.size();

   // Multiple polylines
   if (nbound>1)
    {
     // Does the last node of the polyline connect to the first one
     // of the next one (only up the last but one!)
     for(unsigned i=0;i<nbound-1;i++)
      {
       // Get vector last vertex in current polyline
       unsigned last_vertex = (Boundary_polyline_pt[i]->nvertex())-1;
       Vector<double> v1=Boundary_polyline_pt[i]->
        vertex_coordinate(last_vertex);
       
       // Get vector to first vertex in next polyline
       Vector<double> v2=Boundary_polyline_pt[i+1]->vertex_coordinate(0);
       
       // Work out error
       double error=sqrt(pow(v1[0]-v2[0],2)+pow(v1[1]-v2[1],2));
       
       // Is error accetable?
       if (error>ToleranceForVertexMismatchInPolygons::Tolerable_error)       
        {
         contiguous=false;
         i_offensive=i;
         break;
        }
       // Align
       else
        {
         Boundary_polyline_pt[i+1]->vertex_coordinate(0)=
          Boundary_polyline_pt[i]->vertex_coordinate(last_vertex);
        }
      }
     
     // Does the last one connect to the first one?
     
     // Get vector last vertex last polyline
     unsigned last_vertex = (Boundary_polyline_pt[nbound-1]->nvertex())-1;
     Vector<double> v1=Boundary_polyline_pt[nbound-1]->
      vertex_coordinate(last_vertex);
     
     // Get vector first vertex first polyline
     Vector<double> v2=Boundary_polyline_pt[0]->vertex_coordinate(0);
     double error=sqrt(pow(v1[0]-v2[0],2)+pow(v1[1]-v2[1],2));
     if (error>ToleranceForVertexMismatchInPolygons::Tolerable_error)
      {
       contiguous=false;
       i_offensive=nbound-1;
      } 
     else
      {
       Boundary_polyline_pt[0]->vertex_coordinate(0)=
        Boundary_polyline_pt[nbound-1]->vertex_coordinate(last_vertex);
      }
     
     if (!contiguous)
      {
       std::ostringstream error_stream;
       error_stream
        << "When a Polygon is defined by multiple polylines, the polylines\n"
        << "should define a closed geometry, i.e. the first/last vertex of\n"
        << "adjacent polylines should match.\n\n"
        << "Your polyline number "<< i_offensive 
        <<" has no contiguous neighbour, when judged \nwith the tolerance of "
        << ToleranceForVertexMismatchInPolygons::Tolerable_error
        << " which is specified in the namespace \nvariable "
        << "ToleranceForVertexMismatchInPolygons::Tolerable_error.\n\n"
        << "Feel free to adjust this or to recompile the code without\n"
        << "paranoia if you think this is OK...\n"
        << std::endl;
       throw OomphLibError(error_stream.str(),
                           "TriangleMeshPolygon::TriangleMeshPolygon()",
                           OOMPH_EXCEPTION_LOCATION);
      }

    }
   // Single polyline
   else
    {
     // Should not be closed
     
     // Get vector last vertex in polyline
     unsigned last_vertex = (Boundary_polyline_pt[0]->nvertex())-1;
     Vector<double> v1=Boundary_polyline_pt[0]->
      vertex_coordinate(last_vertex);
     
     // Get vector first vertex first polyline
     Vector<double> v2=Boundary_polyline_pt[0]->vertex_coordinate(0);

     // Gatp
     double gap=sqrt(pow(v1[0]-v2[0],2)+pow(v1[1]-v2[1],2));
     if (gap<ToleranceForVertexMismatchInPolygons::Tolerable_error)
      {
       std::ostringstream error_stream;
       error_stream
        << "When a Polygon is defined by a single polyline, the polyline\n"
        << "should define an open geometry, i.e. the first/last vertex of\n"
        << "the polyline should differ.\n\n"
        << "The first and last vertices of your polyline appear to coincide "
        <<" when judged \nwith the tolerance of "
        << ToleranceForVertexMismatchInPolygons::Tolerable_error
        << " which is specified in the namespace \nvariable "
        << "ToleranceForVertexMismatchInPolygons::Tolerable_error.\n\n"
        << "Feel free to adjust this (e.g. set it to a negative value) or\n"
        << "to recompile the code withoutparanoia if you think this is OK...\n"
        << std::endl;
       throw OomphLibError(error_stream.str(),
                           "TriangleMeshPolygon::TriangleMeshPolygon()",
                           OOMPH_EXCEPTION_LOCATION);
      }
    }
#endif
     
    }
   


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


//=======================================================================
/// \short Constructor: Specify coordinates of a point inside the hole
/// and a vector of pointers to TriangleMeshPolyLines
/// that define the boundary segments of the polygon.
/// Each TriangleMeshPolyLine has its own boundary ID and can contain
/// multiple (straight-line) segments. The optional final argument
/// is a pointer to a Data object whose three values represent 
/// the two displacements of and the rotation angle about the polygon's 
/// centre of mass.
//=======================================================================
TriangleMeshHolePolygon::TriangleMeshHolePolygon(
 const Vector<double>& hole_center,
 const Vector<TriangleMeshPolyLine*>& 
 boundary_polyline_pt,
 Data* centre_displacement_data_pt) :
 TriangleMeshPolygon(boundary_polyline_pt), 
 Hole_coordinate(hole_center),
 Centre_displacement_data_pt(centre_displacement_data_pt),
 External_force_fct_pt(0), External_torque_fct_pt(0), Drag_mesh_pt(0)
{  
 
 // Provide Data for centre-of-mass displacement internally
 if (Centre_displacement_data_pt==0)
  {
   Newmark<2>* timestepper_pt=new Newmark<2>;
   Centre_displacement_data_pt=new Data(timestepper_pt,3);
     
   // I've created it so I have to tidy up too!
   Must_clean_up=true;
  }
 // Data created externally, so somebody else will clean up
 else
  {
   Must_clean_up=false;
  }

 // Centre displacement is internal Data for this element
 add_internal_data(Centre_displacement_data_pt);

 // Build the polyline geom objects
 create_boundary_geom_objects(boundary_polyline_pt); 

 // Original rotation angle is zero
 Phi_c_orig=0.0;

 // Compute coordinates of centre of gravity etc
 Vector<double> r_left(2);
 Vector<double> r_right(2);
 Mass=0.0;
 X_c_orig=0.0;
 Y_c_orig=0.0;
 double inertia_x=0.0;
 double inertia_y=0.0;

 // Loop over polylines
 unsigned nboundary=boundary_polyline_pt.size();
 for (unsigned i=0;i<nboundary;i++)
  {
   // Loop over the segments to get the vertex coordinates
   unsigned nseg=boundary_polyline_pt[i]->nsegment();
   for(unsigned j=0;j<nseg;j++)
    {
     // Get the vertex coordinates
     r_left =boundary_polyline_pt[i]->vertex_coordinate(j);
     r_right=boundary_polyline_pt[i]->vertex_coordinate(j+1);
   
     // Mass (area)
     Mass+=0.5*(r_left[0]*r_right[1]-r_right[0]*r_left[1]);

     // Centroid
     X_c_orig+=(r_left[0]+r_right[0])*
      (r_left[0]*r_right[1]-r_right[0]*r_left[1]);
     Y_c_orig+=(r_left[1]+r_right[1])*
      (r_left[0]*r_right[1]-r_right[0]*r_left[1]);
    }
   if (nboundary==1)
    {
     // Get the vertex coordinates
     r_left =boundary_polyline_pt[0]->vertex_coordinate(nseg);
     r_right=boundary_polyline_pt[0]->vertex_coordinate(0);
   
     // Mass (area)
     Mass+=0.5*(r_left[0]*r_right[1]-r_right[0]*r_left[1]);

     // Centroid
     X_c_orig+=(r_left[0]+r_right[0])*
      (r_left[0]*r_right[1]-r_right[0]*r_left[1]);
     Y_c_orig+=(r_left[1]+r_right[1])*
      (r_left[0]*r_right[1]-r_right[0]*r_left[1]);
    }
  }
   
 // Normalise
 X_c_orig/=(6.0*Mass);
 Y_c_orig/=(6.0*Mass);
   
 // Another loop over polylines for moment of inertia
 for (unsigned i=0;i<nboundary;i++)
  {
   // Loop over the segments to get the vertex coordinates
   unsigned nseg=boundary_polyline_pt[i]->nsegment();
   for(unsigned j=0;j<nseg;j++)
    {
     // Get the vertex coordinates
     r_left =boundary_polyline_pt[i]->vertex_coordinate(j);
     r_right=boundary_polyline_pt[i]->vertex_coordinate(j+1);
       
     // Get moment about centroid
     r_left[0]-=X_c_orig;
     r_left[1]-=Y_c_orig;
     r_right[0]-=X_c_orig;
     r_right[1]-=Y_c_orig;
       
     // Moment of inertia
     inertia_x+=1.0/12.0*(r_left[1]*r_left[1]+
                          r_left[1]*r_right[1]+
                          r_right[1]*r_right[1])*
      (r_left[0]*r_right[1]-r_right[0]*r_left[1]);
       
     inertia_y+=1.0/12.0*(r_left[0]*r_left[0]+
                          r_left[0]*r_right[0]+
                          r_right[0]*r_right[0])*
      (r_left[0]*r_right[1]-r_right[0]*r_left[1]);       
    }
     
   if (nboundary==1)
    {
     // Get the vertex coordinates
     r_left =boundary_polyline_pt[0]->vertex_coordinate(nseg);
     r_right=boundary_polyline_pt[0]->vertex_coordinate(0);
       
     // Get moment about centroid
     r_left[0]-=X_c_orig;
     r_left[1]-=Y_c_orig;
     r_right[0]-=X_c_orig;
     r_right[1]-=Y_c_orig;
       
     // Moment of inertia
     inertia_x+=1.0/12.0*(r_left[1]*r_left[1]+
                          r_left[1]*r_right[1]+
                          r_right[1]*r_right[1])*
      (r_left[0]*r_right[1]-r_right[0]*r_left[1]);
       
     inertia_y+=1.0/12.0*(r_left[0]*r_left[0]+
                          r_left[0]*r_right[0]+
                          r_right[0]*r_right[0])*
      (r_left[0]*r_right[1]-r_right[0]*r_left[1]);    
    }
  }
   
 // Polar moment of inertia is sum of two orthogonal planar moments
 Moment_of_inertia=inertia_x+inertia_y;
   
//    // Tested for circular and elliptical cross section
//    cout << "Mass             : " << Mass << std::endl;
//    cout << "Moment of inertia: " << Moment_of_inertia << std::endl;
//    cout << "X_c              : " << X_c_orig << std::endl;
//    cout << "Y_c              : " << Y_c_orig << std::endl;
//    pause("done");
      
}
 


//===============================================================
/// \short Update the reference configuration by re-setting the original
/// position of the vertices to their current ones, re-set the 
/// original position of the centre of mass, and the displacements 
/// and rotations relative to it
//===============================================================
void TriangleMeshHolePolygon::reset_reference_configuration()
{
 Vector<double> x_orig(2);
 Vector<double> r(2);
 
 // Loop over the polylines and update their vertex positions
 unsigned npoly=Boundary_polyline_pt.size();
 for (unsigned i=0;i<npoly;i++)
  {
   TriangleMeshPolyLine* poly_line_pt=Boundary_polyline_pt[i];
   unsigned nvertex=poly_line_pt->nvertex();
   for (unsigned j=0;j<nvertex;j++)
    {
     x_orig=poly_line_pt->vertex_coordinate(j);
     get_updated_position(x_orig,r);
     poly_line_pt->vertex_coordinate(j)=r;
    }
  }

 // Loop over the boundary geom objects to update their lower left
 // vertices for unique boundary coordinate
 unsigned nbound=Boundary_geom_obj_pt.size();
 for(unsigned j=0;j<nbound;j++)
  {
   // Loop over the geom object
   unsigned ngeom_obj=Boundary_geom_obj_pt[j].size();
   for(unsigned igeom_obj=0;igeom_obj<ngeom_obj;igeom_obj++)
    {
     // Update lower left vertex
     Boundary_geom_obj_pt[j][igeom_obj]->update_lower_left_vertex();
    }     
  }

 // Update coordinates of hole
 Vector<double> orig_hole_coord(Hole_coordinate);
 get_updated_position(orig_hole_coord,Hole_coordinate);
     

 // Update centre of gravity
 double x_displ=Centre_displacement_data_pt->value(0);
 double y_displ=Centre_displacement_data_pt->value(1);
 double phi_displ=Centre_displacement_data_pt->value(2);
 X_c_orig+=x_displ;
 Y_c_orig+=y_displ;
 Phi_c_orig+=phi_displ;

 // Reset displacement and rotation ("non-previous-value" 
 // history values stay)
 TimeStepper* timestepper_pt=Centre_displacement_data_pt->time_stepper_pt();
 unsigned nprev=timestepper_pt->nprev_values();
 for (unsigned t=0;t<nprev;t++)
  {
   Centre_displacement_data_pt->
    set_value(t,0,Centre_displacement_data_pt->value(t,0)-x_displ);

   Centre_displacement_data_pt->
    set_value(t,1,Centre_displacement_data_pt->value(t,1)-y_displ);

   Centre_displacement_data_pt->
    set_value(t,2,Centre_displacement_data_pt->value(t,2)-phi_displ);
  }
}



/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


//===============================================================
/// Update lower left vertex for the current rotation of
/// the associated  TriangleMeshHolePolygon 
//===============================================================
void PolyLineSegmentGeomObject::update_lower_left_vertex()
{
 Vector<double> r_left(2);
 Vector<double> r_right(2);   
 
 // Get current position of left vertex when subjected to
 // the Centre_displacement_data_pt
 Polygon_pt->get_updated_position(R_left,r_left);
 Polygon_pt->get_updated_position(R_right,r_right);
 
 // Swap?
 R_left=r_left;
 R_right=r_right;
 if((r_left[1]>r_right[1])||
    ( (r_left[1]==r_right[1]) && (r_left[0]>r_right[0]) ) )
  {
   Vector<double> tmp_coord(2);
   tmp_coord=R_right;
   R_right=R_left;
   R_left=tmp_coord;
  }
}


//==========================================================
/// Position Vector at Lagrangian coordinate zeta
//==========================================================
void PolyLineSegmentGeomObject::position(const Vector<double>& zeta, 
                                         Vector<double>& r) const
{   
 // The position of the point in the original configuration
 Vector<double> x_orig(2);
 double zeta_scalar=zeta[0];
 interpolated_x_orig(zeta_scalar,x_orig);
 
 // Get new position
 Polygon_pt->get_updated_position(x_orig,r);
}



//===============================================================
/// Return pointer to the j-th Data item that the object's 
/// shape depends on 
//===============================================================
Data* PolyLineSegmentGeomObject::geom_data_pt(const unsigned& j) 
{
#ifdef RANGE_CHECKING
 if (j!=0) 
  { 
   std::ostringstream error_stream;
   error_stream
    << "PolyLineSegmentGeomObject contains just one geom_data_pt\n"
    << "geom_data_pt index has to be 0 and not "<<j
    << std::endl;
   throw OomphLibError(error_stream.str(),
                       "PolyLineSegmentGeomObject::geom_data_pt()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 return Polygon_pt->centre_displacement_data_pt();
}


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////



 //=================================================
 ///  Helper namespace for BCInfo object used
 /// in the identification of boundary elements.
 //=================================================
namespace TriangleBoundaryHelper
{

 /// Structure for Boundary Informations
 struct BCInfo
 {

  /// Face ID
  unsigned Face_id;
 
  /// Boundary ID
  unsigned Boundary;
 
  /// Pointer to bulk finite element
  FiniteElement* FE_pt;
  
 };
 
}

//================================================================
/// Setup lookup schemes which establish which elements are located
/// next to which boundaries (Doc to outfile if it's open).
//================================================================
void TriangleMeshBase::setup_boundary_element_info(std::ostream &outfile)
{

 //Should we document the output here
 bool doc = false;

 if(outfile) doc = true;

 // Number of boundaries
 unsigned nbound=nboundary();
 
 // Wipe/allocate storage for arrays
 Boundary_element_pt.clear();
 Face_index_at_boundary.clear();
 Boundary_element_pt.resize(nbound);
 Face_index_at_boundary.resize(nbound);
 
 // Temporary vector of vectors of pointers to elements on the boundaries: 
 // This is a vector to ensure that order is strictly preserved
 Vector<Vector<FiniteElement*> > vector_of_boundary_element_pt;
 vector_of_boundary_element_pt.resize(nbound);
 
 // Matrix map for working out the fixed face for elements on boundary
 MapMatrixMixed<unsigned,FiniteElement*, int > 
  face_identifier;
 
 // Loop over elements
 //-------------------
 unsigned nel=nelement();

      
 // Get pointer to vector of boundaries that the
 // node lives on
 Vector<std::set<unsigned>*> boundaries_pt(3,0);


 // Data needed to deal with edges through the
 // interior of the domain
 std::map<Edge,unsigned> edge_count;
 std::map<Edge,TriangleBoundaryHelper::BCInfo> 
  edge_bcinfo;
 std::map<Edge,TriangleBoundaryHelper::BCInfo>
  face_info;
 MapMatrixMixed<unsigned,FiniteElement*, int > face_count;
 Vector<unsigned> bonus(nbound);
     
 for (unsigned e=0;e<nel;e++)
  {
   // Get pointer to element
   FiniteElement* fe_pt=finite_element_pt(e);
   
   if (doc) outfile << "Element: " << e << " " << fe_pt << std::endl;
   
   // Only include 2D elements! Some meshes contain interface elements too.
   if (fe_pt->dim()==2)
    {
     // Loop over the element's nodes and find out which boundaries they're on
     // ----------------------------------------------------------------------

     //We need only loop over the corner nodes

     for(unsigned i=0;i<3;i++)
      {
       fe_pt->node_pt(i)->get_boundaries_pt(boundaries_pt[i]);
      }

     //Find the common boundaries of each edge
     Vector<std::set<unsigned> > edge_boundary(3);


     //Edge 0 connects points 1 and 2
     //-----------------------------

     if(boundaries_pt[1] && boundaries_pt[2])
      {
       // Create the corresponding edge
       Edge edge0(fe_pt->node_pt(1),fe_pt->node_pt(2));

       // Update infos about this edge
       TriangleBoundaryHelper::BCInfo info;      
       info.Face_id=0;
       info.FE_pt = fe_pt;
	
       std::set_intersection(boundaries_pt[1]->begin(),boundaries_pt[1]->end(),
                             boundaries_pt[2]->begin(),boundaries_pt[2]->end(),
                             std::insert_iterator<std::set<unsigned> >(
                              edge_boundary[0],edge_boundary[0].begin()));
       std::set<unsigned>::iterator it0=edge_boundary[0].begin(); 

       // Edge does exist:
       if ( edge_boundary[0].size() > 0 )
        {
         info.Boundary=*it0;
         // How many times this edge has been visited
         edge_count[edge0]++;

         // Update edge_bcinfo
         edge_bcinfo.insert(std::make_pair(edge0,info)); 
        }
      }

     //Edge 1 connects points 0 and 2
     //-----------------------------

     if(boundaries_pt[0] && boundaries_pt[2])
      {
       std::set_intersection(boundaries_pt[0]->begin(),boundaries_pt[0]->end(),
                             boundaries_pt[2]->begin(),boundaries_pt[2]->end(),
                             std::insert_iterator<std::set<unsigned> >(
                              edge_boundary[1],edge_boundary[1].begin()));

       // Create the corresponding edge
       Edge edge1(fe_pt->node_pt(0),fe_pt->node_pt(2));

       // Update infos about this edge
       TriangleBoundaryHelper::BCInfo info;
       info.Face_id=1;
       info.FE_pt = fe_pt;
       std::set<unsigned>::iterator it1=edge_boundary[1].begin();

       // Edge does exist:
       if ( edge_boundary[1].size() > 0)
        {
         info.Boundary=*it1;

         // How many times this edge has been visited
         edge_count[edge1]++;  

         // Update edge_bcinfo              
         edge_bcinfo.insert(std::make_pair(edge1,info));
        }
      }
 
     //Edge 2 connects points 0 and 1
     //-----------------------------

     if(boundaries_pt[0] && boundaries_pt[1])
      {
       std::set_intersection(boundaries_pt[0]->begin(),boundaries_pt[0]->end(),
                             boundaries_pt[1]->begin(),boundaries_pt[1]->end(),
                             std::insert_iterator<std::set<unsigned> >(
                              edge_boundary[2],edge_boundary[2].begin()));

       // Create the corresponding edge
       Edge edge2(fe_pt->node_pt(0),fe_pt->node_pt(1));

       // Update infos about this edge
       TriangleBoundaryHelper::BCInfo info;
       info.Face_id=2;
       info.FE_pt = fe_pt;
       std::set<unsigned>::iterator it2=edge_boundary[2].begin();

       // Edge does exist:
       if ( edge_boundary[2].size() > 0)
        {
         info.Boundary=*it2;

         // How many times this edge has been visited
         edge_count[edge2]++;  

         // Update edge_bcinfo
         edge_bcinfo.insert(std::make_pair(edge2,info));
        }
      }


#ifdef PARANOID
     
     // Check if edge is associated with multiple boundaries

     //We now know whether any edges lay on the boundaries
     for(unsigned i=0;i<3;i++)
      {
       //How many boundaries are there
       unsigned count = 0;

       //Loop over all the members of the set and add to the count
       //and set the boundary
       for(std::set<unsigned>::iterator it=edge_boundary[i].begin();
           it!=edge_boundary[i].end();++it)
        {
         ++count;
        }

       //If we're on more than one boundary, this is weird, so die
       if(count > 1)
        {
         std::ostringstream error_stream;
         error_stream << "Edge " << i << " is located on " << 
          count << " boundaries.\n";
         error_stream << "This is rather strange, so I'm going to die\n";
         throw OomphLibError(
          error_stream.str(),
          "TriangleMeshBase::setup_boundary_element_info()",
          OOMPH_EXCEPTION_LOCATION);
        }
      }
     
#endif

     // Now we set the pointers to the boundary sets to zero
     for(unsigned i=0;i<3;i++) {boundaries_pt[i] = 0;}
     
    }
  } //end of loop over all elements




 // Loop over all edges that are located on a boundary
 typedef std::map<Edge,TriangleBoundaryHelper::BCInfo>::iterator ITE;
 for (ITE it=edge_bcinfo.begin();
      it!=edge_bcinfo.end();
      it++)
  {
   Edge current_edge = it->first;
   unsigned  bound=it->second.Boundary; 
   
   //If the edge has been visited only once
   if(edge_count[current_edge]==1)
    {
     // Count the edges that are on the same element and on the same boundary
     face_count(static_cast<unsigned>(bound),it->second.FE_pt)=  
      face_count(static_cast<unsigned>(bound),it->second.FE_pt) + 1;
     
     //If such edges exist, let store the corresponding element
     if( face_count(bound,it->second.FE_pt) > 1)
      {
       // Update edge's infos
       TriangleBoundaryHelper::BCInfo info;
       info.Face_id=it->second.Face_id;
       info.FE_pt = it->second.FE_pt;
       info.Boundary=it->second.Boundary;
       
       // Add it to FIinfo, that stores infos of problematic elements
       face_info.insert(std::make_pair(current_edge,info));
       
       //How many edges on which boundary have to be added
       bonus[bound]++;
      }
     else
      {
       //Add element and face to the appropriate vectors

       // Does the pointer already exits in the vector
       Vector<FiniteElement*>::iterator b_el_it =
        std::find(vector_of_boundary_element_pt[
                   static_cast<unsigned>(bound)].begin(),
                  vector_of_boundary_element_pt[
                   static_cast<unsigned>(bound)].end(),
                  it->second.FE_pt);
       
       //Only insert if we have not found it (i.e. got to the end)
       if(b_el_it == vector_of_boundary_element_pt[
           static_cast<unsigned>(bound)].end())
        {
         vector_of_boundary_element_pt[static_cast<unsigned>(bound)].
          push_back(it->second.FE_pt);
        }
       
       //set_of_boundary_element_pt[static_cast<unsigned>(bound)].insert(
       // it->second.FE_pt);
       face_identifier(static_cast<unsigned>(bound),it->second.FE_pt) = 
        it->second.Face_id;
      }
     
    }
   
  }  //End of "adding-boundaries"-loop
 
 
  
 // Now copy everything across into permanent arrays
 //-------------------------------------------------

 // Loop over boundaries
 for (unsigned i=0;i<nbound;i++)
  {
   // Number of elements on this boundary that have to be added 
   // in addition to other elements
   unsigned bonus1=bonus[i];
   
   // Number of elements on this boundary
   unsigned nel=vector_of_boundary_element_pt[i].size() + bonus1;

   // Allocate storage for the coordinate identifiers
   Face_index_at_boundary[i].resize(nel);

   unsigned e_count=0;
   typedef Vector<FiniteElement*>::iterator IT;
   for (IT it=vector_of_boundary_element_pt[i].begin();
        it!=vector_of_boundary_element_pt[i].end();
        it++)
    {    
     // Recover pointer to element
     FiniteElement* fe_pt=*it;

     // Add to permanent storage
     Boundary_element_pt[i].push_back(fe_pt);

     Face_index_at_boundary[i][e_count] = face_identifier(i,fe_pt);

     // Increment counter
     e_count++;

    }
   // We add the elements that have two or more edges on this boundary
   for (ITE itt= face_info.begin(); itt!= face_info.end(); itt++)
    {
     if (itt->second.Boundary==i)
      {
       // Add to permanent storage
       Boundary_element_pt[i].push_back(itt->second.FE_pt);

       Face_index_at_boundary[i][e_count] = itt->second.Face_id;

       e_count++;
      }

    }

  } //End of loop over boundaries
 


 // Doc?
 //-----
 if (doc)
  {
   // Loop over boundaries
   for (unsigned i=0;i<nbound;i++)
    {
     unsigned nel=Boundary_element_pt[i].size();
     outfile << "Boundary: " << i
             << " is adjacent to " << nel
             << " elements" << std::endl;
     
     // Loop over elements on given boundary
     for (unsigned e=0;e<nel;e++)
      {
       FiniteElement* fe_pt=Boundary_element_pt[i][e];
       outfile << "Boundary element:" <<  fe_pt
               << " Face index of boundary is " 
               <<  Face_index_at_boundary[i][e] << std::endl;
      }
    }
  }
 

 // Lookup scheme has now been setup yet
 Lookup_for_elements_next_boundary_is_setup=true;

}

}
