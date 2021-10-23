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
// Driver code to compare I/O speed for C++ and C-style output


// Various headers
#include "generic.h"
#include "navier_stokes.h"
#include "poisson.h"
#include "unsteady_heat.h"
#include "solid.h"
#include "beam.h"

// 1D Lagrangian mesh 
#include "meshes/one_d_lagrangian_mesh.h"

// 2D quarter circle mesh
#include "meshes/quarter_circle_sector_mesh.h"

// Another 2D mesh
#include "meshes/simple_rectangular_quadmesh.h"

// 3D quarter tube mesh
#include "meshes/quarter_tube_mesh.h"

using namespace std;

using namespace oomph;


////////////////////////////////////////////////////////////////////
// Upgrade quarter circle mesh to become usable with solid mechanics
// elements
////////////////////////////////////////////////////////////////////


//================================================================
/// Elastic quarter circle sector mesh with functionality to
/// attach traction elements to the curved surface. We "upgrade"
/// the RefineableQuarterCircleSectorMesh to become an
/// SolidMesh and equate the Eulerian and Lagrangian coordinates,
/// thus making the domain represented by the mesh the stress-free 
/// configuration. 
/// \n\n
/// The member function \c make_traction_element_mesh() creates
/// a separate mesh of SolidTractionElements that are attached to the
/// mesh's curved boundary (boundary 1). 
//================================================================
template <class ELEMENT>
class ElasticRefineableQuarterCircleSectorMesh :
 public virtual RefineableQuarterCircleSectorMesh<ELEMENT>,
 public virtual SolidMesh
{


public:

 /// Constructor: Build mesh and copy Eulerian coords to Lagrangian
 /// ones so that the initial configuration is the stress-free one.
 ElasticRefineableQuarterCircleSectorMesh<ELEMENT>(GeomObject* wall_pt,
                                         const double& xi_lo,
                                         const double& fract_mid,
                                         const double& xi_hi,
                                         TimeStepper* time_stepper_pt=
                                         &Mesh::Default_TimeStepper) :
  RefineableQuarterCircleSectorMesh<ELEMENT>(wall_pt,xi_lo,fract_mid,xi_hi,
                                             time_stepper_pt)
  {
   /// Make the current configuration the undeformed one by
   /// setting the nodal Lagrangian coordinates to their current
   /// Eulerian ones
   set_lagrangian_nodal_coordinates();
  }


 /// Function to create mesh made of traction elements
 void make_traction_element_mesh(SolidMesh*& traction_mesh_pt)
  {

   // Make new mesh
   traction_mesh_pt = new SolidMesh;

   // Loop over all elements on boundary 1:
   unsigned b=1;
   unsigned n_element = this->nboundary_element(b);
   for (unsigned e=0;e<n_element;e++)
    {
     // The element itself:
     FiniteElement* fe_pt = this->boundary_element_pt(b,e);
     
     // Find the index of the face of element e along boundary b
     int face_index = this->face_index_at_boundary(b,e);
     
     // Create new element
     traction_mesh_pt->add_element_pt(new SolidTractionElement<ELEMENT>
                                      (fe_pt,face_index));
    }
  }

};


////////////////////////////////////////////////////////////////////
// Upgrade simple rectangular quad mesh to become refineable
////////////////////////////////////////////////////////////////////



//==============================start_of_mesh======================
/// Refineable equivalent of the SimpleRectangularQuadMesh.
/// Refinement is performed by the QuadTree-based procedures
/// implemented in the RefineableQuadMesh base class.
//=================================================================
template<class ELEMENT>
class SimpleRefineableRectangularQuadMesh : 
 public virtual SimpleRectangularQuadMesh<ELEMENT>,  
 public RefineableQuadMesh<ELEMENT>
{ 

public: 

 ///  Pass number of elements in the horizontal 
 /// and vertical directions, and the corresponding dimensions.
 /// Timestepper defaults to Static.
 SimpleRefineableRectangularQuadMesh(const unsigned &Nx,
                                     const unsigned &Ny, 
                                     const double &Lx, const double &Ly,
                                     TimeStepper* time_stepper_pt=
                                     &Mesh::Default_TimeStepper) :
  SimpleRectangularQuadMesh<ELEMENT>(Nx,Ny,Lx,Ly,time_stepper_pt)
  {
   // Nodal positions etc. were created in constructor for
   // SimpleRectangularQuadMesh<...> --> We only need to set up 
   // adaptivity information: Associate finite elements with their 
   // QuadTrees and plant them in a QuadTreeForest:
   this->setup_quadtree_forest();

  } // end of constructor
 

 /// Destructor: Empty
 virtual ~SimpleRefineableRectangularQuadMesh() {}

}; // end of mesh



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



//=====================================================================
/// Global output function to compare C++ and C-style output routines
//=====================================================================
void output_both_versions(const string& file_root, Mesh* mesh_pt,
                          const unsigned& npts)
{


 cout << std::endl << "==============================================" 
      << std::endl << std::endl;

 cout << file_root << std::endl << std::endl;
 DocInfo doc_info; 
 char filename[100];

 ofstream some_file;
 
 double total_time_old, total_time_new;

 // C style output
 //-----------------
 {
  cout << "C style" << std::endl;
  doc_info.set_directory("RESLT_C_style");
  
  // Output solution 
  sprintf(filename,"%s/%s.dat",doc_info.directory().c_str(),
          file_root.c_str());
  FILE* file_pt = fopen(filename,"w");
  clock_t t_start = clock();  
  mesh_pt->output(file_pt,npts);
  clock_t t_end = clock();
  total_time_new=double(t_end-t_start)/CLOCKS_PER_SEC;
  cout << "CPU time [sec]: " 
       << total_time_new << std::endl << std::endl;
  fclose(file_pt);


  
 }
 // C++ style output
 //-----------------
 {
  cout << "C++ style" << std::endl;
  doc_info.set_directory("RESLT_Cpp_style");
  
  // Output solution 
  sprintf(filename,"%s/%s.dat",doc_info.directory().c_str(),
          file_root.c_str());
  some_file.open(filename);
  clock_t t_start = clock();  
  mesh_pt->output(some_file,npts);
  clock_t t_end = clock();
  total_time_old=double(t_end-t_start)/CLOCKS_PER_SEC;
  cout << "CPU time [sec]: " 
       << total_time_old << std::endl << std::endl;
  some_file.close();
 }

 cout << "                    RATIO: " << total_time_old/total_time_new
      << std::endl;

}
 








////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////








//=====================================================================
/// Ugly driver code -- brutally loop over lots of mesh/element
/// combinations and produce C/C++-style output. 
//=====================================================================
int main(int argc, char *argv[])
{

 unsigned nplot=5;
 unsigned nelem=1000;


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

 {  
  // Set the undeformed beam to be a straight line at y=0
  GeomObject* undef_beam_pt=new StraightLine(0.0); 
  
  // Build and assign mesh
  double L=1.0;
  Mesh* mesh_pt = 
   new OneDLagrangianMesh<HermiteBeamElement>(nelem,L,undef_beam_pt);
  
  //Assign undeformed geometry
  for(unsigned i=0;i<nelem;i++)
   {      
    //Cast to proper element type
    HermiteBeamElement* elem_pt = 
     dynamic_cast<HermiteBeamElement*>(mesh_pt->element_pt(i));
    
    //Assign the undeformed beam shape
    elem_pt->undeformed_beam_pt() = undef_beam_pt;
   }
  
  // Filename for output
  string file_root="beam";
    
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }
   
   
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
   
 {
  
  // Build the geometric object that describes the outer wall
  GeomObject* curved_boundary_pt = new Ellipse(1.0,1.0);
  
  // The curved boundary of the mesh is defined by the geometric object
  // What follows are the start and end coordinates on the geometric object:
  double xi_lo=0.0;
  double xi_hi=2.0*atan(1.0);
  
  // Fraction along geometric object at which the radial dividing line
  // is placed
  double fract_mid=0.5;
 
 
  //--------------------------------------------------
  {
   // Build and assign mesh
   Mesh* mesh_pt = 
    new ElasticRefineableQuarterCircleSectorMesh<RefineableQPVDElement<2,2> >
    (curved_boundary_pt,xi_lo,fract_mid,xi_hi);
  
   // Filename for output
   string file_root="adapt_solid_pvd_2";
  
   // Compare output
   output_both_versions(file_root,mesh_pt,nplot);
  
   delete mesh_pt;
  }
 
  //--------------------------------------------------
  {
   // Build and assign mesh
   Mesh* mesh_pt = 
    new ElasticRefineableQuarterCircleSectorMesh<RefineableQPVDElement<2,3> >
    (curved_boundary_pt,xi_lo,fract_mid,xi_hi);
  
   // Filename for output
   string file_root="adapt_solid_pvd_3";
  
   // Compare output
   output_both_versions(file_root,mesh_pt,nplot);
  
   delete mesh_pt;
  }
 
 
  //--------------------------------------------------
  {
   // Build and assign mesh
   Mesh* mesh_pt = 
    new ElasticRefineableQuarterCircleSectorMesh<RefineableQPVDElement<2,4> >
    (curved_boundary_pt,xi_lo,fract_mid,xi_hi);
  
   // Filename for output
   string file_root="adapt_solid_pvd_4";
    
   // Compare output
   output_both_versions(file_root,mesh_pt,nplot);
  
   delete mesh_pt;
  }
 
  //--------------------------------------------------
  {
   // Build and assign mesh
   Mesh* mesh_pt = 
    new ElasticRefineableQuarterCircleSectorMesh
    <RefineableQPVDElementWithContinuousPressure<2> >(
     curved_boundary_pt,xi_lo,fract_mid,xi_hi);
  
   // Filename for output
   string file_root="adapt_solid_pvd_cp";
  
   // Compare output
   output_both_versions(file_root,mesh_pt,nplot);
  
   delete mesh_pt;
  }
 
 
  //--------------------------------------------------
  {
   // Build and assign mesh
   Mesh* mesh_pt = 
    new ElasticRefineableQuarterCircleSectorMesh
    <RefineableQPVDElementWithPressure<2> >(
     curved_boundary_pt,xi_lo,fract_mid,xi_hi);
  
   // Filename for output
   string file_root="adapt_solid_pvd_p";
  
   // Compare output
   output_both_versions(file_root,mesh_pt,nplot);
  
   delete mesh_pt;
  }
 
 }


//////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////// 



 // # of elements in x-direction
 unsigned n_x=10;
 
 // # of elements in y-direction
 unsigned n_y=10;
 
 // Domain length in x-direction
 double l_x=1.0;
 
 // Domain length in y-direction
 double l_y=2.0;
 
 
 // Create geom object
 double L_up=1.0;
 double L=5.0;
 double L_down=5.0;
 TimeStepper* time_stepper_pt=new Steady<0>;
 
  GeomObject* Wall_pt=new 
   EllipticalTube(1.0,1.0);


 // Boundaries on object: quarter tube
 Vector<double> xi_lo(2);
 xi_lo[0]=0.0; 
 xi_lo[1]=0.0;

 Vector<double> xi_hi(2);
 xi_hi[0]=L_up+L+L_down;
 xi_hi[1]=0.5*4.0*atan(1.0);
 double frac_mid=0.5;
 unsigned nlayer=11; 

 
 //-------------------------------------------------------------
 {
  
  // Build and assign mesh
  Mesh* mesh_pt=new RefineableQuarterTubeMesh<RefineableQTaylorHoodElement<3> >
   (Wall_pt,xi_lo,frac_mid,xi_hi,nlayer,time_stepper_pt);
  
  
  // Filename for output
  string file_root="adapt_taylor_hood_3";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
  
 }
 
 //-------------------------------------------------------------
 {
  
  // Build and assign mesh
  Mesh* mesh_pt=
   new RefineableQuarterTubeMesh<RefineableQCrouzeixRaviartElement<3> >(
    Wall_pt,xi_lo,frac_mid,xi_hi,nlayer,time_stepper_pt);
  
  
  // Filename for output
  string file_root="adapt_crozier_raviart_3";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
  
 }
 
 
 //--------------------------------------------------
 {
  // Build and assign mesh
  Mesh* mesh_pt =
   new SimpleRectangularQuadMesh<QTaylorHoodElement<2> >(n_x,n_y,l_x,l_y);
  
  // Filename for output
  string file_root="taylor_hood_2";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }
 
 
 //--------------------------------------------------
 {
  // Build and assign mesh
  Mesh* mesh_pt =
   new SimpleRectangularQuadMesh<QCrouzeixRaviartElement<2> >(n_x,n_y,l_x,l_y);
  
  // Filename for output
  string file_root="crozier_raviart_2";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }
 
 //--------------------------------------------------
 {
  // Build and assign mesh
  Mesh* mesh_pt = 
   new SimpleRefineableRectangularQuadMesh
   <RefineableQCrouzeixRaviartElement<2> >(n_x,n_y,l_x,l_y);
  
  // Filename for output
  string file_root="adapt_crozier_raviart_2";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }

 
 //--------------------------------------------------
 {
  // Build and assign mesh
  Mesh* mesh_pt = 
   new SimpleRefineableRectangularQuadMesh
   <RefineableQTaylorHoodElement<2> >(n_x,n_y,l_x,l_y);
  
  // Filename for output
  string file_root="adapt_taylor_hood_2";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }
 
 
 //-------------------------------------------------------------
 {
  
  // Build and assign mesh
  Mesh* mesh_pt=new QuarterTubeMesh<QPoissonElement<3,2> >(
   Wall_pt,xi_lo,frac_mid,xi_hi,nlayer,time_stepper_pt);
  
  
  // Filename for output
  string file_root="poisson_3_2";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
  
 }
 
 //-------------------------------------------------------------
 {
  
  // Build and assign mesh
  Mesh* mesh_pt=new QuarterTubeMesh<QPoissonElement<3,3> >(
   Wall_pt,xi_lo,frac_mid,xi_hi,nlayer,time_stepper_pt);
  
  
  // Filename for output
  string file_root="poisson_3_3";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }
 
 //-------------------------------------------------------------
 {
  
  // Build and assign mesh
  Mesh* mesh_pt=new QuarterTubeMesh<QPoissonElement<3,4> >(
   Wall_pt,xi_lo,frac_mid,xi_hi,nlayer,time_stepper_pt);
  
  
  // Filename for output
  string file_root="poisson_3_4";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
  
 }
 
 
 //-------------------------------------------------------------
 {
  
  // Build and assign mesh
  Mesh* mesh_pt=new RefineableQuarterTubeMesh<RefineableQPoissonElement<3,2> >(
   Wall_pt,xi_lo,frac_mid,xi_hi,nlayer,time_stepper_pt);
  
  
  // Filename for output
  string file_root="adapt_poisson_3_2";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
  
 }
 
 //-------------------------------------------------------------
 {
  
  // Build and assign mesh
  Mesh* mesh_pt=new RefineableQuarterTubeMesh<RefineableQPoissonElement<3,3> >(
   Wall_pt,xi_lo,frac_mid,xi_hi,nlayer,time_stepper_pt);
  
  
  // Filename for output
  string file_root="adapt_poisson_3_3";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
  
 }
 
 //-------------------------------------------------------------
 {
  
  // Build and assign mesh
  Mesh* mesh_pt=new RefineableQuarterTubeMesh<RefineableQPoissonElement<3,4> >(
   Wall_pt,xi_lo,frac_mid,xi_hi,nlayer,time_stepper_pt);
  
  
  // Filename for output
  string file_root="adapt_poisson_3_4";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
  
 }
 
 
 //--------------------------------------------------
 {
  // Build mesh and store pointer 
  Mesh* mesh_pt=new OneDMesh<QPoissonElement<1,2> >(nelem,1.0);
  
  // Filename for output
  string file_root="poisson_1_2";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }
 
 
 //--------------------------------------------------
 {
  // Build mesh and store pointer 
  Mesh* mesh_pt=new OneDMesh<QPoissonElement<1,3> >(nelem,1.0);
  
  // Filename for output
  string file_root="poisson_1_3";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }


 //--------------------------------------------------
 {
  // Build mesh and store pointer 
  Mesh* mesh_pt=new OneDMesh<QPoissonElement<1,4> >(nelem,1.0);
  
  // Filename for output
  string file_root="poisson_1_4";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);

  delete mesh_pt;
 }


 //--------------------------------------------------
 {
  // Build and assign mesh
  Mesh* mesh_pt =
   new SimpleRectangularQuadMesh<QPoissonElement<2,2> >(n_x,n_y,l_x,l_y);

  // Filename for output
  string file_root="poisson_2_2";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);

  delete mesh_pt;
 }



 //--------------------------------------------------
 {
  // Build and assign mesh
  Mesh* mesh_pt = 
   new SimpleRectangularQuadMesh<QPoissonElement<2,3> >(n_x,n_y,l_x,l_y);

  // Filename for output
  string file_root="poisson_2_3";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);

  delete mesh_pt;
 }



 //--------------------------------------------------
 {
  // Build and assign mesh
  Mesh* mesh_pt = 
   new SimpleRectangularQuadMesh<QPoissonElement<2,4> >(n_x,n_y,l_x,l_y);

  // Filename for output
  string file_root="poisson_2_4";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);

  delete mesh_pt;
 }


//--------------------------------------------------
 {
  // Build and assign mesh
  Mesh* mesh_pt =
   new SimpleRefineableRectangularQuadMesh<RefineableQPoissonElement<2,2> >
   (n_x,n_y,l_x,l_y);
  
  // Filename for output
  string file_root="adapt_poisson_2_2";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }
 


 //--------------------------------------------------
 {
  // Build and assign mesh
  Mesh* mesh_pt = 
   new SimpleRefineableRectangularQuadMesh<RefineableQPoissonElement<2,3> >
   (n_x,n_y,l_x,l_y);

  // Filename for output
  string file_root="adapt_poisson_2_3";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);

  delete mesh_pt;
 }



 //--------------------------------------------------
 {
  // Build and assign mesh
  Mesh* mesh_pt = 
   new SimpleRefineableRectangularQuadMesh<RefineableQPoissonElement<2,4> >
   (n_x,n_y,l_x,l_y);

  // Filename for output
  string file_root="adapt_poisson_2_4";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);

  delete mesh_pt;
 }



 
//////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////// 


 //-------------------------------------------------------------
 {
  
  // Build and assign mesh
  Mesh* mesh_pt=new QuarterTubeMesh<QUnsteadyHeatElement<3,2> >(
   Wall_pt,xi_lo,frac_mid,xi_hi,nlayer,time_stepper_pt);
  
  
  // Filename for output
  string file_root="unsteady_heat_3_2";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
  
 }
 
 //-------------------------------------------------------------
 {
  
  // Build and assign mesh
  Mesh* mesh_pt=new QuarterTubeMesh<QUnsteadyHeatElement<3,3> >(
   Wall_pt,xi_lo,frac_mid,xi_hi,nlayer,time_stepper_pt);
  
  
  // Filename for output
  string file_root="unsteady_heat_3_3";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
  
 }
 
 //-------------------------------------------------------------
 {
  
  // Build and assign mesh
  Mesh* mesh_pt=new QuarterTubeMesh<QUnsteadyHeatElement<3,4> >(
   Wall_pt,xi_lo,frac_mid,xi_hi,nlayer,time_stepper_pt);
  
  
  // Filename for output
  string file_root="unsteady_heat_3_4";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
  
 }
 
////////////////////////////////////////////////////////////////
 
 
 //--------------------------------------------------
 {
  // Build mesh and store pointer 
  Mesh* mesh_pt=new OneDMesh<QUnsteadyHeatElement<1,2> >(nelem,1.0);
  
  // Filename for output
  string file_root="unsteady_heat_1_2";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }
 
 
 //--------------------------------------------------
 {
  // Build mesh and store pointer 
  Mesh* mesh_pt=new OneDMesh<QUnsteadyHeatElement<1,3> >(nelem,1.0);
  
  // Filename for output
  string file_root="unsteady_heat_1_3";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }
 
 
 //--------------------------------------------------
 {
  // Build mesh and store pointer 
  Mesh* mesh_pt=new OneDMesh<QUnsteadyHeatElement<1,4> >(nelem,1.0);
  
  // Filename for output
  string file_root="unsteady_heat_1_4";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }
 
 
 
 //--------------------------------------------------
 {
  // Build and assign mesh
  Mesh* mesh_pt =
   new SimpleRectangularQuadMesh<QUnsteadyHeatElement<2,2> >(n_x,n_y,l_x,l_y);
  
  // Filename for output
  string file_root="unsteady_heat_2_2";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }
 
 //-------------------------------------------------- 
 {
  // Build and assign mesh
  Mesh* mesh_pt = 
   new SimpleRectangularQuadMesh<QUnsteadyHeatElement<2,3> >(n_x,n_y,l_x,l_y);
  
  // Filename for output
  string file_root="unsteady_heat_2_3";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }
 
 
 
 //--------------------------------------------------
 {
  // Build and assign mesh
  Mesh* mesh_pt = 
   new SimpleRectangularQuadMesh<QUnsteadyHeatElement<2,4> >(n_x,n_y,l_x,l_y);
  
  // Filename for output
  string file_root="unsteady_heat_2_4";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }
 
 
 //--------------------------------------------------
 {
  // Build and assign mesh
  Mesh* mesh_pt =
   new SimpleRefineableRectangularQuadMesh
   <RefineableQUnsteadyHeatElement<2,2> >(n_x,n_y,l_x,l_y);
  
  // Filename for output
  string file_root="adapt_unsteady_heat_2_2";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }
 
 
 
 //--------------------------------------------------
 {
  // Build and assign mesh
  Mesh* mesh_pt = 
   new SimpleRefineableRectangularQuadMesh
   <RefineableQUnsteadyHeatElement<2,3> >(n_x,n_y,l_x,l_y);
  
  // Filename for output
  string file_root="adapt_unsteady_heat_2_3";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }
 
 
 
 //--------------------------------------------------
 {
  // Build and assign mesh
  Mesh* mesh_pt = 
   new SimpleRefineableRectangularQuadMesh
   <RefineableQUnsteadyHeatElement<2,4> >(n_x,n_y,l_x,l_y);
  
  // Filename for output
  string file_root="adapt_unsteady_heat_2_4";
  
  // Compare output
  output_both_versions(file_root,mesh_pt,nplot);
  
  delete mesh_pt;
 }
 
}


