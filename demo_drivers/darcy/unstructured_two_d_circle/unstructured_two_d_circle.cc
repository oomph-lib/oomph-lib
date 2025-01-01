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

#include "generic.h"
#include "darcy.h"
#include "meshes/triangle_mesh.h"

using namespace oomph;
using namespace std;


//==================================================================
/// Namespace for exact solution and problem parameters
//==================================================================
namespace TestProblem
{

 /// Flag to indicate that we want the forced solution a "point source"
 /// in the mass flux
 bool Use_point_source_solution=false;

 /// Source function
 void source(const Vector<double> &x, Vector<double> &f)
  {
   if (Use_point_source_solution)
    {
     f[0]=0.0;
     f[1]=0.0;
    }
   else
    {
     f[0]=x[0]*x[1]-x[1]*x[1]+2;
     f[1]=x[0]+x[0]*x[0]-0.5*x[1]*x[1]+3;
    }
  }


 /// Mass source strength
 double Mass_source_strength=1000.0;

 /// Mass source strength
 double Mass_source_decay_exponent=50.0;

 /// Mass source
 void mass_source(const Vector<double> &x, double &m)
  {
   if (Use_point_source_solution)
    {
     double r=sqrt(x[0]*x[0]+x[1]*x[1]);
     m=Mass_source_strength*exp(-r*Mass_source_decay_exponent);
    }
   else
    {
     m=0.0;
    }
  }

 /// Pressure around the boundary of the domain
 void boundary_pressure(const double &time,
                        const Vector<double> &x,
                        const Vector<double> &n,
                        double &result)
  {
   if (Use_point_source_solution)
    {
     result=0.0;
    }
   else
    {
     result=2*x[0]+3*x[1]-3.0/2.0;
    }
  }

 /// Exact solution: q1,q2,div_q,p
 void exact_soln(const Vector<double> &x, Vector<double> &soln)
  {

   if (Use_point_source_solution)
    {
     // q[0]
     soln[0]=0.0;
     // q[1]
     soln[1]=0.0;
     // div q
     soln[2]=0.0;
     // p
     soln[3]=0.0;
    }
   else
    {
     // q[0]
     soln[0]=x[0]*x[1]-x[1]*x[1];
     // q[1]
     soln[1]=x[0]+x[0]*x[0]-0.5*x[1]*x[1];
     // div q
     soln[2]=0;
     // p
     soln[3]=2*x[0]+3*x[1]-3.0/2.0;
    }
  }

 /// Target element area for Triangle
 double Element_area=0.01;

 /// The directory in which the solution is output
 std::string Directory="RESLT";

 /// Global function that completes the edge sign setup --
 /// has to be called before projection in unstructured
 /// adaptation
 template<class ELEMENT>
 void edge_sign_setup(Mesh* mesh_pt)
 {
  // The dictionary keeping track of edge signs
  std::map<Edge,unsigned> assignments;
  
  // Loop over all elements
  unsigned n_element = mesh_pt->nelement();
  for(unsigned e=0;e<n_element;e++)
   {
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt->element_pt(e));
  
    // Assign edge signs: Loop over the vertex nodes (always 
    // first 3 nodes for triangles)
    for(unsigned i=0;i<3;i++)
     {
      Node *nod_pt1, *nod_pt2;
      nod_pt1 = el_pt->node_pt(i);
      nod_pt2 = el_pt->node_pt((i+1)%3);
      Edge edge(nod_pt1,nod_pt2);
      unsigned status = assignments[edge];
      
      // This relies on the default value for an int being 0 in a map
      switch(status)
       {
        // If not assigned on either side, give it a + on current side
       case 0:
        assignments[edge]=1;
        break;
        // If assigned + on other side, give it a - on current side
       case 1:
        assignments[edge]=2;
        el_pt->sign_edge(i)=-1;
        break;
        // If assigned - on other side, give it a + on current side
       case 2:
        assignments[edge]=1;
        break;
       }
     } // end of loop over vertex nodes
   } // end of loop over elements


 }
}


//==================================================================
/// Darcy problem
//==================================================================
template<class ELEMENT>
class DarcyProblem : public Problem
{

public:

 /// Constructor
 DarcyProblem();

 /// Post-process results, using unsigned to label files
 void doc_solution(const unsigned &label);

 /// Actions before Newton solve (empty)
 void actions_before_newton_solve() {}

#ifdef ADAPTIVE

 /// Actions before adapt
 void actions_before_adapt()
  {
   doc_solution(20);
   oomph_info << "Done output in actions before adapt\n";

   // Loop over the surface elements
   unsigned n_element = Surface_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Surface_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Surface_mesh_pt->flush_element_and_node_storage();
   
   //  Rebuild the global mesh 
   rebuild_global_mesh();
  }
 
 /// Actions after adapt
 void actions_after_adapt()
  {
   // Re-create pressure bc elements
   create_pressure_elements();

   //  Rebuild the global mesh 
   rebuild_global_mesh();

   // Complete problem setup
   complete_problem_setup();

   doc_solution(21);
   oomph_info << "Done output in actions after adapt\n";
  }

#endif

 /// Doc shape functions
 void doc_shape_functions();

private:

 /// Allocate pressure elements on the top surface
 void create_pressure_elements();

 /// Complete problem setup
 void complete_problem_setup();

#ifdef ADAPTIVE

 /// Pointer to the refineable "bulk" mesh
 RefineableTriangleMesh<ELEMENT>* Bulk_mesh_pt;

#else

 /// Pointer to the "bulk" mesh
 TriangleMesh<ELEMENT>* Bulk_mesh_pt;

#endif

 /// Mesh of boundary condition elements
 Mesh *Surface_mesh_pt;

 /// Trace file for results
 std::ofstream Trace_file;

};





//==================================================================
/// Constructor
//==================================================================
template<class ELEMENT>
DarcyProblem<ELEMENT>::DarcyProblem()
{
 // Mesh
 //-----

 // Closed curve definiting outer boundary
 TriangleMeshClosedCurve *outer_boundary_pt=0;

 // Internal boundary to prevent problem with flux B.C.s
 Vector<TriangleMeshOpenCurve*> inner_open_boundary_pt;
 
 Vector<double> zeta(1);
 Vector<double> posn(2);
 
 Ellipse* outer_boundary_ellipse_pt = new Ellipse(1,1);
 
 Vector<TriangleMeshCurveSection*> outer_curvilinear_boundary_pt(2);
 
 // First bit
 double zeta_start = 0.0;
 double zeta_end = MathematicalConstants::Pi;
 unsigned nsegment =
  (unsigned)max(
   3.0,MathematicalConstants::Pi/std::sqrt(TestProblem::Element_area));
 outer_curvilinear_boundary_pt[0] =
  new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
                            zeta_end, nsegment, 0);
 
 // Second bit
 zeta_start = MathematicalConstants::Pi;
 zeta_end = 2.0*MathematicalConstants::Pi;
 nsegment =
  (unsigned)max(
   3.0,MathematicalConstants::Pi/std::sqrt(TestProblem::Element_area));
 outer_curvilinear_boundary_pt[1] =
  new TriangleMeshCurviLine(outer_boundary_ellipse_pt, zeta_start,
                            zeta_end, nsegment, 1);
 
 outer_boundary_pt =
  new TriangleMeshClosedCurve(outer_curvilinear_boundary_pt);

 // Target element area for Triangle
 double uniform_element_area = TestProblem::Element_area;

 // Use the TriangleMeshParameters object for gathering all
 // the necessary arguments for the TriangleMesh object
 TriangleMeshParameters triangle_mesh_parameters(
   outer_boundary_pt);

 // Take the maximum element area
 triangle_mesh_parameters.element_area() =
  uniform_element_area;


#ifdef ADAPTIVE
 
 // Build adaptive "bulk" mesh
 Bulk_mesh_pt=new RefineableTriangleMesh<ELEMENT>(triangle_mesh_parameters);

 // Create/set error estimator
 Bulk_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 
 // Choose error tolerances 
 if (TestProblem::Use_point_source_solution)
  {
   Bulk_mesh_pt->min_permitted_error()=0.0001; 
   Bulk_mesh_pt->max_permitted_error()=0.001; 
  }
 else
  {
   Bulk_mesh_pt->min_permitted_error()=1.0e-8;
   Bulk_mesh_pt->max_permitted_error()=1.0e-7;
  }

 // Actions before projection
 Bulk_mesh_pt->mesh_update_fct_pt()=&TestProblem::edge_sign_setup<ELEMENT>;

#else

 // Create the buk mesh
 Bulk_mesh_pt = new TriangleMesh<ELEMENT>(triangle_mesh_parameters);

#endif

 // Create mesh of pressure bc elements
 Surface_mesh_pt = new Mesh;
 create_pressure_elements();

 // Combine into global mesh
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_pt);
 build_global_mesh();

 // Complete problem setup
 complete_problem_setup();

 // Assign equation numbers
 oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;


 // Open trace file
 char filename[50];
 sprintf(filename,"%s/trace.dat",TestProblem::Directory.c_str());
 Trace_file.open(filename);

} // end of problem constructor



//===start_of_create_pressure_elements===============================
/// Creat traction elements
//===================================================================
template<class ELEMENT>
void DarcyProblem<ELEMENT>::create_pressure_elements()
{

 // Upper part of boundary
 unsigned bound=0;
  
 // How many bulk elements are next to that boundary 
 unsigned n_neigh = Bulk_mesh_pt->nboundary_element(bound);
 
 // Now loop over bulk elements and create the face elements
 for(unsigned n=0;n<n_neigh;n++)
  {
   // Create the face element
   DarcyFaceElement<ELEMENT>* pressure_element_pt
    = new DarcyFaceElement<ELEMENT>
    (Bulk_mesh_pt->boundary_element_pt(bound,n),
     Bulk_mesh_pt->face_index_at_boundary(bound,n));
   
   // Add to mesh
   Surface_mesh_pt->add_element_pt(pressure_element_pt);

   // Set function pointer
   pressure_element_pt->pressure_fct_pt()=&TestProblem::boundary_pressure;

  }

} // end of assign_traction_elements



//====================================================================
/// Complete problem setup
//====================================================================
template<class ELEMENT>
void DarcyProblem<ELEMENT>::complete_problem_setup()
{
 // Setup the signs for the fluxes
 TestProblem::edge_sign_setup<ELEMENT>(Bulk_mesh_pt);

 // Loop over all elements
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
   
   // Set source fcts
   el_pt->source_fct_pt()=TestProblem::source;
   el_pt->mass_source_fct_pt()=TestProblem::mass_source;
  
#ifdef ADAPTIVE

   // Pin dofs at vertex nodes (because they're only used for projection)
   el_pt->pin_superfluous_darcy_dofs();

#endif

  } // end of loop over elements


 // Apply boundary conditions
 //--------------------------

 // Get the nodal index at which values representing edge fluxes
 // at flux interpolation points are stored
 ELEMENT *el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(0));
 
 // Where are the values stored?
 unsigned n=el_pt->nedge_flux_interpolation_point();
 Vector<unsigned> q_index(n);
 for (unsigned j=0;j<n;j++)
  {
   q_index[j]=el_pt->q_edge_index(j);
  }
 
 // Coordinate vector
 Vector<double> x(2);

 // Impose normal component of flux on lower part of outer boundary
 unsigned ibound=1;
 {
  // Get the number of elements along boundary ibound
  unsigned n_boundary_element=Bulk_mesh_pt->nboundary_element(ibound);
  
  // Loop over the elements along boundary ibound
  for(unsigned e=0;e<n_boundary_element;e++)
   {
    // Upcast the current element to the actual type
    ELEMENT *el_pt=
     dynamic_cast<ELEMENT*>(Bulk_mesh_pt->boundary_element_pt(ibound,e));
    
    // Loop over the edges
    for(unsigned edge=0;edge<3;edge++)
     {
      // Get pointer to node that stores the edge flux dofs for this edge
      Node* nod_pt=el_pt->edge_flux_node_pt(edge);
      
      // Pin/set values for the flux degrees of freedom
      if (nod_pt->is_on_boundary(ibound))
       {
        for (unsigned j=0;j<n;j++)
         {
          nod_pt->pin(q_index[j]);
         }

        // Get face index of face associated with edge
        unsigned f=el_pt->face_index_of_edge(edge);
        
        // Build a temporary face element from which we'll extract
        // the outer unit normal
        DarcyFaceElement<ELEMENT>* face_el_pt=
         new DarcyFaceElement<ELEMENT>(el_pt,f);   
        
        // Loop over the flux interpolation points
        unsigned n_flux_interpolation_points=
         el_pt->nedge_flux_interpolation_point();
        for(unsigned g=0;g<n_flux_interpolation_points;g++)
         {
          // Get the global coords of the flux_interpolation point
          el_pt->edge_flux_interpolation_point_global(edge,g,x);
          
          // Get the exact solution
          Vector<double> exact_soln(4,0.0);
          TestProblem::exact_soln(x,exact_soln);

          // Get unit normal at this flux interpolation point
          Vector<double> s(1);
          el_pt->face_local_coordinate_of_flux_interpolation_point(edge,g,s);
          Vector<double> unit_normal(2);
          face_el_pt->outer_unit_normal(s,unit_normal);
          
#ifdef PARANOID          

          // Sanity check
          Vector<double> x_face(2);
          face_el_pt->interpolated_x(s,x_face);
          if ((x_face[0]-x[0])*(x_face[0]-x[0])+
              (x_face[1]-x[1])*(x_face[1]-x[1])>1.0e-3)
           {
            std::stringstream error;
            error << "Discrepancy in coordinate of flux interpolation point\n"
                  << "(computed by bulk and face elements) for edge " << e 
                  << " and flux int pt " << g << "\n"
                  << "Face thinks node is at: "
                  << x_face[0] << " " << x_face[1] << "\n"
                  << "Bulk thinks node is at: "
                  << x[0] << " " << x[1] << "\n";
            throw OomphLibError(
             error.str(),
             OOMPH_CURRENT_FUNCTION,
             OOMPH_EXCEPTION_LOCATION);
           }
#endif
          
          // Set the boundary flux
          nod_pt->set_value(q_index[g],
                            exact_soln[0]*unit_normal[0]+
                            exact_soln[1]*unit_normal[1]);
          
         } // End of loop over flux interpolation points

        // Don't need face element on that edge any more
        delete face_el_pt;
        
       } // End if for edge on required boundary
     } // End of loop over edges

   } // End of loop over boundary elements
 } // End of loop over boundaries
}

//========================================================================
/// Doc shape functions
//========================================================================
template<class ELEMENT>
void DarcyProblem<ELEMENT>::doc_shape_functions()
{
 ofstream some_file;
 char filename[100];

 // Coarse solution
 unsigned npts_coarse=2;
 sprintf(filename,"%s/coarse_soln.dat",TestProblem::Directory.c_str());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts_coarse);
 some_file.close();

 // Backup for all values (double counting plenty but I don't care...)
 //-------------------------------------------------------------------
 Vector<double> backed_up_value;
 unsigned nel=Bulk_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   // Backup edge-based q values
   unsigned n_q_basis_edge = el_pt->nq_basis_edge();
   for (unsigned j=0;j<n_q_basis_edge;j++)
    {
     backed_up_value.push_back(el_pt->q_edge(j));
    }

   // Backup internal q values
   unsigned n_q_basis = el_pt->nq_basis();
   for (unsigned j=n_q_basis_edge;j<n_q_basis;j++)
    {
     backed_up_value.push_back(el_pt->q_internal(j-n_q_basis_edge));
    }

   // Backup pressures
   unsigned n_p_basis = el_pt->np_basis();
   for (unsigned j=0;j<n_p_basis;j++)
    {
     backed_up_value.push_back(el_pt->p_value(j));
    }
  }
 


 // Loop over all elements to wipe
 //-------------------------------
 for (unsigned e=0;e<nel;e++)
  {
   ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   // Kill edge-based q values
   unsigned n_q_basis_edge = el_pt->nq_basis_edge();
   for (unsigned j=0;j<n_q_basis_edge;j++)
    {
     el_pt->set_q_edge(j,0.0);
    }

   // Kill internal q values
   unsigned n_q_basis = el_pt->nq_basis();
   for (unsigned j=n_q_basis_edge;j<n_q_basis;j++)
    {
     el_pt->set_q_internal(j-n_q_basis_edge,0.0);
    }


   // Kill pressures
   unsigned n_p_basis = el_pt->np_basis();
   for (unsigned j=0;j<n_p_basis;j++)
    {
     el_pt->set_p_value(j,0.0);
    }
  }

 
 // Choose a specific element for which we'll switch on all the dofs
 //-----------------------------------------------------------------
 // one by one
 //-----------
 MeshAsGeomObject* mesh_geom_obj_pt=new MeshAsGeomObject(Bulk_mesh_pt); 
 Vector<double> x(2,0.0);
 Vector<double> s(2,0.0);

 // Get pointer to GeomObject that contains this point and cast
 GeomObject* geom_obj_pt=0;
 mesh_geom_obj_pt->locate_zeta(x,geom_obj_pt,s); 
 ELEMENT* el_pt=dynamic_cast<ELEMENT*>(geom_obj_pt);
 

 // Number of plot points
 unsigned npts=10;
   
 // Doc edge-based q values
 unsigned count_shape_doc=0;
 unsigned n_q_basis_edge = el_pt->nq_basis_edge();
 for (unsigned j=0;j<n_q_basis_edge;j++)
  {
   // Set this dof to one
   el_pt->set_q_edge(j,1.0);

   // Get face index of face associated with this flux basis fct
   // and erect face element
   unsigned f=el_pt->face_index_of_q_edge_basis_fct(j);
   DarcyFaceElement<ELEMENT>* face_el_pt=
    new DarcyFaceElement<ELEMENT>(el_pt,f);   
   
   // Get outer unit normal halfway along edge
   Vector<double> s(1,0.5);
   Vector<double> unit_normal(2);
   Vector<double> x(2);
   face_el_pt->outer_unit_normal(s,unit_normal);
   face_el_pt->interpolated_x(s,x);
   delete face_el_pt;

   // Plot
   sprintf(filename,"%s/q_shape_fct%i.dat",TestProblem::Directory.c_str(),
           count_shape_doc++);
   some_file.open(filename);
   unsigned nel=Bulk_mesh_pt->nelement();
   for (unsigned e=0;e<nel;e++)
    {
     dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
      output(some_file,npts);
    }
   some_file.close();


   sprintf(filename,"%s/outer_unit_normal%i.dat",
           TestProblem::Directory.c_str(),j);
   some_file.open(filename);
   some_file << x[0] << " " 
             << x[1] << " " 
             << unit_normal[0] << " " 
             << unit_normal[1] << " " 
             << std::endl;
   some_file.close();


   sprintf(filename,"%s/flux_interpolation_point%i.dat",
           TestProblem::Directory.c_str(),j);
   some_file.open(filename);
   Vector<double> flux_interpolation_point(2);
   el_pt->edge_flux_interpolation_point_global(j,flux_interpolation_point);
   some_file << flux_interpolation_point[0] << " " 
             << flux_interpolation_point[1] << " " 
             << 1.0 << " " 
             << std::endl;
   some_file.close();

   sprintf(filename,"%s/q_shape_fct_with_project%i.dat",
           TestProblem::Directory.c_str(),j);
   some_file.open(filename);
   for (unsigned e=0;e<nel;e++)
    {
     dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
      output_with_projected_flux(some_file,npts,unit_normal);
    }
   some_file.close();

   // Reset
   el_pt->set_q_edge(j,0.0);
  }
 
 // Doc internal q values
 unsigned n_q_basis = el_pt->nq_basis();
 for (unsigned j=n_q_basis_edge;j<n_q_basis;j++)
  {
   el_pt->set_q_internal(j-n_q_basis_edge,1.0);
   sprintf(filename,"%s/q_shape_fct%i.dat",
           TestProblem::Directory.c_str(),
           count_shape_doc++);
   some_file.open(filename);
   unsigned nel=Bulk_mesh_pt->nelement();
   for (unsigned e=0;e<nel;e++)
    {
     dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
      output(some_file,npts);
    }
   some_file.close();
   el_pt->set_q_internal(j-n_q_basis_edge,0.0);
  }

 // Doc pressures
 count_shape_doc=0;
 unsigned n_p_basis = el_pt->np_basis();
 for (unsigned j=0;j<n_p_basis;j++)
  {
   el_pt->set_p_value(j,1.0);
   sprintf(filename,"%s/p_shape_fct%i.dat",
           TestProblem::Directory.c_str(),
           count_shape_doc++);
   some_file.open(filename);
   unsigned nel=Bulk_mesh_pt->nelement();
   for (unsigned e=0;e<nel;e++)
    {
     dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
      output(some_file,npts);
    }
   some_file.close();
   el_pt->set_p_value(j,0.0);
  }
 

 
 // Restore dofs
 unsigned count=0;
 for (unsigned e=0;e<nel;e++)
  {
   ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   // Reset edge-based q values
   unsigned n_q_basis_edge = el_pt->nq_basis_edge();
   for (unsigned j=0;j<n_q_basis_edge;j++)
    {
     el_pt->set_q_edge(j,backed_up_value[count]);
     count++;
    }

   // Reset internal q values
   unsigned n_q_basis = el_pt->nq_basis();
   for (unsigned j=n_q_basis_edge;j<n_q_basis;j++)
    {
     el_pt->set_q_internal(j-n_q_basis_edge,backed_up_value[count]);
     count++;
    }


   // Reset pressures
   unsigned n_p_basis = el_pt->np_basis();
   for (unsigned j=0;j<n_p_basis;j++)
    {
     el_pt->set_p_value(j,backed_up_value[count]);
     count++;
    }
  }
} 

//========================================================================
/// Write the solution and exact solution to file, and calculate the error
//========================================================================
template<class ELEMENT>
void DarcyProblem<ELEMENT>::doc_solution(const unsigned &label)
{
 ofstream some_file;
 char filename[100];


 // Output computed solution
 //-------------------------
 unsigned npts=5;
 sprintf(filename,"%s/soln%i.dat",TestProblem::Directory.c_str(),label);
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();


 // Output coarse solution
 //-----------------------
 unsigned npts_coarse=2;
 sprintf(filename,"%s/coarse_soln%i.dat",TestProblem::Directory.c_str(),label);
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts_coarse);
 some_file.close();

 // Output exact solution
 //----------------------
 sprintf(filename,"%s/exact_soln%i.dat",TestProblem::Directory.c_str(),label);
 some_file.open(filename);
 Bulk_mesh_pt->output_fct(some_file,
                            npts,
                            TestProblem::exact_soln);
 some_file.close();

 // Output boundary condition elements
 //-----------------------------------
 sprintf(filename,"%s/bc_elements%i.dat",TestProblem::Directory.c_str(),label);
 some_file.open(filename);
 Surface_mesh_pt->output(some_file,npts);
 some_file.close();

 // Doc error
 //----------
 Vector<double> norm(2,0.0);
 Vector<double> error(2,0.0);
 sprintf(filename,"%s/error%i.dat",TestProblem::Directory.c_str(),label);
  some_file.open(filename);
 Bulk_mesh_pt->compute_error(some_file,
                             TestProblem::exact_soln,
                             error,norm);
 some_file.close();

 // Doc error norm:
 cout << std::endl;
 cout << "Norm of exact q : " << sqrt(norm[0]) << std::endl;
 cout << "Norm of exact p : " << sqrt(norm[1]) << std::endl << std::endl;

 cout << "Norm of q error : " << sqrt(error[0]) << std::endl;
 cout << "Norm of p error : " << sqrt(error[1]) << std::endl;
 cout << std::endl << std::endl;

 Trace_file << ndof() << " "
            << sqrt(error[0]) << " "
            << sqrt(error[1]) << " "
            << TestProblem::Element_area << " "
            << sqrt(norm[0]) << " "
            << sqrt(norm[1]) << std::endl;
}



//=====================================================================
/// Driver code for Darcy test problem
//=====================================================================
int main(int argc, char **argv)
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified

 // Output directory
 CommandLineArgs::specify_command_line_flag("--dir",
                                            &TestProblem::Directory);

 // Specify target element area
 CommandLineArgs::specify_command_line_flag("--element_area",
                                            &TestProblem::Element_area);


// Use "point" mass source
 CommandLineArgs::specify_command_line_flag("--point_mass_source");


// Doc shape functions etc then quit
 CommandLineArgs::specify_command_line_flag("--doc_shape_functions");

 // Parse command line
 CommandLineArgs::parse_and_assign();

 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();


 // Doc the shape functions -- use smaller target area
 if (CommandLineArgs::command_line_flag_has_been_set("--doc_shape_functions"))
  {
   oomph_info << "Overwriting element area for doc of shape functions\n";
   TestProblem::Element_area=0.1;
  }

 // Use point mass source?
 if (CommandLineArgs::command_line_flag_has_been_set("--point_mass_source"))
  {
   oomph_info << "Using point mass source\n";
   TestProblem::Use_point_source_solution=true;
  }

#ifdef ADAPTIVE

 // Create problem
 DarcyProblem<ProjectableDarcyElement<TRaviartThomasDarcyElement<1> > > problem;

#else

 // Create problem
 DarcyProblem<TRaviartThomasDarcyElement<1> > problem;

#endif

 // Output initial guess for solution
 problem.doc_solution(0);

 // Doc the shape functions
 if (CommandLineArgs::command_line_flag_has_been_set("--doc_shape_functions"))
  {
   problem.doc_shape_functions();
   exit(0);
  }

#ifdef ADAPTIVE
 
 // Max. number of adaptations
 unsigned max_adapt=1;
 if (TestProblem::Use_point_source_solution)
  {
   max_adapt=3; 
  }

 // Solve the problem with the adaptive Newton's method, allowing
 // up to max_adapt mesh adaptations after every solve.
 problem.newton_solve(max_adapt);
 
#else

 // Solve 
 problem.newton_solve();

#endif

 // Doc result
 problem.doc_solution(1);

 return 0;
}

