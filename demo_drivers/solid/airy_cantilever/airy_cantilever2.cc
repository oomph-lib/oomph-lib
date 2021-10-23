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
//Driver for Airy cantilever beam problem

//Oomph-lib includes
#include "generic.h"
#include "solid.h"
#include "constitutive.h"

//The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

//#define REFINE

namespace oomph
{

//=================start_wrapper==================================
/// Wrapper class for solid elements to modify their output 
/// functions.
//================================================================
template <class ELEMENT>
class MySolidElement : public virtual ELEMENT
{

public:

 /// Constructor: Call constructor of underlying element
 MySolidElement() : ELEMENT() {};

 /// Overload output function:
 void output(std::ostream &outfile, const unsigned &n_plot)
  {

   // Element dimension
   unsigned el_dim = this->dim();

   Vector<double> s(el_dim);
   Vector<double> x(el_dim);
   DenseMatrix<double> sigma(el_dim,el_dim);
   
   switch(el_dim)
    {
     
    case 2:

     //Tecplot header info 
     outfile << "ZONE I=" << n_plot << ", J=" << n_plot << std::endl;
     
     //Loop over element nodes
     for(unsigned l2=0;l2<n_plot;l2++)
      {
       s[1] = -1.0 + l2*2.0/(n_plot-1);
       for(unsigned l1=0;l1<n_plot;l1++)
        {
         s[0] = -1.0 + l1*2.0/(n_plot-1);
         
         // Get Eulerian coordinates and stress
         this->interpolated_x(s,x);
         this->get_stress(s,sigma);

         //Output the x,y,..
         for(unsigned i=0;i<el_dim;i++) 
          {outfile << x[i] << " ";}

         // Output stress
         outfile << sigma(0,0) << " "
                 << sigma(1,0) << " "
                 << sigma(1,1) << " "
                 << std::endl;
        }
      }

     break;
     
    default:

     std::ostringstream error_message;
     error_message << "Output for dim !=2 not implemented" << std::endl;
     throw OomphLibError(error_message.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
  
  }

};



//===========start_face_geometry==============================================
/// FaceGeometry of wrapped element is the same as the underlying element
//============================================================================
template<class ELEMENT>
class FaceGeometry<MySolidElement<ELEMENT> > :
 public virtual FaceGeometry<ELEMENT>
{

public:

 /// Constructor [this was only required explicitly
 /// from gcc 4.5.2 onwards...]
 FaceGeometry() : FaceGeometry<ELEMENT>() {}

};


}






///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Physical_Variables
{

 /// Pointer to strain energy function
 StrainEnergyFunction*Strain_energy_function_pt;

 /// "Mooney Rivlin" coefficient for generalised Mooney Rivlin law
 double C1=1.3;

 /// "Mooney Rivlin" coefficient for generalised Mooney Rivlin law
 double C2=1.3;

 /// Half height of beam
 double H=0.5;

 /// Length of beam
 double L=10.0;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt;

 /// Elastic modulus
 double E=1.0;

 /// Poisson's ratio
 double Nu=0.3;

 /// Uniform pressure
 double P = 0.0;

 /// Constant pressure load. The arguments to this function are imposed
 /// on us by the SolidTractionElements which allow the traction to 
 /// depend on the Lagrangian and Eulerian coordinates x and xi, and on the 
 /// outer unit normal to the surface. Here we only need the outer unit
 /// normal.
 void constant_pressure(const Vector<double> &xi, const Vector<double> &x,
                        const Vector<double> &n, Vector<double> &traction)
 {
  unsigned dim = traction.size();
  for(unsigned i=0;i<dim;i++)
   {
    traction[i] = -P*n[i];
   }
 } // end traction


 /// Non-dim gravity
 double Gravity=0.0;

 /// Non-dimensional gravity as body force
 void gravity(const double& time, 
              const Vector<double> &xi, 
              Vector<double> &b)
 {
  b[0]=0.0;
  b[1]=-Gravity;
 }
 
} //end namespace



//=============begin_problem============================================ 
/// Problem class for the cantilever "beam" structure.
//====================================================================== 
template<class ELEMENT>
class CantileverProblem : public Problem
{

public:

 /// Constructor:
 CantileverProblem(const bool& incompress, const bool& use_fd);
 
 /// Update function (empty)
 void actions_after_newton_solve() {}

 /// Update function (empty)
 void actions_before_newton_solve() {}

#ifdef REFINE

 /// Access function for the solid mesh
 ElasticRefineableRectangularQuadMesh<ELEMENT>*& solid_mesh_pt() 
  {return Solid_mesh_pt;} 

#else

 /// Access function for the solid mesh
 ElasticRectangularQuadMesh<ELEMENT>*& solid_mesh_pt() 
  {return Solid_mesh_pt;} 

#endif


 /// Access function to the mesh of surface traction elements
 SolidMesh*& traction_mesh_pt(){return Traction_mesh_pt;} 

 /// Actions before adapt: Wipe the mesh of traction elements
 void actions_before_adapt();

 /// Actions after adapt: Rebuild the mesh of traction elements
 void actions_after_adapt();

 /// Doc the solution
 void doc_solution();

 /// Run the job -- doc in RESLTi_case
 void run_it(const unsigned& i_case);

private:

 /// Pass pointer to traction function to the
 /// elements in the traction mesh
 void set_traction_pt();

 /// Create traction elements
 void create_traction_elements();

 /// Delete traction elements
 void delete_traction_elements();

 /// Trace file
 ofstream Trace_file;
 
 /// Pointers to node whose position we're tracing
 Node* Trace_node_pt;


#ifdef REFINE

 /// Pointer to solid mesh
 ElasticRefineableRectangularQuadMesh<ELEMENT>* Solid_mesh_pt;

#else

 /// Pointer to solid mesh
 ElasticRectangularQuadMesh<ELEMENT>* Solid_mesh_pt;

#endif


 /// Pointers to meshes of traction elements
 SolidMesh* Traction_mesh_pt;

 /// DocInfo object for output
 DocInfo Doc_info;
 
};


//===========start_of_constructor======================================= 
/// Constructor: 
//====================================================================== 
template<class ELEMENT>
CantileverProblem<ELEMENT>::CantileverProblem(const bool& incompress,
                                              const bool& use_fd) 
{

 // Create the mesh

 // # of elements in x-direction
 unsigned n_x=20;

 // # of elements in y-direction
 unsigned n_y=2;

 // Domain length in x-direction
 double l_x= Global_Physical_Variables::L;

 // Domain length in y-direction
 double l_y=2.0*Global_Physical_Variables::H;

 // Shift mesh downwards so that centreline is at y=0:
 Vector<double> origin(2);
 origin[0]=0.0;
 origin[1]=-0.5*l_y;

#ifdef REFINE

 //Now create the mesh 
 solid_mesh_pt() = new ElasticRefineableRectangularQuadMesh<ELEMENT>(
  n_x,n_y,l_x,l_y,origin);
 
 // Set error estimator
 dynamic_cast<ElasticRefineableRectangularQuadMesh<ELEMENT>* >(
  solid_mesh_pt())->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 
#else
 
 //Now create the mesh 
 solid_mesh_pt() = new ElasticRectangularQuadMesh<ELEMENT>(
  n_x,n_y,l_x,l_y,origin);

#endif
 

 //Assign the physical properties to the elements before any refinement
 //Loop over the elements in the main mesh
 unsigned n_element =solid_mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(solid_mesh_pt()->element_pt(i));
   
   // Set the constitutive law
   el_pt->constitutive_law_pt() =
    Global_Physical_Variables::Constitutive_law_pt;

   //Set the body force
   el_pt->body_force_fct_pt() = Global_Physical_Variables::gravity;
   
   // Get Jacobian by FD?
   if(use_fd)
    {
     el_pt->enable_evaluate_jacobian_by_fd();
    }
   else
    {
     el_pt->disable_evaluate_jacobian_by_fd();
    }

   // Is it incompressible 
   if (incompress)
    {     
     PVDEquationsWithPressure<2>* test_pt = 
      dynamic_cast<PVDEquationsWithPressure<2>*>(
       solid_mesh_pt()->element_pt(i));
     if (test_pt!=0)
      {
       test_pt->set_incompressible();
      }
    }
  }


 // Choose a control node: The last node in the solid mesh
 unsigned nnod=solid_mesh_pt()->nnode();
 Trace_node_pt=solid_mesh_pt()->node_pt(nnod-1);

#ifdef REFINE
 
 // Refine the mesh uniformly
 dynamic_cast<ElasticRefineableRectangularQuadMesh<ELEMENT>* >(
  solid_mesh_pt())->refine_uniformly();
  
#endif
 
 // Construct the traction element mesh
 Traction_mesh_pt=new SolidMesh;
 create_traction_elements();
 
 // Pass pointer to traction function to the elements
 // in the traction mesh
 set_traction_pt();
 
 // Solid mesh is first sub-mesh
 add_sub_mesh(solid_mesh_pt());

 // Add traction sub-mesh
 add_sub_mesh(traction_mesh_pt());

 // Build combined "global" mesh
 build_global_mesh();
 
 // Pin the left boundary (boundary 3) in both directions
 unsigned n_side = mesh_pt()->nboundary_node(3);
 
 //Loop over the nodes
 for(unsigned i=0;i<n_side;i++)
  {
   solid_mesh_pt()->boundary_node_pt(3,i)->pin_position(0);
   solid_mesh_pt()->boundary_node_pt(3,i)->pin_position(1);
  }


 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
  solid_mesh_pt()->element_pt());

 //Attach the boundary conditions to the mesh
 cout << assign_eqn_numbers() << std::endl; 

 
} //end of constructor


//=====================start_of_actions_before_adapt======================
/// Actions before adapt: Wipe the mesh of traction elements
//========================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::actions_before_adapt()
{
 // Kill the traction elements and wipe surface mesh
 delete_traction_elements();
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();

}// end of actions_before_adapt



//=====================start_of_actions_after_adapt=======================
///  Actions after adapt: Rebuild the mesh of traction elements
//========================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::actions_after_adapt()
{
 // Create traction elements from all elements that are 
 // adjacent to boundary 2 and add them to surface meshes
 create_traction_elements();
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();
 
 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<2>::pin_redundant_nodal_solid_pressures(
  solid_mesh_pt()->element_pt());

 // Set pointer to prescribed traction function for traction elements
 set_traction_pt();
 
}// end of actions_after_adapt



//==================start_of_set_traction_pt==============================
/// Set pointer to traction function for the relevant
/// elements in the traction mesh
//========================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::set_traction_pt()
{
 // Loop over the elements in the traction element mesh
 // for elements on the top boundary (boundary 2)
 unsigned n_element=traction_mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid traction element
   SolidTractionElement<ELEMENT> *el_pt = 
    dynamic_cast<SolidTractionElement<ELEMENT>*>
    (traction_mesh_pt()->element_pt(i));

   //Set the traction function
   el_pt->traction_fct_pt() = Global_Physical_Variables::constant_pressure;
  }
 
}// end of set traction pt


 
//============start_of_create_traction_elements==============================
/// Create traction elements 
//=======================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::create_traction_elements()
{
 // Traction elements are located on boundary 2:
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
   Traction_mesh_pt->add_element_pt(new SolidTractionElement<ELEMENT>
                                    (bulk_elem_pt,face_index));   
  }  

 // Pass the pointer to the traction function to the traction elements
 set_traction_pt();
 
} // end of create_traction_elements




//============start_of_delete_traction_elements==============================
/// Delete traction elements and wipe the  traction meshes
//=======================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::delete_traction_elements()
{
 // How many surface elements are in the surface mesh
 unsigned n_element = Traction_mesh_pt->nelement();
 
 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill surface element
   delete Traction_mesh_pt->element_pt(e);
  }
 
 // Wipe the mesh
 Traction_mesh_pt->flush_element_and_node_storage();

} // end of delete_traction_elements



//==============start_doc===========================================
/// Doc the solution
//==================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::doc_solution()
{

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned n_plot = 5; 

 // Output shape of and stress in deformed body
 //--------------------------------------------
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 solid_mesh_pt()->output(some_file,n_plot);
 some_file.close();


 // Output St. Venant solution
 //---------------------------
 sprintf(filename,"%s/exact_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);

 // Element dimension
 unsigned el_dim=2;
 
 Vector<double> s(el_dim);
 Vector<double> x(el_dim);
 Vector<double> xi(el_dim);
 DenseMatrix<double> sigma(el_dim,el_dim);
 
 // Constants for exact (St. Venant) solution
 double a=-1.0/4.0*Global_Physical_Variables::P;
 double b=-3.0/8.0*Global_Physical_Variables::P/Global_Physical_Variables::H;
 double c=1.0/8.0*Global_Physical_Variables::P/
  pow(Global_Physical_Variables::H,3);
 double d=1.0/20.0*Global_Physical_Variables::P/Global_Physical_Variables::H;

 // Loop over all elements to plot exact solution for stresses
 unsigned nel=solid_mesh_pt()->nelement();
 for (unsigned e=0;e<nel;e++)
  {   
   // Get pointer to element
   ELEMENT* el_pt=dynamic_cast<ELEMENT*>(
    solid_mesh_pt()->element_pt(e));
   
   //Tecplot header info 
   some_file << "ZONE I=" << n_plot << ", J=" << n_plot << std::endl;
   
   //Loop over plot points
   for(unsigned l2=0;l2<n_plot;l2++)
    {
     s[1] = -1.0 + l2*2.0/(n_plot-1);
     for(unsigned l1=0;l1<n_plot;l1++)
      {
       s[0] = -1.0 + l1*2.0/(n_plot-1);
       
       // Get Eulerian and Lagrangian coordinates
       el_pt->interpolated_x(s,x);
       el_pt->interpolated_xi(s,xi);
       
       //Output the x,y,..
       for(unsigned i=0;i<el_dim;i++) 
        {some_file << x[i] << " ";}
       
       // Change orientation of coordinate system relative
       // to solution in lecture notes
       double xx=Global_Physical_Variables::L-xi[0];
       double yy=xi[1];

         // Approximate analytical (St. Venant) solution
         sigma(0,0)=c*(6.0*xx*xx*yy-4.0*yy*yy*yy)+
          6.0*d*yy;
         sigma(1,1)=2.0*(a+b*yy+c*yy*yy*yy);
         sigma(1,0)=2.0*(b*xx+3.0*c*xx*yy*yy);
         sigma(0,1)=sigma(1,0);

         // Output stress
         some_file << sigma(0,0) << " "
                   << sigma(1,0) << " "
                   << sigma(1,1) << " "
                   << std::endl;
        }
      }
  }
 some_file.close();

 // Write trace file: Load/displacement characteristics
 Trace_file << Global_Physical_Variables::P  << " " 
            << Trace_node_pt->x(0) << " " 
            << Trace_node_pt->x(1) << " " 
            << std::endl;

 // Increment label for output files
 Doc_info.number()++;

} //end doc





//==============start_run_it========================================
/// Run it
//==================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::run_it(const unsigned& i_case)
{

#ifdef TIME_SOLID_JAC
  PVDEquationsBase<2>::Solid_timer.reset();
#endif

 // Set output directory
 char dirname[100];   

#ifdef REFINE
 sprintf(dirname,"RESLT_refine%i",i_case);
#else
 sprintf(dirname,"RESLT_norefine%i",i_case);
#endif

 Doc_info.set_directory(dirname);

 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",Doc_info.directory().c_str());
 Trace_file.open(filename);


 // Doc solution
 doc_solution();

 // Initial values for parameter values
 Global_Physical_Variables::P=0.0; 
 Global_Physical_Variables::Gravity=0.0;
 
 //Parameter incrementation
 unsigned nstep=1; 
 double p_increment=1.0e-5;   
 for(unsigned i=0;i<nstep;i++)
  {
   // Increment pressure load
   Global_Physical_Variables::P+=p_increment;

   // Solve the problem with Newton's method, allowing
   // up to max_adapt mesh adaptations after every solve.

#ifdef REFINE

   // Max. number of adaptations per solve
   unsigned max_adapt=1;
   
   newton_solve(max_adapt);

#else
 
   newton_solve();

#endif

   // Doc solution
   doc_solution();
   
  }

#ifdef TIME_SOLID_JAC

 std::cout << "Total fill_in... : " 
           <<   PVDEquationsBase<2>::Solid_timer.cumulative_time(0) 
           << std::endl;
 
 std::cout << "Total d_stress_dG: " 
           <<   PVDEquationsBase<2>::Solid_timer.cumulative_time(1) 
           << std::endl;
 
 std::cout << "Total Jac        : " 
           <<   PVDEquationsBase<2>::Solid_timer.cumulative_time(2) 
           << std::endl;
 
  PVDEquationsBase<2>::Solid_timer.reset();
 
#endif

}


//=======start_of_main==================================================
/// Driver for cantilever beam loaded by surface traction and/or
/// gravity
//======================================================================
int main()
{
 
 bool use_fd=false;

 // Is the material incomressible
 bool incompress=false;

 // Number of cases per implementation
 unsigned ncase=5;


 // Loop over fd and analytical implementation
 for (unsigned i=0;i<2;i++)
  {

   // Generalised Hookean constitutive equations
   //-------------------------------------------
   {
    Global_Physical_Variables::Constitutive_law_pt = 
     new GeneralisedHookean(&Global_Physical_Variables::Nu,
                            &Global_Physical_Variables::E);

    incompress=false;
    
#ifdef REFINE
    {
     //Set up the problem with pure displacement based elements
     CantileverProblem<MySolidElement<RefineableQPVDElement<2,3> > > 
      problem(incompress,use_fd);
     
     problem.run_it(0+i*ncase);
    }
#else   
    {
     //Set up the problem with pure displacement based elements
     CantileverProblem<MySolidElement<QPVDElement<2,3> > > 
      problem(incompress,use_fd);
     
     problem.run_it(0+i*ncase);
    }
#endif
    

#ifdef REFINE
    {
     //Set up the problem with continous pressure/displacement
     CantileverProblem<MySolidElement<
      RefineableQPVDElementWithContinuousPressure<2> > > 
      problem(incompress,use_fd); 
     
     problem.run_it(1+i*ncase);
    }
#else
    {
     //Set up the problem with continous pressure/displacement
     CantileverProblem<MySolidElement<
      QPVDElementWithContinuousPressure<2> > > 
      problem(incompress,use_fd); 
     
     problem.run_it(1+i*ncase);
    }
#endif


#ifdef REFINE    
    {
     //Set up the problem with discontinous pressure/displacement
     CantileverProblem<MySolidElement<
      RefineableQPVDElementWithPressure<2> > > problem(incompress,use_fd); 
     
     problem.run_it(2+i*ncase);
    }
#else
    {
     //Set up the problem with discontinous pressure/displacement
     CantileverProblem<MySolidElement<
      QPVDElementWithPressure<2> > > problem(incompress,use_fd); 
     
     problem.run_it(2+i*ncase);
    }
#endif
    
    delete Global_Physical_Variables::Constitutive_law_pt;
    Global_Physical_Variables::Constitutive_law_pt=0;
   }
   

   
   // Incompressible Mooney Rivlin
   //-----------------------------
   {
    Global_Physical_Variables::Strain_energy_function_pt = 
     new MooneyRivlin(&Global_Physical_Variables::C1,
                      &Global_Physical_Variables::C2);
    
    // Define a constitutive law (based on strain energy function)
    Global_Physical_Variables::Constitutive_law_pt = 
     new IsotropicStrainEnergyFunctionConstitutiveLaw(
      Global_Physical_Variables::Strain_energy_function_pt);
    
    incompress=true;


#ifdef REFINE
    {
     //Set up the problem with continous pressure/displacement
     CantileverProblem<MySolidElement<
      RefineableQPVDElementWithContinuousPressure<2> > > 
      problem(incompress,use_fd); 
     
     problem.run_it(3+i*ncase);
    }
#else
    {
     //Set up the problem with continous pressure/displacement
     CantileverProblem<MySolidElement<
      QPVDElementWithContinuousPressure<2> > > 
      problem(incompress,use_fd); 
     
     problem.run_it(3+i*ncase);
    }
#endif


    
#ifdef REFINE
    {
     //Set up the problem with discontinous pressure/displacement
     CantileverProblem<MySolidElement<
      RefineableQPVDElementWithPressure<2> > > 
      problem(incompress,use_fd); 
     
     problem.run_it(4+i*ncase);
    }
#else
    {
     //Set up the problem with discontinous pressure/displacement
     CantileverProblem<MySolidElement<
      QPVDElementWithPressure<2> > > 
      problem(incompress,use_fd); 
     
     problem.run_it(4+i*ncase);
    }
#endif
    
    delete  Global_Physical_Variables::Strain_energy_function_pt;
    Global_Physical_Variables::Strain_energy_function_pt=0;
    
    delete Global_Physical_Variables::Constitutive_law_pt;
    Global_Physical_Variables::Constitutive_law_pt=0;
   }

 
   use_fd=true;
   std::cout << "\n\n\n CR Total fill_in... : bla \n\n\n " << std::endl;
   
  }

} //end of main








