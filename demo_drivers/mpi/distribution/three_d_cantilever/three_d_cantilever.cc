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
//Driver for 3D Airy cantilever beam problem [Note: actually this does not
//use any traction elements, so perhaps it should be renamed?  As with the
//2D example, gravity is the parameter that is incremented.]

//Oomph-lib includes
#include "generic.h"
#include "solid.h"
#include "constitutive.h"

// The mesh
#include "meshes/simple_cubic_mesh.template.h"

// The mesh
#include "meshes/quarter_tube_mesh.h"

using namespace std;

using namespace oomph;

//#define REFINE




/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////

//=========================================================================
/// Simple quarter tube mesh upgraded to become a solid mesh
//=========================================================================
template<class ELEMENT>
class RefineableElasticQuarterTubeMesh :
 public virtual RefineableQuarterTubeMesh<ELEMENT>,
 public virtual SolidMesh
{

public:

 /// Constructor:
 RefineableElasticQuarterTubeMesh(GeomObject* wall_pt,
                                  const Vector<double>& xi_lo,
                                  const double& fract_mid,
                                  const Vector<double>& xi_hi,
                                  const unsigned& nlayer,
                                  TimeStepper* time_stepper_pt=
                                  &Mesh::Default_TimeStepper) :
  QuarterTubeMesh<ELEMENT>(wall_pt,xi_lo,fract_mid,xi_hi,
                           nlayer,time_stepper_pt),
  RefineableQuarterTubeMesh<ELEMENT>(wall_pt,xi_lo,fract_mid,xi_hi,
                                     nlayer,time_stepper_pt),
  SolidMesh()
  {
   //Assign the initial lagrangian coordinates
   set_lagrangian_nodal_coordinates();
  }

 /// Empty Destructor
 virtual ~RefineableElasticQuarterTubeMesh() { }

};

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////

//=========================================================================
/// Simple quarter tube mesh upgraded to become a solid mesh
//=========================================================================
template<class ELEMENT>
class ElasticQuarterTubeMesh : public virtual QuarterTubeMesh<ELEMENT>,
                               public virtual SolidMesh
{

public:

 /// Constructor:
 ElasticQuarterTubeMesh(GeomObject* wall_pt,
                 const Vector<double>& xi_lo,
                 const double& fract_mid,
                 const Vector<double>& xi_hi,
                 const unsigned& nlayer,
                 TimeStepper* time_stepper_pt=
                 &Mesh::Default_TimeStepper) :
  QuarterTubeMesh<ELEMENT>(wall_pt,xi_lo,fract_mid,xi_hi,
                           nlayer,time_stepper_pt),
  SolidMesh()
  {
   //Assign the initial lagrangian coordinates
   set_lagrangian_nodal_coordinates();
  }

 /// Empty Destructor
 virtual ~ElasticQuarterTubeMesh() { }

};

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////


//=========================================================================
/// Simple cubic mesh upgraded to become a solid mesh
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
  SimpleCubicMesh<ELEMENT>(nx,ny,nz,-a,a,-b,b,-c,c,time_stepper_pt),
  RefineableBrickMesh<ELEMENT>(), SolidMesh()
  {

   this->setup_octree_forest();

   //Assign the initial lagrangian coordinates
   set_lagrangian_nodal_coordinates();
  }

 /// Empty Destructor
 virtual ~RefineableElasticCubicMesh() { }

};

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////

//=========================================================================
/// Simple cubic mesh upgraded to become a solid mesh
//=========================================================================
template<class ELEMENT>
class ElasticCubicMesh : public virtual SimpleCubicMesh<ELEMENT>,
                         public virtual SolidMesh
{

public:

 /// Constructor:
 ElasticCubicMesh(const unsigned &nx, const unsigned &ny, const unsigned &nz,
                  const double &a, const double &b, const double &c,
                  TimeStepper* time_stepper_pt=&Mesh::Default_TimeStepper) :
  SimpleCubicMesh<ELEMENT>(nx,ny,nz,-a,a,-b,b,-c,c,time_stepper_pt)
  {
   //Assign the initial lagrangian coordinates
   set_lagrangian_nodal_coordinates();
  }

 /// Empty Destructor
 virtual ~ElasticCubicMesh() { }

};

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



// namespace oomph
// {

// //=================start_wrapper==================================
// /// Wrapper class for solid elements to modify their output
// /// functions.
// //================================================================
// template <class ELEMENT>
// class MySolidElement : public virtual ELEMENT
// {

// public:

//  /// Constructor: Call constructor of underlying element
//  MySolidElement() : ELEMENT() {};

//  /// Overload output function:
//  void output(std::ostream &outfile, const unsigned &n_plot)
//   {

//    // Element dimension
//    unsigned el_dim = this->dim();

//    Vector<double> s(el_dim);
//    Vector<double> x(el_dim);
//    DenseMatrix<double> sigma(el_dim,el_dim);

//    switch(el_dim)
//     {

//     case 2:

//      //Tecplot header info
//      outfile << "ZONE I=" << n_plot << ", J=" << n_plot << std::endl;

//      //Loop over element nodes
//      for(unsigned l2=0;l2<n_plot;l2++)
//       {
//        s[1] = -1.0 + l2*2.0/(n_plot-1);
//        for(unsigned l1=0;l1<n_plot;l1++)
//         {
//          s[0] = -1.0 + l1*2.0/(n_plot-1);

//          // Get Eulerian coordinates and stress
//          this->interpolated_x(s,x);
//          this->get_stress(s,sigma);

//          //Output the x,y,..
//          for(unsigned i=0;i<el_dim;i++)
//           {outfile << x[i] << " ";}

//          // Output stress
//          outfile << sigma(0,0) << " "
//                  << sigma(1,0) << " "
//                  << sigma(1,1) << " "
//                  << std::endl;
//         }
//       }

//      break;

//     default:

//      std::ostringstream error_message;
//      error_message << "Output for dim !=2 not implemented" << std::endl;
//      throw OomphLibError(error_message.str(),
//                          OOMPH_CURRENT_FUNCTION,
//                          OOMPH_EXCEPTION_LOCATION);
//     }

//   }

// };



// //===========start_face_geometry==============================================
// /// FaceGeometry of wrapped element is the same as the underlying element
// //============================================================================
// template<class ELEMENT>
// class FaceGeometry<MySolidElement<ELEMENT> > :
//  public virtual FaceGeometry<ELEMENT>
// {
// };


// }






/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



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

//  // # of elements in x-direction
//  unsigned n_x=20;

//  // # of elements in y-direction
//  unsigned n_y=2;

//  // # of elements in z-direction
//  unsigned n_z=2;

 // Domain length in x-direction
 double l_x= Global_Physical_Variables::L;

//  // Domain length in y-direction
//  double l_y=2.0*Global_Physical_Variables::H;

//  // Domain length in z-direction
//  double l_z=2.0*Global_Physical_Variables::H;


 // Create geometric objects: Elliptical tube with half axes = radius = 1.0
 double radius=1.0;
 GeomObject* wall_pt=new EllipticalTube(radius,radius);

 // Boundaries on object
 Vector<double> xi_lo(2);

 // height of inflow
 xi_lo[0]=0.0;

 // start of Wall_pt
 xi_lo[1]=0.0;

 Vector<double> xi_hi(2);

 // height of outflow
 xi_hi[0]=l_x;

 // end of Wall_pt
 xi_hi[1]=2.0*atan(1.0);

 // # of layers
 unsigned nlayer=6;

 //Radial divider is located half-way along the circumference
 double frac_mid=0.5;


#ifdef REFINE

 //Now create the mesh
 Problem::mesh_pt() = new RefineableElasticQuarterTubeMesh<ELEMENT>
  (wall_pt,xi_lo,frac_mid,xi_hi,nlayer);

 // Set error estimator
 dynamic_cast<RefineableElasticQuarterTubeMesh<ELEMENT>* >(
  mesh_pt())->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 // Error targets for adaptive refinement
 dynamic_cast<RefineableElasticQuarterTubeMesh<ELEMENT>* >(
  mesh_pt())->max_permitted_error()=0.05;
 dynamic_cast<RefineableElasticQuarterTubeMesh<ELEMENT>* >(
  mesh_pt())->min_permitted_error()=0.005;

 // 0.115319 0.00212537

#else

 //Now create the mesh
 Problem::mesh_pt() = new ElasticQuarterTubeMesh<ELEMENT>
  (wall_pt,xi_lo,frac_mid,xi_hi,nlayer);

#endif


 //Assign the physical properties to the elements before any refinement
 //Loop over the elements in the main mesh
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   //Cast to a solid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

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
     PVDEquationsWithPressure<3>* test_pt =
      dynamic_cast<PVDEquationsWithPressure<3>*>(
       mesh_pt()->element_pt(i));
     if (test_pt!=0)
      {
       test_pt->set_incompressible();
      }
    }

  }

 // Choose a control node: The last node in the solid mesh
 unsigned nnod=mesh_pt()->nnode();
 Trace_node_pt=mesh_pt()->node_pt(nnod-1);

 // Pin the left boundary (boundary 0) in all directions
 unsigned b=0;
 unsigned n_side = mesh_pt()->nboundary_node(b);

#ifdef REFINE

 RefineableElasticQuarterTubeMesh<ELEMENT>* solid_mesh_pt=
  dynamic_cast<RefineableElasticQuarterTubeMesh<ELEMENT>* >(mesh_pt());

#else

 ElasticQuarterTubeMesh<ELEMENT>* solid_mesh_pt=
  dynamic_cast<ElasticQuarterTubeMesh<ELEMENT>* >(mesh_pt());

#endif

 oomph_info << "solid_mesh_pt=" << solid_mesh_pt << std::endl;

 //Loop over the nodes
 for(unsigned i=0;i<n_side;i++)
  {
   solid_mesh_pt->boundary_node_pt(0,i)->pin_position(0);
   solid_mesh_pt->boundary_node_pt(0,i)->pin_position(1);
   solid_mesh_pt->boundary_node_pt(0,i)->pin_position(2);
  }


 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(
  mesh_pt()->element_pt());

 //Attach the boundary conditions to the mesh
 cout << assign_eqn_numbers() << std::endl;

 //Reset max residuals
 oomph_info << "Reset max residuals" << std::endl;
 Problem::Max_residuals=1.0e10;
 Problem::Max_newton_iterations=50;

} //end of constructor


//=====================start_of_actions_before_adapt======================
/// Actions before adapt: empty
//========================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::actions_before_adapt()
{

}// end of actions_before_adapt



//=====================start_of_actions_after_adapt=======================
///  Actions after adapt: pin redundant pressures
//========================================================================
template<class ELEMENT>
void CantileverProblem<ELEMENT>::actions_after_adapt()
{
 // Pin the redundant solid pressures (if any)
 PVDEquationsBase<3>::pin_redundant_nodal_solid_pressures(
  mesh_pt()->element_pt());

}// end of actions_after_adapt




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
 sprintf(filename,"%s/soln%i_on_proc%i.dat",Doc_info.directory().c_str(),
         Doc_info.number(),this->communicator_pt()->my_rank());
 some_file.open(filename);
 mesh_pt()->output(some_file,n_plot);
 some_file.close();

 // Select a trace node to be the last available node on each process
 // (to get round the problem of the original trace node possibly not
 //  being available to certain processors)
 unsigned nnod=mesh_pt()->nnode();
 Trace_node_pt=mesh_pt()->node_pt(nnod-1);

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
 PVDEquationsBase<3>::Solid_timer.reset();
#endif

 mesh_pt()->output_boundaries("boundaries.dat");

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
 sprintf(filename,"%s/trace_on_proc%i.dat",Doc_info.directory().c_str(),
         this->communicator_pt()->my_rank());
 Trace_file.open(filename);


 // Doc solution
 doc_solution();

 // Initial values for parameter values
 Global_Physical_Variables::P=0.0;
 Global_Physical_Variables::Gravity=0.0;

 //Parameter incrementation
 unsigned nstep=1;
 double g_increment=1.0e-5;

 for(unsigned i=0;i<nstep;i++)
  {
   // Increment gravity
   Global_Physical_Variables::Gravity+=g_increment;

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

 PVDEquationsBase<3>::doc_timings();

 PVDEquationsBase<3>::Solid_timer.reset();

#endif

}


//=======start_of_main==================================================
/// Driver for cantilever beam loaded by surface traction and/or
/// gravity
//======================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif

 DocInfo mesh_doc_info;
 bool report_stats=true;
 mesh_doc_info.number()=10;
 mesh_doc_info.set_directory("RESLT_MESH");

 // Use finite differencing
 bool use_fd=false;

 // Number of cases per implementation
 unsigned ncase=5;

 // Is the material incomressible
 bool incompress=true;

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
     CantileverProblem<RefineableQPVDElement<3,3> >
      problem(incompress,use_fd);

#ifdef OOMPH_HAS_MPI
     //Distribute it
     std::ifstream input_file;
     std::ofstream output_file;
     char filename[100];

     // Get partition from file
     unsigned n_partition=problem.mesh_pt()->nelement();
     Vector<unsigned> element_partition(n_partition);
     sprintf(filename,"three_d_cantilever_%i_partition.dat",0+i*ncase);
     input_file.open(filename);
     std::string input_string;
     for (unsigned e=0;e<n_partition;e++)
      {
       getline(input_file,input_string,'\n');
       element_partition[e]=atoi(input_string.c_str());
      }

     problem.distribute(element_partition,mesh_doc_info,report_stats);

     problem.check_halo_schemes(mesh_doc_info);
     mesh_doc_info.number()++;
#endif

     problem.run_it(0+i*ncase);
    }
#else
    {
     //Set up the problem with pure displacement based elements
     CantileverProblem<QPVDElement<3,3> >
      problem(incompress,use_fd);

#ifdef OOMPH_HAS_MPI
     //Distribute it
     std::ifstream input_file;
     std::ofstream output_file;
     char filename[100];

     // Get partition from file
     unsigned n_partition=problem.mesh_pt()->nelement();
     Vector<unsigned> element_partition(n_partition);
     sprintf(filename,"three_d_cantilever_%i_partition.dat",0+i*ncase);
     input_file.open(filename);
     std::string input_string;
     for (unsigned e=0;e<n_partition;e++)
      {
       getline(input_file,input_string,'\n');
       element_partition[e]=atoi(input_string.c_str());
      }

     problem.distribute(element_partition,mesh_doc_info,report_stats);

     problem.check_halo_schemes(mesh_doc_info);
     mesh_doc_info.number()++;
#endif

     problem.run_it(0+i*ncase);
    }
#endif


#ifdef REFINE
    {
     //Set up the problem with continous pressure/displacement
     CantileverProblem<RefineableQPVDElementWithContinuousPressure<3> >
      problem(incompress,use_fd);

#ifdef OOMPH_HAS_MPI
     //Distribute it
     std::ifstream input_file;
     std::ofstream output_file;
     char filename[100];

     // Get partition from file
     unsigned n_partition=problem.mesh_pt()->nelement();
     Vector<unsigned> element_partition(n_partition);
     sprintf(filename,"three_d_cantilever_%i_partition.dat",1+i*ncase);
     input_file.open(filename);
     std::string input_string;
     for (unsigned e=0;e<n_partition;e++)
      {
       getline(input_file,input_string,'\n');
       element_partition[e]=atoi(input_string.c_str());
      }

     problem.distribute(element_partition,mesh_doc_info,report_stats);

     problem.check_halo_schemes(mesh_doc_info);
     mesh_doc_info.number()++;
#endif

     problem.run_it(1+i*ncase);
    }
#else
    {
     //Set up the problem with continous pressure/displacement
     CantileverProblem<QPVDElementWithContinuousPressure<3> >
      problem(incompress,use_fd);

#ifdef OOMPH_HAS_MPI
     //Distribute it
     std::ifstream input_file;
     std::ofstream output_file;
     char filename[100];

     // Get partition from file
     unsigned n_partition=problem.mesh_pt()->nelement();
     Vector<unsigned> element_partition(n_partition);
     sprintf(filename,"three_d_cantilever_%i_partition.dat",0+i*ncase);
     input_file.open(filename);
     std::string input_string;
     for (unsigned e=0;e<n_partition;e++)
      {
       getline(input_file,input_string,'\n');
       element_partition[e]=atoi(input_string.c_str());
      }

     problem.distribute(element_partition,mesh_doc_info,report_stats);

     problem.check_halo_schemes(mesh_doc_info);
     mesh_doc_info.number()++;
#endif

     problem.run_it(1+i*ncase);
    }
#endif



#ifdef REFINE
    {
     //Set up the problem with discontinous pressure/displacement
     CantileverProblem<RefineableQPVDElementWithPressure<3> >
      problem(incompress,use_fd);

#ifdef OOMPH_HAS_MPI
     //Distribute it
     std::ifstream input_file;
     std::ofstream output_file;
     char filename[100];

     // Get partition from file
     unsigned n_partition=problem.mesh_pt()->nelement();
     Vector<unsigned> element_partition(n_partition);
     sprintf(filename,"three_d_cantilever_%i_partition.dat",1+i*ncase);
     input_file.open(filename);
     std::string input_string;
     for (unsigned e=0;e<n_partition;e++)
      {
       getline(input_file,input_string,'\n');
       element_partition[e]=atoi(input_string.c_str());
      }

     problem.distribute(element_partition,mesh_doc_info,report_stats);

     problem.check_halo_schemes(mesh_doc_info);
     mesh_doc_info.number()++;
#endif

     problem.run_it(2+i*ncase);
    }
#else
    {
     //Set up the problem with discontinous pressure/displacement
     CantileverProblem<QPVDElementWithPressure<3> >
      problem(incompress,use_fd);

#ifdef OOMPH_HAS_MPI
     //Distribute it
     std::ifstream input_file;
     std::ofstream output_file;
     char filename[100];

     // Get partition from file
     unsigned n_partition=problem.mesh_pt()->nelement();
     Vector<unsigned> element_partition(n_partition);
     sprintf(filename,"three_d_cantilever_%i_partition.dat",0+i*ncase);
     input_file.open(filename);
     std::string input_string;
     for (unsigned e=0;e<n_partition;e++)
      {
       getline(input_file,input_string,'\n');
       element_partition[e]=atoi(input_string.c_str());
      }

     problem.distribute(element_partition,mesh_doc_info,report_stats);

     problem.check_halo_schemes(mesh_doc_info);
     mesh_doc_info.number()++;
#endif

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
     CantileverProblem<RefineableQPVDElementWithContinuousPressure<3> >
      problem(incompress,use_fd);

#ifdef OOMPH_HAS_MPI
     //Distribute it
     std::ifstream input_file;
     std::ofstream output_file;
     char filename[100];

     // Get partition from file
     unsigned n_partition=problem.mesh_pt()->nelement();
     Vector<unsigned> element_partition(n_partition);
     sprintf(filename,"three_d_cantilever_%i_partition.dat",3+i*ncase);
     input_file.open(filename);
     std::string input_string;
     for (unsigned e=0;e<n_partition;e++)
      {
       getline(input_file,input_string,'\n');
       element_partition[e]=atoi(input_string.c_str());
      }

     problem.distribute(element_partition,mesh_doc_info,report_stats);

     problem.check_halo_schemes(mesh_doc_info);
     mesh_doc_info.number()++;
#endif

     problem.run_it(3+i*ncase);
    }
#else
    {
     //Set up the problem with continous pressure/displacement
     CantileverProblem<QPVDElementWithContinuousPressure<3> >
      problem(incompress,use_fd);

#ifdef OOMPH_HAS_MPI
     //Distribute it
     std::ifstream input_file;
     std::ofstream output_file;
     char filename[100];

     // Get partition from file
     unsigned n_partition=problem.mesh_pt()->nelement();
     Vector<unsigned> element_partition(n_partition);
     sprintf(filename,"three_d_cantilever_%i_partition.dat",3+i*ncase);
     input_file.open(filename);
     std::string input_string;
     for (unsigned e=0;e<n_partition;e++)
      {
       getline(input_file,input_string,'\n');
       element_partition[e]=atoi(input_string.c_str());
      }

     problem.distribute(element_partition,mesh_doc_info,report_stats);

     problem.check_halo_schemes(mesh_doc_info);
     mesh_doc_info.number()++;
#endif

     problem.run_it(3+i*ncase);
    }
#endif



#ifdef REFINE
    {
     //Set up the problem with discontinous pressure/displacement
     CantileverProblem<RefineableQPVDElementWithPressure<3> >
      problem(incompress,use_fd);

#ifdef OOMPH_HAS_MPI
     //Distribute it
     std::ifstream input_file;
     std::ofstream output_file;
     char filename[100];

     // Get partition from file
     unsigned n_partition=problem.mesh_pt()->nelement();
     Vector<unsigned> element_partition(n_partition);
     sprintf(filename,"three_d_cantilever_%i_partition.dat",4+i*ncase);
     input_file.open(filename);
     std::string input_string;
     for (unsigned e=0;e<n_partition;e++)
      {
       getline(input_file,input_string,'\n');
       element_partition[e]=atoi(input_string.c_str());
      }

     problem.distribute(element_partition,mesh_doc_info,report_stats);

     problem.check_halo_schemes(mesh_doc_info);
     mesh_doc_info.number()++;
#endif

     problem.run_it(4+i*ncase);
    }
#else
    {
     //Set up the problem with discontinous pressure/displacement
     CantileverProblem<QPVDElementWithPressure<3> >
      problem(incompress,use_fd);

#ifdef OOMPH_HAS_MPI
     //Distribute it
     std::ifstream input_file;
     std::ofstream output_file;
     char filename[100];

     // Get partition from file
     unsigned n_partition=problem.mesh_pt()->nelement();
     Vector<unsigned> element_partition(n_partition);
     sprintf(filename,"three_d_cantilever_%i_partition.dat",4+i*ncase);
     input_file.open(filename);
     std::string input_string;
     for (unsigned e=0;e<n_partition;e++)
      {
       getline(input_file,input_string,'\n');
       element_partition[e]=atoi(input_string.c_str());
      }

     problem.distribute(element_partition,mesh_doc_info,report_stats);

     problem.check_halo_schemes(mesh_doc_info);
     mesh_doc_info.number()++;
#endif

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

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif

} //end of main
