//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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
// Driver for Axisymmetic spin up in a cylinder using black
// box adaptation, using Axisymmetric Taylor Hood and Crouzeix 
// Raviart elements.

// Generic oomph-lib header
#include "generic.h"

// Navier Stokes headers
#include "navier_stokes.h"

// Axisym Navier Stokes headers
#include "axisym_navier_stokes.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

using namespace MathematicalConstants;



//=start_of_namespace====================================================
/// Namespace for physical parameters
//=======================================================================
namespace ValidateRateOfStrain
{

 double A_r=1.0;
 double B_r=2.0;
 double C_r=3.0;


 double A_z=1.0;
 double B_z=2.0;
 double C_z=3.0;

 double A_phi=1.0;
 double B_phi=2.0;
 double C_phi=3.0;

 void get_veloc(const Vector<double>& x, Vector<double>& veloc)
 {
  veloc[0]=(A_r+B_r*x[0]+C_r*x[0])*(B_r*x[1]+C_r*x[1]);
  veloc[1]=(A_z+B_z*x[0]+C_z*x[0])*(B_z*x[1]+C_z*x[1]);
  veloc[2]=(A_phi+B_phi*x[0]+C_phi*x[0])*(B_phi*x[1]+C_phi*x[1]);  
 }

 double fd_strain(const unsigned& i, const unsigned& j, 
                  const Vector<double>& x)
 {
  // copy to keep x constant
  Vector<double> x_local(x);
  
  // Get base velocity component
  Vector<double> veloc(3);
  get_veloc(x_local,veloc);
  double base=veloc[i];

  // Take fd step
  double epsilon=1.0e-8;
  x_local[j]+=epsilon;
  get_veloc(x_local,veloc);
  double adv=veloc[i];

  // Return FD approx.
  return (adv-base)/epsilon;
 }

 // Get exact strain
 void exact_strain(const Vector<double>& x, DenseMatrix<double>& strain_rate)
 {
  Vector<double> veloc(3);
  get_veloc(x,veloc);
  double ur=veloc[0];
  //double uz=veloc[1];
  double uphi=veloc[2];

  double durdr  =fd_strain(0,0,x);
  double durdz  =fd_strain(0,1,x);
  //double durdphi=fd_strain(0,2,x);

  double duzdr  =fd_strain(1,0,x);
  double duzdz  =fd_strain(1,1,x);
  //double duzdphi=fd_strain(1,2,x); 

  double duphidr  =fd_strain(2,0,x);
  double duphidz  =fd_strain(2,1,x);
  //double duphidphi=fd_strain(2,2,x);

  // Assign strain rates without negative powers of the radius
  // and zero those with:
  strain_rate(0,0)=durdr;
  strain_rate(0,1)=0.5*(durdz+duzdr);
  strain_rate(1,0)=strain_rate(0,1);
  strain_rate(0,2)=0.0;
  strain_rate(2,0)=strain_rate(0,2);
  strain_rate(1,1)=duzdz;
  strain_rate(1,2)=0.5*duphidz;
  strain_rate(2,1)=strain_rate(1,2);
  strain_rate(2,2)=0.0;

  
  // Overwrite the strain rates with negative powers of the radius
  // unless we're at the origin
  if (abs(x[0])>1.0e-16)
   {
    double inverse_radius=1.0/x[0];
    strain_rate(0,2)=0.5*(duphidr-inverse_radius*uphi);
    strain_rate(2,0)=strain_rate(0,2);
    strain_rate(2,2)=inverse_radius*ur;
   }

  // Deliberately introduce error to check things
  double error_factor=1.0;
  for (unsigned i=0;i<3;i++)
   {
    for (unsigned j=0;j<3;j++)
     {
      strain_rate(i,j)*=error_factor;
     }
   }

 }

} // end of namespace



//=start_of_namespace====================================================
/// Namespace for physical parameters
//=======================================================================
namespace Global_Physical_Variables
{
 /// Reynolds number
 double Re=5.0;

 /// Womersley number
 double ReSt=5.0;

} // end of namespace



//==start_of_problem_class============================================
/// \short Refineable Rotating cylinder problem in rectangular 
/// axisymmetric domain.
//====================================================================
template<class ELEMENT, class TIMESTEPPER>
class RefineableRotatingCylinderProblem : public Problem
{

public:

 /// Constructor
 RefineableRotatingCylinderProblem(
  const unsigned &nx, const unsigned &ny,
  const double &lx, const double &ly);

 /// Destructor: Empty
 ~RefineableRotatingCylinderProblem(){}

 /// Update the after solve (empty)
 void actions_after_newton_solve() {}

 /// \short Update the problem specs before solve. 
 /// (Re-)set velocity boundary conditions. 
 void actions_before_newton_solve()
  { 
   // Overwrite with no flow along the solid boundaries (0,1,2)
   unsigned num_bound = mesh_pt()->nboundary();
   for(unsigned ibound=0;ibound<num_bound;ibound++)
    {
     unsigned num_nod= mesh_pt()->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       for (unsigned i=0;i<3;i++)
	{
         if (ibound<3)// For the solid walls only boundaries 0,1,2
          {
           //Get the radial values
           double r = mesh_pt()->boundary_node_pt(ibound,inod)->x(0);
           switch (i)
            {
            case 2:// azimuthal velocity
             mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,i,r);
             break;
            case 1:// axial velocity
             mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,i,0.0);
             break;
            case 0:// radial velocity
             mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,i,0.0);
             break;
            }
          }
         // Overwrite with no radial or azimuthal flow on symmetry boundary (3)
         if (ibound==3)
          {
           if (i!=1)
            {
             mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,i,0.0);
            }
          }
        }
      }
    }
  } // end of actions before solve

 /// After adaptation: Pin pressure again (the previously pinned
 /// value might have disappeared) and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableAxisymmetricNavierStokesEquations::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());
   
   // Pin redudant pressure dofs
   RefineableAxisymmetricNavierStokesEquations::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   // Now set the pressure in first element at 'node' 0 to 0.0
   fix_pressure(0,0,0.0);
  }
 
 /// \short Access function for the mesh
 RefineableRectangularQuadMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(
    Problem::mesh_pt());
  }

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Do unsteady run up to max. time with given timestep
 void unsteady_run(const double& t_max, const double& dt, string dir_name, 
                   const unsigned &maximum_ref_level,
                   const unsigned &minimum_ref_level);

 /// Validate rate of strain
 void validate_strain();

private:

 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Fix the pressure at that element
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  }

 /// Trace file
 ofstream Trace_file;

}; // end of problem class



//==start_of_constructor==================================================
/// Constructor for RefineableDrivenCavity problem 
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
RefineableRotatingCylinderProblem<ELEMENT,TIMESTEPPER>::
RefineableRotatingCylinderProblem(
 const unsigned &nx, const unsigned &ny, 
 const double &lx, const double &ly)
{ 
 // Allocate the timestepper -- This constructs the time object as well
 add_time_stepper_pt(new TIMESTEPPER());

 // Setup mesh

 // Build and assign mesh
 Problem::mesh_pt() = 
  new RefineableRectangularQuadMesh<ELEMENT>(nx,ny,lx,ly,time_stepper_pt());

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;
 
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here all nodes on all solid boundaries (0,1,2).
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Pin values for radial velocity on all boundaries.
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(0); 
     // Pin values for axial velocity on all solid boundaries
     if (ibound!=3)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(1);
      }
     // Pin values for azimuthal velocity on all boundaries
     mesh_pt()->boundary_node_pt(ibound,inod)->pin(2);
    }
  } // end loop over boundaries

 //Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to Reynolds
 // number
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
   //Assign time 
   el_pt->time_pt() = time_pt();
   //The mesh remains fixed
   el_pt->disable_ALE();
   //Set the Reynolds number, etc
   el_pt->re_pt() = &Global_Physical_Variables::Re;
   el_pt->re_st_pt() = &Global_Physical_Variables::ReSt;
  }
 
 // Pin redudant pressure dofs
 RefineableAxisymmetricNavierStokesEquations::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 
 // Now set the pressure in first element at 'node' 0 to 0.0
 fix_pressure(0,0,0.0);
 
 //Attach the boundary conditions to the mesh
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} // end of constructor

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void RefineableRotatingCylinderProblem<ELEMENT,TIMESTEPPER>::
doc_solution(DocInfo& doc_info)
{ 
 // Doc
 Trace_file << time_pt()->time() << " " 
            << mesh_pt()->max_permitted_error() << " "
            << mesh_pt()->min_permitted_error() << " "
            << mesh_pt()->max_error() << " "
            << mesh_pt()->min_error() << " " << std::endl;


 // Doc elemental errors
 if (true)
 {
  Mesh* tmp_mesh_pt=mesh_pt();
  unsigned nel=mesh_pt()->nelement();
  Vector<double> elemental_error(nel);
  mesh_pt()->spatial_error_estimator_pt()->get_element_errors(tmp_mesh_pt,
                                                              elemental_error, 
                                                              doc_info);
 }

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);

 // Write file as a tecplot text object
 some_file << "TEXT X=2.5,Y=93.6,F=HELV,HU=POINT,C=BLUE,H=26,T=\"time = " 
           << time_pt()->time() << "\"";
 // ...and draw a horizontal line whose length is proportional
 // to the elapsed time
 some_file << "GEOMETRY X=2.5,Y=98,T=LINE,C=BLUE,LT=0.4" << std::endl;
 some_file << "1" << std::endl;
 some_file << "2" << std::endl;
 some_file << " 0 0" << std::endl;
 some_file << time_pt()->time()*20.0 << " 0" << std::endl;

 // Write dummy zones
 some_file << "ZONE I=2,J=2" << std::endl;
 some_file << "-0.05 -0.05 -0.007 -0.007 -0.05 -0.05" << std::endl;
 some_file << "1.05 -0.05 -0.007 -0.007 -0.05 -0.05" << std::endl;
 some_file << "-0.05 1.05 -0.007 -0.007 -0.05 -0.05" << std::endl;
 some_file << "1.05 1.05 -0.007 -0.007 -0.05 -0.05" << std::endl;
 some_file << "ZONE I=2,J=2" << std::endl;
 some_file << "-0.05 -0.05 0.007 0.007 1.05 2.55" << std::endl;
 some_file << "1.05 -0.05 0.007 0.007 1.05 2.55" << std::endl;
 some_file << "-0.05 1.05 0.007 0.007 1.05 2.55" << std::endl;
 some_file << "1.05 1.05 0.007 0.007 1.05 2.55" << std::endl;

 some_file.close();


 // Doc strain 
 sprintf(filename,"%s/strain_rate%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);

 some_file 
  << "VARIABLES=\"r\",\"z\",\"e<SUB>rr</SUB>\",\"e<SUB>rz</SUB>\","
  << "\"e<SUB>r<GREEK>q</GREEK></SUB>\","
  << "\"e<SUB>zz</SUB>\",\"e<SUB>z<GREEK>q</GREEK></SUB>\","
  << "\"e<SUB><GREEK>qq</GREEK></SUB>\"" << std::endl;

 //Vector of local and global coordinates
 Vector<double> s(2);
 Vector<double> x(3);

 // Rate of strain
 DenseMatrix<double> strain_rate(3,3);

 // Loop over elements to plot
 unsigned nelem=mesh_pt()->nelement();
 for (unsigned e=0;e<nelem;e++)
  {
   // Get pointer to element
   ELEMENT* el_pt=dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   // Tecplot header info
   some_file << el_pt->tecplot_zone_string(npts);
 
   // Loop over plot points
   unsigned num_plot_points=el_pt->nplot_points(npts);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     
     // Get local coordinates of plot point
     el_pt->get_s_plot(iplot,npts,s);
     
     // Get global coordinates of plot point
     el_pt->interpolated_x(s,x);
     
     // Get strain
     el_pt->strain_rate(s,strain_rate);

     some_file << x[0] << " " << x[1] << " " 
               << strain_rate(0,0) << " " 
               << strain_rate(0,1) << " " 
               << strain_rate(0,2) << " " 
               << strain_rate(1,1) << " " 
               << strain_rate(1,2) << " " 
               << strain_rate(2,2) << " "
               << std::endl;
    }
  }
 some_file.close();
 
} // end of doc_solution



//==start_of_unsteady_run=================================================
/// Perform run up to specified time
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void RefineableRotatingCylinderProblem<ELEMENT,TIMESTEPPER>::validate_strain()
{

 /// Assign validation velocity field
 Vector<double> x(2);
 Vector<double> veloc(3);
 unsigned nnod=mesh_pt()->nnode();
 for (unsigned j=0;j<nnod;j++)
  {
   x[0]=mesh_pt()->node_pt(j)->x(0);
   x[1]=mesh_pt()->node_pt(j)->x(1);
   ValidateRateOfStrain::get_veloc(x,veloc);
   mesh_pt()->node_pt(j)->set_value(0,veloc[0]);
   mesh_pt()->node_pt(j)->set_value(1,veloc[1]);
   mesh_pt()->node_pt(j)->set_value(2,veloc[2]);
  }

 // Number of plot points
 unsigned nplot=5; 

 // Output solution 
 ofstream outfile1("fe_strain.dat");
 ofstream outfile2("ex_strain.dat");

 //Vector of local coordinates
 Vector<double> s(2);
 
 // Rate of strain
 DenseMatrix<double> strain_rate(3,3);

 // Loop over elements to plot
 unsigned nelem=mesh_pt()->nelement();
 for (unsigned e=0;e<nelem;e++)
  {
   // Get pointer to element
   ELEMENT* el_pt=dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   // Tecplot header info
   outfile1 << el_pt->tecplot_zone_string(nplot);
   outfile2 << el_pt->tecplot_zone_string(nplot);
 
   // Loop over plot points
   unsigned num_plot_points=el_pt->nplot_points(nplot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     
     // Get local coordinates of plot point
     el_pt->get_s_plot(iplot,nplot,s);
     
     // Get global coordinates of plot point
     el_pt->interpolated_x(s,x);
     
     // Get strain
     el_pt->strain_rate(s,strain_rate);
     outfile1 << x[0] << " " << x[1] << " " 
              << strain_rate(0,0) << " " 
              << strain_rate(0,1) << " " 
              << strain_rate(0,2) << " " 
              << strain_rate(1,1) << " " 
              << strain_rate(1,2) << " " 
              << strain_rate(2,2) << " "
              << std::endl;

     // Get exact strain
     ValidateRateOfStrain::exact_strain(x,strain_rate);
     outfile2 << x[0] << " " << x[1] << " " 
              << strain_rate(0,0) << " " 
              << strain_rate(0,1) << " " 
              << strain_rate(0,2) << " " 
              << strain_rate(1,1) << " " 
              << strain_rate(1,2) << " " 
              << strain_rate(2,2) << " "
              << std::endl;
    }
  }
 outfile1.close();
 outfile2.close();
}



//==start_of_unsteady_run=================================================
/// Perform run up to specified time
//========================================================================
template<class ELEMENT, class TIMESTEPPER>
void RefineableRotatingCylinderProblem<ELEMENT,TIMESTEPPER>::unsteady_run(
 const double& t_max, const double& dt, string dir_name, 
 const unsigned &maximum_ref_level,const unsigned &minimum_ref_level)
{
 // Setup labels for output
 //-------------------------
 DocInfo doc_info;

 // Output directory
 doc_info.set_directory(dir_name); 

 // Step number
 doc_info.number()=0;

 // Open trace file
 char filename[100];   
 sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
 Trace_file.open(filename);

 // Initialise Trace file
 Trace_file << "time" << ", " 
            << "max permitted error" << ", "
            << "min permitted error" << ", "
            << "max error" << ", "
            << "min error" << ", " << std::endl;

 // Set IC
 assign_initial_values_impulsive(dt);

 // Over-ride the maximum and minimum permitted errors
 mesh_pt()->max_permitted_error() = 1.0e-2; //Default = 1.0e-3
 mesh_pt()->min_permitted_error() = 1.0e-3; //Default = 1.0e-5
 
 // Over-ride the maximum and minimum permitted refinement levels
 mesh_pt()->max_refinement_level() = maximum_ref_level;
 mesh_pt()->min_refinement_level() = minimum_ref_level;

// Initial refinement level
 for (unsigned count=0;count<minimum_ref_level;count++)
  {
   refine_uniformly();
  }
 
 // Max. number of spatial adaptations per timestep
 unsigned max_adapt=maximum_ref_level-minimum_ref_level;



//  validate_strain();
//  pause("called validate_strain");


 // Number of steps
 unsigned ntsteps=unsigned(t_max/dt);
 
 // First timestep?
 bool first=true;
 
 //Doc initial solution
 doc_solution(doc_info);

 //Increment counter for solutions 
 doc_info.number()++;
 
 // Specify normalising factor explicitly
 Z2ErrorEstimator* error_pt=dynamic_cast<Z2ErrorEstimator*>(
  mesh_pt()->spatial_error_estimator_pt());
 error_pt->reference_flux_norm() = 0.01;


 // Timestepping loop
 for (unsigned istep=0;istep<ntsteps;istep++)
  {
   // Take fixed timestep with spatial adaptivity
   unsteady_newton_solve(dt,max_adapt,first);
   
   // set first to false, since no longer first step
   first=false; 
   // set maximum adaptations to 1
   max_adapt=1;
   
   // Output solution
   doc_solution(doc_info);
   // Increment counter for solutions 
   doc_info.number()++;
   
   cout << istep <<std::endl;
  }
 
} // end of unsteady run

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
 


//==start_of_main======================================================
/// Driver for RotatingCylinderProblem
//=====================================================================
int main(int argc, char* argv[])
{
 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 ///Maximum Time
 double t_max=1.0;
 ///Duration of Time-step
 double dt=0.01;
 
 // number of elements in x direction
 unsigned nx=2;
 // number of elements in y direction
 unsigned ny=2;
 // length in x direction
 double lx=1.0;
 // length in y direction
 double ly=1.3;

 /// maximum refinement level
 unsigned maximum_ref_level=4;

 /// minimum refinement level
 unsigned minimum_ref_level=1;

 // If validation run use smaller number of timesteps
 if (CommandLineArgs::Argc>1)
  {
   ///Maximum Time
   t_max=0.02;
  }

 // Do Taylor Hood 
 {
  RefineableRotatingCylinderProblem<
   RefineableAxisymmetricQTaylorHoodElement,BDF<2> > 
   problem(nx,ny,lx,ly);

  cout << "Doing RefineableAxisymmetricQTaylorHoodElement" << std::endl;
  problem.unsteady_run(t_max,dt,"RESLT_TH",
                       maximum_ref_level,minimum_ref_level);
 }
  
 // Do Crouzeix Raviart
 {
  RefineableRotatingCylinderProblem<
   RefineableAxisymmetricQCrouzeixRaviartElement, BDF<2> > 
   problem(nx,ny,lx,ly);

  cout << "Doing RefineableAxisymmetricQCrouzeixRaviartElement" << std::endl;

  problem.unsteady_run(t_max,dt,"RESLT_CR",
                       maximum_ref_level,minimum_ref_level);
 }

} // end of main






