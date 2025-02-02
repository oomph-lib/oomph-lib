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
// Driver for 2-D Channel with changing width, using Taylor Hood 
// and Crouzeix Raviart elements.

// Generic oomph-lib header
#include "generic.h"

// Navier Stokes headers
#include "navier_stokes.h"

// The mesh
#include "meshes/channel_spine_mesh.h"

using namespace std;

using namespace oomph;

//==start_of_namespace===================================================
/// Namespace for physical parameters
//=======================================================================
namespace Global_Physical_Variables
{

 /// Reynolds number
 double Re=100;

} // end_of_namespace





/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
// Sinusoidal wall as geometric object
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



//=========================================================================
/// Geometric object representing a sinusoidal wall, parametrised by
///  \f[ x = \zeta \f]
///  \f[ y = H + A \sin\left(\frac{2\pi \left(\zeta-\zeta_{\mbox{min}}\right)}
///                         {\zeta_{\mbox{max}}-\zeta_{\mbox{min}}}
///                          \right)\f]
//=========================================================================
class SinusoidalWall : public GeomObject
{

public:

 /// Constructor:  Pass height, amplitude, zeta min and zeta max
 /// (all are pinned by default)
 SinusoidalWall(const double& height, const double& amplitude,
               const double& zeta_min, const double& zeta_max)
  : GeomObject(1,2)
  {
   // Make space for the geometric data
   Geom_data_pt.resize(1);
   
   // Create data: Four value, no time dependence, free by default
   Geom_data_pt[0] = new Data(4);

   // I've created the data, I have to clean up
   Must_clean_up=true;
   
   // Pin the data
   for (unsigned i=0;i<4;i++) {Geom_data_pt[0]->pin(i);}
   
   // Give it a value: Initial height
   Geom_data_pt[0]->set_value(0,height);
   Geom_data_pt[0]->set_value(1,amplitude);
   Geom_data_pt[0]->set_value(2,zeta_min);
   Geom_data_pt[0]->set_value(3,zeta_max);
  }


 /// Constructor:  One item of geometric data, containing
 /// four values:
 /// \code
 ///  Geom_data_pt[0]->value(0) = height
 ///  Geom_data_pt[0]->value(1) = amplitude
 ///  Geom_data_pt[0]->value(2) = zeta_min
 ///  Geom_data_pt[0]->value(3) = zeta_max
 /// \endcode
 SinusoidalWall(const Vector<Data*>& geom_data_pt) : GeomObject(1,2)
  {
#ifdef PARANOID
   if(geom_data_pt.size()!=1)
    {
     std::ostringstream error_stream;
     error_stream 
      << "Wrong size for geom_data_pt " << geom_data_pt.size() << std::endl;
     if (geom_data_pt[0]->nvalue()!=1)
      {
       error_stream << "Wrong geom_data_pt[0]->nvalue() " 
                    << geom_data_pt[0]->nvalue() << std::endl;
      }
     
     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif
   Geom_data_pt.resize(1);
   Geom_data_pt[0]=geom_data_pt[0];

   // Data has been created externally: Must not clean up
   Must_clean_up=false;
  }



 /// Destructor:  Clean up if necessary
 ~SinusoidalWall()
  {
   // Do I need to clean up?
   if (Must_clean_up)
    {
     delete Geom_data_pt[0];
     Geom_data_pt[0]=0;
    }
  }

 /// Access function for horizontal half axis
 double& height(){return *Geom_data_pt[0]->value_pt(0);}

 /// Access function for vertical half axis
 double& amplitude(){return *Geom_data_pt[0]->value_pt(1);}


 /// Position vector at Lagrangian coordinate zeta 
 void position(const Vector<double>& zeta, Vector<double>& r) const
  {
#ifdef PARANOID
   if (r.size()!=Ndim)
    {
     throw OomphLibError("The vector r has the wrong dimension\n",
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif
   // Get parameters
   double H = Geom_data_pt[0]->value(0);
   double A = Geom_data_pt[0]->value(1);
   double zeta_min = Geom_data_pt[0]->value(2);
   double zeta_max = Geom_data_pt[0]->value(3);
   double Lz = zeta_max-zeta_min;
   
   // Position Vector
   r[0] = zeta[0];
   r[1] = H + A*sin(MathematicalConstants::Pi*(zeta[0]-zeta_min)/Lz);
  }


 /// Parametrised position on object: r(zeta). Evaluated at
 /// previous timestep. t=0: current time; t>0: previous
 /// timestep. 
 void position(const unsigned& t, const Vector<double>& zeta,
               Vector<double>& r) const
  {
#ifdef PARANOID
   if (t>Geom_data_pt[0]->time_stepper_pt()->nprev_values())
    {
     std::ostringstream error_stream;
     error_stream 
      << "t > nprev_values() in SpikedLine:  " << t << " " 
      << Geom_data_pt[0]->time_stepper_pt()->nprev_values() << std::endl;

     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif

   // Get parameters
   double H = Geom_data_pt[0]->value(t,0);
   double A = Geom_data_pt[0]->value(t,1);
   double zeta_min = Geom_data_pt[0]->value(t,2);
   double zeta_max = Geom_data_pt[0]->value(t,3);
   double Lz = zeta_max-zeta_min;
   
   // Position Vector at time level t
   r[0] = zeta[0];
   r[1] = H + A*sin(MathematicalConstants::Pi*(zeta[0]-zeta_min)/Lz);
  }



 /// Derivative of position Vector w.r.t. to coordinates: 
 /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i). 
 /// Evaluated at current time.
 virtual void dposition(const Vector<double>& zeta, 
                        DenseMatrix<double> &drdzeta) const
  {
   // Get parametres
   double A = Geom_data_pt[0]->value(1);
   double zeta_min = Geom_data_pt[0]->value(2);
   double zeta_max = Geom_data_pt[0]->value(3);
   double Lz = zeta_max-zeta_min;

   // Tangent vector
   drdzeta(0,0)=1.0;
   drdzeta(0,1)=MathematicalConstants::Pi*A*
    cos(MathematicalConstants::Pi*(zeta[0]-zeta_min)/Lz)/Lz;
  }


 /// 2nd derivative of position Vector w.r.t. to coordinates: 
 /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ = ddrdzeta(alpha,beta,i). 
 /// Evaluated at current time.
 virtual void d2position(const Vector<double>& zeta, 
                         RankThreeTensor<double> &ddrdzeta) const
  {
   // Get parametres
   double A = Geom_data_pt[0]->value(1);
   double zeta_min = Geom_data_pt[0]->value(2);
   double zeta_max = Geom_data_pt[0]->value(3);
   double Lz = zeta_max-zeta_min;

   // Derivative of tangent vector
   ddrdzeta(0,0,0)=0.0;
   ddrdzeta(0,0,1)=-MathematicalConstants::Pi*MathematicalConstants::Pi*
    A*sin(MathematicalConstants::Pi*(zeta[0]-zeta_min)/Lz)/(Lz*Lz);
  }


 /// Posn Vector and its  1st & 2nd derivatives
 /// w.r.t. to coordinates:
 /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i). 
 /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ = ddrdzeta(alpha,beta,i). 
 /// Evaluated at current time.
 virtual void d2position(const Vector<double>& zeta, Vector<double>& r,
                         DenseMatrix<double> &drdzeta,
                         RankThreeTensor<double> &ddrdzeta) const
  {
   // Get parametres
   double H = Geom_data_pt[0]->value(0);
   double A = Geom_data_pt[0]->value(1);
   double zeta_min = Geom_data_pt[0]->value(2);
   double zeta_max = Geom_data_pt[0]->value(3);
   double Lz = zeta_max-zeta_min;

   // Position Vector
   r[0] = zeta[0];
   r[1] = H + A*sin(MathematicalConstants::Pi*(zeta[0]-zeta_min)/Lz);

   // Tangent vector
   drdzeta(0,0)=1.0;
   drdzeta(0,1)=MathematicalConstants::Pi*A*
    cos(MathematicalConstants::Pi*(zeta[0]-zeta_min)/Lz)/Lz;

   // Derivative of tangent vector
   ddrdzeta(0,0,0)=0.0;
   ddrdzeta(0,0,1)=-MathematicalConstants::Pi*MathematicalConstants::Pi*
    A*sin(MathematicalConstants::Pi*(zeta[0]-zeta_min)/Lz)/(Lz*Lz);
  }

 /// How many items of Data does the shape of the object depend on?
 unsigned ngeom_data() const {return Geom_data_pt.size();}
 
 /// Return pointer to the j-th Data item that the object's 
 /// shape depends on 
 Data* geom_data_pt(const unsigned& j) {return Geom_data_pt[j];}
 

private:

 /// Vector of pointers to Data items that affects the object's shape
 Vector<Data*> Geom_data_pt;

 /// Do I need to clean up?
 bool Must_clean_up;

};

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Channel flow through a non-uniform channel whose geometry is defined 
/// by a spine mesh.
//====================================================================
template<class ELEMENT>
class ChannelSpineFlowProblem : public Problem
{

public:

 /// Constructor
 ChannelSpineFlowProblem();
 
 /// Destructor: (empty)
 ~ChannelSpineFlowProblem(){}

 /// Update the problem specs before solve. 
 /// Set velocity boundary conditions just to be on the safe side...
 void actions_before_newton_solve()
  { 
   // Update the mesh
   mesh_pt()->node_update();

  } // end_of_actions_before_newton_solve


 /// Update the after solve (empty)
 void actions_after_newton_solve(){}

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

private:
 
 /// Width of channel
 double Ly;

 
}; // end_of_problem_class



//==start_of_constructor==================================================
/// Constructor for ChannelSpineFlow problem 
//========================================================================
template<class ELEMENT>
ChannelSpineFlowProblem<ELEMENT>::ChannelSpineFlowProblem()
{ 

 // Setup mesh

 // # of elements in x-direction in left region 
 unsigned Nx0=3;
 // # of elements in x-direction in centre region
 unsigned Nx1=12;
 // # of elements in x-direction in right region
 unsigned Nx2=8;

 // # of elements in y-direction
 unsigned Ny=10;

 // Domain length in x-direction in left region
 double Lx0=0.5;
 // Domain length in x-direction in centre region
 double Lx1=0.7;
 // Domain length in x-direction in right region
 double Lx2=1.5;
 
 // Domain length in y-direction
 Ly=1.0;
 
 // Build geometric object that represents the sinusoidal bump on
 // the upper wall:

 // 40% indendentation
 double amplitude_upper = -0.4*Ly;
 // Minimum and maximum coordinates of bump
 double zeta_min=Lx0;
 double zeta_max=Lx0+Lx1;
 GeomObject* UpperWall = 
  new SinusoidalWall(Ly,amplitude_upper,zeta_min,zeta_max);
 
 // Build and assign mesh -- pass pointer to geometric object
 // that represents the sinusoidal bump on the upper wall
 Problem::mesh_pt() = new ChannelSpineMesh<ELEMENT>(Nx0,Nx1,Nx2,Ny,
                                                    Lx0,Lx1,Lx2,Ly,
                                                    UpperWall);


 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here: All boundaries are Dirichlet boundaries, except on boundary 1
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     if (ibound!=1)
      {
       // Loop over values (u and v velocities)
       for (unsigned i=0;i<2;i++)
        {
         mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
        }
      }
     else
      {
       // Parallel outflow ==> no-slip
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(1); 
      }
    }
  } // end loop over boundaries

 
 // No slip on stationary upper and lower walls (boundaries 0 and 2) 
 // and parallel outflow (boundary 1) 
 for (unsigned ibound=0;ibound<num_bound-1;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     if (ibound!=1)
      {
       for (unsigned i=0;i<2;i++)
        {
         mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
        }
      }
     else
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
      }
    }
  }


 // Setup parabolic inflow along boundary 3:
 unsigned ibound=3; 
 unsigned num_nod= mesh_pt()->nboundary_node(ibound);
 for (unsigned inod=0;inod<num_nod;inod++)
  {
   double y=mesh_pt()->boundary_node_pt(ibound,inod)->x(1);
   // Parallel, parabolic inflow
   mesh_pt()->boundary_node_pt(ibound,inod)->set_value(0,y*(Ly-y));
   mesh_pt()->boundary_node_pt(ibound,inod)->set_value(1,0.0);
  }

 // Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor: Pass 
 // pointer to Reynolds number
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Physical_Variables::Re;
  } 

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
} // end_of_constructor




//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void ChannelSpineFlowProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();
 
} // end_of_doc_solution



//==start_of_main======================================================
/// Driver for channel flow problem with spine mesh.
//=====================================================================
int main()
{

 // Set output directory
 DocInfo doc_info;
 doc_info.set_directory("RESLT");
 doc_info.number()=0;
 
 // Solve problem with Taylor Hood elements
 //---------------------------------------
 {
  //Build problem
  ChannelSpineFlowProblem<SpineElement<QTaylorHoodElement<2> > >
   problem;
  
  // Solve the problem with automatic adaptation
  problem.newton_solve();
  
  //Output solution
  problem.doc_solution(doc_info);
  // Step number
  doc_info.number()++;

 } // end of Taylor Hood elements
 
 
 // Solve problem with Crouzeix Raviart elements
 //--------------------------------------------
 {
  // Build problem
  ChannelSpineFlowProblem<SpineElement<QCrouzeixRaviartElement<2> > >
   problem;
  
  // Solve the problem with automatic adaptation
  problem.newton_solve();
  
  //Output solution
  problem.doc_solution(doc_info);
  // Step number
  doc_info.number()++;
  
 } // end of Crouzeix Raviart elements
      

     
} // end_of_main




//====================================================================
/// Create the files to illustrate the sparse node update operation
//====================================================================
void doc_sparse_node_update()
{

 // Set output directory
 DocInfo doc_info;
 doc_info.set_directory("RESLT");

 // Setup mesh
 
 // # of elements in x-direction in left region 
 unsigned Nx0=1;
 // # of elements in x-direction in centre region
 unsigned Nx1=5;
 // # of elements in x-direction in right region
 unsigned Nx2=1;
 
 // # of elements in y-direction
 unsigned Ny=5;
 
 // Domain length in x-direction in left region
 double Lx0=0.5;
 // Domain length in x-direction in centre region
 double Lx1=1.0;
 // Domain length in x-direction in right region
 double Lx2=0.5;
 
 // Domain length in y-direction
 double Ly=1.0;
 
 double amplitude_upper = -0.4*Ly;
 double zeta_min=Lx0;
 double zeta_max=Lx0+Lx1;
 GeomObject* UpperWall = 
  new SinusoidalWall(Ly,amplitude_upper,zeta_min,zeta_max);
 
 // Build and assign mesh
 ChannelSpineMesh<SpineElement<QTaylorHoodElement<2> > >* mesh_pt =
  new ChannelSpineMesh<SpineElement<QTaylorHoodElement<2> > >(Nx0,Nx1,Nx2,Ny,
                                                             Lx0,Lx1,Lx2,Ly,
                                                             UpperWall);
 
 // Update *all* nodal positions
 mesh_pt->node_update();
 
 // Change the amplitude
 dynamic_cast<SinusoidalWall*>(mesh_pt->wall_pt())->amplitude()=0.5;
 
 unsigned count=0;
 ofstream some_file;
 char filename[100];
 
 // Number of plot points
 unsigned npts=5; 
 
 // Output solution 
 sprintf(filename,"%s/mesh_update%i.dat",doc_info.directory().c_str(),
         count);
 count++;
 
 // Loop over spines in centre
 unsigned n_node = mesh_pt->nnode();
 for (unsigned inode=0;inode<n_node;inode++)
  {
   SpineNode* node_pt=dynamic_cast<SpineNode*>(
    mesh_pt->node_pt(inode));
   if (node_pt->node_update_fct_id()==1)
    {
     node_pt->node_update();
     // Output solution 
     some_file.open(filename);
     sprintf(filename,"%s/mesh_update%i.dat",doc_info.directory().c_str(),
             count);
     count++;
     mesh_pt->output(some_file,npts);
     some_file.close();
    }
  }
 
}

     












