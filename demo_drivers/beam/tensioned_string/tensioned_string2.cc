//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2024 Matthias Heil and Andrew Hazel
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
//Driver function for a simple beam proble

//OOMPH-LIB includes
#include "generic.h"
#include "beam.h"
#include "meshes/one_d_lagrangian_mesh.h"

using namespace std;

using namespace oomph;

//=========================================================================
/// Steady, straight 1D line in 2D space with stretch_ratio
///  \f[ x = \zeta stretch_ratio \f]
///  \f[ y = 0  \f]
//=========================================================================
class NewStraightLine: public GeomObject
{
public:
 
  /// Constructor derives from GeomObject(1, 2)
  /// Pass stretch_ratio to this class (defaults to 1 in which case
  /// the Lagrangian coordinate is the arclength)
  NewStraightLine(const double& stretch_ratio=1.0) : GeomObject(1, 2)
  {
    Stretch_ratio = stretch_ratio;
  }

  /// Broken copy constructor
  NewStraightLine(const NewStraightLine& dummy) = delete;

  /// Broken assignment operator
  void operator=(const NewStraightLine&) = delete;

  /// Position Vector at Lagrangian coordinate zeta
  void position(const Vector<double>& zeta, Vector<double>& r) const
  {
    r[0] = zeta[0] * Stretch_ratio;
    r[1] = 0.0;
  }


  /// Derivative of position Vector w.r.t. to coordinates:
  /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
  /// Evaluated at current time.
  virtual void dposition(const Vector<double>& zeta,
                         DenseMatrix<double>& drdzeta) const
  {
    // Tangent vector
    drdzeta(0, 0) = Stretch_ratio;
    drdzeta(0, 1) = 0.0;
  }


  /// 2nd derivative of position Vector w.r.t. to coordinates:
  /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
  /// ddrdzeta(alpha,beta,i). Evaluated at current time.
  virtual void d2position(const Vector<double>& zeta,
                          RankThreeTensor<double>& ddrdzeta) const
  {
    // Derivative of tangent vector
    ddrdzeta(0, 0, 0) = 0.0;
    ddrdzeta(0, 0, 1) = 0.0;
  }


  /// Posn Vector and its  1st & 2nd derivatives
  /// w.r.t. to coordinates:
  /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
  /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
  /// ddrdzeta(alpha,beta,i).
  /// Evaluated at current time.
  virtual void d2position(const Vector<double>& zeta,
                          Vector<double>& r,
                          DenseMatrix<double>& drdzeta,
                          RankThreeTensor<double>& ddrdzeta) const
  {
    // Position Vector
    r[0] = zeta[0] * Stretch_ratio;
    r[1] = 0.0;

    // Tangent vector
    drdzeta(0, 0) = Stretch_ratio;
    drdzeta(0, 1) = 0.0;

    // Derivative of tangent vector
    ddrdzeta(0, 0, 0) = 0.0;
    ddrdzeta(0, 0, 1) = 0.0;
  }

private:

 /// Define the length of the beam
 double Stretch_ratio;
 
};




//========start_of_namespace========================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{
 /// Non-dimensional thickness
 double H;

 /// 2nd Piola Kirchhoff pre-stress
 double Sigma0;

 /// Pressure load
 double P_ext;

 /// Load function: Apply a constant external pressure to the beam
 void load(const Vector<double>& xi, const Vector<double> &x,
           const Vector<double>& N, Vector<double>& load)
 {
  for(unsigned i=0;i<2;i++) {load[i] = -P_ext*N[i];}
 }

} // end of namespace

//======start_of_problem_class==========================================
/// Beam problem object
//======================================================================
class ElasticBeamProblem : public Problem
{
public:
 
 /// Constructor: The arguments are the number of elements, 
 /// the length of domain
 ElasticBeamProblem(const unsigned &n_elem, const double &length,
                    const double& stretch_ratio);
 
 /// Conduct a parameter study
 void parameter_study();
 
 /// Return pointer to the mesh
 OneDLagrangianMesh<HermiteBeamElement>* mesh_pt() 
  {return dynamic_cast<OneDLagrangianMesh<HermiteBeamElement>*>
    (Problem::mesh_pt());}

 /// No actions need to be performed after a solve
 void actions_after_newton_solve() {}

 /// No actions need to be performed before a solve
 void actions_before_newton_solve() {}

private:

 /// Pointer to the node whose displacement is documented
 Node* Doc_node_pt;

 /// Length of domain (in terms of the Lagrangian coordinates)
 double Length;

 /// Pointer to geometric object that represents the beam's undeformed shape
 GeomObject* Undef_beam_pt;

}; // end of problem class


//=============start_of_constructor=====================================
/// Constructor for elastic beam problem
//======================================================================
ElasticBeamProblem::ElasticBeamProblem(const unsigned &n_elem,
                                       const double &length,
                                       const double& stretch_ratio) : Length(length)
{
 // Set the undeformed beam to be a straight line at y=0
 Undef_beam_pt=new NewStraightLine(stretch_ratio); 

 // Create the (Lagrangian!) mesh, using the geometric object
 // Undef_beam_pt to specify the initial (Eulerian) position of the
 // nodes.

 // Compensate for stretched coordinate by adjusting nominal length (in terms
 // of the coordinate)
 double compensated_length=length/stretch_ratio;
 Problem::mesh_pt() = 
  new OneDLagrangianMesh<HermiteBeamElement>(n_elem,compensated_length,Undef_beam_pt);

 // Set the boundary conditions: Each end of the beam is fixed in space
 // Loop over the boundaries (ends of the beam)
 for(unsigned b=0;b<2;b++)
  {
   // Pin displacements in both x and y directions
   // [Note: The mesh_pt() function has been overloaded
   //  to return a pointer to the actual mesh, rather than
   //  a pointer to the Mesh base class. The current mesh is derived
   //  from the SolidMesh class. In such meshes, all access functions
   //  to the nodes, such as boundary_node_pt(...), are overloaded
   //  to return pointers to SolidNodes (whose position can be
   //  pinned) rather than "normal" Nodes.]
   mesh_pt()->boundary_node_pt(b,0)->pin_position(0); 
   mesh_pt()->boundary_node_pt(b,0)->pin_position(1); 
  }
 
 //Find number of elements in the mesh
 unsigned n_element = mesh_pt()->nelement();
 
 //Loop over the elements to set physical parameters etc.
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast to the specific element type
   HermiteBeamElement *elem_pt = 
    dynamic_cast<HermiteBeamElement*>(mesh_pt()->element_pt(e));
   
   // Set physical parameters for each element:
   elem_pt->sigma0_pt() = &Global_Physical_Variables::Sigma0;
   elem_pt->h_pt() = &Global_Physical_Variables::H;

   // Set the load Vector for each element
   elem_pt->load_vector_fct_pt() = &Global_Physical_Variables::load;

   // Set the undeformed shape for each element
   elem_pt->undeformed_beam_pt() = Undef_beam_pt;
  } // end of loop over elements

 // Choose node at which displacement is documented (halfway along -- provided
 // we have an odd number of nodes; complain if this is not the
 // case because the comparison with the exact solution will be wrong 
 // otherwise!)
 unsigned n_nod=mesh_pt()->nnode();
 if (n_nod%2!=1)
  {
   cout << "Warning: Even number of nodes " << n_nod << std::endl;
   cout << "Comparison with exact solution will be misleading..." << std::endl;
  }
 Doc_node_pt=mesh_pt()->node_pt((n_nod+1)/2-1);
 
 // Assign the global and local equation numbers
 cout << "# of dofs " << assign_eqn_numbers() << std::endl;

} // end of constructor


//=======start_of_parameter_study==========================================
/// Solver loop to perform parameter study
//=========================================================================
void ElasticBeamProblem::parameter_study()
{
 // Over-ride the default maximum value for the residuals
 Problem::Max_residuals = 1.0e10;
 
 // Set the increments in control parameters
 double pext_increment = 0.001;
 
 // Set initial values for control parameters 
 Global_Physical_Variables::P_ext = 0.0 - pext_increment;
 
 // Create label for output
 DocInfo doc_info;
 
 // Set output directory -- this function checks if the output
 // directory exists and issues a warning if it doesn't.
 doc_info.set_directory("RESLT");
 
 // Open a trace file
 ofstream trace("RESLT/trace_beam.dat");
 
 // Write a header for the trace file
 trace << 
  "VARIABLES=\"p_e_x_t\",\"d\"" << 
  ", \"p_e_x_t_(_e_x_a_c_t_)\"" << std::endl;
 
 // Output file stream used for writing results
 ofstream file;
 // String used for the filename
 char filename[100]; 

 // Loop over parameter increments
 unsigned nstep=10;
 for(unsigned i=1;i<=nstep;i++)
  {
   // Increment pressure
   Global_Physical_Variables::P_ext += pext_increment;

   // Solve the system
   newton_solve();
    
   // Calculate exact solution for `string under tension' (applicable for
   // small wall thickness and pinned ends)

   // The tangent of the angle beta
   double tanbeta =-2.0*Doc_node_pt->x(1)/Length;

   double exact_pressure = 0.0;
   //If the beam has deformed, calculate the pressure required
   if(tanbeta!=0)
    {
      
      //Calculate the opening angle alpha
      double alpha = 2.0*atan(2.0*tanbeta/(1.0-tanbeta*tanbeta));

      // Jump back onto the main branch if alpha>180 degrees
      if (alpha<0) alpha+=2.0*MathematicalConstants::Pi;

     // Green strain:
     double gamma=0.5*(0.25*alpha*alpha/(sin(0.5*alpha)*sin(0.5*alpha))-1.0);

     //Calculate the exact pressure
     exact_pressure=Global_Physical_Variables::H*
      (Global_Physical_Variables::Sigma0+gamma)*alpha/Length;
    } 
   
   // Document the solution
   sprintf(filename,"RESLT/beam%i.dat",i);
   file.open(filename);
   mesh_pt()->output(file,5);
   file.close();
   
   // Write trace file: Pressure, displacement and exact solution
   // (for string under tension)
   trace << Global_Physical_Variables::P_ext  << " " 
         << abs(Doc_node_pt->x(1))
         << " " << exact_pressure 
         << std::endl;
  }
 
} // end of parameter study

//========start_of_main================================================
/// Driver for beam (string under tension) test problem.
/// Redo with non-arclength parametrisation.
//=====================================================================
int main(int argc, char** argv)
{
 
 // Store command line arguments
 CommandLineArgs::setup(argc, argv);
 
 // Stretch ratio
 double stretch_ratio=1.0;
 CommandLineArgs::specify_command_line_flag(
  "--stretch_ratio",&stretch_ratio);

 // Set the 2nd Piola Kirchhoff prestress
 Global_Physical_Variables::Sigma0=0.1; 
 CommandLineArgs::specify_command_line_flag(
  "--sigma0",&Global_Physical_Variables::Sigma0);

  // Parse command line
  CommandLineArgs::parse_and_assign();

  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();

 // Set the non-dimensional thickness 
 Global_Physical_Variables::H=0.01; 
  
 // Set the length of domain
 double L = 10.0;

 // Number of elements (choose an even number if you want the control point 
 // to be located at the centre of the beam)
 unsigned n_element = 10;

 // Construst the problem
 ElasticBeamProblem problem(n_element,L,stretch_ratio);

 // Conduct parameter study
 problem.parameter_study();

} // end of main

