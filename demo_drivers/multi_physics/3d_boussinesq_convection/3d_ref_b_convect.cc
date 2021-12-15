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
//Driver for a multi-physics problem that couples the Navier--Stokes
//equations to the advection diffusion equations,
//giving Boussinesq convection

//Oomph-lib headers, we require the generic, advection-diffusion
//and navier-stokes elements.
#include "generic.h"
#include "advection_diffusion.h"
#include "navier_stokes.h"

// The mesh is our standard rectangular quadmesh
#include "meshes/simple_cubic_mesh.h"

// Use the oomph and std namespaces
using namespace oomph;
using namespace std;

//================================================================
/// Function-type-object to perform comparison of elements
//================================================================
class ElementCmp
{
public:

 /// Comparison. Are the values identical or not?
 bool operator()(GeneralisedElement* const &x, GeneralisedElement* const &y)
  const
  {
   FiniteElement* cast_x = dynamic_cast<FiniteElement*>(x);
   FiniteElement* cast_y = dynamic_cast<FiniteElement*>(y);

   if((cast_x ==0) || (cast_y==0)) {return 0;}
   else
    {return cast_x->node_pt(0)->x(2) < cast_y->node_pt(0)->x(2);}
  }
};


//============================================================================
/// Boussinesq preconditioner. This extracts upper/lower triangular
/// blocks in the 3x3 overall block matrix structure arising from
/// the monolithic discretisation of Boussinesq problems.
/// Dofs are decomposed into fluid velocity, pressure
/// and temperature unknowns. NavierStokesSchurComplementPreconditioner is used
/// as the inexact solver for the fluid block; SuperLU (in
/// its incarnation as an "exact" preconditioner) is used for
/// the temperature block. By default we retain the fluid on temperature off
/// diagonal blocks.
//=============================================================================
class BoussinesqPreconditioner :
 public BlockPreconditioner<CRDoubleMatrix>
{

public :

 /// Constructor: By default use block triangular form with retained
 /// fluid on temperature terms. Requires a non-constant problem pointer
 /// because of the underlying NavierStokesSchurComplementPreconditioner.
 BoussinesqPreconditioner(Problem* problem_pt)
  {
   // Default setting: Fluid onto temperature
   Retain_temperature_onto_fluid_terms=false;
   Retain_fluid_onto_temperature_terms=true;

   // Create the Navier Stokes preconditioner
   Navier_stokes_preconditioner_pt =
    new NavierStokesSchurComplementPreconditioner(problem_pt);
   Navier_stokes_preconditioner_pt->disable_doc_time();

   //Set the temperature preconditioner
   Temperature_preconditioner_pt = new SuperLUPreconditioner;

   //Initialise the P and F block preconditioners
   P_preconditioner_pt = 0;
   F_preconditioner_pt = 0;

#ifdef OOMPH_HAS_HYPRE
//If we have compiled with MPI,
//only use HYPRE if it's been initialised
#ifdef OOMPH_HAS_MPI
   if(MPI_Helpers::mpi_has_been_initialised())
#endif
    {
     //Set up the internal preconditioners
     P_preconditioner_pt = new HyprePreconditioner;

     // Set parameters for use as preconditioner on Poisson-type problem
     Hypre_default_settings::set_defaults_for_3D_poisson_problem(
      static_cast<HyprePreconditioner*>(P_preconditioner_pt));

     // Use Hypre for the Schur complement block
     Navier_stokes_preconditioner_pt->
      set_p_preconditioner(P_preconditioner_pt);

     F_preconditioner_pt = new HyprePreconditioner;

     // Set parameters for use as preconditioner in for momentum
     // block in Navier-Stokes problem
     Hypre_default_settings::set_defaults_for_navier_stokes_momentum_block(
      static_cast<HyprePreconditioner*>(F_preconditioner_pt));

     // Use Hypre for momentum block
     Navier_stokes_preconditioner_pt->
      set_f_preconditioner(F_preconditioner_pt);

     //Set the Temperature preconditioner to also use AMG
     delete Temperature_preconditioner_pt; Temperature_preconditioner_pt = 0;

     Temperature_preconditioner_pt = new HyprePreconditioner;
     Hypre_default_settings::set_defaults_for_navier_stokes_momentum_block(
      static_cast<HyprePreconditioner*>(Temperature_preconditioner_pt));
    }
#endif

   // Preconditioner hasn't been set up yet.
   Preconditioner_has_been_setup=false;

   // Null out pointers to additional block matrices not already
   // stored in Navier Stokes block preconditioner
   Block_matrix_0_1_pt=0;
   Block_matrix_1_0_pt=0;

   // set Doc_time to false
   Doc_time = false;

   // set the mesh to be the main mesh of the problem
   this->set_nmesh(1);
   this->set_mesh(0,problem_pt->mesh_pt());

   // Set up the (non-const) mesh pointer for the schur complement
   // preconditioner.
   Navier_stokes_preconditioner_pt->set_navier_stokes_mesh
    (problem_pt->mesh_pt());
  }


 /// Destructor: Clean up.
 ~BoussinesqPreconditioner()
  {
   // Do what it says
   clean_up_memory();

   //Delete the P block preconditioners
   delete P_preconditioner_pt;  P_preconditioner_pt = 0;

   //Delete the F block preconditioner
   delete F_preconditioner_pt; F_preconditioner_pt = 0;

   //Delete the Navier-Stokes preconditioner (inexact solver)
   delete Navier_stokes_preconditioner_pt; Navier_stokes_preconditioner_pt = 0;

   //Delete the temperature preconditioner (inexact solver)
   delete Temperature_preconditioner_pt; Temperature_preconditioner_pt = 0;
  }


 /// Broken copy constructor
 BoussinesqPreconditioner(const BoussinesqPreconditioner&)
  {
   BrokenCopy::broken_copy("BoussinesqPreconditioner");
  }


 /// Switch to block-diagonal preconditioner
 void use_block_diagonal_version()
  {
   Retain_temperature_onto_fluid_terms=false;
   Retain_fluid_onto_temperature_terms=false;
  }

 /// Switch to block-triangular preconditioner in which
 /// action of fluid dofs onto temperature equations is retained
 void use_block_triangular_version_with_fluid_on_temperature()
  {
   Retain_temperature_onto_fluid_terms=false;
   Retain_fluid_onto_temperature_terms=true;
  }

 /// Switch to block-triangular preconditioner in which
 /// action of temperature dofs onto fluid equations is retained
 void use_block_triangular_version_with_temperature_on_fluid()
  {
   Retain_temperature_onto_fluid_terms=true;
   Retain_fluid_onto_temperature_terms=false;
  }

 /// Setup the preconditioner
 void setup();

 /// Apply preconditioner to r
 void preconditioner_solve(const DoubleVector &r,
                           DoubleVector &z);

 /// Access function to the Navier Stokes preconditioner (inexact solver)
 NavierStokesSchurComplementPreconditioner*
 navier_stokes_preconditioner_pt() const
  {
   return Navier_stokes_preconditioner_pt;
  }

 /// Enable documentation of timings
 void enable_doc_time() {Doc_time = true;}

 /// Disable documentation of timings
 void disable_doc_time() {Doc_time = false;}


private:

 /// Clean up memory, delete the blocks allocated in setup
 void clean_up_memory()
  {
   if (Preconditioner_has_been_setup)
    {
     if (Retain_temperature_onto_fluid_terms)
      {
       delete Block_matrix_0_1_pt;
       Block_matrix_0_1_pt=0;
      }
     if (Retain_fluid_onto_temperature_terms)
      {
       delete Block_matrix_1_0_pt;
       Block_matrix_1_0_pt=0;
      }
    }
  }

 /// Pointer the Navier Stokes preconditioner (inexact solver)
 NavierStokesSchurComplementPreconditioner* Navier_stokes_preconditioner_pt;

 /// Pointer to the temperature preconditioner  (inexact solver)
 Preconditioner* Temperature_preconditioner_pt;

 /// Inexact solver for P block
 Preconditioner* P_preconditioner_pt;

 /// Inexact solver for F block
 Preconditioner* F_preconditioner_pt;

 /// Pointer to fluid/temperature interaction matrix
 CRDoubleMatrix* Block_matrix_0_1_pt;
// MatrixVectorProduct* Mat_vec_operator_0_1_pt;

 /// Pointer to temperature/fluid interaction matrix
 CRDoubleMatrix* Block_matrix_1_0_pt;
// MatrixVectorProduct* Mat_vec_operator_1_0_pt;

 /// Boolean indicating the preconditioner has been set up
 bool Preconditioner_has_been_setup;

 /// Boolean flag used to indicate that the temperature onto fluid
 /// interaction terms are to be retained
 bool Retain_temperature_onto_fluid_terms;

 /// Boolean flag used to indicate that the fluid onto temperature
 /// interaction terms are to be retained
 bool Retain_fluid_onto_temperature_terms;

 /// Set Doc_time to true for outputting results of timings
 bool Doc_time;
 };


/// ///////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////
// Boussinesq preconditioner member functions
/// ///////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////




//=============================================================================
/// Setup the preconditioner. Note: Matrix must be a CRDoubleMatrix.
//=============================================================================
void BoussinesqPreconditioner::setup()
 {

  // Clean up any existing data
  this->clean_up_memory();

  // Get total number of dofs.
  unsigned n_dof=mesh_pt(0)->finite_element_pt(0)->dim()+2;

  // Number of fluid dofs: one less
  unsigned n_fluid_dof = n_dof-1;

  // This preconditioner has two types of DOF: Fluid dofs and temperature dofs

  // Initialise block type to zero (correct for fluid dofs)
  Vector<unsigned> dof_to_block_map(n_dof,0);

  // Temperature dofs are "type one"
  dof_to_block_map[n_dof-1] = 1;

  // Call block setup for this preconditioner
  this->block_setup(dof_to_block_map);

  // Block mapping for the subsidiary Navier Stokes preconditioner:
  // blocks 0 and 1 in the Boussinesq preconditioner are also blocks 0 and 1
  // in the subsidiary Navier Stokes one.
  Vector<unsigned> ns_dof_lookup(n_fluid_dof);
  for (unsigned i = 0; i < n_fluid_dof; i++)
   {
    ns_dof_lookup[i] = i;
   }

  // Turn the Navier Stokes Schur complement preconditioner into a
  // subsidiary preconditioner of this preconditioner
  Navier_stokes_preconditioner_pt->
   turn_into_subsidiary_block_preconditioner(this,ns_dof_lookup);

  // Setup the Navier Stokes preconditioner:
  Navier_stokes_preconditioner_pt->setup(matrix_pt());

  // Extract the additional blocks we need for Boussinesq:

  // Temperature matrix
  CRDoubleMatrix* block_matrix_1_1_pt = new CRDoubleMatrix;
  this->get_block(1,1,*block_matrix_1_1_pt);

  // Temperature on fluid terms (if needed)

  if (Retain_temperature_onto_fluid_terms)
   {
    delete Block_matrix_0_1_pt; Block_matrix_0_1_pt = 0;
    Block_matrix_0_1_pt = new CRDoubleMatrix;
    this->get_block(0,1,*Block_matrix_0_1_pt);
   }

  // Fluid on temperature terms (if needed)
  if (Retain_fluid_onto_temperature_terms)
   {
    delete Block_matrix_1_0_pt; Block_matrix_1_0_pt = 0;
    Block_matrix_1_0_pt = new CRDoubleMatrix;
    this->get_block(1,0,*Block_matrix_1_0_pt);
   }


  // Setup the temperature preconditioner (inexact solver)
  double t_start = TimingHelpers::timer();
  Temperature_preconditioner_pt->setup(block_matrix_1_1_pt);
  double t_end = TimingHelpers::timer();
  double setup_time= t_end-t_start;
  delete block_matrix_1_1_pt; block_matrix_1_1_pt = 0;


 // Output times
 if(Doc_time)
  {
   oomph_info << "Temperature sub-preconditioner setup time [sec]: "
              << setup_time << "\n";

   // Doc density of retained interaction block
   if (Retain_temperature_onto_fluid_terms)
    {
     oomph_info
      << "Fill level of temperature on fluid blocks (C_ut and C_pt): " <<
      double(Block_matrix_0_1_pt->nnz())/
      double(Block_matrix_0_1_pt->nrow()*
             Block_matrix_0_1_pt->ncol())*100.0
      << "%  " << std::endl;
    }

   // Doc density of retained interaction block
   if (Retain_fluid_onto_temperature_terms)
    {
     oomph_info
      << "Fill level of fluid on temperature blocks (C_tu and C_tp): " <<
      double(Block_matrix_1_0_pt->nnz())/
      double(Block_matrix_1_0_pt->nrow()*
             Block_matrix_1_0_pt->ncol())*100.0
      << "%  " << std::endl;
    }
  }

  // We're done (and we stored some data)
  Preconditioner_has_been_setup=true;


 }


//======================================================================
/// Apply preconditioner to Vector r
//======================================================================
void BoussinesqPreconditioner::preconditioner_solve(const DoubleVector &r,
                                             DoubleVector &z)
{
 // if z is not setup then give it the same distribution
 if (!z.built())
  {
   z.build(r.distribution_pt(),0.0);
  }

 // Make copy of residual vector (to overcome const-ness
 DoubleVector res(r);


 // Retain off-diagonals that represent effect of temperature on fluid
 //-------------------------------------------------------------------
 if (Retain_temperature_onto_fluid_terms)
  {

   // Working vectors
   DoubleVector temp_temperature_vec;
   DoubleVector temp_fluid_vec;

   // Copy temperature values from residual to temp_vec:
   // Loop over all entries in the global vector (this one
   // includes tempearture, velocity and pressure dofs in some random fashion)
   get_block_vector(1,res,temp_temperature_vec);

   // Solve temperature system by back-substitution
   // with LU-decomposed temperature-temperature Jacobian
   Temperature_preconditioner_pt->preconditioner_solve(temp_temperature_vec,
                                                       temp_temperature_vec);
   this->return_block_vector(1,temp_temperature_vec,z);

   // NOTE: temp_temperature_vec now contains z_t = S^{-1} r_t

   // Multiply C_{ut} by z_t
   Block_matrix_0_1_pt->multiply(temp_temperature_vec,temp_fluid_vec);
//   Mat_vec_operator_0_1_pt->multiply(temp_temperature_vec,temp_fluid_vec);
   temp_temperature_vec.clear();

   // Subtract from fluid residual vector for fluid solve
   DoubleVector another_temp_vec;
   this->get_block_vector(0,res,another_temp_vec);
   another_temp_vec -= temp_fluid_vec;
   this->return_block_vector(0,another_temp_vec,res);

   // now apply the navier stokes lsc preconditioner
   Navier_stokes_preconditioner_pt->preconditioner_solve(res,z);
  }


 // Retain off-diagonals that represent effect of fluid on temperature
 //-------------------------------------------------------------------
 // (or diagonal preconditioner)
 //-----------------------------
 else
  {

   // Call fluid preconditioner for fluid block
   Navier_stokes_preconditioner_pt->preconditioner_solve(res,z);

   // Working vectors
   DoubleVector temp_temperature_vec;

   // get the temperature vector
   get_block_vector(1,res,temp_temperature_vec);

   // Do matrix vector products with fluid onto temperature coupling matrices:
   if (Retain_fluid_onto_temperature_terms)
    {
     DoubleVector temp_fluid_vec;
     get_block_vector(0,z,temp_fluid_vec);

     // Auxiliary vector to hold the matrix vector product of the
     // fluid-onto-temperature coupling matrices with the fluid solutions:
     DoubleVector aux_vec;

     // Multiply C_{tu} by z_u
     Block_matrix_1_0_pt->multiply(temp_fluid_vec, aux_vec);

     // ...and subtract from r_t:
     temp_temperature_vec-=aux_vec;
    }

   // Solve temperature system by back-substitution
   // with LU-decomposed temperature-temperature Jacobian
   Temperature_preconditioner_pt->preconditioner_solve(temp_temperature_vec,
                                                 temp_temperature_vec);

   // Now copy result_vec (i.e. z_t) back into the global vector z.
   // Loop over all entries in the global results vector z:
   return_block_vector(1,temp_temperature_vec,z);
  }
}



//============start_element_class============================================
/// A RefineableElement class that solves the
/// Boussinesq approximation of the Navier--Stokes
/// and energy equations by coupling two pre-existing classes.
/// The RefineableQAdvectionDiffusionElement
/// with bi-quadratic interpolation for the
/// scalar variable (temperature) and
/// RefineableQCrouzeixRaviartElement which solves the Navier--Stokes equations
/// using bi-quadratic interpolation for the velocities and a discontinuous
/// bi-linear interpolation for the pressure. Note that we are free to
/// choose the order in which we store the variables at the nodes. In this
/// case we choose to store the variables in the order fluid velocities
/// followed by temperature. We must, therefore, overload the function
/// AdvectionDiffusionEquations<DIM>::u_index_adv_diff() to indicate that
/// the temperature is stored at the DIM-th position not the 0-th. We do not
/// need to overload the corresponding function in the
/// NavierStokesEquations<DIM> class because the velocities are stored
/// first. Finally, we choose to use the flux-recovery calculation from the
/// fluid velocities to provide the error used in the mesh adaptation.
//==========================================================================
template<unsigned DIM>
class RefineableBuoyantQCrouzeixRaviartElement:
public virtual RefineableQAdvectionDiffusionElement<DIM,3>,
public virtual RefineableQCrouzeixRaviartElement<DIM>
{

private:

 /// Pointer to a new physical variable, the Rayleigh number
 double* Ra_pt;

 /// The static default value of the Rayleigh number
 static double Default_Physical_Constant_Value;

public:
 /// Constructor: call the underlying constructors and
 /// initialise the pointer to the Rayleigh number to address the default
 /// value of 0.0
 RefineableBuoyantQCrouzeixRaviartElement() :
  RefineableQAdvectionDiffusionElement<DIM,3>(),
  RefineableQCrouzeixRaviartElement<DIM>()
  {
   Ra_pt = &Default_Physical_Constant_Value;
  }

 /// The required number of values stored at the nodes is
 /// the sum of the required values of the two single-physics elements. This
 /// step is generic for any composed element of this type.
 inline unsigned required_nvalue(const unsigned &n) const
  {return (RefineableQAdvectionDiffusionElement<DIM,3>::required_nvalue(n) +
           RefineableQCrouzeixRaviartElement<DIM>::required_nvalue(n));}

 /// Access function for the Rayleigh number (const version)
 const double &ra() const {return *Ra_pt;}

 /// Access function for the pointer to the Rayleigh number
 double* &ra_pt() {return Ra_pt;}


 /// Final override for disable ALE
 void disable_ALE()
  {
   //Disable ALE in both sets of equations
   RefineableNavierStokesEquations<DIM>::disable_ALE();
   RefineableAdvectionDiffusionEquations<DIM>::disable_ALE();
  }

 /// Final override for enable ALE
 void enable_ALE()
  {
   //Enable ALE in both sets of equations
   RefineableNavierStokesEquations<DIM>::enable_ALE();
   RefineableAdvectionDiffusionEquations<DIM>::enable_ALE();
  }



 /// Number of scalars/fields output by this element. Broken
 /// virtual. Needs to be implemented for each new specific element type.
 /// Temporary dummy
 unsigned nscalar_paraview() const
 {
  throw OomphLibError(
   "This function hasn't been implemented for this element",
   OOMPH_CURRENT_FUNCTION,
   OOMPH_EXCEPTION_LOCATION);
  
  // Dummy unsigned
  return 0;
 }
 
 /// Write values of the i-th scalar field at the plot points. Broken
 /// virtual. Needs to be implemented for each new specific element type.
 /// Temporary dummy
 void scalar_value_paraview(std::ofstream& file_out,
                            const unsigned& i,
                            const unsigned& nplot)  const
 {
  throw OomphLibError(
   "This function hasn't been implemented for this element",
   OOMPH_CURRENT_FUNCTION,
   OOMPH_EXCEPTION_LOCATION);
 }
 
 
 /// Name of the i-th scalar field. Default implementation
 /// returns V1 for the first one, V2 for the second etc. Can (should!) be
 /// overloaded with more meaningful names.
 std::string scalar_name_paraview(const unsigned& i) const
  {
   return "V"+StringConversion::to_string(i);
  }
 
 ///  Overload the standard output function with the broken default
 void output(ostream &outfile)
  {FiniteElement::output(outfile);}

 /// Output function:
 ///  x,y,u   or    x,y,z,u at Nplot^DIM plot points
 void output(ostream &outfile, const unsigned &nplot)
  {
   //vector of local coordinates
   Vector<double> s(DIM);

   // Tecplot header info
   outfile << this->tecplot_zone_string(nplot);

   // Loop over plot points
   unsigned num_plot_points=this->nplot_points(nplot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     // Get local coordinates of plot point
     this->get_s_plot(iplot,nplot,s);

     // Output the position of the plot point
     for(unsigned i=0;i<DIM;i++)
      {outfile << this->interpolated_x(s,i) << " ";}

     // Output the fluid velocities at the plot point
     for(unsigned i=0;i<DIM;i++)
      {outfile << this->interpolated_u_nst(s,i) << " ";}

     // Output the fluid pressure at the plot point
     outfile << this->interpolated_p_nst(s)  << " ";

     // Output the temperature at the plot point
     outfile << this->interpolated_u_adv_diff(s) << " " << std::endl;
    }
   outfile << std::endl;

   // Write tecplot footer (e.g. FE connectivity lists)
   this->write_tecplot_zone_footer(outfile,nplot);

  }

 /// C-style output function:  Broken default
 void output(FILE* file_pt)
  {FiniteElement::output(file_pt);}

 ///  C-style output function: Broken default
 void output(FILE* file_pt, const unsigned &n_plot)
  {FiniteElement::output(file_pt,n_plot);}

 /// Output function for an exact solution: Broken default
 void output_fct(ostream &outfile, const unsigned &Nplot,
                 FiniteElement::SteadyExactSolutionFctPt
                 exact_soln_pt)
  {FiniteElement::output_fct(outfile,Nplot,exact_soln_pt);}


 /// Output function for a time-dependent exact solution.
 /// Broken default
 void output_fct(ostream &outfile, const unsigned &Nplot,
                 const double& time,
                 FiniteElement::UnsteadyExactSolutionFctPt
                 exact_soln_pt)
  {
   FiniteElement::
    output_fct(outfile,Nplot,time,exact_soln_pt);
  }

 /// Overload the index at which the temperature
 /// variable is stored. We choose to store is after the fluid velocities.
 unsigned u_index_adv_diff() const
  {return DIM;}

 /// Number of vertex nodes in the element is obtained from the
 /// geometric element.
 unsigned nvertex_node() const
  {return QElement<DIM,3>::nvertex_node();}

 /// Pointer to the j-th vertex node in the element,
 /// Call the geometric element's function.
 Node* vertex_node_pt(const unsigned& j) const
  {return QElement<DIM,3>::vertex_node_pt(j);}

 /// The total number of continously interpolated values is
 /// DIM+1 (DIM fluid velocities and one temperature).
 unsigned ncont_interpolated_values() const
  {return DIM+1;}


 /// Get the continuously interpolated values at the local coordinate s.
 /// We choose to put the fluid velocities first, followed by the
 /// temperature.
 void get_interpolated_values(const Vector<double>&s,  Vector<double>& values)
  {
   //Storage for the fluid velocities
   Vector<double> nst_values;

   //Get the fluid velocities from the fluid element
   RefineableQCrouzeixRaviartElement<DIM>::
    get_interpolated_values(s,nst_values);

   //Storage for the temperature
   Vector<double> advection_values;

   //Get the temperature from the advection-diffusion element
   RefineableQAdvectionDiffusionElement<DIM,3>::
    get_interpolated_values(s,advection_values);

   //Add the fluid velocities to the values vector
   for(unsigned i=0;i<DIM;i++) {values.push_back(nst_values[i]);}

   //Add the concentration to the end
   values.push_back(advection_values[0]);
  }



 /// Get all continuously interpolated values at the local
 /// coordinate s at time level t (t=0: present; t>0: previous).
 /// We choose to put the fluid velocities first, followed by the
 /// temperature
 void get_interpolated_values(const unsigned& t, const Vector<double>&s,
                              Vector<double>& values)
  {
   //Storage for the fluid velocities
   Vector<double> nst_values;

   //Get the fluid velocities from the fluid element
   RefineableQCrouzeixRaviartElement<DIM>::
    get_interpolated_values(t,s,nst_values);

   //Storage for the temperature
   Vector<double> advection_values;

   //Get the temperature from the advection-diffusion element
   RefineableQAdvectionDiffusionElement<DIM,3>::
    get_interpolated_values(s,advection_values);

   //Add the fluid velocities to the values vector
   for(unsigned i=0;i<DIM;i++) {values.push_back(nst_values[i]);}

   //Add the concentration to the end
   values.push_back(advection_values[0]);

  } // end of get_interpolated_values



 /// The additional hanging node information must be set up
 /// for both single-physics elements.
 void further_setup_hanging_nodes()
  {
   RefineableQCrouzeixRaviartElement<DIM>::further_setup_hanging_nodes();
   RefineableQAdvectionDiffusionElement<DIM,3>::further_setup_hanging_nodes();
  }



 /// Call the rebuild_from_sons functions for each of the
 /// constituent multi-physics elements.
 void rebuild_from_sons(Mesh* &mesh_pt)
  {
   RefineableQAdvectionDiffusionElement<DIM,3>::rebuild_from_sons(mesh_pt);
   RefineableQCrouzeixRaviartElement<DIM>::rebuild_from_sons(mesh_pt);
  }



 /// Call the underlying single-physics element's further_build()
 /// functions and make sure that the pointer to the Rayleigh number
 /// is passed to the sons
 void further_build()
  {
   RefineableQCrouzeixRaviartElement<DIM>::further_build();
   RefineableQAdvectionDiffusionElement<DIM,3>::further_build();

   //Cast the pointer to the father element to the specific
   //element type
   RefineableBuoyantQCrouzeixRaviartElement<DIM>* cast_father_element_pt
    = dynamic_cast<RefineableBuoyantQCrouzeixRaviartElement<DIM>*>(
     this->father_element_pt());

   //Set the pointer to the Rayleigh number to be the same as that in
   //the father
   this->Ra_pt = cast_father_element_pt->ra_pt();
  } //end of further build



 /// The recovery order is that of the NavierStokes elements.
 unsigned nrecovery_order()
  {return RefineableQCrouzeixRaviartElement<DIM>::nrecovery_order();}

 /// The number of compound fluxes is two (one for the fluid and
 /// one for the temperature)
 unsigned ncompound_fluxes() {return 2;}

 /// The number of Z2 flux terms is the same as that in
 /// the fluid element plus that in the advection-diffusion element
 unsigned num_Z2_flux_terms()
  {
   return (RefineableQCrouzeixRaviartElement<DIM>::num_Z2_flux_terms() +
           RefineableQAdvectionDiffusionElement<DIM,3>::num_Z2_flux_terms());
  }



 /// Get the Z2 flux from the fluid element
 void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
  {
   //Find the number of fluid fluxes
   unsigned n_fluid_flux =
    RefineableQCrouzeixRaviartElement<DIM>::num_Z2_flux_terms();
   //Fill in the first flux entries as the velocity entries
   RefineableQCrouzeixRaviartElement<DIM>::get_Z2_flux(s,flux);

   //Find the number of temperature fluxes
   unsigned n_temp_flux =
    RefineableQAdvectionDiffusionElement<DIM,3>::num_Z2_flux_terms();
   Vector<double> temp_flux(n_temp_flux);
   //Get the temperature flux
   RefineableQAdvectionDiffusionElement<DIM,3>::get_Z2_flux(s,temp_flux);

   //Add the temperature flux to the end of the flux vector
   for(unsigned i=0;i<n_temp_flux;i++)
    {
     flux[n_fluid_flux+i] = temp_flux[i];
    }
  } //end of get_Z2_flux

 /// Fill in which flux components are associated with the fluid
 /// measure and which are associated with the temperature measure
 void get_Z2_compound_flux_indices(Vector<unsigned> &flux_index)
  {
   //Find the number of fluid fluxes
   unsigned n_fluid_flux =
    RefineableQCrouzeixRaviartElement<DIM>::num_Z2_flux_terms();
   //Find the number of temperature fluxes
   unsigned n_temp_flux =
    RefineableQAdvectionDiffusionElement<DIM,3>::num_Z2_flux_terms();

   //The fluid fluxes are first
   //The values of the flux_index vector are zero on entry, so we
   //could omit this line
   for(unsigned i=0;i<n_fluid_flux;i++) {flux_index[i] = 0;}
   //Set the temperature fluxes (the last set of fluxes
   for(unsigned i=0;i<n_temp_flux;i++) {flux_index[n_fluid_flux + i] = 1;}
  }


 /// Validate against exact solution at given time
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points and compute L2 error
 /// and L2 norm of velocity solution over element
 /// Overload to broken default
 void compute_error(ostream &outfile,
                    FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt,
                    const double& time,
                    double& error, double& norm)
  {FiniteElement::compute_error(outfile,exact_soln_pt,
                                time,error,norm);}

 /// Validate against exact solution.
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points and compute L2 error
 /// and L2 norm of velocity solution over element
 /// Overload to broken default.
void compute_error(ostream &outfile,
                    FiniteElement::SteadyExactSolutionFctPt exact_soln_pt,
                    double& error, double& norm)
  {FiniteElement::
   compute_error(outfile,exact_soln_pt,error,norm);}

 /// Overload the wind function in the advection-diffusion equations.
 /// This provides the coupling from the Navier--Stokes equations to the
 /// advection-diffusion equations because the wind is the fluid velocity.
 void get_wind_adv_diff(const unsigned& ipt,
                        const Vector<double> &s, const Vector<double>& x,
                        Vector<double>& wind) const
  {
   //The wind function is simply the velocity at the points
   this->interpolated_u_nst(s,wind);
  }


 /// Overload the body force in the navier-stokes equations
 /// This provides the coupling from the advection-diffusion equations
 /// to the Navier--Stokes equations, the body force is the
 /// temperature multiplied by the Rayleigh number acting in the
 /// direction opposite to gravity.
 void get_body_force_nst(const double& time, const unsigned& ipt,
                         const Vector<double> &s, const Vector<double> &x,
                         Vector<double> &result)
  {
   // Get vector that indicates the direction of gravity from
   // the Navier-Stokes equations
   Vector<double> gravity(NavierStokesEquations<DIM>::g());

   // Temperature-dependent body force:
   for (unsigned i=0;i<DIM;i++)
    {
     result[i] = -gravity[i]*this->interpolated_u_adv_diff(s)*ra();
    }
  }

 /// Fill in the constituent elements' contribution to the residual vector.
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the residuals of the Navier-Stokes equations
   RefineableNavierStokesEquations<DIM>::fill_in_contribution_to_residuals(
    residuals);

   //Call the residuals of the advection-diffusion equations
   RefineableAdvectionDiffusionEquations<DIM>::
    fill_in_contribution_to_residuals(residuals);
  }


 /// Compute the element's residual Vector and the jacobian matrix
 /// using full finite differences, the default implementation
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
#ifdef USE_FD_JACOBIAN_FOR_REFINEABLE_BUOYANT_Q_ELEMENT
   FiniteElement::fill_in_contribution_to_jacobian(residuals,jacobian);
#else
   //Calculate the Navier-Stokes contributions (diagonal block and residuals)
   RefineableNavierStokesEquations<DIM>::
    fill_in_contribution_to_jacobian(residuals,jacobian);

   //Calculate the advection-diffusion contributions
   //(diagonal block and residuals)
   RefineableAdvectionDiffusionEquations<DIM>::
    fill_in_contribution_to_jacobian(residuals,jacobian);

   //We now fill in the off-diagonal (interaction) blocks analytically
   this->fill_in_off_diagonal_jacobian_blocks_analytic(residuals,jacobian);
#endif
  } //End of jacobian calculation

 /// Add the element's contribution to its residuals vector,
 /// jacobian matrix and mass matrix
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals, DenseMatrix<double> &jacobian,
  DenseMatrix<double> &mass_matrix)
  {
   //Call the (broken) version in the base class
   FiniteElement::fill_in_contribution_to_jacobian_and_mass_matrix(
     residuals,jacobian,mass_matrix);
  }

 /// Compute the contribution of the off-diagonal blocks
 /// analytically.
 void fill_in_off_diagonal_jacobian_blocks_analytic(
  Vector<double> &residuals, DenseMatrix<double> &jacobian)
  {
   //Perform another loop over the integration loops using the information
   //from the original elements' residual assembly loops to determine
   //the conributions to the jacobian

   // Local storage for pointers to hang_info objects
   HangInfo *hang_info_pt=0, *hang_info2_pt=0;

   //Local storage for the index in the nodes at which the
   //Navier-Stokes velocities are stored (we know that this should be 0,1,2)
   unsigned u_nodal_nst[DIM];
   for(unsigned i=0;i<DIM;i++) {u_nodal_nst[i] = this->u_index_nst(i);}

   //Local storage for the  index at which the temperature is stored
   const unsigned u_nodal_adv_diff = this->u_index_adv_diff();

   //Find out how many nodes there are
   const unsigned n_node = this->nnode();

   //Set up memory for the shape and test functions and their derivatives
   Shape psif(n_node), testf(n_node);
   DShape dpsifdx(n_node,DIM), dtestfdx(n_node,DIM);

   //Number of integration points
   const unsigned n_intpt = this->integral_pt()->nweight();

   //Get Physical Variables from Element
   double Ra = this->ra();
   double Pe = this->pe();
   Vector<double> gravity = this->g();

   //Integers to store the local equations and unknowns
   int local_eqn=0, local_unknown=0;

   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Get the integral weight
     double w = this->integral_pt()->weight(ipt);

     //Call the derivatives of the shape and test functions
     double J =
      this->dshape_and_dtest_eulerian_at_knot_nst(ipt,psif,dpsifdx,
                                                  testf,dtestfdx);

     //Premultiply the weights and the Jacobian
     double W = w*J;

     //Calculate local values of temperature derivatives
     //Allocate
     Vector<double> interpolated_du_adv_diff_dx(DIM,0.0);

     // Loop over nodes
     for(unsigned l=0;l<n_node;l++)
      {
       //Get the nodal value
       double u_value = this->nodal_value(l,u_nodal_adv_diff);
       //Loop over the derivative directions
       for(unsigned j=0;j<DIM;j++)
        {
         interpolated_du_adv_diff_dx[j] += u_value*dpsifdx(l,j);
        }
      }

     //Assemble the Jacobian terms
     //---------------------------

     //Loop over the test functions/eqns
     for(unsigned l=0;l<n_node;l++)
      {
       //Local variables to store the number of master nodes and
       //the weight associated with the shape function if the node is hanging
       unsigned n_master=1;
       double hang_weight=1.0;

       //Local bool (is the node hanging)
       bool is_node_hanging = this->node_pt(l)->is_hanging();

       //If the node is hanging, get the number of master nodes
       if(is_node_hanging)
        {
         hang_info_pt = this->node_pt(l)->hanging_pt();
         n_master = hang_info_pt->nmaster();
        }
       //Otherwise there is just one master node, the node itself
       else
        {
         n_master = 1;
        }

       //Loop over the master nodes
       for(unsigned m=0;m<n_master;m++)
        {
         //If the node is hanging get weight from master node
         if(is_node_hanging)
          {
           //Get the hang weight from the master node
           hang_weight = hang_info_pt->master_weight(m);
          }
         else
          {
           // Node contributes with full weight
           hang_weight = 1.0;
          }


         //Assemble derivatives of Navier Stokes momentum w.r.t. temperature
         //-----------------------------------------------------------------

         // Loop over velocity components for equations
         for(unsigned i=0;i<DIM;i++)
          {

           //Get the equation number
           if(is_node_hanging)
            {
             //Get the equation number from the master node
             local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                              u_nodal_nst[i]);
            }
           else
            {
             // Local equation number
             local_eqn = this->nodal_local_eqn(l,u_nodal_nst[i]);
            }

           if(local_eqn >= 0)
            {
             //Local variables to store the number of master nodes
             //and the weights associated with each hanging node
             unsigned n_master2=1;
             double hang_weight2=1.0;

             //Loop over the nodes for the unknowns
             for(unsigned l2=0;l2<n_node;l2++)
              {
               //Local bool (is the node hanging)
               bool is_node2_hanging = this->node_pt(l2)->is_hanging();

               //If the node is hanging, get the number of master nodes
               if(is_node2_hanging)
                {
                 hang_info2_pt = this->node_pt(l2)->hanging_pt();
                 n_master2 = hang_info2_pt->nmaster();
                }
               //Otherwise there is one master node, the node itself
               else
                {
                 n_master2 = 1;
                }

               //Loop over the master nodes
               for(unsigned m2=0;m2<n_master2;m2++)
                {
                 if(is_node2_hanging)
                  {
                   //Read out the local unknown from the master node
                   local_unknown =
                    this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                         u_nodal_adv_diff);
                   //Read out the hanging weight from the master node
                   hang_weight2 = hang_info2_pt->master_weight(m2);
                  }
                 else
                  {
                   //The local unknown number comes from the node
                   local_unknown = this->nodal_local_eqn(l2,u_nodal_adv_diff);
                   //The hang weight is one
                   hang_weight2 = 1.0;
                  }

                 if(local_unknown >= 0)
                  {
                   //Add contribution to jacobian matrix
                   jacobian(local_eqn,local_unknown)
                    += -gravity[i]*psif(l2)*Ra*testf(l)*
                    W*hang_weight*hang_weight2;
                  }
                }
              }
            }
          }


         //Assemble derivative of adv diff eqn w.r.t. fluid veloc
         //------------------------------------------------------
         {
          //Get the equation number
          if(is_node_hanging)
           {
            //Get the equation number from the master node
            local_eqn = this->local_hang_eqn(hang_info_pt->master_node_pt(m),
                                             u_nodal_adv_diff);
           }
          else
           {
            // Local equation number
            local_eqn = this->nodal_local_eqn(l,u_nodal_adv_diff);
           }

          //If it's not pinned
          if(local_eqn >= 0)
           {
            //Local variables to store the number of master nodes
            //and the weights associated with each hanging node
            unsigned n_master2=1;
            double hang_weight2=1.0;

            //Loop over the nodes for the unknowns
            for(unsigned l2=0;l2<n_node;l2++)
             {
              //Local bool (is the node hanging)
              bool is_node2_hanging = this->node_pt(l2)->is_hanging();

              //If the node is hanging, get the number of master nodes
              if(is_node2_hanging)
               {
                hang_info2_pt = this->node_pt(l2)->hanging_pt();
                n_master2 = hang_info2_pt->nmaster();
               }
              //Otherwise there is one master node, the node itself
              else
               {
                n_master2 = 1;
               }

              //Loop over the master nodes
              for(unsigned m2=0;m2<n_master2;m2++)
               {
                //If the node is hanging
                if(is_node2_hanging)
                 {
                  //Read out the hanging weight from the master node
                  hang_weight2 = hang_info2_pt->master_weight(m2);
                 }
                //If the node is not hanging
                else
                 {
                  //The hang weight is one
                  hang_weight2 = 1.0;
                 }

                //Loop over the velocity degrees of freedom
                for(unsigned i2=0;i2<DIM;i2++)
                 {
                  //If the node is hanging
                  if(is_node2_hanging)
                   {
                    //Read out the local unknown from the master node
                    local_unknown =
                     this->local_hang_eqn(hang_info2_pt->master_node_pt(m2),
                                          u_nodal_nst[i2]);
                   }
                  else
                   {
                    //The local unknown number comes from the node
                    local_unknown=this->nodal_local_eqn(l2, u_nodal_nst[i2]);
                   }

                  //If it's not pinned
                  if(local_unknown >= 0)
                   {
                    //Add the contribution to the jacobian matrix
                    jacobian(local_eqn,local_unknown)
                     -= Pe*psif(l2)*interpolated_du_adv_diff_dx[i2]*testf(l)
                     *W*hang_weight*hang_weight2;
                   }
                 }
               }
             }
           }
         }

        }
      }
    }
  } //End of function

 ///
 /// Perhaps this deserves a comment rather than three slashes?
 /// I did not write this function...
 unsigned ndof_types() const
  {
   return DIM+2;
  }


 /// Create a list of pairs for all unknowns in this element,
 /// so that the first entry in each pair contains the global equation
 /// number of the unknown, while the second one contains the number
 /// of the DOF that this unknown is associated with.
 /// (Function can obviously only be called if the equation numbering
 /// scheme has been set up.)
 void get_dof_numbers_for_unknowns(
  std::list<std::pair<unsigned long,unsigned> >& dof_lookup_list) const
  {
   // number of nodes
   unsigned n_node = this->nnode();

   // number of pressure values
   unsigned n_press = this->npres_nst();

   // temporary pair (used to store dof lookup prior to being added to list)
   std::pair<unsigned,unsigned> dof_lookup;

   // pressure dof number
   unsigned pressure_dof_number = DIM;

   // loop over the pressure values
   for (unsigned n = 0; n < n_press; n++)
    {
     // determine local eqn number
     int local_eqn_number = this->p_local_eqn(n);

     // ignore pinned values - far away degrees of freedom resulting
     // from hanging nodes can be ignored since these are be dealt
     // with by the element containing their master nodes
     if (local_eqn_number >= 0)
      {
       // store dof lookup in temporary pair: First entry in pair
       // is global equation number; second entry is dof type
       dof_lookup.first = this->eqn_number(local_eqn_number);
       dof_lookup.second = pressure_dof_number;

       // add to list
       dof_lookup_list.push_front(dof_lookup);
      }
    }

   // loop over the nodes
   for (unsigned n = 0; n < n_node; n++)
    {
     //loop over velocity dofs
     for (unsigned v = 0; v < DIM; v++)
      {
       // determine local eqn number
       int local_eqn_number = this->nodal_local_eqn(n, v);

       // ignore pinned values
       if (local_eqn_number >= 0)
        {
         // store dof lookup in temporary pair: First entry in pair
         // is global equation number; second entry is dof type
         dof_lookup.first = this->eqn_number(local_eqn_number);
         dof_lookup.second = v;

         // add to list
         dof_lookup_list.push_front(dof_lookup);
        }
      }

     // Temperature is stored at the end
     unsigned v=DIM;
     {
      // determine local eqn number
      int local_eqn_number = this->nodal_local_eqn(n, v);

      // ignore pinned values
      if (local_eqn_number >= 0)
       {
        // store dof lookup in temporary pair: First entry in pair
        // is global equation number; second entry is dof type
        dof_lookup.first = this->eqn_number(local_eqn_number);

        // Temperture is the last type
        dof_lookup.second = DIM+1;

        // add to list
        dof_lookup_list.push_front(dof_lookup);

       }
     }
    }
  }



};


//===================================================================
/// Set the default value of the Rayleigh number to be zero
//===================================================================
template<>
double RefineableBuoyantQCrouzeixRaviartElement<3>::
Default_Physical_Constant_Value=0.0;



/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////



//======start_of_namespace============================================
/// Namespace for the physical parameters in the problem
//====================================================================
namespace Global_Physical_Variables
{

 /// Peclet number
 double Peclet=10.0;

 /// Reynolds number
 double Reynolds = 10.0;

 /// Rayleigh number
 double Rayleigh = 100.0;

 /// Gravity vector
 Vector<double> Direction_of_gravity(3);

} // end_of_namespace


/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////



//=======start_of_problem_class=======================================
/// 3D Convection  problem on rectangular domain, discretised
/// with refineable elements. The specific type
/// of element is specified via the template parameter.
//====================================================================
template<class ELEMENT>
class RefineableConvectionProblem : public Problem
{

public:

 /// Constructor
 RefineableConvectionProblem();

 /// Destructor. Clean up all allocated memory
 ~RefineableConvectionProblem()
  {
   //Delete the mesh's error estimator
   delete mesh_pt()->spatial_error_estimator_pt();
   //Delete the mesh
   delete Problem::mesh_pt();

   //Can we cast the solver to an iterative linear solver
   GMRES<CRDoubleMatrix>* iterative_linear_solver_pt =
    dynamic_cast<GMRES<CRDoubleMatrix>*>(this->linear_solver_pt());
   //If so delete the preconditioner and the solver
   if(iterative_linear_solver_pt)
    {
     delete iterative_linear_solver_pt->preconditioner_pt();
     delete iterative_linear_solver_pt;
    }
  }

 /// Update the problem specs before solve:
 void actions_before_newton_solve();

 /// Update the problem after solve (empty)
 void actions_after_newton_solve(){}

 /// Overloaded version of the problem's access function to
 /// the mesh. Recasts the pointer to the base Mesh object to
 /// the actual mesh type.
 RefineableSimpleCubicMesh<ELEMENT>* mesh_pt()
  {
   return dynamic_cast<RefineableSimpleCubicMesh<ELEMENT>*>(
    Problem::mesh_pt());
  } //end of access function to specic mesh

 /// Actions before adapt:(empty)
 void actions_before_adapt() {}

 /// Actions after adaptation (empty)
 void actions_after_adapt() {}

 /// Doc the solution.
 void doc_solution();

private:

 /// DocInfo object
 DocInfo Doc_info;

 //Storage for length
 double Lz;


}; // end of problem class


//=======start_of_constructor=============================================
/// Constructor for adaptive thermal convection problem
//========================================================================
template<class ELEMENT>
RefineableConvectionProblem<ELEMENT>::
RefineableConvectionProblem()
{
 // Set output directory
 Doc_info.set_directory("RESLT");

 // # of elements in x-direction
 unsigned n_x=4;

 // # of elements in y-direction
 unsigned n_y=4;

 // # of elements in the z-direction
 unsigned n_z=12;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;

 // Domain length in the z-direction
 double l_z = 3.0;

 //Set the internal value
 Lz = l_z;

 // Build the mesh
 RefineableSimpleCubicMesh<ELEMENT>* cast_mesh_pt =
  new RefineableSimpleCubicMesh<ELEMENT>(n_x,n_y,n_z,l_x,l_y,l_z);

 //Set the problem's mesh pointer
 Problem::mesh_pt() = cast_mesh_pt;

 //Build iterative linear solver
 GMRES<CRDoubleMatrix>* iterative_linear_solver_pt = new
  GMRES<CRDoubleMatrix>;

 // Set maximum number of iterations
 iterative_linear_solver_pt->max_iter() = 100;

 // Set tolerance
 iterative_linear_solver_pt->tolerance() = 1.0e-8;

 BoussinesqPreconditioner* prec_pt = new BoussinesqPreconditioner(this);

 //Set the preconditioner
 iterative_linear_solver_pt->preconditioner_pt() = prec_pt;

 //Set the linear solver
 linear_solver_pt() = iterative_linear_solver_pt;

 // Create/set error estimator
 cast_mesh_pt->spatial_error_estimator_pt()=new Z2ErrorEstimator;

 // Set error targets for adaptive refinement
 cast_mesh_pt->max_permitted_error()=0.5e-3;
 cast_mesh_pt->min_permitted_error()=0.5e-4;


 // Set the boundary conditions for this problem: All nodes are
 // free by default -- only need to pin the ones that have Dirichlet
 // conditions here

 //Loop over the boundaries
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   //Set the maximum index to be pinned (all values by default)
   unsigned val_max=4;

   //If we are on the inlet or side walls pin everything
   //If we are on the outlet pin x and y velocities
   if(ibound==5) {val_max = 2;}

   //If we are on the side walls pin the x-velocity only
   if((ibound==2) || (ibound==4)) {val_max=1;}

   //Loop over the number of nodes on the boundry
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     //Loop over the desired values stored at the nodes and pin
     for(unsigned j=0;j<val_max;j++)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(j);
      }
    }
  }

 // Complete the build of all elements so they are fully functional

 // Loop over the elements to set up element-specific
 // things that cannot be handled by the (argument-free!) ELEMENT
 // constructor.
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralsedElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   // Set the Peclet number
   el_pt->pe_pt() = &Global_Physical_Variables::Peclet;

   // Set the Reynolds number (1/Pr in our non-dimensionalisation)
   el_pt->re_pt() = &Global_Physical_Variables::Reynolds;

   // Set ReSt (also 1/Pr in our non-dimensionalisation)
   el_pt->re_st_pt() = &Global_Physical_Variables::Reynolds;

   // Set the Rayleigh number
   el_pt->ra_pt() = &Global_Physical_Variables::Rayleigh;

   //Set Gravity vector
   el_pt->g_pt() = &Global_Physical_Variables::Direction_of_gravity;
  }

 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << endl;

} // end of constructor


//===================start_actions_before_newton_solve===========================
/// Update the problem specs before solve: (Re-)set boundary conditions
//========================================================================
template<class ELEMENT>
void RefineableConvectionProblem<ELEMENT>::actions_before_newton_solve()
{
 //Set the halfwith of the jet
 double x_halfwidth = 0.5;
 double y_halfwidth = 0.1875;
 //Set the peak flow of the jet
 double u_jet = 1.0;
 //Set the coordinates of the centre of the jet
 double x_jet = 0.5;  double y_jet = 0.375;
 //Bluntness of jet
 int alpha = 10;


 // Loop over the boundaries
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary
   unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   for(unsigned inod=0;inod<num_nod;inod++)
    {
     // Get pointer to node
     Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);

     //If we are on the inlet
     if(ibound==0)
      {
       //No x-velocity
       nod_pt->set_value(0,0.0);
       //No y-velocity
       nod_pt->set_value(1,0.0);

       //Default z-velocity
       double zvel=0.0;

       //get x and y values
       double x = nod_pt->x(0);
       double y = nod_pt->x(1);

       //Jetlike z-velocity in square/rectangle
       if((std::abs(x-x_jet) <= x_halfwidth) &&
          (std::abs(y-y_jet) <= y_halfwidth))
        {
         zvel =   u_jet*(1.0 - pow(std::abs(x-x_jet)/x_halfwidth,alpha))*
          (1.0 - pow(std::abs(y-y_jet)/y_halfwidth,alpha));
        }

       //Set the velocity
       nod_pt->set_value(2,zvel);

       //Set the temperature (hot)
       nod_pt->set_value(3,0.5);
      }

     //If we are on the upper and lower walls
     if((ibound==1) || (ibound==3))
      {
       //No slip on all velocity components
       for(unsigned j=0;j<3;j++) {nod_pt->set_value(j,0.0);}
       //The temperature has a discontinuity halfway along
       if(nod_pt->x(2) < 0.5*Lz)
        {
         nod_pt->set_value(3,0.5);
        }
       else
        {
         nod_pt->set_value(3,-0.5);
        }
      }

     //If we are on the side walls only pin the x-velocity
     if((ibound==2) || (ibound==4))
      {
       //No penetration on x-component only
       nod_pt->set_value(0,0.0);
      }

     //If we are on the end then set transverse velocities to zero
     if(ibound==5)
      {
       //No slip on transverse velocity components
       for(unsigned j=0;j<2;j++) {nod_pt->set_value(j,0.0);}
      }
    }
  }

}  // end of actions before solve



//====================start_of_doc_solution===============================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void RefineableConvectionProblem<ELEMENT>::doc_solution()
{
 //Declare an output stream and filename
 ofstream some_file;
 char filename[100];

 // Number of plot points: npts x npts
 unsigned npts=5;

 // Output solution
 //-----------------
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();

 Doc_info.number()++;
} // end of doc


//===============start_of_main========================================
/// Driver code for 3D Boussinesq convection problem with
/// adaptivity.
//====================================================================
int main(int argc, char **argv)
{
 //Uncomment for a parallel version
//#ifdef OOMPH_HAS_MPI
// // Set up MPI_Helpers
// MPI_Helpers::init(argc,argv);
//#endif

 // Set the direction of gravity
 Global_Physical_Variables::Direction_of_gravity[0] = 0.0;
 Global_Physical_Variables::Direction_of_gravity[1] = -1.0;
 Global_Physical_Variables::Direction_of_gravity[2] = 0.0;

 // Create the problem with 3D twenty-seven-node refineable elements.
 RefineableConvectionProblem<
  RefineableBuoyantQCrouzeixRaviartElement<3> > problem;

 //Solve the problem with (up to) two levels of adaptation
 problem.newton_solve();

 //Document the solution
 problem.doc_solution();

//#ifdef OOMPH_HAS_MPI
// // finalize MPI
// MPI_Helpers::finalize();
//#endif
} // end of main
