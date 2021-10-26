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
// Driver for space-time pulsatile Poiseuille flow

// Extra files necessary for std::strcpy on other machines
#include <cstring>
#include <limits>

// Generic routines
#include "generic.h"

// The Navier-Stokes equations
#include "navier_stokes.h"

// The space-time Navier-Stokes equations
#include "space_time_navier_stokes.h"

// The oscillating cylinder 2D spatial mesh
#include "meshes/rectangle_with_moving_cylinder_mesh.h"
#include "meshes/extruded_cube_mesh_from_quad_mesh_with_macro_elements.h"

// Add in the block preconditioning machinery
#include "space_time_block_preconditioner.h"

// Include oomph namespace
using namespace oomph;

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

// The specialisation of the PMLLayerElement and FaceGeometry element
// has to be inside the oomph namespace
namespace oomph
{
  //======start_of_MyRefineableQTaylorHoodElement=============================
  /// Overloaded element that allows projection of use as PML element
  //==========================================================================
  class MyRefineableQTaylorHoodElement :
    public virtual PMLElementBase<2>,
    public virtual RefineableQTaylorHoodElement<2>
  {
  public:
    /// Constructor
    MyRefineableQTaylorHoodElement() {}

    /// Pure virtual function in which we specify the
    /// values to be pinned (and set to zero) on the outer edge of
    /// the "pml" layer. None since we're not using this as pml functionality
    void values_to_be_pinned_on_outer_pml_boundary(Vector<unsigned>& values_to_pin)
    {}
  };

  //======start_of_PMLLayerElement============================================
  /// Policy class defining the elements to be used in the PML layers. Same!
  //==========================================================================
  template<>
  class PMLLayerElement<MyRefineableQTaylorHoodElement> :
    public virtual MyRefineableQTaylorHoodElement
  {
  public:
    /// Constructor: Call the constructor for the
    /// appropriate QElement
    PMLLayerElement() : MyRefineableQTaylorHoodElement()
    {}
  };

  //======start_of_FaceGeometry===============================================
  /// Face geometry for element is the same as that for the underlying
  /// wrapped element
  //==========================================================================
  template<>
  class FaceGeometry<MyRefineableQTaylorHoodElement> :
    public virtual QElement<1,3>
  {
  public:
    /// Constructor: Call the constructor for the 1D quadratic element
    FaceGeometry() : QElement<1,3>() {}
  };
} // End of namespace oomph

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//======start_of_OscillatingCylinder_class============================
/// Oscillating cylinder class
//====================================================================
class OscillatingCylinder : public GeomObject
{
public:
  /// Constructor: Pass in the radius, the amplitude of the cylinder
  /// motion, the simulation Strouhal number and a pointer to time object.
  OscillatingCylinder(double* radius_pt,
                      double* amplitude_pt,
                      Time* time_pt) :
    GeomObject(1,2),
    Radius_pt(radius_pt),
    Amplitude_pt(amplitude_pt),
    Time_pt(time_pt)
  {}

  /// Destructor: Empty
  virtual ~OscillatingCylinder() {}

  /// Access function for the amplitude (lvalue)
  double& amplitude()
  {
    // Return the value of the amplitude
    return *Amplitude_pt;
  } // End of amplitude

  /// Access function for the Time pointer
  Time* time_pt()
  {
    // Return the Time pointer
    return Time_pt;
  } // End of time_pt

  /// Current position vector to material point at Lagrangian
  /// coordinate xi (steady version)
  void position(const Vector<double>& xi,
                Vector<double>& r) const
  {
    // X-coordinate
    r[0]=(*Radius_pt)*cos(xi[0]);

    // Y-coordinate
    r[1]=(*Radius_pt)*sin(xi[0]);
  } // End of position

  /// Current position vector to material point at Lagrangian
  /// coordinate xi (unsteady version). Implementation includes a
  /// transition phase where the cylinder oscillates to a smaller
  /// amplitude than the target value. Used to ensure that the solution
  /// isn't drastically different to that at the next time step. This
  /// can be disabled by setting Use_transition_phase to false.
  void position(const unsigned& t,
                const Vector<double>& xi,
                Vector<double>& r) const
  {
    // We can only use this if t=0
    if (t!=0)
    {
      // Create an output stream
      std::ostringstream error_message_stream;

      // Create an error message
      error_message_stream << "Trying to access the position at time history "
                           << "value, t: " << t << ".\nThis doesn't make sense "
                           << "in the space-time mesh!" << std::endl;

      // Throw the error message
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // Calculate the coordinate before translation
    position(xi,r);

    // Get current time
    double time=Time_pt->time(t);

    // Scaling factor
    double arg=2.0*MathematicalConstants::Pi;

    // Calculate the translation
    double translation=(*Amplitude_pt)*sin(arg*time);

    // Update the y-coordinate
    r[1]+=translation;
  } // End of position

  /// Parametrised position on object: r(zeta). Evaluated at
  /// the continuous time value, t.
  virtual void position(const double& t,
                        const Vector<double>& xi,
                        Vector<double>& r) const
  {
    // Calculate the coordinate before translation
    position(xi,r);

    // Scaling factor
    double arg=2.0*MathematicalConstants::Pi;

    // Calculate the translation
    double translation=(*Amplitude_pt)*sin(arg*t);

    // Update the y-coordinate
    r[1]+=translation;
  } // End of position

  /// Velocity at any given point on the rigid cylinder at time, t
  virtual void velocity(const double& t,Vector<double>& u) const
  {
    // Scaling factor
    double arg=2.0*MathematicalConstants::Pi;

    // Zero velocity component in the x-direction
    u[0]=0.0;

    // Calculate the (non-zero) velocity component in the y-direction
    u[1]=arg*(*Amplitude_pt)*cos(arg*t);
  } // End of velocity

private:

  /// Radius of the cylinder
  double* Radius_pt;

  /// Non-dimensionalised amplitude of the cylinder motion
  double* Amplitude_pt;

  /// Pointer to the current time in the problem
  Time* Time_pt;
}; // End of OscillatingCylinder class

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//======start_of_NodeReordering_namespace=============================
/// Contains helper function to reorganise nodes
//====================================================================
namespace NodeReordering
{
  /// Function for ordering nodes. Return true if first node's position
  /// is "before" second nodes. Dimension 0 checked first, then... until they
  /// are different (by more than tol=1e-10). If they are both in exactly
  /// the same place an error is thrown.
  inline bool node_global_position_comparison(Node* nd1_pt, Node* nd2_pt)
  {
    // If we're comparing a node with itself
    if (nd1_pt==nd2_pt)
    {
      // Don't do anything (passing true makes std::sort break...)
      return false;
    }

    // Get the number of dimensions stored by the node
    unsigned n_dim=nd1_pt->ndim();

    // A vector containing the indices in the order to check them. Order by
    // time slices first, then order the nodes in each time-slice by their x
    // position and finally order the nodes by their y-position (3D version)
    unsigned chosen_index[]= {2,0,1};

    // The coordinate indices in the order to check them (2D version)
    unsigned chosen_index_2d[]= {0,1};

    // Allocate space for the index we're going to order by
    unsigned j=0;

    // Loop over the spatial coordinates
    for (unsigned i=0; i<n_dim; i++)
    {
      // If we're in 3D (space-time)
      if (n_dim==3)
      {
        // The index we're going to order by
        j=chosen_index[i];
      }
      else
      {
        // The index we're going to order by
        j=chosen_index_2d[i];
      }

      // Check to see if the points occupy the same position (in index j)
      if (std::abs(nd1_pt->x(j)-nd2_pt->x(j))>1e-10)
      {
        // Check to see if node 1 is before node 2
        if (nd1_pt->x(j)<nd2_pt->x(j))
        {
          // Node 1 is before node 2 so return true
          return true;
        }
        // If node 2 is after node 1
        else
        {
          // Node 1 isn't before node 2 so return false
          return false;
        }
      } // if (std::abs(nd1_pt->x(j)-nd2_pt->x(j))>1e-10)
    } // for (unsigned i=0;i<n_dim;i++)

    // Create an output stream
    std::ostringstream error_message_stream;

    // Construct the error message
    error_message_stream << "Nodes are at the same point to ~1e-10! "
                         << "The difference is "
                         << std::abs(nd1_pt->x(j)-nd2_pt->x(j)) << std::endl;

    // Throw the error message
    throw OomphLibError(error_message_stream.str(),
                        OOMPH_EXCEPTION_LOCATION,
                        OOMPH_CURRENT_FUNCTION);
  } // End of node_global_position_comparison


  /// Get a vector of the nodes in the order in which they are
  /// encountered when stepping through the elements (similar to
  /// reorder_nodes() but without changing the mesh's node vector).
  void get_node_reordering(Mesh* mesh_pt,
                           Vector<Node*>& reordering,
                           const bool& use_old_ordering)
  {
    // If the user wants to use the original ordering
    if (use_old_ordering)
    {
      // Setup map to check if nodes have been done yet
      std::map<Node*,bool> done;

      // Loop over all nodes
      unsigned nnod=mesh_pt->nnode();

      // Initialise the vector
      reordering.assign(nnod,0);

      // Return immediately if there are no nodes: Note assumption:
      // Either all the elements' nodes stored here or none. If only a subset
      // is stored in the Node_pt vector we'll get a range checking error below
      // (only if run with range checking, of course).
      if (nnod==0)
      {
        // Return immediately
        return;
      }

      // Loop over the nodes in the mesh
      for (unsigned j=0; j<nnod; j++)
      {
        // Indicate whether or not the node has been swapped
        done[mesh_pt->node_pt(j)]=false;
      }

      // Initialise counter for number of nodes
      unsigned long count=0;

      // Get the number of elements in the mesh
      unsigned nel=mesh_pt->nelement();

      // Loop over all elements
      for (unsigned e=0; e<nel; e++)
      {
        // Upcase FiniteElement (or derived) class object
        FiniteElement* el_pt=mesh_pt->finite_element_pt(e);

        // If the upcast was successful
        if (el_pt!=0)
        {
          // Get the number of nodes in this element
          unsigned nnod=el_pt->nnode();

          // Loop over nodes in element
          for (unsigned j=0; j<nnod; j++)
          {
            // Get a pointer to the j-th node in the element
            Node* nod_pt=el_pt->node_pt(j);

            // Has node been done yet?
            if (!done[nod_pt])
            {
              // Insert into node vector. NOTE: If you get a seg fault/range
              // checking error here then you probably haven't added all the
              // elements' nodes to the Node_pt vector -- this is most likely
              // to arise in the case of meshes of face elements (though they
              // usually don't store the nodes at all so if you have any
              // problems here there's something unusual/not quite right in
              // any case... For this reason we don't range check here by
              // default (not even under paranoia) but force you turn on proper
              // (costly) range checking to track this down...
              reordering[count]=nod_pt;

              // Indicate that the node has been done
              done[nod_pt]=true;

              // Increase counter
              count++;
            }
          } // for (unsigned j=0;j<nnod;j++)
        } // if (el_pt!=0)
      } // for (unsigned e=0;e<nel;e++)

      // Sanity check
      if (count!=nnod)
      {
        // Create an error message
        std::string error_message="Trouble: Number of nodes hasn't stayed ";

        // Finish off the message
        error_message+="constant during reordering!\n";

        // Throw an error
        throw OomphLibError(error_message,
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    }
    else
    {
      // Copy node vector out
      unsigned n_node=mesh_pt->nnode();

      // Resize the node ordering vector
      reordering.resize(n_node);

      // Loop over the nodes
      for (unsigned i=0; i<n_node; i++)
      {
        // Assign the i-th node pointer entry
        reordering[i]=mesh_pt->node_pt(i);
      }

      // Now sort the nodes lexicographically
      std::sort(reordering.begin(),reordering.end(),
                &node_global_position_comparison);
    } // if (use_old_ordering)
  } // End of get_node_reordering


  /// Reorder nodes in the order in which they are encountered when
  /// stepping through the elements
  void reorder_nodes(Mesh* mesh_pt,const bool& use_old_ordering)
  {
    // Create storage for the reordered nodes
    Vector<Node*> reordering;

    // Get the reordered nodes (without altering the mesh's node vector)
    get_node_reordering(mesh_pt,reordering,use_old_ordering);

    // Get the number of nodes in the mesh
    unsigned n_node=mesh_pt->nnode();

    // Loop over all of the nodes
    for (unsigned i=0; i<n_node; i++)
    {
      // Replace the Mesh's i-th node pointer with the reordered node pointer
      mesh_pt->node_pt(i)=reordering[i];
    }
  } // End of reorder_nodes
} // End of NodeReordering

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//======start_of_GlobalParameters_namespace===========================
/// Global parameters for the problem
//====================================================================
namespace GlobalParameters
{
  /// --------------------------Cylinder Motion----------------------------
  /// Amplitude of the cylinder motion used by Williamson & Roshko
  /// (1988). As a side note, since the (simulation) cylinder has unit
  /// diameter (i.e. D=1) this is actually the dimensionless amplitude.
  double Amplitude=0.25;

  /// The target amplitude; used if we're going to do a parameter
  /// sweep through the amplitude-wavelength plane. This is generally the
  /// second parameter for which parameter continuation is used (if, of
  /// course, the value of Amplitude is different to Amplitude_target).
  double Amplitude_target=0.50;

  /// The number of steps used to reach the target amplitude
  unsigned N_amplitude_step=5;

  /// The ratio between the cylinder excitation period, T_e, and
  /// the stationary cylinder vortex-shedding period, T_s. Explicitly,
  /// we have, Period_ratio=T_e/T_s. If the value of this parameter ever
  /// changes then the function
  ///                 update_physical_parameter(),
  /// (defined in this namespace) MUST be called immediately after.
  ///
  /// NOTE: We use the ratio T_e/T_s (instead of T_s/T_e; used by
  /// Leontini et al. 2006) to match the x-axis in the experimental
  /// phase diagram of Williamson & Roshko (1988) (with the y-axis given
  /// by the non-dimensionalised amplitude parameter).
  double Period_ratio=1.0;

  /// The target Period_ratio value; used if we're going to do a
  /// parameter sweep through the amplitude-wavelength plane. This will
  /// normally be the last parameter for which parameter continuation is
  /// used (again, if the value of Period_ratio is different to
  /// Period_ratio_target).
  double Period_ratio_target=1.0;

  /// The number of steps used to reach the target Period_ratio value
  unsigned N_period_ratio_step=1;
  /// --------------------------Cylinder Motion----------------------------


  /// ---------------------Navier-Stokes Parameters------------------------
  /// Set the (current) Reynolds number. A pointer to this variable
  /// is provided to elements to make them fully functional. As this is
  /// used to calculate the Womersley number (=Re*St), the function
  ///                 update_physical_parameter(),
  /// (defined in this namespace) MUST be called immediately after editing
  /// the value of this variable.
  double Re=2.0;

  /// The target Reynolds number used in simulations. If this is
  /// different to the current Reynolds number then the first part of this
  /// simulation will (or at least, should) work to reach the target
  /// Reynolds number. Also, the value this number takes dictates the
  /// choice of Strouhal number (=St).
  double Re_target=10.0;

  /// The number of steps used to reach the target Reynolds number.
  /// In general the target Reynolds number will be reached using natural
  /// parameter continuation as pseudo-arc-length continuation doesn't
  /// seem to be a reliable method for it (or at least that's what was
  /// observed in coarse-grid simulations).
  unsigned N_re_step=2;

  /// The Strouhal number, St, is normally defined as
  ///                             St=L/UT,
  /// where L, U and T are the defining length-scale, flow speed and
  /// defining time scale. Physically, this parameter represents the
  /// ratio of convection time-scale to that of our artificial time-
  /// scale. The natural length-scale to use here is the cylinder
  /// diameter and U is defined to be the speed of the flow at the
  /// inlet. To capture time-periodic solutions (which have the same
  /// period as a single cylinder oscillation) we non-dimensionalise
  /// time on the cylinder oscillation time-scale, i.e. T_{e}. This
  /// allows us to use a mesh with fixed (unit) length in the time-
  /// direction; avoiding any re-meshing. Thus, we define
  ///                       St=St_{oomph}=D/UT_{e}.
  /// However, this is not the only time-scale of the flow. There is
  /// also the time-scale of the vortex-shedding past a stationary
  /// cylinder. The Strouhal number associated with this (which is
  /// dependent on the choice of Reynolds number) is
  ///                    St(Re)=St_{vort}(Re)=D/UT_{s},
  /// where T_{s} is the period of vortex-shedding. This is also the
  /// Strouhal number computed by:
  ///                 calculate_strouhal_number(Re),
  /// which is defined further below.
  /// It then follows:
  ///               St_{oomph}=(D/UT_{e})
  ///                         =(D/UT_{s})*(T_{s}/T_{e})
  ///                         =St_{vort}/Period_ratio.
  /// In general, we always choose St_{vort} to be that associated with
  /// the value of Re_target (which will remain the same throughout a
  /// simulation). As such, the only time the value of
  ///                St=St_{oomph}=St_{vort}/Period_ratio,
  /// may change is when the value of Period_ratio changes.
  ///
  /// NOTE: This value is also used to scale the cylinder velocity since,
  /// in dimensional units, nodes on the cylinder boundary satisfy:
  ///                u_{node}^{*}=dr^{*}_{cyl}/dt^{*},
  /// where the asterisk indicates the quantity is in dimensional units.
  /// After non-dimensionalising and rearranging we find
  ///            u_{node}=(D/UT_{e})*dr_{cyl}/dt
  ///                    =St_{oomph}*dr_{cyl}/dt
  ///                    =(St_{vort}/Period_ratio)*dr_{cyl}/dt.
  double St=1.0;

  /// The Womersley number (=Re*St), otherwise denoted as ReSt,
  /// is dependent on the value of the Reynolds number and the Strouhal
  /// number. If either value is changed then the function
  ///                 update_physical_parameter(),
  /// will be called which, in turn, updates the Womersley number value.
  double ReSt=Re*St;

  /// Function to calculate the Strouhal number appropriate for
  /// this simulation. The value chosen corresponds to the Strouhal
  /// number associated with the vortex-shedding frequency of the flow
  /// past a stationary cylinder.
  ///
  /// NOTE (1): The Re-St values for 46<Re<180 can be found in:
  ///         Williamson, C.H.K, (1988)."Defining a universal and
  ///         continuous Strouhal–Reynolds number relationship for
  ///         the laminar vortex shedding of a circular cylinder".
  /// Noteworthy Re-St value(s):
  ///                 (i) Re=100 --> St=0.1643.
  ///
  /// NOTE (2): This should only be used for Reynolds numbers above 46
  /// (where the Hopf bifurcation occurs for the flow past a stationary
  /// cylinder) and below (roughly) 180 (where a secondary bifurcation
  /// occurs at which point the flow becomes 3D).
  ///
  /// NOTE (3): This should only ever be called through the function
  ///                 update_physical_parameter();
  /// to compute the Strouhal number. It should never be used on its
  /// own (as it will only be appropriate for a unit Period_ratio value).
  double calculate_strouhal_number(const double& re)
  {
    // The min. Reynolds number
    double min_re=46.0;

    // The max. Reynolds number
    double max_re=180.0;

    // The first coefficient of the Re-St polynomial
    double a=-3.3265;

    // The second coefficient of the Re-St polynomial
    double b=0.1816;

    // The third coefficient of the Re-St polynomial
    double c=0.00016;

    // If we're above the maximum Reynolds number
    if (re>max_re)
    {
      // Throw an error
      throw OomphLibError("Don't know what to do for this Reynolds number!",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // If we're below the minimum Reynolds number
    else if (re<min_re)
    {
      // Just return the Strouhal value at the minimum Reynolds number
      return a/min_re+b+c*min_re;
    }
    // Otherwise, use the relationship in the Williamson paper
    else
    {
      // Return the Strouhal value at this Reynolds number
      return a/re+b+c*re;
    }
  } // End of calculate_strouhal_number

  /// Update physical parameters. This updates:
  ///               (1) The Reynolds number, and;
  ///               (2) The Strouhal number,
  /// and should ALWAYS be called after either the Reynolds number or
  /// Period_ratio value changes.
  void update_physical_parameters()
  {
    // Update the Strouhal number
    St=calculate_strouhal_number(Re)/Period_ratio;

    // Update the Womersley number
    ReSt=Re*St;
  } // End of update_physical_parameters
  /// ---------------------Navier-Stokes Parameters------------------------

  /// -------------------------Domain Properties---------------------------
  /// Radius of the cylinder
  double Radius=0.5;

  /// Pointer to the cylinder
  OscillatingCylinder* Cylinder_pt=0;

  /// The radius of the annular region surrounding the cylinder
  double Annular_region_radius=1.0;

  /// Height of domain
  double Height=20.0;

  /// X-coordinate of upstream end of domain
  double X_left=-10.0;

  /// X-coordinate of downstream end of domain
  double X_right=40.0;

  /// Length of square central box domain
  double Length_of_central_box=10.0;

  /// Number of uniform refinements before the mesh extrusion
  unsigned N_uniform_refinement_before_solve=1;

  /// The length of the mesh in the time direction
  double L_t=1.0;

  /// The number of elements in the time direction
  unsigned N_t=25;

  /// Update mesh parameters. This is (and only needs to be)
  /// once per simulation, during the setup of the mesh. This decides
  /// how thick the "fine resolution" layer of elements around the
  /// cylinder is. This layer is used to accurately resolve the boundary
  /// layer close to the cylinder.
  void update_mesh_parameters()
  {
    // Sanity check: N_t has to be odd!
    if ((N_t%2)==0)
    {
      // Throw an error
      OomphLibWarning("Method requires an odd number of time slices!\n",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
    }

    // Update the radius of the annular region from the updated parameter
    // values. NOTE: The annular rings are used to resolve the boundary
    // layers so they should not be made too large (hence the use of std::min)
    Annular_region_radius=
      Radius+std::min(2.0*Radius,0.5*((0.5*Length_of_central_box)-Radius));
  } // End of update_mesh_parameters
  /// -------------------------Domain Properties---------------------------


  /// ----------------------------Solver Info------------------------------
  /// Variable to choose which preconditioner to use. The actual
  /// preconditioner we choose to use is defined by the enumeration class
  /// implemented in the problem
  unsigned Preconditioner=0;

  /// Storage for the number of dof types in the mesh. Will be
  /// assigned in the function assign_time_slice_id()
  unsigned N_dof_type=0;

  /// Helper function which sets up the mapping between DOF types
  /// and which block they should be assigned to. This relies on the concept
  /// of "time slices" in the space-time formulation. All dofs in a given
  /// time slice will be aggregrated together
  void set_up_dof_to_block_mapping(Vector<unsigned>& dof_to_block_map)
  {
    // Resize the vector
    dof_to_block_map.resize(N_dof_type);

    // Loop over the dofs
    for (unsigned i=0; i<N_dof_type; i++)
    {
      // How many unique dof types are there per element? i.e. The velocities
      // and pressure in the first time slice constitute the first 3 then the
      // velocities in the middle time slice (in the element noting that the
      // NS elements use quadratic interpolation in time). The velocities and
      // pressure in the final time slice of the element are not included
      // because they are stored as the first 3 dof types in the next element
      // (in the time direction).
      unsigned n_unique_dof_per_element=3;

      // How many unique dof types are there per element after we've aggregated
      // the velocity components together? i.e. the velocities and pressure in
      // the first time slice in the element constitute one (aggregated) dof
      unsigned n_unique_aggregated_dof_per_element=1;

      // What (local) elemental dof does this (global) dof correspond to?
      unsigned i_local_dof=i%n_unique_dof_per_element;

      // Which elemental time slice does this dof correspond to?
      unsigned i_temporal=(i-i_local_dof)/n_unique_dof_per_element;

      // The first time slice in the element (u,v and p)
      if ((i_local_dof==0)||(i_local_dof==1)||(i_local_dof==2))
      {
        // Calculate the i-th entry
        dof_to_block_map[i]=i_temporal*n_unique_aggregated_dof_per_element;
      }
      else
      {
        // Create an output stream
        std::ostringstream error_message_stream;

        // Create an error message
        error_message_stream << "There should only be 3 unique dofs per element. "
                             << "Instead, you have " << n_unique_dof_per_element
                             << " unique dofs per element." << std::endl;

        // Throw the error message
        throw OomphLibError(error_message_stream.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    } // for (unsigned i=0;i<n_dof_types;i++)
  } // End of set_up_dof_to_block_mapping
  /// ----------------------------Solver Info------------------------------


  /// -----------------------Documentation Helpers-------------------------
  // DocInfo object for documentation
  DocInfo Doc_info;

  /// Trace file to doc. asymmetry norm data in
  std::ofstream Trace_file;

  /// Helper variable to indicate whether or not to document the solution
  bool Document_solution=true;

  /// Number of plot points (in each direction)
  unsigned N_plot_point=2;

  /// Document the maximum deformation inside the central box
  void doc_maximum_central_box_deformation()
  {
    // Calculate the distance from the edge of the annular ring to the
    // box boundary. NOTE: We check from the annular ring because the
    // region between the cylinder and annular ring is made rigid so
    // no compression occurs there.
    double compression_region_width=(Length_of_central_box/2.0-
                                     Annular_region_radius);

    // Calculate the current compression of the inner box
    double compression_ratio=((compression_region_width-Amplitude)/
                              compression_region_width);

    // If the compression is too large for the mesh
    if (compression_ratio<0.0)
    {
      // Used to create an error message
      std::ostringstream error_message_stream;

      // Create an error message
      error_message_stream << "The cylinder amplitude exceeds the size of "
                           << "the central box! Make the box larger!"
                           << std::endl;

      // Throw an error to the user
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    // If the compression is large then warn the user
    else if (compression_ratio<0.5)
    {
      // Used to create a warning message
      std::ostringstream warning_message_stream;

      // Create a warning message
      warning_message_stream << "Maximal mesh compression results in elements "
                             << "being reduced\nto "
                             << compression_ratio*100.0
                             << "% of their original width. It is therefore\n"
                             << "recommended the central box be made larger."
                             << std::endl;

      // Throw a warning to the user
      OomphLibWarning(warning_message_stream.str(),
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      // Document the compression (special ASCII characters used to output
      // the text in bold red)
      oomph_info << "\033[1;31m"
                 << "\nMaximum element compression ratio inside central box: "
                 << "\033[0m"
                 << compression_ratio*100.0 << "%" << std::endl;
    }
  } // End of doc_maximum_central_box_deformation

  /// Helper function that takes a pointer to one of the problem
  /// parameters and returns a string to denote it. Used in the (generic)
  /// parameter continuation functions to tell the user which parameter
  /// we're changing and what value it has at each iteration
  std::string parameter_to_string(const double* const parameter_pt)
  {
    // Are we dealing with the Reynolds number?
    if (parameter_pt==&Re)
    {
      // Return a string to denote the Reynolds number
      return "Re";
    }
    else if (parameter_pt==&Amplitude)
    {
      // Return a string to denote the amplitude
      return "A";
    }
    else if (parameter_pt==&Period_ratio)
    {
      // Return a string to denote the period ratio
      return "T_e/T_s";
    }
    else
    {
      // Create an output stream
      std::ostringstream warning_message_stream;

      // Create an error message
      warning_message_stream << "Don't know what to denote this parameter with "
                             << "so I'm just going to return an empty string..."
                             << std::endl;

      // Provide a warning
      OomphLibWarning(warning_message_stream.str(),
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);

      // Return an empty string (this way is faster than simply returning "")
      return std::string();
    } // if (parameter_pt==&Re)
  } // End of parameter_to_string
  /// -----------------------Documentation Helpers-------------------------


  /// ---------------------------Miscellaneous-----------------------------
  /// Find a node on the centerline
  /// N.B. We are modifying the *pointer* el_centerline_pt not the actual
  /// data. If we just pass a pointer to the element in then (from outside)
  /// we are only given a copy of the pointer which is discarded after the
  /// function call (but we need it to stay alive!). To edit the pointer
  /// itself we have to pass a pointer to the pointer (hence the **).
  void find_node_on_centerline(Mesh* mesh_pt,
                               FiniteElement** el_centerline_pt,
                               unsigned& node_index)
  {
    // Number of elements in the mesh
    unsigned n_element=mesh_pt->nelement();

    // Number of nodes in an element
    unsigned n_el_node=mesh_pt->finite_element_pt(0)->nnode();

    // Number of spatial dimensions
    unsigned n_spatial_dim=2;

    // Eulerian position
    Vector<double> x(n_spatial_dim+1,0.0);

    // Loop over the nodes in the mesh
    for (unsigned i=0; i<n_element; i++)
    {
      // Get a pointer to this element
      FiniteElement* el_pt=mesh_pt->finite_element_pt(i);

      // Loop over the nodes in the i-th element
      for (unsigned j=0; j<n_el_node; j++)
      {
        // Get a pointer to the j-th node in the i-th element
        Node* nod_pt=el_pt->node_pt(j);

        // Might want pressure information later too so make sure it's not pinned
        if (!(nod_pt->is_pinned(n_spatial_dim)))
        {
          // Loop over the coordinates
          for (unsigned k=0; k<n_spatial_dim+1; k++)
          {
            // Get the i-th coordinate
            x[k]=nod_pt->x(k);
          }

          // Check if the node lies on the centerline, outside the central box
          // and on the initial time-boundary
          if ((x[0]>(0.5*Length_of_central_box))&&
              (std::abs(x[1])<1.0e-10)&&
              (std::abs(x[2])<1.0e-10))
          {
            // Store a pointer to the chosen element
            (*el_centerline_pt)=el_pt;

            // Store the local nodal number
            node_index=j;

            // We're done; we only need one node
            return;
          }
        } // if (!(el_pt->node_pt(i)->is_pinned(n_spatial_dim)))
      } // for (unsigned j=0;j<n_el_node;j++)
    } // for (unsigned i=0;i<n_element;i++)
  } // End of find_node_on_centerline

  /// Function to round a double to the nearest integral value
  double round(const double& d)
  {
    // Round it
    return std::floor(d+0.5);
  } // End of round
  /// ---------------------------Miscellaneous-----------------------------
} // End of GlobalParameters

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//======start_of_Subsidiary_Preconditioner_Helper=======================
/// The function get_new_preconditioner() returns an instance of
/// the chosen preconditioner to be used as a subsidiary preconditioner
/// in a GeneralPurposeBlockPreconditioner
//======================================================================
namespace SubsidiaryPreconditionerHelper
{
  /// Function used to generate instances of the chosen subsidiary
  /// preconditioner to solve the (subsidiary) block systems
  Preconditioner* get_new_preconditioner()
  {
    // Return a new instance of SuperLU
    return new SuperLUPreconditioner;
  } // End of get_new_preconditioner
} // End of namespace SubsidiaryPreconditionerHelper


//=====start_of_problem_class=========================================
/// NavierStokes problem
//====================================================================
template<class ELEMENT>
class NavierStokesProblem : public Problem
{
public:

  /// Constructor
  NavierStokesProblem();

  /// Destructor
  ~NavierStokesProblem();

  /// Update the problem specs before solve (empty)
  void actions_before_newton_solve() {}

  /// Update the problem specs before solve (empty)
  void actions_after_newton_solve() {}

  /// Update the problem specs after an increase in a parameter
  void actions_after_parameter_increase(double* const& parameter_pt);

  /// Create the space-time mesh with the chosen number of elements in
  /// the time direction and the chosen spatial resolution
  void create_spacetime_mesh();

  /// Assign the appropriate boundary conditions and enforce periodicity
  /// in the time direction
  void apply_boundary_conditions();

  /// The mixed order elements use linear interpolation in time so the
  /// only nodes which contribute to the unknowns in the system are those that
  /// lie on the temporal boundaries of the elements. Thus, all nodes that do
  /// not lie on these boundaries need to be pinned (otherwise we'd get zero
  /// rows in the system matrix making it singular...).
  void pin_redundant_temporal_nodes();

  /// Assign the time slice number to each element
  void assign_time_slice_id();

  /// Assign the chosen solver to this Problem (and preconditioner if
  /// so desired)
  void set_up_spacetime_solver();

  /// Complete problem setup: do anything else that's needed to make the
  /// elements fully functional (e.g. pass pointers to problem parameters)
  void complete_problem_setup();

  /// Use continuation in a particular parameter. This should make
  /// things sufficiently generic so *natural* continuation can be used
  /// interchangeably with ease for any given input parameter.
  /// NOTE: It's important that the function:
  ///
  ///            actions_after_parameter_increase(...)
  ///
  /// be overloaded so that anything that needs to be updated after a change
  /// in the input parameter is indeed updated.
  void run_natural_continuation(double& parameter,
                                const double& parameter_target,
                                const unsigned& max_n_parameter_step);

  /// Doc the solution
  void doc_solution(const bool& doc_spacetime_soln=false);

private:

  /// Oomph-lib iterative linear solver
  IterativeLinearSolver* Solver_pt;

  /// Preconditioner
  Preconditioner* Prec_pt;

  // Enumeration of the preconditioners
  enum
  {
    Diagonal_preconditioner=0,
    Lower_triangular_preconditioner=1
  };

  /// Pointer to the space-time mesh
  ExtrudedCubeMeshFromQuadMesh<ELEMENT>* Bulk_mesh_pt;

  /// Pointer to the spatial mesh. NOTE: We keep this so that we
  /// retain a pointer to the Domain which will tell the ExtrudedDomain
  /// how to move the nodes when we want to do a parameter sweep in the
  /// cylinder oscillation amplitude...
  RefineableQuadMeshWithMovingCylinder<
  MyRefineableQTaylorHoodElement>* Spatial_mesh_pt;

  // Boolean variable to check if periodicity has been set up
  bool Periodicity_has_been_enforced;

  // Enumeration of the boundaries of the space-time mesh
  enum
  {
    Lower_wall_boundary_id=0,
    Outflow_boundary_id=1,
    Upper_wall_boundary_id=2,
    Inflow_boundary_id=3,
    Cylinder_surface_boundary_id=4,
    Initial_time_boundary_id=5,
    Final_time_boundary_id=6
  };
}; // End of NavierStokesProblem class


//========start_of_constructor============================================
/// Constructor for NavierStokes problem in cubic domain
//========================================================================
template<class ELEMENT>
NavierStokesProblem<ELEMENT>::NavierStokesProblem() :
  Solver_pt(0),
  Prec_pt(0),
  Bulk_mesh_pt(0),
  Spatial_mesh_pt(0),
  Periodicity_has_been_enforced(false)
{
  // Update the mesh-related parameters
  GlobalParameters::update_mesh_parameters();

  // Update the physical parameters
  GlobalParameters::update_physical_parameters();

  // Generate the space-time mesh (stored as Bulk_mesh_pt)
  create_spacetime_mesh();

  // Assign the appropriate boundary conditions
  apply_boundary_conditions();

  // Assign the time-slice ID to each element
  // NOTE: This must be done after the boundary conditions have been assigned
  // otherwise it won't know about the periodicity.)
  assign_time_slice_id();

  // Complete the problem setup to make the elements fully functional
  complete_problem_setup();

  // Make it pseudo-traction-free
  SpaceTimeNavierStokesMixedOrderEquations<2>::Gamma[0]=0.0;
  SpaceTimeNavierStokesMixedOrderEquations<2>::Gamma[1]=0.0;

  // Use the block lower triangular preconditioner
  GlobalParameters::Preconditioner=Lower_triangular_preconditioner;

  // Call the auxiliary solver setup function
  set_up_spacetime_solver();

  // Assign equation numbers
  oomph_info << "\nNumber of equations: " << assign_eqn_numbers() << std::endl;
} // End of NavierStokesProblem


//========start_of_destructor=============================================
/// Destructor for NavierStokes problem in cubic domain
//========================================================================
template<class ELEMENT>
NavierStokesProblem<ELEMENT>::~NavierStokesProblem()
{
  // If the user is NOT using the default linear solver
  if (dynamic_cast<SuperLUSolver*>(linear_solver_pt())==0)
  {
    // Assign the linear solver
    delete linear_solver_pt();

    // Make it a null pointer
    linear_solver_pt()=0;

    // Make the (privately stored) solver pointer a null pointer
    Solver_pt=0;

    // We would only use a preconditioner if we're not using SuperLU so
    // check here if the preconditioner pointer has been set
    if (Prec_pt!=0)
    {
      // Delete the preconditioner
      delete Prec_pt;

      // Make it a null pointer
      Prec_pt=0;
    }
  } // if (dynamic_cast<SuperLUSolver*>(linear_solver_pt())==0)

  // Delete the space-time mesh
  delete Bulk_mesh_pt;

  // Make the pointer null
  Bulk_mesh_pt=0;

  // Delete the spatial mesh
  delete Spatial_mesh_pt;

  // Make the pointer null
  Spatial_mesh_pt=0;

  // Delete the cylinder
  delete GlobalParameters::Cylinder_pt;

  // Make it a null pointer
  GlobalParameters::Cylinder_pt=0;
} // End of ~NavierStokesProblem


//=====start_of_actions_after_parameter_increase============================
/// Update the problem specs after an increase in a parameter
//==========================================================================
template<class ELEMENT>
void NavierStokesProblem<ELEMENT>::
actions_after_parameter_increase(double* const& parameter_pt)
{
  // If we have been passed a pointer to the Reynolds number or period ratio
  if ((parameter_pt==&GlobalParameters::Re)||
      (parameter_pt==&GlobalParameters::Period_ratio))
  {
    // Call the physical parameters update function
    GlobalParameters::update_physical_parameters();

    // Update the boundary conditions (for the update of the velocities on
    // the cylinder boundary)
    apply_boundary_conditions();
  }
  else if (parameter_pt==&GlobalParameters::Amplitude)
  {
    // Update the nodal positions
    Bulk_mesh_pt->node_update();

    // Update the boundary conditions (for the update of the velocities on
    // the cylinder boundary)
    apply_boundary_conditions();
  }
  else
  {
    // Tell the user
    oomph_info << "\nCalled actions_after_parameter_increase(...) but I"
               << "\ndon't know what to do for this parameter. I'm going to "
               << "\nassume I'm not meant to do anything here. I hope you know"
               << "\nwhat you're doing...\n" << std::endl;
  }
} // End of actions_after_parameter_increase

//=====start_of_create_spacetime_mesh=======================================
/// Helper function to create the space-time mesh (to be assigned to
/// Bulk_mesh_pt) with the chosen number of elements in the time direction
/// and an appropriate spatial resolution (to capture the time-periodic
/// solution properly).
//==========================================================================
template<class ELEMENT>
void NavierStokesProblem<ELEMENT>::create_spacetime_mesh()
{
  // Storage for the start time
  double start_t=0.0;

  //--------------------------
  // Generate 2D spatial mesh:
  //--------------------------
  // Record the start time
  start_t=TimingHelpers::timer();

  // Use BDF2
  add_time_stepper_pt(new Steady<0>);

  // Make a new cylinder
  GlobalParameters::Cylinder_pt=
    new OscillatingCylinder(&GlobalParameters::Radius,
                            &GlobalParameters::Amplitude,
                            time_pt());

  // Make a new mesh and assign its pointer
  Spatial_mesh_pt=
    new RefineableQuadMeshWithMovingCylinder<MyRefineableQTaylorHoodElement>(
    GlobalParameters::Cylinder_pt,
    GlobalParameters::Annular_region_radius,
    GlobalParameters::Length_of_central_box,
    GlobalParameters::X_left,
    GlobalParameters::X_right,
    GlobalParameters::Height);

  // Loop over the refinements
  for (unsigned i=0; i<GlobalParameters::N_uniform_refinement_before_solve; i++)
  {
    // Refine the mesh
    Spatial_mesh_pt->refine_uniformly();
  }

  // Quick statistics about the mesh and it's setup
  oomph_info << "\nNumber of nodes in spatial mesh: "
             << Spatial_mesh_pt->nnode()
             << "\nNumber of elements in spatial mesh: "
             << Spatial_mesh_pt->nelement()
             << "\nTime taken to generate refined spatial mesh [sec]: "
             << TimingHelpers::timer()-start_t << std::endl;

  //---------------------------------
  // Generate the 3D space-time mesh:
  //---------------------------------
  // Indicate that we want the mesh extrusion time to be doc-ed
  MeshExtrusionHelpers::Mesh_extrusion_helper.enable_doc_mesh_setup_time();

  // Create the extruded mesh
  Bulk_mesh_pt=new ExtrudedCubeMeshFromQuadMesh<ELEMENT>
  (Spatial_mesh_pt,GlobalParameters::N_t,GlobalParameters::L_t);

  // The created space-time mesh is the only mesh so assign it
  Problem::mesh_pt()=Bulk_mesh_pt;

  // Record the start time
  start_t=TimingHelpers::timer();

  // Pin the redundant temporal nodes
  pin_redundant_temporal_nodes();

  // Reorder the nodes (don't use the regular setup ordering)
  NodeReordering::reorder_nodes(Bulk_mesh_pt,false);

  // Document the setup time
  oomph_info << "\nTime taken for redundant node pinning/node reordering [sec]: "
             << TimingHelpers::timer()-start_t << std::endl;

  // Check the maximum deformation of the mesh due to the cylinder motion
  GlobalParameters::doc_maximum_central_box_deformation();

  // Tell the user we've finished
  oomph_info << "\nCompleted mesh generation and documentation!" << std::endl;
} // End of create_spacetime_mesh


//=====start_of_pin_redundant_temporal_nodes================================
/// The mixed order elements use linear interpolation in time so the
/// only nodes which contribute to the unknowns in the system are those that
/// lie on the temporal boundaries of the elements. Thus, all nodes that do
/// not lie on these boundaries need to be pinned (otherwise we'd get zero
/// rows in the system matrix making it singular...).
//==========================================================================
template<class ELEMENT>
void NavierStokesProblem<ELEMENT>::pin_redundant_temporal_nodes()
{
  // Number of nodes in each direction
  unsigned n_node_1d=Bulk_mesh_pt->finite_element_pt(0)->nnode_1d();

  // Number of nodes in a space-time element
  unsigned n_el_node=Bulk_mesh_pt->finite_element_pt(0)->nnode();

  // Sanity check: only works for 3D space-time elements (2D space + 1D time)
  if (n_el_node!=std::pow(n_node_1d,3))
  {
    // Throw an error
    throw OomphLibError("Can currently only deal with 3D space-time elements!",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }

  // Get the number of elements in the mesh
  unsigned n_element=Bulk_mesh_pt->nelement();

  // Loop over the elements
  for (unsigned i=0; i<n_element; i++)
  {
    // Loop over the nodes
    for (unsigned j=0; j<n_el_node; j++)
    {
      // Storage for the local time slice ID (0<=i_temporal<=NNODE_1D-1)
      unsigned j_temporal=0;

      // The spatial node number
      unsigned j_spatial=j%(n_node_1d*n_node_1d);

      // Which local time slice are we in?
      j_temporal=(j-j_spatial)/(n_node_1d*n_node_1d);

      // If we're not on first/final elemental time slice
      if ((j_temporal!=0)&&(j_temporal!=n_node_1d-1))
      {
        // Get a pointer to the j-th node in the i-th element
        Node* node_pt=Bulk_mesh_pt->finite_element_pt(i)->node_pt(j);

        // Get the number of unknowns at this node
        unsigned n_value=node_pt->nvalue();

        // Loop over the unknowns
        for (unsigned k=0; k<n_value; k++)
        {
          // Pin the k-th unknown at this node
          node_pt->pin(k);
        }
      } // if ((j_temporal!=0)&&(j_temporal!=n_node_1d-1))
    } // for (unsigned j=0;j<n_el_node;j++)
  } // for (unsigned i=0;i<n_element;i++)
} // End of pin_redundant_temporal_nodes


//========================================================================
/// Set up the solver for this Problem
//========================================================================
template<class ELEMENT>
void NavierStokesProblem<ELEMENT>::set_up_spacetime_solver()
{
  // Create oomph-lib iterative linear solver
  Solver_pt=new GMRES<CRDoubleMatrix>;

  // Use RHS preconditioning
  //dynamic_cast<GMRES<CRDoubleMatrix>*>(Solver_pt)->set_preconditioner_RHS();
  dynamic_cast<GMRES<CRDoubleMatrix>*>(Solver_pt)->set_preconditioner_LHS();

  // Set the tolerance
  Solver_pt->tolerance()=1.0e-10;

  // Maximum number of iterations
  Solver_pt->max_iter()=200;

  // Set linear solver
  linear_solver_pt()=Solver_pt;

  //-----------------------------------------
  // Create the master-level preconditioners:
  //-----------------------------------------
  // Do we want to document the memory usage?
  bool document_memory_usage=true;

  // Solve the diagonal blocks associated with each time-slice separately
  if (GlobalParameters::Preconditioner==Diagonal_preconditioner)
  {
    // Create a new instance of the space-time preconditioner
    Prec_pt=new BlockDiagonalPreconditioner<CRDoubleMatrix>;
  }
  // Solve the block lower-triangular part of the system matrix
  else if (GlobalParameters::Preconditioner==Lower_triangular_preconditioner)
  {
    // Create a new instance of the space-time preconditioner
    Prec_pt=new BandedBlockTriangularPreconditioner<CRDoubleMatrix>;

    // Indicate that we're using a block lower triangular solve
    dynamic_cast<BandedBlockTriangularPreconditioner<CRDoubleMatrix>*>
    (Prec_pt)->lower_triangular();

    // The order of the interpolation in the time direction (linear)
    unsigned temporal_order=1;

    // Provide the bandwidth; only subdiagonal block entries
    dynamic_cast<BandedBlockTriangularPreconditioner<CRDoubleMatrix>*>
    (Prec_pt)->set_block_bandwidth(temporal_order);

    // If we want to document the memory usage
    if (document_memory_usage)
    {
      dynamic_cast<BandedBlockTriangularPreconditioner<CRDoubleMatrix>*>
      (Prec_pt)->enable_doc_memory_usage();
    }
  }
  // If the user provided an invalid input
  else
  {
    // Throw an error
    throw OomphLibError("Invalid choice of preconditioner.",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
  }

  // Allocate space for the DOF to block map; this tells the master block
  // preconditioner which DOF types to aggregate
  Vector<unsigned> dof_to_block_map;

  // Call the auxiliary function which sets up the mapping
  GlobalParameters::set_up_dof_to_block_mapping(dof_to_block_map);

  // Create an upcasted pointer to the master preconditioner
  GeneralPurposeBlockPreconditioner<CRDoubleMatrix>* upcasted_master_prec_pt=
    dynamic_cast<GeneralPurposeBlockPreconditioner<CRDoubleMatrix>*>(Prec_pt);

  // Build silently!
  upcasted_master_prec_pt->enable_silent_preconditioner_setup();

  // Pass the DOF-to-block map to the preconditioner
  upcasted_master_prec_pt->set_dof_to_block_map(dof_to_block_map);

  // Specify the subsidiary block preconditioner
  upcasted_master_prec_pt->set_subsidiary_preconditioner_function(
    SubsidiaryPreconditionerHelper::get_new_preconditioner);

  // Pass a pointer to the (space-time) mesh
  upcasted_master_prec_pt->add_mesh(Bulk_mesh_pt);

  // Now assign the preconditioner to the linear solver
  Solver_pt->preconditioner_pt()=Prec_pt;
} // End of set_up_spacetime_solver


//========================================================================
/// Complete problem setup: pass pointers to physical variables.
//========================================================================
template<class ELEMENT>
void NavierStokesProblem<ELEMENT>::complete_problem_setup()
{
  // Get the number of elements in the bulk mesh
  unsigned n_bulk_element=Bulk_mesh_pt->nelement();

  // Loop over the bulk elements
  for (unsigned e=0; e<n_bulk_element; e++)
  {
    // Upcast to a fluid element
    ELEMENT* el_pt=dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

    // Set the Reynolds number (set equal to the Womersley number)
    el_pt->re_pt()=&GlobalParameters::Re;

    // Set the Strouhal number
    el_pt->re_st_pt()=&GlobalParameters::ReSt;
  }
} // End of complete_problem_setup


//=========start_of_element_to_ijk_map_setup==============================
/// Set up the map which maps a given element e to it's (i,j,k)
/// coordinates in the mesh where i,j and k are integers (such that
/// 0<=i<N_x, 0<=j<N_y, 0<=k<N_t).
//========================================================================
template<class ELEMENT>
void NavierStokesProblem<ELEMENT>::assign_time_slice_id()
{
  // Number of dimensions
  unsigned n_dim=3;

  // Space for the local coordinates (at the center of the element)
  Vector<double> s(n_dim,0.0);

  // The ID of the coordinate we want, i.e. the ID of the time-direction
  unsigned time_index=n_dim-1;

  // Storage for the time coordinate
  double t=0.0;

  // Get the number of elements in the mesh
  unsigned n_element=Bulk_mesh_pt->nelement();

  // Loop over the elements
  for (unsigned i=0; i<n_element; i++)
  {
    // Get the Eulerian coordinates at the center of the element
    t=Bulk_mesh_pt->finite_element_pt(i)->interpolated_x(s,time_index);

    // The length of an element in the time direction
    double el_length=GlobalParameters::L_t/double(GlobalParameters::N_t);

    // Subtract the length of half an element off (since we're at the centre
    // of the element using s=(0,0,0))
    t-=el_length/2.0;

    // Divide by the length of an element
    t/=el_length;

    // Store it (have to round first otherwise 0.999 would become 0)
    unsigned id=unsigned(GlobalParameters::round(t));

    // Upcast the element
    ELEMENT* el_pt=dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));

    // Assign the time slice ID
    el_pt->set_time_slab_id(id);

    // There are 5 unique dofs per element 3 in common with neighbouring
    // time slices (i.e. u, v and p) so 3 unique ones per element. We can
    // either use an impulsive start (i.e. pin dofs in the first time slice)
    // or periodic boundary conditions (the dofs in the final time slice
    // are the dofs in the first time slice) so the actual number of dof
    // types is always the same
    GlobalParameters::N_dof_type=3*GlobalParameters::N_t;

    // Finally, tell it how many dof types there are in the mesh
    el_pt->set_ndof_types(GlobalParameters::N_dof_type);
  } // for (unsigned i=0;i<n_element;i++)
} // End of assign_time_slice_id


//=========start_of_apply_boundary_conditions====================
/// Assign the appropriate boundary conditions, i.e. periodicity
/// in the t-direction. In the x and y-direction apply Dirichlet boundary
/// conditions. In summary:
///             Boundary 0 (t=0) -- Periodic in time (with boundary 5)
///             Boundary 1 (y=0) -- Dirichlet
///             Boundary 2 (x=1) -- Dirichlet
///             Boundary 3 (y=1) -- Dirichlet
///             Boundary 4 (x=0) -- Dirichlet
///             Boundary 5 (t=1) -- Periodic in time (with boundary 0)
//========================================================================
template<class ELEMENT>
void NavierStokesProblem<ELEMENT>::apply_boundary_conditions()
{
  // Get the start time
  double start_t=TimingHelpers::timer();

  // Number of SPATIAL dimensions
  unsigned n_dim=3;

  // Storage for the spatial coordinates
  Vector<double> spatial_coordinates(n_dim-1,0.0);

  // Storage for the solution (also passes the pressure field)
  Vector<double> u_exact(n_dim,0.0);

  // Number of nodes on t=0 boundary
  unsigned n_boundary0_node=Bulk_mesh_pt->
                            nboundary_node(Initial_time_boundary_id);

  // Number of nodes on t=1 boundary
  unsigned n_boundary1_node=Bulk_mesh_pt->
                            nboundary_node(Final_time_boundary_id);

  // Get the number of boundaries in the mesh
  unsigned n_boundary=Bulk_mesh_pt->nboundary();

  // Loop over the boundaries
  for (unsigned b=0; b<n_boundary; b++)
  {
    // Get the number of nodes on the b-th boundary
    unsigned n_node=Bulk_mesh_pt->nboundary_node(b);

    // If we're on either of the time boundaries
    if ((b==Initial_time_boundary_id)||(b==Final_time_boundary_id))
    {
      // If we need to set up periodicity
      if (!Periodicity_has_been_enforced)
      {
        //----------------------------
        // Establish nodal periodicity
        //----------------------------
        // Make sure there are as many nodes on boundary 0 as there are on
        // boundary (n_boundary-1)
        if (n_boundary0_node!=n_boundary1_node)
        {
          // Create an output stream
          std::ofstream output_file;

          // Open a file
          output_file.open("RESLT/nodes_b0.csv");

          // Loop over the nodes on the t=0 boundary
          for (unsigned i=0; i<n_boundary0_node; i++)
          {
            // Output the coordinates of the i-th node on the t=0 boundary
            Bulk_mesh_pt->boundary_node_pt(Initial_time_boundary_id,i)->
            output(output_file);
          }

          // Close the file
          output_file.close();

          // Open a file
          output_file.open("RESLT/nodes_b1.csv");

          // Loop over the nodes on the t=0 boundary
          for (unsigned i=0; i<n_boundary1_node; i++)
          {
            // Output the coordinates of the i-th node on the t=1 boundary
            Bulk_mesh_pt->boundary_node_pt(Final_time_boundary_id,i)->
            output(output_file);
          }

          // Close the file
          output_file.close();

          // Create an output stream
          std::ostringstream error_message_stream;

          // Create an error message
          error_message_stream << "Different number of nodes on t=0 and t=1 "
                               << "boundary!\nThere are " << n_boundary0_node
                               << " nodes on boundary " << Initial_time_boundary_id
                               << " and " << n_boundary1_node
                               << " nodes on boundary " << Final_time_boundary_id
                               << "!" << std::endl;

          // Throw an error
          throw OomphLibError(error_message_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        }

        // Loop over the nodes on the t=0 boundary
        for (unsigned i=0; i<n_boundary0_node; i++)
        {
          // Get the pointer to the associated node
          Node* node0_pt=Bulk_mesh_pt->boundary_node_pt(Initial_time_boundary_id,i);

          // Boolean to indicate whether or not the neighbour has been found
          bool has_neighbour_node_been_found=false;

          // Loop over the nodes on the t=1 boundary
          for (unsigned j=0; i<n_boundary1_node; j++)
          {
            // Get the pointer to the associated node
            Node* node1_pt=Bulk_mesh_pt->boundary_node_pt(Final_time_boundary_id,j);

            // Distance value
            double distance=0.0;

            // Loop over the entries of x
            for (unsigned k=0; k<n_dim-1; k++)
            {
              // Update the distance (2 norm)
              distance+=pow(((node0_pt->x(k))-(node1_pt->x(k))),2.0);
            }

            // Square root it
            distance=std::sqrt(distance);

            // Check if it matches to within a reasonable tolerance
            if (std::fabs(distance)<Tree::max_neighbour_finding_tolerance())
            {
              // Make the nodes periodic; the node on the t=1 boundary now points
              // to the node on the t=0 boundary.
              node1_pt->make_periodic(node0_pt);

              // We've found the neighbouring node
              has_neighbour_node_been_found=true;

              // We're done; break out!
              break;
            }
          } // for (unsigned i=0;i<n_boundary0_node;i++)

          // If we get here and we haven't found the neighbouring node, something's
          // wrong so throw an error
          if (!has_neighbour_node_been_found)
          {
            // Throw an error
            throw OomphLibError("Couldn't find neighbouring node!",
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);
          }
        } // for (unsigned i=0;i<n_boundary0_node;i++)

        // If we've got here then the periodicity has been set up
        Periodicity_has_been_enforced=true;
      } // if (!Periodicity_has_been_enforced)
    }
    // If we're not on the time boundaries
    else
    {
      // If we're on the lower/upper wall
      if ((b==Lower_wall_boundary_id)||(b==Upper_wall_boundary_id))
      {
        // Loop over the nodes on the b-th boundary
        for (unsigned n=0; n<n_node; n++)
        {
          // Get a pointer to the n-th node on the b-th boundary
          Node* node_pt=Bulk_mesh_pt->boundary_node_pt(b,n);

          // Loop over the velocity components
          for (unsigned i=0; i<n_dim-1; i++)
          {
            // Pin the i-th velocity component
            node_pt->pin(i);
          }

          // Apply tow-tank boundary conditions (nonzero horizontal velocity)
          node_pt->set_value(0,1.0);

          // Apply tow-tank boundary conditions (no vertical velocity)
          node_pt->set_value(1,0.0);
        } // for (unsigned n=0;n<n_node;n++)
      }
      // If we're at the outflow boundary
      else if (b==Outflow_boundary_id)
      {
        // Don't actually do anything; it's better than forcing parallel flow
        // at the outlet (no numerical oscillations/boundary layer near exit)
      }
      // If we're at the inflow boundary
      else if (b==Inflow_boundary_id)
      {
        // Loop over the nodes
        for (unsigned n=0; n<n_node; n++)
        {
          // Pointer to the boundary node
          Node* boundary_node_pt=Bulk_mesh_pt->boundary_node_pt(b,n);

          // Pin the horizontal velocity
          boundary_node_pt->pin(0);

          // Pin the vertical velocity
          boundary_node_pt->pin(1);

          // Uniform flow
          boundary_node_pt->set_value(0,1.0);

          // Parallel flow means no vertical velocity component
          boundary_node_pt->set_value(1,0.0);
        }
      }
      // If we're dealing with the cylinder surface (apply no-slip conditions)
      else if (b==Cylinder_surface_boundary_id)
      {
        // Loop over the nodes
        for (unsigned i_nod=0; i_nod<n_node; i_nod++)
        {
          // Pointer to the boundary node
          Node* boundary_node_pt=Bulk_mesh_pt->boundary_node_pt(b,i_nod);

          // Pin the horizontal velocity at every node
          boundary_node_pt->pin(0);

          // Pin the vertical velocity at every node
          boundary_node_pt->pin(1);

          // Get the time value associated with this node
          double time=boundary_node_pt->x(n_dim-1);

          // Storage for the cylinder velocity
          Vector<double> u(2,0.0);

          // Get the velocity of the cylinder at this point in time
          GlobalParameters::Cylinder_pt->velocity(time,u);

          //--------------------------------------------------------------------
          // Apply no-slip condition for NS on a moving wall node noting:
          //        Uu = (D/T_e) dR/dt => u = St dR/dt = St dR_cyl/dt,
          // where the Strouhal number is defined as St = D/(UT_e).
          //--------------------------------------------------------------------
          // Set the horizontal velocity value
          boundary_node_pt->set_value(0,GlobalParameters::St*u[0]);

          // Set the vertical velocity value
          boundary_node_pt->set_value(1,GlobalParameters::St*u[1]);
        }
      } // if ((b==Lower_wall_boundary_id)||(b==Upper_wall_boundary_id))
    } // if ((b==Initial_time_boundary_id)||(b==Final_time_boundary_id))
  } // for (unsigned b=0;b<n_bound;b++)

  // Record the end time
  double end_t=TimingHelpers::timer();

  // Compute the time taken
  double boundary_condition_application_time=end_t-start_t;

  // Output the setup time to the screen
  oomph_info << "Time taken for application of boundary conditions [sec]: "
             << boundary_condition_application_time << std::endl;
} // End of apply_boundary_conditions


//=======start_of_run_natural_continuation=================================
/// Use continuation in a particular parameter. This should make
/// things sufficiently generic so natural/arc-length continuation can
/// be used interchangeably with ease for any given input parameter.
/// NOTE: It's important that the function:
///
///            actions_after_parameter_increase(...)
///
/// be overloaded so that anything that needs to be updated after a change
/// in the input parameter is indeed updated.
///
/// Inputs:
/// -------
///        (1) The parameter that we're running continuation in;
///        (2) The target parameter value;
///        (4) The max. number of steps used to reach parameter_target;
///        (5) (Optional) Whether or not to document the solution at
///            each step.
//=========================================================================
template<class ELEMENT>
void NavierStokesProblem<ELEMENT>::
run_natural_continuation(double& parameter,
                         const double& parameter_target,
                         const unsigned& max_n_parameter_step)
{
  // If we can actually do anything
  if ((std::fabs(parameter_target-parameter)>=1.0e-14)&&
      (max_n_parameter_step!=0))
  {
    // Store the string used to denote the input parameter
    std::string parameter_string=
      GlobalParameters::parameter_to_string(&parameter);

    // Copy the initial parameter value
    double initial_parameter=parameter;

    // The size of the amplitude value increment
    double parameter_increment=((parameter_target-parameter)/
                                double(max_n_parameter_step));

    // Tell the user what we're doing
    oomph_info << ANSIEscapeCode::Red
               << "\nStarting natural continuation process for "
               << parameter_string << "!" << ANSIEscapeCode::Reset
               << "\n\nInitial " << parameter_string << " value: " << parameter
               << "\nTarget " << parameter_string << " value: "
               << parameter_target << "\nMax. number of parameter steps: "
               << max_n_parameter_step << "\n\n"
               << ANSIEscapeCode::Red
               << "Parameter value cases (" << parameter_string << "): "
               << ANSIEscapeCode::Reset << std::endl;

    // Loop over the increments
    for (unsigned i=0; i<max_n_parameter_step+1; i++)
    {
      // Tell the user
      oomph_info << " (" << i << ") "
                 << initial_parameter+i*parameter_increment
                 << std::endl;
    } // for (unsigned i=0;i<max_n_parameter_step+1;i++)

    // Vector to contain the Problem dofs
    DoubleVector dofs_backup;

    // The maximum number of times tono halve the parameter increment
    unsigned max_n_reattempt=30;

    // The number of times we've halved the parameter increment
    unsigned n_reattempt=0;

    // Loop over the increments
    while (std::fabs(parameter_target-parameter)>1.0e-14)
    {
      // Get the dofs
      get_dofs(dofs_backup);

      // Increment the parameter value
      parameter+=parameter_increment;

      // Tell the user
      oomph_info << "\n" << ANSIEscapeCode::Red
                 << "Solving for " << parameter_string << " value: "
                 << ANSIEscapeCode::Reset << parameter << std::endl;

      // Update anything that needs to be changed after a change in this parameter
      actions_after_parameter_increase(&parameter);

      // Try doing a solve
      try
      {
        // Solve for this parameter value
        newton_solve();
      }
      // If the simulation threw up an error
      catch (NewtonSolverError& error)
      {
        // Make sure we haven't had to try too many times
        if (n_reattempt<max_n_reattempt)
        {
          // Decrement the parameter value
          parameter-=parameter_increment;

          // Make all the necessary updates after a change in this parameter
          actions_after_parameter_increase(&parameter);

          // Reset the dofs
          set_dofs(dofs_backup);

          // Half the parameter increment
          parameter_increment*=0.5;

          // Tell the user
          oomph_info << "\n" << ANSIEscapeCode::Red
                     << "Solve failed! Re-attempt " << n_reattempt+1
                     << " -- halving parameter increment (for "
                     << parameter_string << ") to: "
                     << ANSIEscapeCode::Reset;

          // The parameter increment is actually just this increment
          oomph_info << parameter_increment << std::endl;

          // Indicate that we're trying again
          n_reattempt++;

          // Try again...
          continue;
        }
        // If we've had to try too many times
        else
        {
          // Create an output stream
          std::ostringstream error_message_stream;

          // Create the error message
          error_message_stream << "Re-attempted a solve too many times. "
                               << "Exiting here!" << std::endl;

          // Throw an error
          throw OomphLibError(error_message_stream.str(),
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
        } // if (n_reattempt<max_n_reattempt)
      } // try
    } // while (std::fabs(parameter_target-parameter)<1.0e-08)
  } // if ((std::fabs(parameter_target-parameter)>1.0e-10)&&(max...
} // End of run_natural_continuation


//=======start_of_doc_solution============================================
/// Document the solution
//========================================================================
template<class ELEMENT>
void NavierStokesProblem<ELEMENT>::doc_solution(const bool& doc_spacetime_soln)
{
  // Start the clock
  double timer_s=TimingHelpers::timer();

  // Make an ofstream object to output the solution
  std::ofstream some_file;

  // Storage for the filename
  char filename[100];

  // Number of plot points to use for the big space-time solution
  unsigned n_plot_point=3;

  //-----------------
  // Output solution:
  //-----------------
  // Create the filename suffix
  sprintf(filename,"%s/soln%i",
          GlobalParameters::Doc_info.directory().c_str(),
          GlobalParameters::Doc_info.number());

  // Convert it to a string to allow ease of replacing dots
  std::string filename_as_string(filename);

  // Replace all dots with the string "pt"
  std::replace(filename_as_string.begin(),filename_as_string.end(),'.','p');

  // Now append the filename extension
  filename_as_string+=".dat";

  // Open a file with the constructed filename
  some_file.open(filename_as_string.c_str());

  // Set the precision of the outputted data
  some_file.precision(20);

  // Output the (numerically) approximated solution
  Bulk_mesh_pt->output(some_file,n_plot_point);

  // We're done; close the file
  some_file.close();

  // Finally, output the time taken
  oomph_info << "Total time for documentation [sec]: "
             << TimingHelpers::timer()-timer_s << std::endl;

  // Increment counter for solutions
  GlobalParameters::Doc_info.number()++;
} // End of doc_solution

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////

//=======start_of_main====================================================
/// Driver code for unsteady heat equation
//========================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
  // Initialise MPI
  MPI_Helpers::init(argc,argv);

  // Switch off output modifier
  oomph_info.output_modifier_pt()=&default_output_modifier;

  // Switch off oomph_info output for all processors but rank 0
  if (MPI_Helpers::communicator_pt()->my_rank()!=0)
  {
    oomph_info.stream_pt()=&oomph_nullstream;
    OomphLibWarning::set_stream_pt(&oomph_nullstream);
    OomphLibError::set_stream_pt(&oomph_nullstream);
  }
  else
  {
    oomph_info << "\n=====================================================\n"
               << "Number of processors: "
               << MPI_Helpers::communicator_pt()->nproc()
               << "\n=====================================================\n"
               << std::endl;
  }
#endif

  // Output directory
  GlobalParameters::Doc_info.set_directory("RESLT");

  //----------------------
  // Problem setup & solve
  //----------------------
  // Typedef the block preconditionable element type
  typedef BlockPrecQTaylorHoodMixedOrderSpaceTimeElement ELEMENT;

  // Start the clock
  double timer_s=TimingHelpers::timer();

  // Build problem
  NavierStokesProblem<ELEMENT> problem;

  // Tell the user
  oomph_info << ANSIEscapeCode::Red
             << "\nSolving with the following problem parameters:"
             << ANSIEscapeCode::Reset
             << ANSIEscapeCode::Red << "\n - Re: " << ANSIEscapeCode::Reset
             << GlobalParameters::Re
             << ANSIEscapeCode::Red << "\n - ReSt: " << ANSIEscapeCode::Reset
             << GlobalParameters::ReSt
             << ANSIEscapeCode::Red << "\n - A: " << ANSIEscapeCode::Reset
             << GlobalParameters::Amplitude
             << ANSIEscapeCode::Red << "\n - T_e/T_s: " << ANSIEscapeCode::Reset
             << GlobalParameters::Period_ratio
             << std::endl;

  // Solve the problem again
  problem.newton_solve();

  // Doc the solution
  problem.doc_solution();

  // Use parameter continuation in the Reynolds number
  problem.run_natural_continuation(GlobalParameters::Re,
                                   GlobalParameters::Re_target,
                                   GlobalParameters::N_re_step);

  // Doc the solution
  problem.doc_solution();

  // Use parameter continuation for the amplitude value
  problem.run_natural_continuation(GlobalParameters::Amplitude,
                                   GlobalParameters::Amplitude_target,
                                   GlobalParameters::N_amplitude_step);

  // Doc the solution
  problem.doc_solution();

  // Tell the user we're done
  oomph_info << "\nSimulation complete!\nTotal time for simulation [sec]: "
             << TimingHelpers::timer()-timer_s << std::endl;

#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif
} // End of main
