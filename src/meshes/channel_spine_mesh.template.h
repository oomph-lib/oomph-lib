// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
#ifndef OOMPH_CHANNEL_SPINE_MESH_HEADER
#define OOMPH_CHANNEL_SPINE_MESH_HEADER

// oomph-lib includes
#include "../generic/spines.h"
#include "rectangular_quadmesh.template.h"

namespace oomph
{
  //======================================================================
  /// Spine mesh class derived from standard 2D mesh.
  /// The mesh contains a StraightLine GeomObject which defines the height
  /// of the left and right regions (0,2) and another GeomObject is passed
  /// to the constructor to define the height in the central region.
  //======================================================================
  template<class ELEMENT>
  class ChannelSpineMesh : public RectangularQuadMesh<ELEMENT>, public SpineMesh
  {
  public:
    /// Constructor: Pass number of elements in x-direction in regions
    /// 0,1 and 2, number of elements in y-direction, length in x direction in
    /// regions 0,1 and 2, height mesh, pointer to the GeomObject defining the
    /// heightof the central region and pointer to timestepper (defaults to
    /// Steady timestepper)
    ChannelSpineMesh(const unsigned& nx0,
                     const unsigned& nx1,
                     const unsigned& nx2,
                     const unsigned& ny,
                     const double& lx0,
                     const double& lx1,
                     const double& lx2,
                     const double& h,
                     GeomObject* wall_pt,
                     TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// Constructor: Pass number of elements in x-direction in regions
    /// 0,1 and 2, number of elements in y-direction, length in x direction in
    /// regions 0,1 and 2, height mesh, pointer to the GeomObject defining the
    /// heightof the central region, a boolean flag to indicate whether or not
    /// the mesh is periodic and pointer to timestepper (defaults to Steady
    /// timestepper)
    ChannelSpineMesh(const unsigned& nx0,
                     const unsigned& nx1,
                     const unsigned& nx2,
                     const unsigned& ny,
                     const double& lx0,
                     const double& lx1,
                     const double& lx2,
                     const double& h,
                     GeomObject* wall_pt,
                     const bool& periodic_in_x,
                     TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper);

    /// Access functions for pointers to the \f$ i \f$ -th element in
    /// the left region.
    FiniteElement*& left_element_pt(const unsigned long& i)
    {
      return Left_element_pt[i];
    }

    /// Access functions for pointers to the \f$ i \f$ -th element in
    /// the centre region.
    FiniteElement*& centre_element_pt(const unsigned long& i)
    {
      return Centre_element_pt[i];
    }

    /// Access functions for pointers to the \f$ i \f$ -th element in
    /// the right region.
    FiniteElement*& right_element_pt(const unsigned long& i)
    {
      return Right_element_pt[i];
    }

    /// Number of elements in left region
    unsigned long nleft() const
    {
      return Left_element_pt.size();
    }

    /// Number of elements in centre region
    unsigned long ncentre() const
    {
      return Centre_element_pt.size();
    }

    /// Number of elements in right region
    unsigned long nright() const
    {
      return Right_element_pt.size();
    }

    /// Number of elements in bulk
    unsigned long nbulk() const
    {
      unsigned long Nbulk = Left_element_pt.size() + Centre_element_pt.size() +
                            Right_element_pt.size();
      return Nbulk;
    }

    /// Reorder the elements so we loop over them vertically first
    /// (advantageous in "wide" domains if a frontal solver is used).
    void element_reorder();

    /// General node update function implements pure virtual function
    /// defined in SpineMesh base class and performs specific node update
    /// actions: along vertical spines
    virtual void spine_node_update(SpineNode* spine_node_pt)
    {
      // Get spine node's fraction along the spine
      double W = spine_node_pt->fraction();

      // Get local coordinates
      Vector<double> s_wall(1);
      s_wall[0] = spine_node_pt->spine_pt()->geom_parameter(0);

      // Get position vector to wall
      Vector<double> position(2);
      spine_node_pt->spine_pt()->geom_object_pt(0)->position(s_wall, position);

      // Set the value of y
      spine_node_pt->x(1) = this->Ymin + W * position[1];
    }

    /// Return the value of the x-coordinate at the node given by the
    /// local node number (xnode, ynode) in the element (xelement,yelement).
    /// The description is in a "psudeo" two-dimensional coordinate system,
    /// so the range of xelement is [0,Nx-1], yelement is [0,Ny-1], and
    /// that of xnode and ynode is [0,Np-1]. The default is to return
    /// nodes that are equally spaced in the x coodinate.
    virtual double x_spacing_function(unsigned xelement,
                                      unsigned xnode,
                                      unsigned yelement,
                                      unsigned ynode)
    {
      // Calculate the values of equal increments in nodal values in left region
      double xstep1 = Lx0 / ((this->Np - 1) * Nx0);
      // Calculate the values of equal increments in nodal values in centre
      // region
      double xstep2 = Lx1 / ((this->Np - 1) * Nx1);
      // Calculate the values of equal increments in nodal values in right
      // region
      double xstep3 = Lx2 / ((this->Np - 1) * Nx2);

      // left region
      if (xelement < Nx0)
      {
        // Return the appropriate value
        return (this->Xmin + xstep1 * ((this->Np - 1) * xelement + xnode));
      }
      // centre region
      else if (xelement < Nx0 + Nx1)
      {
        // Return the appropriate value
        return (Lx0 + xstep2 * ((this->Np - 1) * (xelement - Nx0) + xnode));
      }
      // right region
      else if (xelement < Nx0 + Nx1 + Nx2)
      {
        // Return the appropriate value
        return (Lx0 + Lx1 +
                xstep3 * ((this->Np - 1) * (xelement - Nx0 - Nx1) + xnode));
      }
      else
      {
        throw OomphLibError("Should not have got here",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      // Dummy return to keep compiler from barking
      return 0.0;
    }

    /// Access function for spines in left region
    Spine*& left_spine_pt(const unsigned long& i)
    {
#ifdef PARNOID
      if (i > Nleft_spine)
      {
        throw OomphLibError("Arguemnt out of range",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
#endif
      return Spine_pt[i];
    }

    /// Access function for spines in centre region
    Spine*& centre_spine_pt(const unsigned long& i)
    {
      if (i > Ncentre_spine)
      {
        throw OomphLibError("Arguemnt out of range",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        return Spine_pt[Nleft_spine + i];
      }
    }

    /// Access function for spines in right region
    Spine*& right_spine_pt(const unsigned long& i)
    {
      if (i > Nright_spine)
      {
        throw OomphLibError("Arguemnt out of range",
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        return Spine_pt[Nleft_spine + Ncentre_spine - 1 + i];
      }
    }

    /// Access function for the number of spines in the left region
    unsigned nleft_spine()
    {
      return Nleft_spine;
    }

    /// Access function for the number of spines in the centre region
    unsigned ncentre_spine()
    {
      return Ncentre_spine;
    }

    /// Access function for the number of spines in the right region
    unsigned nright_spine()
    {
      return Nright_spine;
    }

    /// Access function to the GeomObject for upper wall
    GeomObject* wall_pt()
    {
      return Wall_pt;
    }

    /// Access function to the GeomObject for the straight upper wall
    GeomObject* straight_wall_pt()
    {
      return Straight_wall_pt;
    }

  protected:
    /// Vector of pointers to element in the left region
    Vector<FiniteElement*> Left_element_pt;

    /// Vector of pointers to element in the centre region
    Vector<FiniteElement*> Centre_element_pt;

    /// Vector of pointers to element in the right region
    Vector<FiniteElement*> Right_element_pt;

    /// Helper function to actually build the channel-spine mesh
    /// (called from various constructors)
    virtual void build_channel_spine_mesh(TimeStepper* time_stepper_pt);

    /// Number of elements in the left region
    unsigned Nx0;

    /// Number of elements in the centre region
    unsigned Nx1;

    /// Number of elements in the right region
    unsigned Nx2;

    /// Length of left region
    double Lx0;

    /// Length of centre region
    double Lx1;

    /// Length of right region
    double Lx2;

    /// Number of spines in left region
    unsigned Nleft_spine;

    /// Number of spines in centre region
    unsigned Ncentre_spine;

    /// Number of spines in right region
    unsigned Nright_spine;

    /// GeomObject for upper wall
    GeomObject* Wall_pt;

    /// GeomObject for the straight upper wall
    GeomObject* Straight_wall_pt;
  };

} // namespace oomph

#endif
