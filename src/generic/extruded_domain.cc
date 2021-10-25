// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
// oomph-lib headers
#include "extruded_domain.h"
#include "extruded_macro_element.h"

namespace oomph
{
  //=================================================================
  /// Constructor
  //=================================================================
  ExtrudedDomain::ExtrudedDomain(Domain* domain_pt,
                                 const unsigned& n_extruded_element,
                                 const double& extrusion_length)
    : Domain(),
      Domain_pt(domain_pt),
      N_extruded_element(n_extruded_element),
      T_min(0.0),
      T_max(extrusion_length)
  {
    // The number of macro elements to create in the extruded domain
    unsigned n_macro_element =
      N_extruded_element * (Domain_pt->nmacro_element());

    // There are four macro elements
    Macro_element_pt.resize(n_macro_element);

    // Build the 3D extruded macro elements
    for (unsigned i = 0; i < n_macro_element; i++)
    {
      // Create the i-th macro element
      Macro_element_pt[i] = new QExtrudedMacroElement<3>(this, i);
    }
  } // End of ExtrudedDomain


  //=================================================================
  /// Constructor
  //=================================================================
  ExtrudedDomain::ExtrudedDomain(Domain* domain_pt,
                                 const unsigned& n_extruded_element,
                                 const double& t_min,
                                 const double& t_max)
    : Domain(),
      Domain_pt(domain_pt),
      N_extruded_element(n_extruded_element),
      T_min(t_min),
      T_max(t_max)
  {
    // The number of macro elements to create in the extruded domain
    unsigned n_macro_element =
      N_extruded_element * (Domain_pt->nmacro_element());

    // There are four macro elements
    Macro_element_pt.resize(n_macro_element);

    // Build the 3D extruded macro elements
    for (unsigned i = 0; i < n_macro_element; i++)
    {
      // Create the i-th macro element
      Macro_element_pt[i] = new QExtrudedMacroElement<3>(this, i);
    }

    /*
      // DRAIG: Might be worth having this as an output function if it's not
      // already in the base class...
      // Loop over the number of macro elements
      for (unsigned i=0;i<n_macro_element;i++)
      {
      // Create an output stream
      std::ofstream some_file;

      // Allocate space for the filename
      char filename[1000];

      // Number of plot points
      unsigned n_pts=20;

      // Create the file name
      sprintf(filename,"RESLT/extruded_element%i.dat",i);

      // Open a file with the created filename
      some_file.open(filename);

      // Output macro element
      macro_element_pt(i)->output(some_file,n_pts);

      // Close the file
      some_file.close();
      }
      exit(0);
    */
  } // End of ExtrudedDomain


  //=================================================================
  /// Number of macro elements in domain
  //=================================================================
  unsigned ExtrudedDomain::nmacro_element()
  {
    // Return the size of the macro element container
    return Macro_element_pt.size();
  } // End of nmacro_element


  //=================================================================
  /// Access to i-th extruded macro element
  //=================================================================
  ExtrudedMacroElement* ExtrudedDomain::macro_element_pt(const unsigned& i)
  {
    // Return a pointer to the i-th macro element in storage
    return dynamic_cast<ExtrudedMacroElement*>(Macro_element_pt[i]);
  } // End of macro_element_pt


  //=================================================================
  /// Vector representation of the i_macro-th macro element
  /// boundary i_direct (e.g. N/S/W/E in 2D spatial = 3D space-time).
  /// NOTE: Some extra care has to be taken here to translate the
  /// OcTree enumeration to the QuadTree enumeration (in the
  /// appropriate manner) so that the original Domain object can be
  /// used to calculate the global coordinate associated with the
  /// provided local coordinates.
  //=================================================================
  void ExtrudedDomain::macro_element_boundary(const unsigned& time,
                                              const unsigned& i_macro,
                                              const unsigned& i_direct,
                                              const Vector<double>& s,
                                              Vector<double>& x)
  {
    // Make sure that time=0 otherwise this doesn't make sense
    if (time != 0)
    {
      // Create an output stream
      std::ostringstream error_message_stream;

      // Create an error message
      error_message_stream << "This output function outputs a space-time\n"
                           << "representation of the domain. As such, it\n"
                           << "does not make sense to output the domain\n"
                           << "at a previous time level!" << std::endl;

      // Throw an error
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }

    // The number of dimensions
    unsigned n_dim = 3;

    // The number of macro elements in the original spatial Domain
    unsigned n_macro = Domain_pt->nmacro_element();

    // Calculate the shifted macro element number (to correspond to the
    // correct macro element in the original spatial Domain)
    unsigned true_macro_id = i_macro % n_macro;

    // Which layer of the extruded domain are we in? Layer 0 corresponds to
    // elements closest to t=T_min and the layer N_extruded_element-1
    // corresponds to elements closest to t=T_max
    unsigned i_layer = (i_macro - true_macro_id) / n_macro;

    // Calculate the width of an extrusion layer
    double layer_width = (T_max - T_min) / double(N_extruded_element);

    // Storage for the global (spatial) coordinates
    Vector<double> x_project(n_dim - 1, 0.0);

    // Calculate the time value associated with the start of the i-layer-th
    // extrusion layer
    double t_lower = (T_min + i_layer * layer_width);

    // Initialise the time value
    double t = t_lower;

    // If we're on one of the edges (once projected to 2D in the t-direction)
    if ((i_direct == OcTreeNames::L) || (i_direct == OcTreeNames::R) ||
        (i_direct == OcTreeNames::D) || (i_direct == OcTreeNames::U))
    {
      // Update the time value to get the time value associated
      // with the input surface coordinates
      t += 0.5 * (1.0 + s[1]) * layer_width;

      // Get the local coordinate associated with the surface coordinate
      // in the projected Domain. This is the first surface coordinate
      // on the left, right, down and up faces of the 3D macro element.
      Vector<double> s_project(1, s[0]);

      // If we're on the left face (or the West face once projected to 2D)
      if (i_direct == OcTreeNames::L)
      {
        // Call the corresponding function through the provided Domain object
        Domain_pt->macro_element_boundary(
          t, true_macro_id, QuadTreeNames::W, s_project, x_project);
      }
      // If we're on the right face (or the East face once projected to 2D)
      else if (i_direct == OcTreeNames::R)
      {
        // Call the corresponding function through the provided Domain object
        Domain_pt->macro_element_boundary(
          t, true_macro_id, QuadTreeNames::E, s_project, x_project);
      }
      // If we're on the down face (or the South face once projected to 2D)
      else if (i_direct == OcTreeNames::D)
      {
        // Call the corresponding function through the provided Domain object
        Domain_pt->macro_element_boundary(
          t, true_macro_id, QuadTreeNames::S, s_project, x_project);
      }
      // If we're on the up face (or the North face once projected to 2D)
      else if (i_direct == OcTreeNames::U)
      {
        // Call the corresponding function through the provided Domain object
        Domain_pt->macro_element_boundary(
          t, true_macro_id, QuadTreeNames::N, s_project, x_project);
      }
    }
    else if ((i_direct == OcTreeNames::B) || (i_direct == OcTreeNames::F))
    {
      // Grab the i_macro-th macro element from the Domain object
      MacroElement* macro_element_pt =
        Domain_pt->macro_element_pt(true_macro_id);

      // If we're on the back face (or the whole element at t=T_max)
      if (i_direct == OcTreeNames::B)
      {
        // Call the macro_map function to calculate the global (spatial)
        // coordinates associated with the given surface coordinates
        macro_element_pt->macro_map(t, s, x_project);
      }
      // If we're on the front face (or the whole element at t=T_max)
      else if (i_direct == OcTreeNames::F)
      {
        // Update the time value to get the time value associated
        // with the input surface coordinates
        t += layer_width;

        // Call the macro_map function to calculate the global (spatial)
        // coordinates associated with the given surface coordinates
        macro_element_pt->macro_map(t, s, x_project);
      }
    }
    else
    {
      // Create an output stream
      std::ostringstream error_message_stream;

      // Create an error message
      error_message_stream << "Incorrect face enumeration input! Should either "
                           << "be L,R,D,U,B or F. You input "
                           << OcTree::Direct_string[i_direct] << "!"
                           << std::endl;

      // We should never get here so throw an error if we do
      throw OomphLibError(error_message_stream.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    } // if ((i_direct==OcTreeNames::L)||(i_direct==OcTreeNames::R)||...

    // Loop over the global (spatial) coordinates
    for (unsigned i = 0; i < n_dim - 1; i++)
    {
      // Copy the global coordinates over
      x[i] = x_project[i];
    }

    // Finally, add in the final coordinate (the time value)
    x[n_dim - 1] = t;
  } // End of macro_element_boundary
} // End of namespace oomph
