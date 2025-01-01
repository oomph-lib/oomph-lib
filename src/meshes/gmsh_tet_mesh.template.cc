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

#ifndef OOMPH_GMSH_TET_MESH_TEMPLATE_CC
#define OOMPH_GMSH_TET_MESH_TEMPLATE_CC


#include "gmsh_tet_mesh.template.h"


namespace oomph
{
  //======================================================================
  /// Adapt problem based on specified elemental error estimates
  //======================================================================
  template<class ELEMENT>
  void RefineableGmshTetMesh<ELEMENT>::adapt(const Vector<double>& elem_error)
  {
    double t_start = 0.0;

    //###################################
    t_start = TimingHelpers::timer();
    //###################################

    // Get refinement targets
    Vector<double> target_size(elem_error.size());
    double max_edge_ratio =
      this->compute_volume_target(elem_error, target_size);
    // Get maximum target volume
    unsigned n = target_size.size();
    double max_size = 0.0;
    double min_size = DBL_MAX;
    for (unsigned e = 0; e < n; e++)
    {
      if (target_size[e] > max_size) max_size = target_size[e];
      if (target_size[e] < min_size) min_size = target_size[e];
    }

    oomph_info << "Maximum target size: " << max_size << std::endl;
    oomph_info << "Minimum target size: " << min_size << std::endl;
    oomph_info << "Number of elements to be refined " << this->Nrefined
               << std::endl;
    oomph_info << "Number of elements to be unrefined " << this->Nunrefined
               << std::endl;
    oomph_info << "Max edge ratio " << max_edge_ratio << std::endl;

    double orig_max_size, orig_min_size;
    this->max_and_min_element_size(orig_max_size, orig_min_size);
    oomph_info << "Max/min element size in original mesh: " << orig_max_size
               << " " << orig_min_size << std::endl;

    //##################################################################
    oomph_info
      << "adapt: Time for getting volume targets                      : "
      << TimingHelpers::timer() - t_start << " sec " << std::endl;
    //##################################################################

    // Should we bother to adapt?
    if ((Nrefined > 0) || (Nunrefined > this->max_keep_unrefined()) ||
        (max_edge_ratio > this->max_permitted_edge_ratio()))
    {
      if (!((Nrefined > 0) || (Nunrefined > max_keep_unrefined())))
      {
        oomph_info << "Mesh regeneration triggered by edge ratio criterion\n";
      }


      // Are we dealing with a solid mesh?
      SolidMesh* solid_mesh_pt = dynamic_cast<SolidMesh*>(this);


      // If the mesh is a solid mesh then do the mapping based on the
      // Eulerian coordinates
      bool use_eulerian_coords = false;
      if (solid_mesh_pt != 0)
      {
        use_eulerian_coords = true;
      }

      MeshAsGeomObject* mesh_geom_obj_pt = 0;

#ifdef OOMPH_HAS_CGAL

      // Make cgal-based bin
      CGALSamplePointContainerParameters cgal_params(this);
      if (use_eulerian_coords)
      {
        cgal_params.enable_use_eulerian_coordinates_during_setup();
      }
      mesh_geom_obj_pt = new MeshAsGeomObject(&cgal_params);

#else

      std::ostringstream error_message;
      error_message << "Non-CGAL-based target size transfer not implemented.\n";
      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);

      // Do something here...

#endif

      // Set up a map from pointer to element to its number
      // in the mesh
      std::map<GeneralisedElement*, unsigned> element_number;
      unsigned nelem = this->nelement();
      for (unsigned e = 0; e < nelem; e++)
      {
        element_number[this->element_pt(e)] = e;
      }

      // Get min/max coordinates
      Vector<std::pair<double, double>> min_and_max_coordinates =
        mesh_geom_obj_pt->sample_point_container_pt()
          ->min_and_max_coordinates();

      // Setup grid dimensions for size transfer
      double dx =
        (min_and_max_coordinates[0].second - min_and_max_coordinates[0].first);
      double dy =
        (min_and_max_coordinates[1].second - min_and_max_coordinates[1].first);
      double dz =
        (min_and_max_coordinates[2].second - min_and_max_coordinates[2].first);

      double dx_target =
        this->Gmsh_parameters_pt->dx_for_target_size_transfer();
      unsigned nx = unsigned(dx / dx_target);
      unsigned ny = unsigned(dy / dx_target);
      unsigned nz = unsigned(dz / dx_target);

      dx /= double(nx);
      dy /= double(ny);
      dz /= double(nz);


      // Size transfer via hard disk -- yikes...
      std::string target_size_file_name =
        this->Gmsh_parameters_pt->target_size_file_name();

      std::ofstream target_size_file;
      target_size_file.open(target_size_file_name.c_str());
      target_size_file << min_and_max_coordinates[0].first << " "
                       << min_and_max_coordinates[1].first << " "
                       << min_and_max_coordinates[2].first << " " << std::endl;
      target_size_file << dx << " " << dy << " " << dz << std::endl;
      target_size_file << nx + 1 << " " << ny + 1 << " " << nz + 1 << std::endl;


      // Doc target areas
      int counter =
        this->Gmsh_parameters_pt->counter_for_filename_gmsh_size_transfer();
      std::string stem =
        this->Gmsh_parameters_pt->stem_for_filename_gmsh_size_transfer();
      std::ofstream latest_sizes_file;
      bool doc_target_areas = false;
      if ((stem != "") && (counter >= 0))
      {
        doc_target_areas = true;
        std::string filename =
          stem + oomph::StringConversion::to_string(counter) + ".dat";
        latest_sizes_file.open(filename.c_str());
        latest_sizes_file << "ZONE I=" << nx + 1 << ", J=" << ny + 1
                          << ", K=" << nz + 1 << std::endl;
        this->Gmsh_parameters_pt->counter_for_filename_gmsh_size_transfer()++;
      }


      Vector<double> x(3);
      for (unsigned i = 0; i <= nx; i++)
      {
        x[0] = min_and_max_coordinates[0].first + double(i) * dx;
        for (unsigned j = 0; j <= ny; j++)
        {
          x[1] = min_and_max_coordinates[1].first + double(j) * dy;
          for (unsigned k = 0; k <= nz; k++)
          {
            x[2] = min_and_max_coordinates[2].first + double(k) * dz;

            // Try the specified number of nearest sample points for Newton
            // search then just settle on the nearest one
            Vector<double> s(3);
            GeomObject* geom_obj_pt = 0;
            unsigned max_sample_points =
              this->Gmsh_parameters_pt
                ->max_sample_points_for_limited_locate_zeta_during_target_size_transfer();

#ifdef OOMPH_HAS_CGAL

            dynamic_cast<CGALSamplePointContainer*>(
              mesh_geom_obj_pt->sample_point_container_pt())
              ->limited_locate_zeta(x, max_sample_points, geom_obj_pt, s);

#else

            std::ostringstream error_message;
            error_message
              << "Non-CGAL-based target size transfer not implemented.\n";
            throw OomphLibError(error_message.str(),
                                OOMPH_CURRENT_FUNCTION,
                                OOMPH_EXCEPTION_LOCATION);

            // Do something here...

#endif


#ifdef PARANOID
            if (geom_obj_pt == 0)
            {
              std::stringstream error_message;
              error_message << "Limited locate zeta failed for zeta = [ "
                            << x[0] << " " << x[1] << " " << x[2]
                            << " ]. Makes no sense!\n";
              throw OomphLibError(error_message.str(),
                                  OOMPH_CURRENT_FUNCTION,
                                  OOMPH_EXCEPTION_LOCATION);
            }
            else
            {
#endif

              FiniteElement* fe_pt = dynamic_cast<FiniteElement*>(geom_obj_pt);

#ifdef PARANOID
              if (fe_pt == 0)
              {
                std::stringstream error_message;
                error_message
                  << "Cast to FE for GeomObject returned by limited "
                  << "locate zeta failed for zeta = [ " << x[0] << " " << x[1]
                  << " " << x[2] << " ]. Makes no sense!\n";
                throw OomphLibError(error_message.str(),
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
              }
              else
              {
#endif

                // What's the target size of the element that contains this
                // point
                double tg_size =
                  pow(target_size[element_number[fe_pt]], 1.0 / 3.0);
                target_size_file << tg_size << " ";

                // Doc?
                if (doc_target_areas)
                {
                  latest_sizes_file << x[0] << " " << x[1] << " " << x[2] << " "
                                    << tg_size << std::endl;
                }

#ifdef PARANOID
              }
            }
#endif
          }
          target_size_file << std::endl;
        }
      }
      target_size_file.close();

      if (doc_target_areas)
      {
        latest_sizes_file.close();
      }

      // Build new mesh, reading area information from file
      bool use_mesh_grading_from_file = true;
      RefineableGmshTetMesh<ELEMENT>* new_mesh_pt =
        new RefineableGmshTetMesh<ELEMENT>(this->Gmsh_parameters_pt,
                                           use_mesh_grading_from_file,
                                           this->Time_stepper_pt);

      /* // keep around for debugging */
      /* unsigned nplot=5; */
      /* ofstream vtu_file; */
      /* vtu_file.open("new_mesh.vtu"); */
      /* new_mesh_pt->output_paraview(vtu_file,nplot); */
      /* vtu_file.close(); */

      //###################################
      t_start = TimingHelpers::timer();
      //###################################

      ProjectionProblem<ELEMENT>* project_problem_pt = 0;

      // Check that the projection step is not disabled
      if (!this->Gmsh_parameters_pt->projection_is_disabled())
      {
        // Project current solution onto new mesh
        //---------------------------------------
        project_problem_pt = new ProjectionProblem<ELEMENT>;
        project_problem_pt->mesh_pt() = new_mesh_pt;
        project_problem_pt->project(this);

        oomph_info
          << "adapt: Time for project soln onto new mesh                : "
          << TimingHelpers::timer() - t_start << " sec " << std::endl;
      }
      //###################################
      t_start = TimingHelpers::timer();
      //###################################

      // Flush the old mesh
      unsigned nnod = nnode();
      for (unsigned j = nnod; j > 0; j--)
      {
        delete Node_pt[j - 1];
        Node_pt[j - 1] = 0;
      }
      unsigned nel = nelement();
      for (unsigned e = nel; e > 0; e--)
      {
        delete Element_pt[e - 1];
        Element_pt[e - 1] = 0;
      }

      // Now copy back to current mesh
      //------------------------------
      nnod = new_mesh_pt->nnode();
      Node_pt.resize(nnod);
      nel = new_mesh_pt->nelement();
      Element_pt.resize(nel);
      for (unsigned j = 0; j < nnod; j++)
      {
        Node_pt[j] = new_mesh_pt->node_pt(j);
      }
      for (unsigned e = 0; e < nel; e++)
      {
        Element_pt[e] = new_mesh_pt->element_pt(e);
      }

      // Copy the boundary schemes
      unsigned nbound = new_mesh_pt->nboundary();
      Boundary_element_pt.resize(nbound);
      Face_index_at_boundary.resize(nbound);
      Boundary_node_pt.resize(nbound);
      for (unsigned b = 0; b < nbound; b++)
      {
        unsigned nel = new_mesh_pt->nboundary_element(b);
        Boundary_element_pt[b].resize(nel);
        Face_index_at_boundary[b].resize(nel);
        for (unsigned e = 0; e < nel; e++)
        {
          Boundary_element_pt[b][e] = new_mesh_pt->boundary_element_pt(b, e);
          Face_index_at_boundary[b][e] =
            new_mesh_pt->face_index_at_boundary(b, e);
        }
        unsigned nnod = new_mesh_pt->nboundary_node(b);
        Boundary_node_pt[b].resize(nnod);
        for (unsigned j = 0; j < nnod; j++)
        {
          Boundary_node_pt[b][j] = new_mesh_pt->boundary_node_pt(b, j);
        }
      }

      // Also copy over the new boundary and region information
      unsigned n_region = new_mesh_pt->nregion();

      // Only bother if we have regions
      if (n_region > 1)
      {
        // Deal with the region information first
        this->Region_element_pt.resize(n_region);
        this->Region_attribute.resize(n_region);
        for (unsigned i = 0; i < n_region; i++)
        {
          // Copy across region attributes (region ids!)
          this->Region_attribute[i] = new_mesh_pt->region_attribute(i);

          // Find the number of elements in the region
          unsigned r = this->Region_attribute[i];
          unsigned n_region_element = new_mesh_pt->nregion_element(r);
          this->Region_element_pt[i].resize(n_region_element);
          for (unsigned e = 0; e < n_region_element; e++)
          {
            this->Region_element_pt[i][e] =
              new_mesh_pt->region_element_pt(r, e);
          }
        }

        // Now the boundary region information
        this->Boundary_region_element_pt.resize(nbound);
        this->Face_index_region_at_boundary.resize(nbound);

        // Now loop over the boundaries
        for (unsigned b = 0; b < nbound; ++b)
        {
          // Loop over the regions
          for (unsigned i = 0; i < n_region; ++i)
          {
            unsigned r = this->Region_attribute[i];

            unsigned n_boundary_el_in_region =
              new_mesh_pt->nboundary_element_in_region(b, r);
            if (n_boundary_el_in_region > 0)
            {
              // Allocate storage in the map
              this->Boundary_region_element_pt[b][r].resize(
                n_boundary_el_in_region);
              this->Face_index_region_at_boundary[b][r].resize(
                n_boundary_el_in_region);

              // Copy over the information
              for (unsigned e = 0; e < n_boundary_el_in_region; ++e)
              {
                this->Boundary_region_element_pt[b][r][e] =
                  new_mesh_pt->boundary_element_in_region_pt(b, r, e);
                this->Face_index_region_at_boundary[b][r][e] =
                  new_mesh_pt->face_index_at_boundary_in_region(b, r, e);
              }
            }
          }
        } // End of loop over boundaries

      } // End of case when more than one region


      // Flush the mesh
      new_mesh_pt->flush_element_and_node_storage();

      // Delete the mesh and the problem
      delete new_mesh_pt;
      delete project_problem_pt;

      //##################################################################
      oomph_info
        << "adapt: Time for moving nodes etc. to actual mesh          : "
        << TimingHelpers::timer() - t_start << " sec " << std::endl;
      //##################################################################


      // Solid mesh?
      if (solid_mesh_pt != 0)
      {
        // Warning
        std::stringstream error_message;
        error_message
          << "Lagrangian coordinates are currently not projected but are\n"
          << "are re-set during adaptation. This is not appropriate for\n"
          << "real solid mechanics problems!\n";
        OomphLibWarning(error_message.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);

        // Reset Lagrangian coordinates
        dynamic_cast<SolidMesh*>(this)->set_lagrangian_nodal_coordinates();
      }

      double max_area;
      double min_area;
      this->max_and_min_element_size(max_area, min_area);
      oomph_info << "Max/min element size in adapted mesh: " << max_area << " "
                 << min_area << std::endl;
    }
    else
    {
      oomph_info << "Not enough benefit in adaptation.\n";
      Nrefined = 0;
      Nunrefined = 0;
    }
  }
} // namespace oomph


#endif
