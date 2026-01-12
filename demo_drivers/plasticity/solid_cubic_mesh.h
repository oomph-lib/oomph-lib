// This file is part of the MercuryDPM project (https://www.mercurydpm.org).
// Copyright (c), The MercuryDPM Developers Team. All rights reserved.
// License: BSD 3-Clause License; see the LICENSE file in the root directory.

#ifndef SOLID_CUBIC_MESH_H
#define SOLID_CUBIC_MESH_H

// #include "generic.h"
#include "solid.h"
#include "meshes/simple_cubic_mesh.h"

namespace oomph
{
    template<class ELEMENT>
    class SolidCubicMesh : public virtual SimpleCubicMesh<ELEMENT>, public virtual SolidMesh
    {
    public:

        SolidCubicMesh(const unsigned &nx, const unsigned &ny,
                       const unsigned &nz,
                       const double &a,
                       const double &b,
                       const double &c,
                       TimeStepper *time_stepper_pt = &Mesh::Default_TimeStepper) :
                       SimpleCubicMesh<ELEMENT>(nx, ny, nz, -a, a, -b, b, -c, c, time_stepper_pt),
                       SolidMesh()
        {
            // Set the boundary coordinates of the nodes:
            Vector<double> zeta(2, 0.0);
            Vector<double> x(3, 0.0);
            
            // Loop over boundaries
            for(unsigned b : {0,1,2,3,4,5})
            {
                // Loop over the nodes on that boundary 
                const unsigned n_node = this->nboundary_node(b);
                for(unsigned l=0; l<n_node; l++)
                {
                    // Safe cast to SolidNode
                    SolidNode* node_pt = dynamic_cast<SolidNode*>(this->boundary_node_pt(b, l));
                    if (!node_pt)
                    {
                         throw OomphLibError("Node on boundary is not a SolidNode",
                                             OOMPH_CURRENT_FUNCTION,
                                             OOMPH_EXCEPTION_LOCATION);
                    }

                    for(unsigned i=0; i<3; i++)
                    {
                        x[i] = node_pt->x(i);
                    }
                    
                    // calculate boundary coordinate of node
                    calculate_boundary_coordinate_of_node(x, b, -a, a, -b, b, -c, c, zeta);
                    
                    // assign boundary coordinate of node
                    node_pt->add_to_boundary(b);
                    node_pt->set_coordinates_on_boundary(b, zeta);
                }
                // We have set local boundary coordinates
                this->set_boundary_coordinate_exists(b);
            }

            // Assign the initial lagrangian coordinates
            this->set_lagrangian_nodal_coordinates();
        }   


        /// \short Constructor:
        // nx, ny, nz: number of elements in the x, y, and z directions
        // xMax, yMax, zMax: dimensions of the cube (assume the center of the cube is at the origin)
        // timeStepper: defaults to Steady.
        SolidCubicMesh(const unsigned& nx, const unsigned& ny, const unsigned& nz,
                        const double& xMin, const double& xMax, const double& yMin,
                        const double& yMax, const double& zMin, const double& zMax,
                        TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper) :
                        SimpleCubicMesh<ELEMENT>(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax, time_stepper_pt),
                        SolidMesh()
        {
            // Set the boundary coordinates of the nodes:
            Vector<double> zeta(2, 0.0);
            // Loop over boundaries
            for(unsigned b : {0,1,2,3,4,5})
            {
                // Loop over the nodes on that boundary
                const unsigned n_node = this->nboundary_node(b);
                for(unsigned l=0; l<n_node; l++)
                {
                SolidNode* node_pt = this->boundary_node_pt(b, l);
                Vector<double> x(3, 0.0);
                for(unsigned i=0; i<3; i++)
                {
                x[i] = node_pt->x(i);
                }
                // calculate boundary coordinate of node
                calculate_boundary_coordinate_of_node(x, b, xMin, xMax, yMin, yMax, zMin, zMax, zeta);
                // assign boundary coordinate of node
                node_pt->add_to_boundary(b);
                node_pt->set_coordinates_on_boundary(b, zeta);
                }
                // We have set local boundary coordinates
                this->set_boundary_coordinate_exists(b);
            }

            //Assign the initial lagrangian coordinates
            set_lagrangian_nodal_coordinates();
        }

        void setup_boundary_element_info()
        {
            SimpleCubicMesh<ELEMENT>::setup_boundary_element_info();
        }

        void setup_boundary_element_info(std::ostream& outfile)
        {
            SimpleCubicMesh<ELEMENT>::setup_boundary_element_info(outfile);
        }

        static double zeta_linear(const double& min, const double& max, const double& val)
        {
            return (2.0 * val - max - min) / (max - min);
        }

    private:
        void calculate_boundary_coordinate_of_node(const Vector<double>& x, const unsigned& b,
                                                   const double& xMin, const double& xMax, 
                                                   const double& yMin, const double& yMax, 
                                                   const double& zMin, const double& zMax,
                                                   Vector<double>& zeta)
        {
            // Standard Oomph SimpleCubicMesh numbering:
            // 0: x=min, 1: x=max, 2: y=min, 3: y=max, 4: z=min, 5: z=max
            switch (b)
            {
            case 0: // Left (x min) -> use y, z
            case 1: // Right (x max) -> use y, z
                zeta[0] = zeta_linear(yMin, yMax, x[1]);
                zeta[1] = zeta_linear(zMin, zMax, x[2]);
                break;
            case 2: // Bottom (y min) -> use x, z
            case 3: // Top (y max) -> use x, z
                zeta[0] = zeta_linear(xMin, xMax, x[0]);
                zeta[1] = zeta_linear(zMin, zMax, x[2]);
                break;
            case 4: // Front (z min) -> use x, y
            case 5: // Back (z max) -> use x, y
                zeta[0] = zeta_linear(xMin, xMax, x[0]);
                zeta[1] = zeta_linear(yMin, yMax, x[1]);
                break;
            default:
                std::string error_string = std::to_string(b);
                error_string += " is not a valid boundary."; 
                throw OomphLibError(error_string,
                                    OOMPH_CURRENT_FUNCTION,
                                    OOMPH_EXCEPTION_LOCATION);
            }
        }
    };

 // Implement the RefineableSolidCubicMesh. Adds boundary coordinates.
 template<class ELEMENT>
 class RefineableSolidCubicMesh : public virtual RefineableSimpleCubicMesh<ELEMENT>, public virtual SolidMesh
 {
 public:
  
  RefineableSolidCubicMesh(const unsigned &nx, const unsigned &ny,
                 const unsigned &nz,
                 const double &a,
                 const double &b,
                 const double &c,
                 TimeStepper *time_stepper_pt = &Mesh::Default_TimeStepper) :
                 SimpleCubicMesh<ELEMENT>(nx, ny, nz, -a, a, -b, b, -c, c, time_stepper_pt),
                 RefineableSimpleCubicMesh<ELEMENT>(nx, ny, nz, -a, a, -b, b, -c, c, time_stepper_pt),
                 SolidMesh()
  {
   // Set the boundary coordinates of the nodes:
   Vector<double> zeta(2, 0.0);
   // Loop over boundaries
   for(unsigned b : {0,1,2,3,4,5})
   {
    // Loop over the nodes on that boundary
    const unsigned n_node = this->nboundary_node(b);
    for(unsigned l=0; l<n_node; l++)
    {
     SolidNode* node_pt = this->boundary_node_pt(b, l);
     Vector<double> x(3, 0.0);
     for(unsigned i=0; i<3; i++)
     {
      x[i] = node_pt->x(i);
     }
     // calculate boundary coordinate of node
     calculate_boundary_coordinate_of_node(x, b, -a, a, -b, b, -c, c, zeta);
     // assign boundary coordinate of node
     node_pt->add_to_boundary(b);
     node_pt->set_coordinates_on_boundary(b, zeta);
    }
    // We have set local boundary coordinates
    this->set_boundary_coordinate_exists(b);
   }
   
   //Assign the initial lagrangian coordinates
   set_lagrangian_nodal_coordinates();
  }

  /// \short Constructor:
  // nx, ny, nz: number of elements in the x, y, and z directions
  // xMax, yMax, zMax: dimensions of the cube (assume the center of the cube is at the origin)
  // timeStepper: defaults to Steady.
  RefineableSolidCubicMesh(const unsigned& nx, const unsigned& ny, const unsigned& nz,
                 const double& xMin, const double& xMax, const double& yMin,
                 const double& yMax, const double& zMin, const double& zMax,
                 TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper) :
                 SimpleCubicMesh<ELEMENT>(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax, time_stepper_pt),
                 RefineableSimpleCubicMesh<ELEMENT>(nx, ny, nz, xMin, xMax, yMin, yMax, zMin, zMax, time_stepper_pt),
                 SolidMesh()
  {
   // Set the boundary coordinates of the nodes:
   Vector<double> zeta(2, 0.0);
   // Loop over boundaries
   for(unsigned b : {0,1,2,3,4,5})
   {
    // Loop over the nodes on that boundary
    const unsigned n_node = this->nboundary_node(b);
    for(unsigned l=0; l<n_node; l++)
    {
     SolidNode* node_pt = this->boundary_node_pt(b, l);
     Vector<double> x(3, 0.0);
     for(unsigned i=0; i<3; i++)
     {
      x[i] = node_pt->x(i);
     }
     // calculate boundary coordinate of node
     calculate_boundary_coordinate_of_node(x, b, xMin, xMax, yMin, yMax, zMin, zMax, zeta);
     // assign boundary coordinate of node
     node_pt->add_to_boundary(b);
     node_pt->set_coordinates_on_boundary(b, zeta);
    }
    // We have set local boundary coordinates
    this->set_boundary_coordinate_exists(b);
   }

   //Assign the initial lagrangian coordinates
   set_lagrangian_nodal_coordinates();
  }

  void setup_boundary_element_info()
  {
   RefineableSimpleCubicMesh<ELEMENT>::setup_boundary_element_info();
  }

  void setup_boundary_element_info(std::ostream& outfile)
  {
   RefineableSimpleCubicMesh<ELEMENT>::setup_boundary_element_info(outfile);
  }

  static double zeta_linear(const double& min, const double& max, const double& x)
  {
   return (2.0*x - max - min)/(max-min);
  }

 private:
  void calculate_boundary_coordinate_of_node(const Vector<double>& x, const unsigned& b,
                                             const double& xMin, const double& xMax, const double& yMin,
                                             const double& yMax, const double& zMin, const double& zMax,
                                             Vector<double>& zeta)
  {
   switch (b)
   {
   case 0:
    zeta[0] = zeta_linear(xMin, xMax, x[0]);
    zeta[1] = zeta_linear(yMin, yMax, x[1]);
    break;
   case 1:
    zeta[0] = zeta_linear(xMin, xMax, x[0]);
    zeta[1] = zeta_linear(zMin, zMax, x[2]);
    break;
   case 2:
    zeta[0] = zeta_linear(yMin, yMax, x[1]);
    zeta[1] = zeta_linear(zMin, zMax, x[2]);
    break;
   case 3:
    zeta[0] = zeta_linear(xMin, xMax, x[0]);
    zeta[1] = zeta_linear(zMin, zMax, x[2]);
    break;
   case 4:
    zeta[0] = zeta_linear(yMin, yMax, x[1]);
    zeta[1] = zeta_linear(zMin, zMax, x[2]);
    break;
   case 5:
    zeta[0] = zeta_linear(xMin, xMax, x[0]);
    zeta[1] = zeta_linear(yMin, yMax, x[1]);
    break;
   default:
    // Change to oomph error
    std::string error_string = std::to_string(b);
    error_string += " is not a valid boundary."; 
    throw OomphLibError(error_string,
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
    break;
   }
  }
 };

} // End namespace

#endif
