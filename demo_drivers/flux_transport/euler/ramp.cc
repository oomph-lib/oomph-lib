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
//Demonstration code to solve a one-dimensional Poisson equation
//using the OOMPH-LIB base classes only

//Include std C++ IO library
#include <iostream>
#include <iomanip>

//Include the OOMPH-LIB base classes
#include "generic.h"
#include "flux_transport.h"

using namespace oomph;
using namespace std;


namespace Global
{
 Vector<double> x0(2,0.0);
 double Gamma = 1.4;
 double Beta = 5.0;
 
 const double pi = MathematicalConstants::Pi;

 /// Function that determines the initial conditions
 void exact_solution(const double &t, const Vector<double> &x,
                     Vector<double> &u)
 {
  double exp_term  = std::exp(1.0 - 
                              (x[0] - t - x0[0])*
                              (x[0] - t - x0[0]) -
                              (x[1] - x0[1])*
                              (x[1] - x0[1]));
  
  //The density
  u[0] = pow(
   (1.0 - ((Gamma-1.0)/(16*Gamma*pi*pi))*Beta*Beta*exp_term*exp_term),
   (1.0/(Gamma -1.0)));

  //The velocities
  u[2] = 1.0 - Beta*exp_term*(x[1] - x0[1])/(2.0*pi);
  u[3] = Beta*exp_term*(x[0] - x0[0])/(2.0*pi);

  //The pressure
  double pres = pow(u[0],Gamma);

  //Now calculate the energy
  u[1] = pres/(Gamma-1.0) + 0.5*u[0]*(u[2]*u[2] + u[3]*u[3]);
  
  //Finally convery the velocities to momenta
  u[2] *= u[0];
  u[3] *= u[0];
 }
}

//Create a quadmesh of DG elements
template<class ELEMENT>
class TwoDDGMesh : public DGMesh
{
 //Map that will store the neighbours
 std::map<std::pair<FiniteElement*,unsigned>,FiniteElement*> Neighbour_map;

 //map of face_reflection elements
 std::map<std::pair<FiniteElement*,int>,FaceElement*> Face_element_pt;

public:
 //Constructor
 TwoDDGMesh(const unsigned &Nx, const unsigned &Ny, 
          TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
  {
   Vector<double> s_fraction;
   //Lengths of the mesh
   double Lx = 1.0;
   double Ly = 1.0;

   double llx = -0.5;
   double lly = 0.0;
   //Angle (in degrees of ramp)
   double theta = 10.0;
   const double pi = MathematicalConstants::Pi;
   //Work out the length of each element in the x and y directions
   //(Assuming uniform spacing)
   double el_length[2] = {Lx/(double)Nx, Ly/(double)Ny};

   //Vector of booleans for boundary information
   std::vector<bool> boundary_info;

   //loop over the elements in x
   for(unsigned ex=0;ex<Nx;ex++)
    {
     //loop over the element in y
     for(unsigned ey=0;ey<Ny;ey++)
      {
       //Create a new DG element
       ELEMENT* local_element_pt = new ELEMENT;
       //Find the number of nodes
       const unsigned n_node = local_element_pt->nnode();
       //Have we constructed
       bool constructed = false;
       
       //If we are on a boundary
       if(ex==0) 
        {
         const unsigned n_p  = local_element_pt->nnode_1d();
         boundary_info.resize(n_node);
         for(unsigned n=0;n<n_node;n++) {boundary_info[n] = false;}
         //The left hand side is on a boundary
         for(unsigned n=0;n<n_p;n++) 
          {
           boundary_info[n*n_p] = true;
          }
         
         //Lower left
         if(ey==0)
          {
          for(unsigned n=1;n<n_p;n++) {boundary_info[n] = true;}
          }
         
         //Upper left
         if(ey==Ny-1)
          {
           for(unsigned n=1;n<n_p;n++) {boundary_info[n_p*(n_p-1) + n] = true;}
          }

         //Now construct
         local_element_pt->construct_boundary_nodes_and_faces(
          this, boundary_info, time_stepper_pt);

         constructed = true;
        }
       //Right boundary
       else if(ex==Nx-1)
        {
         const unsigned n_p  = local_element_pt->nnode_1d();
         boundary_info.resize(n_node);
         for(unsigned n=0;n<n_node;n++) {boundary_info[n] = false;}
         //The right-hand side is on a boundary
         for(unsigned n=0;n<n_p;n++) 
          {
           boundary_info[n_p-1 + n*n_p] = true;
          }
         
         //Lower right
         if(ey==0)
          {
          for(unsigned n=0;n<n_p-1;n++) {boundary_info[n] = true;}
          }
         
         //Upper right
         if(ey==Ny-1)
          {
           for(unsigned n=0;n<n_p-1;n++) 
            {boundary_info[n_p*(n_p-1) + n] = true;}
          }

         //Now construct
         local_element_pt->construct_boundary_nodes_and_faces(
          this, boundary_info, time_stepper_pt);
         
         constructed = true;
        }
       //Otherwise it's in the middle
       else
        {
         //If on the bottom or top
         if((ey==0) || (ey==Ny-1))
          {
           const unsigned n_p  = local_element_pt->nnode_1d();
           boundary_info.resize(n_node);
           for(unsigned n=0;n<n_node;n++) {boundary_info[n] = false;}

           //If the bottom
           if(ey==0)
            {
             for(unsigned n=0;n<n_p;n++) 
              {
               boundary_info[n] = true;
              }
            }
           
           //If on the top
           if(ey==Ny-1)
            {
             for(unsigned n=0;n<n_p;n++) 
              {boundary_info[n_p*(n_p-1) + n] = true;}
            }
           

           //Now construct
           local_element_pt->construct_boundary_nodes_and_faces(
            this, boundary_info, time_stepper_pt);
           
           constructed = true;
          }
        }
       
       //If we haven't already constructed
       if(!constructed)
        {
         //Construct nodes and faces
         local_element_pt->construct_nodes_and_faces(this,time_stepper_pt);
        }

       //Find the lower left corner of the element
       double ll_corner[2] = {llx + ex*el_length[0],lly + ey*el_length[1]};
       //For each of the nodes set the position
       for(unsigned n=0;n<n_node;n++)
        {
         //Get pointer to the node
         Node* nod_pt = local_element_pt->node_pt(n);
         //Get the relative position in local coordinates
         local_element_pt->local_fraction_of_node(n,s_fraction);
         //Loop over the coordinates and set the position
         for(unsigned i=0;i<2;i++)
          {
           nod_pt->x(i) = ll_corner[i] + el_length[i]*s_fraction[i];
          }
         
         //Add the ramp part
         if(nod_pt->x(0) >= 0.0)
          {
           //Work out position of ramp
           double y_min = nod_pt->x(0)*tan(theta*pi/180.0);
           //Now rescale the height
           double frac = nod_pt->x(1)*(Ly-y_min)/Ly;
           //And add it
           nod_pt->x(1) = y_min + frac;
          }

         //Add each node to the node list
         Node_pt.push_back(nod_pt);
        }
       //Now add the element to the list
       Element_pt.push_back(local_element_pt);
      }
    }

   //Set up the boundary nodes
   this->set_nboundary(4);

   //Now let's add the boundary nodes
   //Left-hand and right-hand sides
   for(unsigned e=0;e<Ny;e++)
    {
     //Left
     FiniteElement* elem_pt = this->finite_element_pt(e);
     unsigned n_p = elem_pt->nnode_1d();
     for(unsigned n=0;n<n_p;n++)
      {
       this->add_boundary_node(3,elem_pt->node_pt(n_p*n));
      }
     //Right
     elem_pt = this->finite_element_pt(Ny*(Nx-1)+e);
     n_p = elem_pt->nnode_1d();
     for(unsigned n=0;n<n_p;n++)
      {
       this->add_boundary_node(1,elem_pt->node_pt(n_p-1 + n_p*n));
      }
    }

   //Top and bottom
   for(unsigned e=0;e<Nx;e++)
    {
     //Bottom
     FiniteElement* elem_pt = this->finite_element_pt(Ny*e);
     unsigned n_p = elem_pt->nnode_1d();
     for(unsigned n=0;n<n_p;n++)
      {
       this->add_boundary_node(0,elem_pt->node_pt(n));
      }
     
     //Now add the reflecting element for the south face
     Face_element_pt[std::make_pair(elem_pt,2)] = 
      new DGEulerFaceReflectionElement<ELEMENT>(elem_pt,-2);

     //Top
     elem_pt = this->finite_element_pt(Ny-1 + Ny*e);
     n_p = elem_pt->nnode_1d();
     for(unsigned n=0;n<n_p;n++)
      {
       this->add_boundary_node(2,elem_pt->node_pt((n_p-1)*n_p + n));
      }
    }

   //Now loop over all the elements and set up the neighbour map
   for(unsigned ex=0;ex<Nx;ex++)
    {
     for(unsigned ey=0;ey<Ny;ey++)
      {
       //Get pointer to the element
       unsigned element_index = ex*Ny + ey;
       FiniteElement* local_el_pt = finite_element_pt(element_index);
       
       //Storage for indices of neighbours
       int index[4];

       //North neighbour (just one up)
       if(ey < (Ny-1)) {index[0] = element_index + 1;}
       //If at top, no neighbour
       else {index[0] = -1;}

       //East neighbour (just one across)
       if(ex < (Nx-1)) {index[1] = element_index + Ny;}
       //If at side no neighbour
       else {index[1] = -1;}

       //South neighbour 
       if(ey > 0) {index[2] = element_index - 1;}
       else {index[2] = -1;}

       //West neighbour
       if(ex > 0) {index[3] = element_index - Ny;}
       else {index[3] = -1;}

       //Now store the details in the mape
       for(unsigned i=0;i<4;i++)
        {
         if(index[i]!=-1)
          {
           Neighbour_map[std::make_pair(local_el_pt,i)] =
            finite_element_pt(index[i]);
          }
        }
       }
    }
   
   //Now set up the neighbour information
   /*const unsigned n_element = this->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     dynamic_cast<ELEMENT*>(this->element_pt(e))->setup_face_neighbour_info();
     }*/
   this->setup_face_neighbour_info();
  }

 //We can just use the map here
 void neighbour_finder(FiniteElement* const &bulk_element_pt,
                       const int &face_index,
                       const Vector<double> &s_bulk, 
                       FaceElement* &face_element_pt,
                       Vector<double> &s_face)
  {
   //We have a single face coordinate of size 1
   s_face.resize(1);
   
   ELEMENT* elem_pt = 0;

   //Now things differ depending upon which face we are located
   switch(face_index)
    {
     //North face
    case 2:    
     //Get the neighbour
     elem_pt =  dynamic_cast<ELEMENT*>(
      Neighbour_map[std::make_pair(bulk_element_pt,0)]);
     
     if(elem_pt!=0)
      {
       //The neighbouring face element is the south face of the neighbour
       face_element_pt = elem_pt->face_element_pt(2);
      }
     //Otherwise free outflow
     else
      {
       face_element_pt = 
        dynamic_cast<ELEMENT*>(bulk_element_pt)->face_element_pt(0);
      }
     
       //Then set the face coordinate
       s_face[0] = s_bulk[0];
       break;

     //East face
    case 1:
     //Get the neighbour
     elem_pt = 
      dynamic_cast<ELEMENT*>(
       Neighbour_map[std::make_pair(bulk_element_pt,1)]);

     if(elem_pt!=0)
      {
       face_element_pt = elem_pt->face_element_pt(3);
      }
     else
      {
       face_element_pt = 
        dynamic_cast<ELEMENT*>(bulk_element_pt)->face_element_pt(1);
      }

     //Then set the face coordinate
     s_face[0] = s_bulk[1];
     break;
     
     //South face
    case -2:
     //Get the neighbour
     elem_pt = 
      dynamic_cast<ELEMENT*>(
       Neighbour_map[std::make_pair(bulk_element_pt,2)]);

     if(elem_pt!=0)
      {
       face_element_pt = elem_pt->face_element_pt(0);
      }
     //Otherwise use reflection element
     else
      {
       face_element_pt = Face_element_pt[std::make_pair(bulk_element_pt,2)];
       if(face_element_pt==0) 
        {
         throw OomphLibError("Something wrong",
                             OOMPH_CURRENT_FUNCTION,
                             OOMPH_EXCEPTION_LOCATION);
                             }
      }

     //Then set the face coordiante
     s_face[0] = s_bulk[0];
     break;

     //West face
    case -1:
     //Get the neighbour
     elem_pt  = 
      dynamic_cast<ELEMENT*>
      (Neighbour_map[std::make_pair(bulk_element_pt,3)]);
     
     if(elem_pt!=0)
      {
       face_element_pt = elem_pt->face_element_pt(1);
      }
     else
      {
       face_element_pt = dynamic_cast<ELEMENT*>(bulk_element_pt)
        ->face_element_pt(3);
      }
     
     s_face[0] = s_bulk[1];
     break;
     
    default:
   throw OomphLibError("Coordinate is on no face or not on a unique face",
                       OOMPH_CURRENT_FUNCTION,
                       OOMPH_EXCEPTION_LOCATION);
   break;
    }
  }
 
};


template<class ELEMENT>
class TwoDDGProblem : public Problem
{
public:

 TwoDDGMesh<ELEMENT>* mesh_pt() 
  {return dynamic_cast<TwoDDGMesh<ELEMENT>*>(Problem::mesh_pt());}
 

 ~TwoDDGProblem()
  {
   delete Problem::mesh_pt();
   delete this->explicit_time_stepper_pt();
  }


 TwoDDGProblem(const unsigned &Nx, const unsigned &Ny)
  {
   this->set_explicit_time_stepper_pt(new RungeKutta<4>);
   this->disable_info_in_newton_solve();

   this->enable_discontinuous_formulation();

   Problem::mesh_pt() = new TwoDDGMesh<ELEMENT>(Nx,Ny);

   unsigned n_element = Problem::mesh_pt()->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     ELEMENT* cast_element_pt =
      dynamic_cast<ELEMENT*>(Problem::mesh_pt()->element_pt(e));
     
     cast_element_pt->gamma_pt() = &Global::Gamma;

    }

   //Pin the inlet boundary
   {
    unsigned b=3;
    {
     const unsigned n_node = mesh_pt()->nboundary_node(b);
     for(unsigned n=0;n<n_node;n++)
      {
       Node* nod_pt = mesh_pt()->boundary_node_pt(b,n);
       const unsigned n_value = nod_pt->nvalue();
       for(unsigned i=0;i<n_value;i++) {nod_pt->pin(i);}
      }
    }
   }

      
   std::cout << "How many " << assign_eqn_numbers() << std::endl;
  }

 void apply_initial_conditions()
  {  
   double M = 2.0;
   const double gamma = Global::Gamma;
   const double E = 1.0/(gamma-1.0) + 0.5*M*M*gamma;
   //Loop over all the nodes in the mesh
   unsigned n_node = mesh_pt()->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     Node* nod_pt = mesh_pt()->node_pt(n);
     //Set the values
     nod_pt->set_value(0,gamma);
     nod_pt->set_value(1,E);
     nod_pt->set_value(2,M*gamma);
     nod_pt->set_value(3,0.0);
    }
  }


 void limit()
  {
   const double eps = 1.0e-14;
   //First stage is to calculate the averages of all elements
   const unsigned n_element = this->mesh_pt()->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     //Calculate the averages
     dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))
      ->allocate_memory_for_averages();
     dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->calculate_averages();
    }

   Vector<double> s(2,0.0), x(2);

   const unsigned n_flux = 4;

   const unsigned n_dim = 2;


   //Store pointers to neighbours
   DenseMatrix<ELEMENT*> bulk_neighbour(n_element,4);
   //Storage for the unknowns
   Vector<double> interpolated_u(n_flux);
   
   //Now we need to actually do the limiting, but let's check things first
   for(unsigned e=0;e<n_element;e++)
    {
     //Store the element
     ELEMENT* cast_element_pt = 
      dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
     //Find the centre of the element
     cast_element_pt->interpolated_x(s,x);

     
     //Storage for the approximation to the central gradient from all four
     //faces
     RankThreeTensor<double> grad(4,n_flux,2);
     //Storage for the neighbours size
     Vector<double> size(4);

     //Get the local coordinate of the nodes
     Vector<double> s_face(1), s_bulk(2);
     //Get the central value
     DenseMatrix<double> prim_centre(2,n_flux);
     //My Centre
     for(unsigned i=0;i<n_flux;i++) 
      {interpolated_u[i] = cast_element_pt->average_value(i);}
     prim_centre(0,0) = interpolated_u[0];
     prim_centre(0,1) = cast_element_pt->pressure(interpolated_u);
     prim_centre(0,2) = interpolated_u[2]/interpolated_u[0];
     prim_centre(0,3) = interpolated_u[3]/interpolated_u[0];

     //Now store the average primitivee values 
     for(unsigned i=0;i<n_flux;i++)
      {cast_element_pt->average_prim_value(i) = prim_centre(0,i);}

     //Loop over all the faces
     for(unsigned f=0;f<4;f++)
      {
       DGFaceElement* face_element_pt = 
        dynamic_cast<DGFaceElement*>(cast_element_pt->face_element_pt(f));
       
       //Get the average of the nodal points on the face 
       //(***NOT GENERAL**)
       const unsigned n_face_node = face_element_pt->nnode();

       //Storage for the average nodal values
       DenseMatrix<double> face_average(n_face_node,n_flux);
       
       //Storage for the neighbouring cell's centre average
       DenseMatrix<double> cell_average(n_face_node,n_flux);

       //Storage for the sizes of the neighbouring bulk elements
       Vector<double> neighbour_size(n_face_node);
   
       //Now the storage for the nodal values and centre values
       DenseMatrix<double> nodal_x(n_face_node,n_dim);
       //Now the storage for the centre values
       DenseMatrix<double> centre_x(n_face_node,n_dim);
                                           
       //Neighbours face
       FaceElement* neighbour_face_pt=0;

       //Now calculate those average values
       for(unsigned n=0;n<n_face_node;n++)
        {
         //Get the local coordinates of the node
         face_element_pt->local_coordinate_of_node(n,s_face);
         //Get the global coordinate of the node
         for(unsigned i=0;i<n_dim;i++)
          {
           nodal_x(n,i) = face_element_pt->interpolated_x(s_face,i);
          }

         //Get the flux at the node in this face
         face_element_pt->interpolated_u(s_face,interpolated_u);

         /*std::cout << "Values at face :";
          for(unsigned i=0;i<n_flux;i++)
           {
            std::cout << interpolated_u[i] << " ";
           }
           std::cout << std::endl;*/

         for(unsigned i=0;i<n_flux;i++)
          {
           face_average(n,i) = interpolated_u[i];
          }

         //Now get the bulk coordinate
         face_element_pt->get_local_coordinate_in_bulk(s_face,s_bulk);
         //Now get the neighbour
         cast_element_pt->
          get_neighbouring_face_and_local_coordinate(
           face_element_pt->face_index(),s_bulk,neighbour_face_pt,s_face);
         
         //Get the fluxes of the neighbour
         dynamic_cast<DGFaceElement*>(neighbour_face_pt)
          ->interpolated_u(s_face,interpolated_u);
         //Add the flux from the neighbour
         ELEMENT* neighbour_bulk_element_pt = dynamic_cast<ELEMENT*>(
          neighbour_face_pt->bulk_element_pt());

         /*std::cout << "Values at face neighbour :";
         for(unsigned i=0;i<n_flux;i++)
          {
           std::cout << interpolated_u[i] << " ";
          }
         std::cout << std::endl;
          
         std::cout << "SANITY CHECK " << face_element_pt << " " 
                   << neighbour_face_pt << " "
                   << cast_element_pt << " " 
                   << neighbour_bulk_element_pt << "\n";*/
         

         //Get the centre of the neighbouring bulk element
         for(unsigned i=0;i<n_dim;i++) {s_bulk[i] = 0.0;}
         //Do it
         for(unsigned i=0;i<n_dim;i++)
          {
           centre_x(n,i) = neighbour_bulk_element_pt
            ->interpolated_x(s_bulk,i);
          }
           
         for(unsigned i=0;i<n_flux;i++)
          {
           face_average(n,i) += interpolated_u[i];
           face_average(n,i) *= 0.5;
           //Get the neighbour
           cell_average(n,i) = neighbour_bulk_element_pt->average_value(i);
          }
         
         //Get the size of the neighbour
         neighbour_size[n] = neighbour_bulk_element_pt->average_value(n_flux);

         //Only bother to store the first one
         if(n==0)
          {
           bulk_neighbour(e,f) = neighbour_bulk_element_pt;
          }
        }
       
       //Now we have the face averages and the cell averages
       //Check that we have a conforming mesh
       for(unsigned n=1;n<n_face_node;n++)
        {
         if((std::abs(centre_x(0,0) - centre_x(n,0)) > eps) || 
            (std::abs(centre_x(0,1) - centre_x(n,1)) > eps))
         {
          throw OomphLibError("Mesh is not conforming in limiter",
                              OOMPH_CURRENT_FUNCTION,
                              OOMPH_EXCEPTION_LOCATION);
         }
        }
       

       //Set centre_x(0,0) to be the current element
       for(unsigned i=0;i<n_dim;i++)
        {
         centre_x(0,i) = x[i];
        }

       //Set the size of the neighbour
       size[f] = neighbour_size[0];

       //Convert our values to primitive variables
       DenseMatrix<double> prim_nodal(2,n_flux);

       //Face averages (corner nodes)
       //std::cout << "Nodal values \n";
       for(unsigned n=0;n<2;n++)
        {
         //Load the fluxes
         for(unsigned i=0;i<n_flux;i++) 
          {interpolated_u[i] = face_average(n,i);}
         
         //Output the fluxes
         /*std::cout << n << " : ";
         for(unsigned i=0;i<n_flux;i++)
          {std::cout << interpolated_u[i] << " ";}
          std::cout << std::endl;*/

         prim_nodal(n,0) = interpolated_u[0];
         prim_nodal(n,1) = cast_element_pt->pressure(interpolated_u);
         prim_nodal(n,2) = interpolated_u[2]/interpolated_u[0];
         prim_nodal(n,3) = interpolated_u[3]/interpolated_u[0];
        }

       //Neighbours centre (use value from first node
       for(unsigned i=0;i<n_flux;i++) 
        {interpolated_u[i] = cell_average(0,i);}

       /*std::cout << "Neighbour centre : ";
       for(unsigned i=0;i<n_flux;i++)
        {std::cout << interpolated_u[i] << " ";}
        std::cout << std::endl;*/

       prim_centre(1,0) = interpolated_u[0];
       //Should really use neighbour here
       prim_centre(1,1) = cast_element_pt->pressure(interpolated_u);
       prim_centre(1,2) = interpolated_u[2]/interpolated_u[0];
       prim_centre(1,3) = interpolated_u[3]/interpolated_u[0];


       //We can now estimate the gradients

      
       //Output first
       /*std::cout << "Centre " << x[0] << " " << x[1] << " : ";
       for(unsigned i=0;i<n_flux;i++)
        {std::cout << prim_centre(0,i) << " ";}
       
       std::cout << "\n Neighbour " << centre_x(0,0) << " " 
                 << centre_x(0,1) << " : ";
       for(unsigned i=0;i<n_flux;i++)
        {std::cout << prim_centre(1,i) << " ";}

       for(unsigned n=0;n<2;n++)
        {
         std::cout << "\n Node " << nodal_x(n,0) << " " 
                   << nodal_x(n,1) << " : ";
         for(unsigned i=0;i<n_flux;i++)
          {std::cout << prim_nodal(n,i) << " ";}
        }
        std::cout << std::endl;*/


       double x_len = centre_x(1,0) - centre_x(0,0);
       double y_len = centre_x(1,1) - centre_x(0,1);

       double X_len = nodal_x(1,0) - nodal_x(0,0);
       double Y_len = nodal_x(1,1) - nodal_x(0,1);

       //If we have ghost cell's handle it!
       if((std::abs(x_len) < eps) && (std::abs(y_len) < eps))
        {
         //The distance is twice the distance from the centre to
         //the midpoint of the edges
         double X_mid = nodal_x(0,0) + 0.5*X_len;
         double Y_mid = nodal_x(0,1) + 0.5*Y_len;

         x_len = 2.0*(X_mid - centre_x(0,0));
         y_len = 2.0*(Y_mid - centre_x(0,1));
        }

       double A = X_len*y_len + x_len*Y_len;

       for(unsigned i=0;i<n_flux;i++)
        {
         grad(f,i,0) = 
          ((prim_centre(1,i) - prim_centre(0,i))*Y_len
           - (prim_nodal(1,i) - prim_nodal(0,i))*y_len)/A;
         grad(f,i,1) =
          ((prim_centre(1,i) - prim_centre(0,i))*X_len
           - (prim_nodal(1,i) - prim_nodal(0,i))*x_len)/A;
        }
      }

     //Let's output the face gradients
     /*for(unsigned f=0;f<4;f++)
      {
       std::cout << "Neighbour " << size[f] << " : ";
       for(unsigned i=0;i<n_flux;i++)
        {
         std::cout << grad(f,i,0) << " " << grad(f,i,1) << " ";
        }
       std::cout << std::endl;
       }*/

     //Zero out the gradients
     for(unsigned i=0;i<n_flux;i++)
      {
       for(unsigned j=0;j<2;j++)
        {
         cast_element_pt->average_gradient(i,j) = 0.0;
        }
      }

     double sum = 0.0;
     for(unsigned f=0;f<4;f++)
      {
       const double A = size[f];
       sum += A;
       for(unsigned i=0;i<n_flux;i++)
        {
         for(unsigned j=0;j<2;j++)
          {
           cast_element_pt->average_gradient(i,j) += A*grad(f,i,j);
          }
        }
      }

     //Divide by the combined sum
     for(unsigned i=0;i<n_flux;i++)
      {
       for(unsigned j=0;j<2;j++)
        {
         cast_element_pt->average_gradient(i,j) /= sum;
        }
      }
    }

   //Let's have a look at it
   /*for(unsigned e=0;e<n_element;e++)
    {
     std::cout << "Approximate gradient in element " << e << " : ";
     for(unsigned i=0;i<n_flux;i++)
      {
       for(unsigned j=0;j<2;j++)
        {
         std::cout << dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
          average_gradient(i,j) << " " ;
        }
      }
     std::cout << "\n";
     }*/

   //OK now we use the weighted gradients to construct our limited gradient
   const double eps2 = 1.0e-10;

   for(unsigned e=0;e<n_element;e++)
    {
     ELEMENT* cast_element_pt = 
      dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
     
     /*std::cout << "Base element \n";
       cast_element_pt->output(std::cout);*/

     DenseMatrix<double> g(n_flux,4,0.0);

     //Get the neighbours
     for(unsigned f=0;f<4;f++)
      {
       ELEMENT* neighbour_bulk_element_pt = bulk_neighbour(e,f);

       /*std::cout << "Neighbour face" << f << "\n";
         neighbour_bulk_element_pt->output(std::cout);*/


       //Now sort it out
       for(unsigned i=0;i<n_flux;i++)
        {
         for(unsigned j=0;j<2;j++)
          {
           g(i,f) += 
            std::pow(neighbour_bulk_element_pt->average_gradient(i,j),2.0);
          }
        }
      }

     DenseMatrix<double> limited_gradient(n_flux,2,0.0);
     
     //Let's now see it
     for(unsigned i=0;i<n_flux;i++)
      {

       //Calcaulte the denominators
       double denom = 4.0*eps2;
       for(unsigned f=0;f<4;f++) {denom += g(i,f)*g(i,f)*g(i,f);}

       //Calculate the weights
       Vector<double> w(4,0.0);
       w[0] = (g(i,1)*g(i,2)*g(i,3) + eps2)/denom;
       w[1] = (g(i,0)*g(i,2)*g(i,3) + eps2)/denom;
       w[2] = (g(i,0)*g(i,1)*g(i,3) + eps2)/denom;
       w[3] = (g(i,0)*g(i,1)*g(i,2) + eps2)/denom;
            
       //Now reconstruct the limited gradient
       for(unsigned f=0;f<4;f++)
        {
         for(unsigned j=0;j<2;j++)
          {
           limited_gradient(i,j) += 
            w[f]*bulk_neighbour(e,f)->average_gradient(i,j);
          }
        }
      }

     //std::cout << "In element " << e << " :\n ";
     //We should now have the limited gradients for all fluxe
     /*for(unsigned i=0;i<n_flux;i++)
      {
       std::cout << cast_element_pt->average_gradient(i,0) << " "
                 << cast_element_pt->average_gradient(i,1) << " : ";
       std::cout << limited_gradient(i,0) << " "
                 << limited_gradient(i,1) << "\n";
      }
      std::cout << "\n";*/

     //Now we have limited the primitive variables we have to evaluate the 
     //values at each node
     const unsigned n_node = cast_element_pt->nnode();
     //Find the centre of the element
     cast_element_pt->interpolated_x(s,x);
     

     //Work out distances to the nodes
     DenseMatrix<double> dx(n_node,2);
     for(unsigned n=0;n<n_node;n++)
      {
       Node* nod_pt = cast_element_pt->node_pt(n);
       for(unsigned i=0;i<n_dim;i++) {dx(n,i) = nod_pt->x(i) - x[i];}
      }
     
     //New gradients of primitive values
     DenseMatrix<double> dprim(n_node,n_flux);
     

     //Reconstruct limited primitive value gradients
     for(unsigned n=0;n<n_node;n++)
      {
       for(unsigned i=0;i<n_flux;i++)
        {
         dprim(n,i) = dx(n,0)*limited_gradient(i,0) 
          + dx(n,1)*limited_gradient(i,1);
        }
      }

     //Limited values
     DenseMatrix<double> limited_value(n_node,n_flux);
     //Density
     for(unsigned n=0;n<n_node;n++)
      {
       limited_value(n,0) = 
        cast_element_pt->average_prim_value(0) + dprim(n,0);
      }
     
     //If any are too small correct
     double min = 1.0e-2;
     bool loop_flag=false;
     do
      {
       bool limit_more=false;
       for(unsigned n=0;n<n_node;n++)
        {
         //If we're too small, reduce the gradient
         if(limited_value(n,0) < min)
          {
           limit_more=true;
           break;
          }
        }

       //Now do we limit more
       if(limit_more)
        {
         std::cout << "Limiting density gradient further \n";
         limited_gradient(0,0) *= 0.5;
         limited_gradient(0,1) *= 0.5;
         //Calculated the densities again
         for(unsigned n=0;n<n_node;n++)
          {
           limited_value(n,0) = cast_element_pt->average_prim_value(0)
            + dx(n,0)*limited_gradient(0,0) 
            + dx(n,1)*limited_gradient(0,1);
          }
         //Loop again
         loop_flag = true;
        }
      } while(loop_flag);

     //Reconstruct the momentum
     for(unsigned n=0;n<n_node;n++)
      {
       for(unsigned i=0;i<2;i++)
        {
         limited_value(n,2+i) = 
          cast_element_pt->average_value(2+i)
          + cast_element_pt->average_prim_value(0)*dprim(n,2+i)
          + cast_element_pt->average_prim_value(2+i)*dprim(n,0);
        }
      }
           
     //Now reconstruct the energy
     const double gamma = cast_element_pt->gamma();
     
     for(unsigned n=0;n<n_node;n++)
      {
       limited_value(n,1) = cast_element_pt->average_value(1)
        + (1.0/(gamma-1.0))*dprim(n,1)
        + 0.5*dprim(n,0)*(cast_element_pt->average_prim_value(2)*
                          cast_element_pt->average_prim_value(2) +
                          cast_element_pt->average_prim_value(3)*
                          cast_element_pt->average_prim_value(3))
        + cast_element_pt->average_prim_value(0)*(
         cast_element_pt->average_prim_value(2)*dprim(2,n) +
         cast_element_pt->average_prim_value(3)*dprim(3,n));
      }

     //Check for negative pressures
     bool negative_pressure = false;
     for(unsigned n=0;n<n_node;n++)
      {
       for(unsigned i=0;i<n_flux;i++) 
        {interpolated_u[i] = limited_value(n,i);}

       //Get the pressure 
       double p = cast_element_pt->pressure(interpolated_u);
       if(p < eps2) {negative_pressure = true; break;}
      }
     
     //Correct negative pressure by limiting energy gradient to zero
     if(negative_pressure)
      {
       std::cout << "Correcting negative pressure\n";
       //If the average value is zero, then set the average value from
       //the average velocities
       double average_E = cast_element_pt->average_value(1);
       double rho_u = cast_element_pt->average_value(2);
       double rho_v = cast_element_pt->average_value(3);
       double rho = cast_element_pt->average_value(0);

       //Find the average kinetic energy
       double kin = 0.5*(rho_u*rho_u + rho_v*rho_v)/rho;

       //If the average energy is less than the average kinetic energy
       //we will have negative pressures (we don't want this!)
       if(average_E < kin)
        {
         if(average_E < 0) 
          {
           throw OomphLibError("Real problem energy negative\n",
                               OOMPH_CURRENT_FUNCTION,
                               OOMPH_EXCEPTION_LOCATION);
          }
         
         //Keep the ratio of the momentum transport the same
         double ratio = rho_u/rho_v;
         
         //Find the sign of rho_v
         int sign = 1;
         if(rho_v < 0.0) {sign *= -1;} 

         double sum_square_momenta = 2.0*average_E*rho;

         //Reset the values preserving the sign of v
         rho_v = sign*sqrt(sum_square_momenta/(1.0 + ratio*ratio));
         rho_u = ratio*rho_v; 
        }
       
       for(unsigned n=0;n<n_node;n++)
        {
         //Set consistent averages
         limited_value(n,0) = rho;
         limited_value(n,1) = average_E;
         limited_value(n,2) = rho_u;
         limited_value(n,3) = rho_v;
        }
      }

     //Will have to think about boundary conditions

     //Finally set the limited values
     for(unsigned n=0;n<n_node;n++)
      {
       //Cache the node
       Node* nod_pt = cast_element_pt->node_pt(n);
       for(unsigned i=0;i<n_flux;i++)
        {
         //Only limit if it's not pinned!
         if(!nod_pt->is_pinned(i))
          {
           /*std::cout << "Limiting " <<
            nod_pt->value(i) << "  to "  <<
            limited_value(n,i) << "\n";*/
           //Do it
           nod_pt->set_value(i,limited_value(n,i));
          }
        }
      }
    }

  }



 void parameter_study(std::ostream &trace, const bool &disc)
  {
   this->enable_mass_matrix_reuse();
   double dt = 0.0001;
   apply_initial_conditions();
   
   Vector<double> error(4,0.0);
   char filename[100];
   
   unsigned count=1;
   for(unsigned t=0;t<5000;t++)
    {
     explicit_timestep(dt);
     //Now limit
     //limit();
     if(count==50)
      {
       if(disc)
        {
         sprintf(filename,"disc_%li_time%g.dat",mesh_pt()->nelement(),time());
        }
       else
        {
         sprintf(filename,"cont_%li_time%g.dat",mesh_pt()->nelement(),time());
        }
       std::ofstream outfile(filename);
       mesh_pt()->output(outfile);
       count=0;
      }

     trace << mesh_pt()->nelement() << " " 
           << this->time() << std::endl;
     ++count;
    }
  }
};


//----------------------------MAIN FUNCTION-------------------------------

int main()
{
 ofstream trace("trace_disc.dat");
  
 unsigned n_element = 20;
 
 //for(unsigned i=0;i<3;i++)
  {
   TwoDDGProblem<DGSpectralEulerElement<2,2> > 
    problem(n_element,n_element);
   problem.parameter_study(trace,true);
   n_element *= 2;
  }
 
 trace.close();

 return 1;
}
