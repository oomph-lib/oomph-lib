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
 
 const double pi = 4.0*atan(1.0);
}

//Create a quadmesh of DG elements
template<class ELEMENT>
class TwoDDGMesh : public DGMesh
{
 //Map that will store the neighbours
 std::map<std::pair<FiniteElement*,int>,FiniteElement*> Neighbour_map;

 //Vector of face_reflection elements
 Vector<FaceElement*> Face_element_pt;

public:
 //Constructor
 TwoDDGMesh(const unsigned &n_x, const unsigned &n_y, 
          TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
  {
   Vector<double> s_fraction;
   //Lengths of the mesh
   double Lx = 2.4;
   double Ly = 0.8;

   unsigned Nx = 4*n_x, Ny = 4*n_y;

   //Work out the length of each element in the x and y directions
   //(Assuming uniform spacing)
   //Create the flow bit
   double el_length[2] = {Lx/(double)Nx, Ly/(double)Ny};

   //Vector of booleans for boundary information
   std::vector<bool> boundary_info;

   unsigned Nx1 = n_x;
   unsigned Ny1 = 5*n_y;
   double llx1 = 0.0;
   double lly1 = 0.0;

   //loop over the elements in x
   for(unsigned ex=0;ex<Nx1;ex++)
    {
     //loop over the element in y
     for(unsigned ey=0;ey<Ny1;ey++)
      {
       //Create a new DG element
       ELEMENT* local_element_pt = new ELEMENT;
       //Find the number of nodes
       const unsigned n_node = local_element_pt->nnode();
       //Have we constructed
       bool constructed = false;
       //Left boundary 
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
         if(ey==Ny1-1)
          {
           for(unsigned n=1;n<n_p;n++) {boundary_info[n_p*(n_p-1) + n] = true;}
          }

         //Now construct
         local_element_pt->construct_boundary_nodes_and_faces(
          this, boundary_info, time_stepper_pt);

         constructed = true;
        }
       //Right boundary
       else if(ex==Nx1-1)
        {
         const unsigned n_p  = local_element_pt->nnode_1d();
         boundary_info.resize(n_node);
         for(unsigned n=0;n<n_node;n++) {boundary_info[n] = false;}
         
         //Bottom right
         if(ey==0)
          {
           //The bottom is a boundary
           for(unsigned n=0;n<n_p-1;n++) {boundary_info[n] = true;}
          }

         
         //The right-hand side is on a boundary in the early bit
         //of the step
         if(ey < n_y)
          {
           for(unsigned n=0;n<n_p;n++) 
            {
             boundary_info[n_p-1 + n*n_p] = true;
            }
          }

         //If top right
         if(ey==Ny1-1)
          {
           //The top is a boundary
           for(unsigned n=0;n<n_p;n++) 
            {boundary_info[n_p*(n_p-1) + n] = true;}
          }
         
         local_element_pt->construct_boundary_nodes_and_faces(
            this, boundary_info, time_stepper_pt);
         
         constructed = true;
        }
       //Otherwise it's in the middle
       else
        {
         //If on the bottom or top
         if((ey==0) || (ey==Ny1-1))
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
           if(ey==Ny1-1)
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
       double ll_corner[2] = {llx1 + ex*el_length[0],lly1 + ey*el_length[1]};
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
         
         //Add each node to the node list
         Node_pt.push_back(nod_pt);
        }
       //Now add the element to the list
       Element_pt.push_back(local_element_pt);
      }
    }
 

   //Now do the main bit of the mesh
   double llx = 0.6;
   double lly = 0.2;
  
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
       //Right boundary
       if(ex==Nx-1)
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
 //Left-hand side
 for(unsigned e=0;e<Ny1;e++)
  {
   FiniteElement* const elem_pt = this->finite_element_pt(e);
   unsigned n_p = elem_pt->nnode_1d();
   for(unsigned n=0;n<n_p;n++)
    {
     this->add_boundary_node(3,elem_pt->node_pt(n_p*n));
    }
  }

 //Right hand side
 for(unsigned e=0;e<Ny;e++)
  {
   FiniteElement* const elem_pt = 
    this->finite_element_pt(Ny1*Nx1 + Ny*(Nx-1)+e);
   unsigned n_p = elem_pt->nnode_1d();
   for(unsigned n=0;n<n_p;n++)
    {
     this->add_boundary_node(1,elem_pt->node_pt(n_p-1 + n_p*n));
    }
  }

 
      
 //Top
  {
   //First bit
   for(unsigned e=0;e<Nx1;e++)
    {
     FiniteElement* const elem_pt = this->finite_element_pt(Ny1-1 + Ny1*e);
     unsigned n_p = elem_pt->nnode_1d();
     for(unsigned n=0;n<n_p;n++)
      {
       this->add_boundary_node(2,elem_pt->node_pt((n_p-1)*n_p + n));
      }
    }
   
   //Second bit
   for(unsigned e=0;e<Nx;e++)
    {
     FiniteElement* const elem_pt = 
      this->finite_element_pt(Nx1*Ny1 + Ny-1 + Ny*e);
     unsigned n_p = elem_pt->nnode_1d();
     for(unsigned n=0;n<n_p;n++)
      {
       this->add_boundary_node(2,elem_pt->node_pt((n_p-1)*n_p + n));
      }
    }
  }
 
 //Bottom
  {
   //First bit
   for(unsigned e=0;e<Nx1;e++)
    {
     //Bottom
     FiniteElement* const elem_pt = this->finite_element_pt(Ny1*e);
     unsigned n_p = elem_pt->nnode_1d();
     for(unsigned n=0;n<n_p;n++)
      {
       this->add_boundary_node(0,elem_pt->node_pt(n));
      }
    }
   
   //Add the right-hand of the lower elements
   for(unsigned e=0;e<n_y;e++)
    {
     FiniteElement* const elem_pt = this->finite_element_pt(Ny1*(Nx1-1) + e);
     unsigned n_p = elem_pt->nnode_1d();
     for(unsigned n=0;n<n_p;n++)
      {
       this->add_boundary_node(0,elem_pt->node_pt(n_p-1 + n_p*n));
      }
    }
   
   //Second bit
   for(unsigned e=0;e<Nx;e++)
    {
     //Bottom
     FiniteElement* const elem_pt = 
      this->finite_element_pt(Nx1*Ny1  + Ny*e);
     unsigned n_p = elem_pt->nnode_1d();
     for(unsigned n=0;n<n_p;n++)
      {
       this->add_boundary_node(0,elem_pt->node_pt(n));
      }
    }
  }

  /* for(unsigned b=0;b<4;b++)
    {
     std::cout << "Boundary : " << b << "\n";
     unsigned nbound = this->nboundary_node(b);
     for(unsigned n=0;n<nbound;n++)
      {
       std::cout << "( " << this->boundary_node_pt(b,n)->x(0) << ", "
                 << this->boundary_node_pt(b,n)->x(1) << ")\n";
      }
      }*/

   //Now loop over all the elements and set up the neighbour map
  for(unsigned ex=0;ex<Nx;ex++)
   {
    for(unsigned ey=0;ey<Ny;ey++)
     {
      //Get pointer to the element
      unsigned element_index = Nx1*Ny1  + ex*Ny + ey;
      FiniteElement* local_el_pt = finite_element_pt(element_index);
      
      //Storage for indices of neighbours
      int index[4];
      
      //North neighbour (just one up)
      if(ey < (Ny-1)) {index[0] = element_index + 1;}
      //If at top (no neighbour)
      else {index[0] = -1;}

       //East neighbour (just one across)
       if(ex < (Nx-1)) {index[1] = element_index + Ny;}
       //If at side (no neighbour)
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
          if(index[i] != -1)
           {
            Neighbour_map[std::make_pair(local_el_pt,i)] =
             finite_element_pt(index[i]);
           }
        }
      }
    }
 
  //Now sort out the maps for the first bit
  for(unsigned ex=0;ex<Nx1;ex++)
   {
    for(unsigned ey=0;ey<Ny1;ey++)
     {
      //Get pointer to the element
      unsigned element_index = ex*Ny1 + ey;
      FiniteElement* local_el_pt = finite_element_pt(element_index);
      
      //Storage for indices of neighbours
      int index[4];
      
      //North neighbour (just one up)
      if(ey < (Ny1-1)) {index[0] = element_index + 1;}
      //If at top (no neighbour)
      else {index[0] = -1;}

       //East neighbour (just one across)
       if(ex < (Nx1-1)) {index[1] = element_index + Ny1;}
       //If at side (no neighbour)
       else 
        {
         if(ey >= n_y) {index[1] = element_index + Ny;}
         else {index[1] = -1;}
        }

       //South neighbour 
       if(ey > 0) {index[2] = element_index - 1;}
       else {index[2] = -1;}

       //West neighbour
       if(ex > 0) {index[3] = element_index - Ny1;}
       else {index[3] = -1;}

       //Now store the details in the mape
       for(unsigned i=0;i<4;i++)
        {
         if(index[i] != -1)
          {
           Neighbour_map[std::make_pair(local_el_pt,i)] =
            finite_element_pt(index[i]);
          }
        }
      }
    }
  

 //Now set up the neighbour information
 const unsigned n_element = this->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   dynamic_cast<ELEMENT*>(this->element_pt(e))->setup_face_neighbour_info();
  }
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

   ELEMENT* elem_pt=0;
   
   //Now things differ depending upon which face we are located
   switch(face_index)
    {
     //North face
    case 2:
     //Get the neighbour
     elem_pt = dynamic_cast<ELEMENT*>(
      Neighbour_map[std::make_pair(bulk_element_pt,0)]);

     if(elem_pt!=0)
      {
       //The neighbouring face element is the south face of the neighbour
       face_element_pt = elem_pt->face_element_pt(2);
      }
     //Otherwise we build in the face reflection element
     else
      {
       face_element_pt = new DGEulerFaceReflectionElement<ELEMENT>(
        bulk_element_pt,2);
       Face_element_pt.push_back(face_element_pt);

       //dynamic_cast<ELEMENT*>(bulk_element_pt)->face_element_pt(0);
      }
     
     //Then set the face coordinate
     s_face[0] = s_bulk[0];
     break;

     //East face
    case 1:
      //Get the neighbour
     elem_pt = dynamic_cast<ELEMENT*>(
      Neighbour_map[std::make_pair(bulk_element_pt,1)]);

     //The neighbouring face element is the west faec of the neighbour
     if(elem_pt!=0)
      {
       face_element_pt = elem_pt->face_element_pt(3);
      }
     else
      {
       ELEMENT* cast_bulk_element_pt = dynamic_cast<ELEMENT*>(bulk_element_pt);
       //If we're on the outlet set the face to be the face
       if(cast_bulk_element_pt->face_element_pt(1)->node_pt(0)->
          is_on_boundary(1))
        {
         face_element_pt = 
          dynamic_cast<ELEMENT*>(bulk_element_pt)->face_element_pt(1);
        }
       //Otherwise we must be on the corner
       else
        {
         face_element_pt = new DGEulerFaceReflectionElement<ELEMENT>(
          bulk_element_pt,1);
         Face_element_pt.push_back(face_element_pt);
        }
      }
     
     //Then set the face coordinate
     s_face[0] = s_bulk[1];
     break;
     
     //South face
    case -2:
      //Get the neighbour
     elem_pt = dynamic_cast<ELEMENT*>(
      Neighbour_map[std::make_pair(bulk_element_pt,2)]);

     //the neighbouring face element is the north face of the neighbour
     if(elem_pt!=0)
      {
       face_element_pt = elem_pt->face_element_pt(0);
      }
     else
      {
       //Build in the face reflection element
       face_element_pt = new DGEulerFaceReflectionElement<ELEMENT>(
        bulk_element_pt,-2);
       Face_element_pt.push_back(face_element_pt);
       //face_element_pt = 
       // dynamic_cast<ELEMENT*>(bulk_element_pt)->face_element_pt(2);
      }
     
     //Then set the face coordiante
     s_face[0] = s_bulk[0];
     break;

     //West face
    case -1:
      //Get the neighbour
     elem_pt = dynamic_cast<ELEMENT*>(
      Neighbour_map[std::make_pair(bulk_element_pt,3)]);
     
     //The neighbouring face element is the east face of the neighbour
     if(elem_pt!=0)
      {
       face_element_pt = elem_pt->face_element_pt(1);
      }
     else
      {
       face_element_pt = 
        dynamic_cast<ELEMENT*>(bulk_element_pt)->face_element_pt(3);
      }
     s_face[0] = s_bulk[1];
     break;
     
    default:
   throw OomphLibError("Coordinate is on no face or not on a unique face",
                       "DGMesh::neighbour_finder()",
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
   this->shut_up_in_newton_solve() = true;

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
   const double gamma = Global::Gamma;
   const double E = 1.0/(gamma-1.0) + 4.5*gamma;
   //Loop over all the nodes in the mesh
   unsigned n_node = mesh_pt()->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     Node* nod_pt = mesh_pt()->node_pt(n);
     //Set the values
     nod_pt->set_value(0,gamma);
     nod_pt->set_value(1,E);
     nod_pt->set_value(2,3*gamma);
     nod_pt->set_value(3,0.0);
    }
  }

 void parameter_study(std::ostream &trace, const bool &disc)
  {
   this->enable_mass_matrix_reuse();
   double dt = 0.001;
   apply_initial_conditions();
   
   Vector<double> error(4,0.0);
   char filename[100];
   
   unsigned count=1;
   for(unsigned t=0;t<500;t++)
    {
     explicit_timestep(dt);
     if(count==1)
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
  
 unsigned n_element = 4;
 
 for(unsigned i=0;i<3;i++)
  {
   TwoDDGProblem<DGSpectralEulerElement<2,2> > 
    problem(n_element,n_element);
   problem.parameter_study(trace,true);
   n_element *= 2;
  }
 
 trace.close();

 return 1;
}
