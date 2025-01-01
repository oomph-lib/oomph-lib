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

public:
 //Constructor
 TwoDDGMesh(const unsigned &Nx, const unsigned &Ny, 
          TimeStepper* time_stepper_pt = &Mesh::Default_TimeStepper)
  {
   Vector<double> s_fraction;
   //Lengths of the mesh
   double Lx = 10.0;
   double Ly = 10.0;

   double llx = -5.0;
   double lly = -5.0;
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
         //Now set the value
         nod_pt->set_value(0,1.0);
         
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
       unsigned index[4];

       //North neighbour (just one up)
       if(ey < (Ny-1)) {index[0] = element_index + 1;}
       //If at top, return index at bottom (periodic conditions)
       else {index[0] = ex*Ny;}

       //East neighbour (just one across)
       if(ex < (Nx-1)) {index[1] = element_index + Ny;}
       //If at side return index at other side (periodic conditions)
       else {index[1] = ey;}

       //South neighbour 
       if(ey > 0) {index[2] = element_index - 1;}
       else {index[2] = element_index + Ny-1;}

       //West neighbour
       if(ex > 0) {index[3] = element_index - Ny;}
       else {index[3] = element_index + (Nx-1)*Ny;}

       //Now store the details in the mape
       for(unsigned i=0;i<4;i++)
        {
         Neighbour_map[std::make_pair(local_el_pt,i)] =
          finite_element_pt(index[i]);
        }
      }
    }
   
   //Now set up the neighbour information
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
   
   //Now things differ depending upon which face we are located
   switch(face_index)
    {
     //North face
    case 2:
     //The neighbouring face element is the south face of the neighbour
     face_element_pt = 
      dynamic_cast<ELEMENT*>(
       Neighbour_map[std::make_pair(bulk_element_pt,0)])->
      face_element_pt(2);
     //Then set the face coordinate
     s_face[0] = s_bulk[0];
     break;

     //East face
    case 1:
     //The neighbouring face element is the west faec of the neighbour
     face_element_pt = 
      dynamic_cast<ELEMENT*>(
       Neighbour_map[std::make_pair(bulk_element_pt,1)])->
      face_element_pt(3);
     //Then set the face coordinate
     s_face[0] = s_bulk[1];
     break;
     
     //South face
    case -2:
     //the neighbouring face element is the north face of the neighbour
     face_element_pt = 
      dynamic_cast<ELEMENT*>(
       Neighbour_map[std::make_pair(bulk_element_pt,2)])
      ->face_element_pt(0);
     //Then set the face coordiante
     s_face[0] = s_bulk[0];
     break;

     //West face
    case -1:
     //The neighbouring face element is the east face of the neighbour
     face_element_pt = 
      dynamic_cast<ELEMENT*>
      (Neighbour_map[std::make_pair(bulk_element_pt,3)])->
      face_element_pt(1);
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

   //Set the boundary conditions
   Vector<double> u_initial(4);
   Vector<double> x(2);
   //Pin all boundaries
   for(unsigned b=0;b<4;b++)
    {
     const unsigned n_node = mesh_pt()->nboundary_node(b);
     for(unsigned n=0;n<n_node;n++)
      {
       Node* nod_pt = mesh_pt()->boundary_node_pt(b,n);
       const unsigned n_value = nod_pt->nvalue();
       for(unsigned i=0;i<n_value;i++) {nod_pt->pin(i);}
       //Now set the values
       for(unsigned i=0;i<2;i++) {x[i] = nod_pt->x(i);}
       Global::exact_solution(0.0,x,u_initial);
       for(unsigned i=0;i<n_value;i++) {nod_pt->set_value(i,u_initial[i]);}
      }
    }

      
   std::cout << "How many " << assign_eqn_numbers() << std::endl;
  }
 
 /// Compute the complete errors in the problem
 void compute_error(const double &t, Vector<double> &error)
  {
   error.initialise(0.0);
   
   Vector<double> local_error(4);
   Vector<double> local_norm(4);
   
   const unsigned n_element = Problem::mesh_pt()->nelement();
   //Do the timestep
   for(unsigned e=0;e<n_element;e++)
    {
     dynamic_cast<ELEMENT*>(
      Problem::mesh_pt()->element_pt(e))
      ->compute_error(std::cout,
                      Global::exact_solution,t,local_error,
                      local_norm);  
     for(unsigned i=0;i<4;i++)
      {
       error[i] += local_error[i];
      }
   }
  }


 void apply_boundary_conditions(const double &t)
  {
   //Set the boundary conditions
   Vector<double> u_initial(4);
   Vector<double> x(2);
   //Pin all boundaries
   for(unsigned b=0;b<4;b++)
    {
     const unsigned n_node = mesh_pt()->nboundary_node(b);
     for(unsigned n=0;n<n_node;n++)
      {
       Node* nod_pt = mesh_pt()->boundary_node_pt(b,n);
       const unsigned n_value = nod_pt->nvalue();
       //Now set the values
       for(unsigned i=0;i<2;i++) {x[i] = nod_pt->x(i);}
       Global::exact_solution(t,x,u_initial);
       for(unsigned i=0;i<n_value;i++) {nod_pt->set_value(i,u_initial[i]);}
      }
    }
  }

 void actions_after_explicit_stage()
  {
   apply_boundary_conditions(this->time());
  }


 void apply_initial_conditions()
  {
   //Storage for the coordinates
   Vector<double> x(2);
   //Storage for the initial condition
   Vector<double> initial_u(4);
   //Loop over all the nodes in the mesh
   unsigned n_node = mesh_pt()->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     Node* nod_pt = mesh_pt()->node_pt(n);
     x[0] = nod_pt->x(0);
     x[1] = nod_pt->x(1);
     //Get the initial conditions
     Global::exact_solution(0.0,x,initial_u);
     //Set them
     for(unsigned i=0;i<4;i++)
      {
       nod_pt->set_value(i,initial_u[i]);
      }
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
   for(unsigned t=0;t<50;t++)
    {
     explicit_timestep(dt);
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

     compute_error(this->time(),error);
     trace << mesh_pt()->nelement() << " " 
           << this->time() << " " << sqrt(error[0]) << 
      " " << sqrt(error[1]) << " " << sqrt(error[2])
           << " " << sqrt(error[3]) << std::endl;
     ++count;
    }
  }
};


//----------------------------MAIN FUNCTION-------------------------------

int main()
{
 ofstream trace("trace_disc.dat");
  
 unsigned n_element = 8;
 
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
