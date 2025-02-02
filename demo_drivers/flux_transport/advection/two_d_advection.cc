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
 //Set the value of pi
 const double pi = MathematicalConstants::Pi;

 //Set a diagonal wind
 void two_d_wind(const oomph::Vector<double> &x,
                 oomph::Vector<double> &wind)
 {
  wind.resize(2);
  wind[0] = 1.0; wind[1] = 1.0;
 }

 /// Function that determines the initial conditions
  void initial_condition(const double &time, const Vector<double> &x,
                        Vector<double> &u)
 {

  u[0] = sin(2*pi*x[0])*sin(2*pi*x[1]);
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
   double Lx = 1.0;
   double Ly = 1.0;

   double llx = -0.5;
   double lly = -0.5;
   //Work out the length of each element in the x and y directions
   //(Assuming uniform spacing)
   double el_length[2] = {Lx/(double)Nx, Ly/(double)Ny};

   //loop over the elements in x
   for(unsigned ex=0;ex<Nx;ex++)
    {
     //loop over the element in y
     for(unsigned ey=0;ey<Ny;ey++)
      {
       //Create a new DG element
       ELEMENT* local_element_pt = new ELEMENT;
       //Construct nodes and faces
       local_element_pt->construct_nodes_and_faces(this,time_stepper_pt);
       //Find the number of nodes
       unsigned n_node = local_element_pt->nnode();
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

   //Setup the connections between the neighbouring faces
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
     
     cast_element_pt->wind_fct_pt() = &Global::two_d_wind;
    }

      
   std::cout << "How many " << assign_eqn_numbers() << std::endl;
  }
 
 /// Compute the complete errors in the problem
 void compute_error(const double &t, Vector<double> &error)
  {
   error[0] = 0.0;
   
   Vector<double> local_error(1);
   Vector<double> local_norm(1);
   
   const unsigned n_element = Problem::mesh_pt()->nelement();
   //Do the timestep
   for(unsigned e=0;e<n_element;e++)
    {
     dynamic_cast<ScalarAdvectionEquations<2>*>(
      Problem::mesh_pt()->element_pt(e))
      ->compute_error(std::cout,
                      Global::initial_condition,t,local_error,
                      local_norm);  
     error[0] += local_error[0];
    }
  }


 void apply_initial_conditions()
  {
   //Storage for the coordinates
   Vector<double> x(2);
   //Storage for the initial condition
   Vector<double> initial_u(1);
   //Loop over all the nodes in the mesh
   unsigned n_node = mesh_pt()->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     x[0] = mesh_pt()->node_pt(n)->x(0);
     x[1] = mesh_pt()->node_pt(n)->x(1);
     //Get the initial conditions
     Global::initial_condition(0.0,x,initial_u);

     mesh_pt()->node_pt(n)->set_value(0,initial_u[0]);
    }
  }

 void parameter_study(std::ostream &trace, const bool &disc)
  {
   this->enable_mass_matrix_reuse();
   double dt = 0.001;
   apply_initial_conditions();
   
   Vector<double> error(1,0.0);
   char filename[100];
   
   unsigned count=1;
   for(unsigned t=0;t<500;t++)
    {
     explicit_timestep(dt);
     if(count==250)
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
           << this->time() << " " << sqrt(error[0]) << std::endl;
     ++count;
    }
  }
};


//----------------------------MAIN FUNCTION-------------------------------

int main()
{
 ofstream trace("trace_disc.dat");
  
 unsigned n_element = 2;
 
 for(unsigned i=0;i<3;i++)
  {
   TwoDDGProblem<DGSpectralScalarAdvectionElement<2,3> > 
    problem(n_element,n_element);
   problem.parameter_study(trace,true);
   n_element *= 2;
  }
 
 trace.close();

 return 1;
}
