//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2022 Matthias Heil and Andrew Hazel
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
//Demonstration code to solve classic shock-tube problems
//(1D euler equations)
//using the OOMPH-LIB base classes only

//Include std C++ IO library
#include <iostream>
#include <iomanip>

//Include the OOMPH-LIB base classes
#include "generic.h"
#include "flux_transport.h"
#include "meshes/one_d_mesh.h"

using namespace oomph;
using namespace std;

namespace Global
{
 double Gamma = 1.4;
 
 void initial_left(Vector<double> &u, const bool &sod)
 {
  //Sod
  double ro, vel, pres;
  if(sod)
   {
    ro = 1.0;
    vel = 0.0;
    pres = 1.0;
   }
  else
   {
    //Lax
    ro = 0.445;
    vel = 0.698;
    pres = 3.528;
   }

  double m = vel*ro;
  double E = pres/(Gamma-1.0) + 0.5*ro*vel*vel;

  u[0] = ro;
  u[1] = E;
  u[2] = m;
 }

 void initial_right(Vector<double> &u, const bool &sod)
 {
  double ro, vel, pres;
  //Sod
  if(sod)
   {
    ro = 0.125;
    vel = 0.0;
    pres = 0.1;
   }
  else
   {
    //Lax
    ro = 0.5;
    vel = 0.0;
    pres = 0.571;
   }

  double m = vel*ro;
  double E = pres/(Gamma-1.0) + 0.5*ro*vel*vel;

  u[0] = ro;
  u[1] = E;
  u[2] = m;
 }
}


//----------------------ONE DIMENSIONAL (LINE) MESH---------------------------

//============================================================================
/// A simple one dimensional mesh: uniformly spaced nodes in the domain x=0,1
//============================================================================
template<class ELEMENT>
class OneDimMesh : public DGMesh
{

public:

 /// Mesh Constructor. The argument is the desired number of elements
 OneDimMesh(const unsigned &n_element, TimeStepper* time_stepper_pt =
  &Mesh::Default_TimeStepper)
 {
  double X_min = -1.0;
  double X_max = 1.0;

  //Resize the vector of pointers to elements: there are n_element elements
  Element_pt.resize(n_element); 
  //Loop over and construct all elements
  for(unsigned e=0;e<n_element;e++)
   {
    Element_pt[e] = new ELEMENT;
   }

  //Construct nodes and faces for the elements
  {
   //First element
   std::vector<bool> boundary_info(finite_element_pt(0)->nnode(),false);
   boundary_info[0] = true;
   dynamic_cast<ELEMENT*>(element_pt(0))
    ->construct_boundary_nodes_and_faces(this,boundary_info,time_stepper_pt);
   
   //Construct normal nodes for middle elements
   for(unsigned e=1;e<n_element-1;e++)
    {
     dynamic_cast<ELEMENT*>(element_pt(e))->
      construct_nodes_and_faces(this,time_stepper_pt);
    }
   
   //Final element   
   unsigned n_node = finite_element_pt(n_element-1)->nnode();
   boundary_info.resize(n_node,false);
   boundary_info[0] = false;
   boundary_info[n_node-1] = true;
   dynamic_cast<ELEMENT*>(element_pt(n_element-1))
    ->construct_boundary_nodes_and_faces(this,boundary_info,time_stepper_pt);
  }


  //We've now created all the nodes -- let's set their positions:
  //Elements are uniformly spaced
  double el_length = (X_max - X_min)/(double)n_element;
  //Storage for fractional position of node
  Vector<double> s_fraction;

  //Loop over the elements again
  for(unsigned e=0;e<n_element;e++)
   {
    //Locally cache the element
    FiniteElement* elem_pt = finite_element_pt(e);
    //Locally cache the pointer to the first node
    Node* nod_pt = elem_pt->node_pt(0);
    //Set it's position
    nod_pt->x(0) = X_min + e*el_length;
    //Add to the Mesh's list of nodes
    Node_pt.push_back(nod_pt);

    unsigned n_node = elem_pt->nnode();
    //Now we must loop over all other nodes
    for(unsigned n=1;n<n_node;n++)
     {
      //Get the pointer to the n-th node
      nod_pt = elem_pt->node_pt(n);

      //Find the fractional position of the node
      elem_pt->local_fraction_of_node(n,s_fraction);

      //Now set the position of the node
      nod_pt->x(0) = X_min + (e + s_fraction[0])*el_length;

      //Add to the list of nodes in the mesh
      Node_pt.push_back(nod_pt);
     }
   }


  //Find the total number of nodes created
  const unsigned n_node = this->nnode();

  //There are two boundaries
  set_nboundary(2);

  //Boundary 0 contains the first node in the mesh:
  add_boundary_node(0,Node_pt[0]);

  //Boundary 1 contains the final node in the mesh:
  add_boundary_node(1,Node_pt[n_node-1]); 

  //Check
  std::cout << "Constructed mesh with " << this->nelement() 
            << " Elements and " << this->nnode() << " Nodes" << std::endl;

  //Now set up the neighbour information
  this->setup_face_neighbour_info();

 } // End of constructor

 /// Output the face element
 void output_faces(std::ostream &outfile)
  {
   //Loop over the elements
   unsigned n_element = this->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     dynamic_cast<ELEMENT*>(element_pt(e))->output_faces(outfile);
    }
  }

 //Find the neighbour 
 void neighbour_finder(FiniteElement* const &bulk_element_pt,
                       const int &face_index,
                       const Vector<double> &s_bulk,
                       FaceElement* &face_element_pt,
                       Vector<double> &s_face)
  {
   //The face coordinate is always a vector of size zero
   s_face.resize(0);

   //Not efficient but will work
   const unsigned n_element = this->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     if(finite_element_pt(e)==bulk_element_pt)
      {
       //Now what's up
       switch(face_index)
        {
         //left face
        case -1:
         //If not the first
         if(e!=0)
          {
           //Right face of previous element
           face_element_pt = 
            dynamic_cast<ELEMENT*>(Element_pt[e-1])->face_element_pt(1);
          }
         //Otherwie just set the face to be the face of the first element
         else
          {
           face_element_pt = 
            //dynamic_cast<ELEMENT*>(Element_pt[n_element-1])->face_element_pt(1);
           dynamic_cast<ELEMENT*>(Element_pt[e])->face_element_pt(0);
          }
         break;

         //Right face
        case 1:
         //If not the last
         if(e!=n_element-1)
          {
           //Left face of next element
           face_element_pt = 
            dynamic_cast<ELEMENT*>(Element_pt[e+1])->face_element_pt(0);
          }
         //Otherwie make periodic
         else
          {
           face_element_pt = 
            //dynamic_cast<ELEMENT*>(Element_pt[0])->face_element_pt(0);
           dynamic_cast<ELEMENT*>(Element_pt[e])->face_element_pt(1);
          }
         break;


         //Otherwise complain
        default:
         throw OomphLibError(
          "Coordinate is not on any face",
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
        }
       
       //Break out of the for loop
       break;
      }
    }

  }

}; // End of OneDimMesh class.

//-----------------------PROBLEM CLASS-------------------------------------

class EulerProblem : public Problem 
{
  public:

  /// Problem constructor: 
 EulerProblem()
  {
   this->set_explicit_time_stepper_pt(new RungeKutta<4>);
   this->add_time_stepper_pt(new BDF<4>);
   this->disable_info_in_newton_solve();
  }
 
 /// Destructor: clean up memory
 ~EulerProblem()
  {
   delete Problem::mesh_pt();
   delete this->time_stepper_pt();
   delete this->explicit_time_stepper_pt();
  }
 
 //Setup the initial conditions
  void set_initial_conditions(const double &dt, const bool &sod)
   {
    this->initialise_dt(dt);
    const unsigned n_element = Problem::mesh_pt()->nelement();
   
    Vector<double> initial_u(3);
    //Set to a given function
    for(unsigned e=0;e<n_element;e++)
     {
      bool left = true;
      //If in the second half
      if(e >= n_element/2) {left = false;}
      
      //Cache the element
      FiniteElement* const elem_pt = mesh_pt()->finite_element_pt(e);
      const unsigned n_node = elem_pt->nnode();
      for(unsigned n=0;n<n_node;n++)
       {  
        //Cache the node
        Node* nod_pt = elem_pt->node_pt(n);
        if(left) {Global::initial_left(initial_u,sod);}
        else {Global::initial_right(initial_u,sod);}
        

        //for(unsigned i=0;i<3;i++)
        // {
        //  nod_pt->set_value(i,initial_u[i]);
        // }

        //Initialise any previous values (impulsive start)
        const unsigned n_prev = nod_pt->time_stepper_pt()->nprev_values();
        for(unsigned t=0;t<=n_prev;t++)
         {
          for(unsigned i=0;i<3;i++)
           {
            nod_pt->set_value(t,i,initial_u[i]);
           }
         }
       }
     }
   }


 
 void parameter_study(const bool &sod)
  {
   this->enable_mass_matrix_reuse();
   this->time() = 0.0;
   double dt = 0.0001;
   this->set_initial_conditions(dt,sod);
   
   Vector<double> error(1,0.0);
   
   char filename[100];
   ofstream outfile;
   
   if(sod)
    {
     sprintf(filename,"sod_%li_time%g.dat",mesh_pt()->nelement(),
             this->time());
    }
   else
    {
     sprintf(filename,"lax_%li_time%g.dat",mesh_pt()->nelement(),
             this->time());
    }
   outfile.open(filename);
   mesh_pt()->output(outfile);
   outfile.close();

    unsigned count=1;
    for(unsigned i=0;i<2000;i++)
     {
      //Take an explicit timestep
      explicit_timestep(dt);

      //Take an implicit timestep
      //unsteady_newton_solve(dt);

      if(count==500)
       {
        if(sod)
         {
          sprintf(filename,"sod_%li_time%g.dat",mesh_pt()->nelement(),
                  this->time());
         }
        else
         {
          sprintf(filename,"lax_%li_time%g.dat",mesh_pt()->nelement(),
                  this->time());
         }
        outfile.open(filename);
        mesh_pt()->output(outfile);
        outfile.close();
        count=0;
       }
      

      ++count;
     }
   }
};


//==========================================================================
/// Define the Problem which creates the mesh, applies the
/// boundary conditions, and assigns equation numbers.
//==========================================================================
template<class ELEMENT>
class DGProblem : public EulerProblem
 {
  /// Pointer to a slope limiter
  SlopeLimiter* Slope_limiter_pt;
  

  public:

  /// Problem constructor: 
  DGProblem(const unsigned &n_element) 
   {
    //Make the formulation discontinuous
    this->enable_discontinuous_formulation();

    //Create a OneDimMesh Mesh object and set it to be the problem's mesh.
    //The element type, TwoNodePoissonElement, is passed  as a template 
    //parameter to the mesh. The argument to the constructor indicates
    //the number of elements in the mesh.
    Problem::mesh_pt() = new OneDimMesh<ELEMENT>(n_element,
                                                 this->time_stepper_pt());

    //Set the slope limiter
    Slope_limiter_pt = new MinModLimiter(5.0,true);
    
    
    for(unsigned i=0;i<3;i++)
     {
      mesh_pt()->boundary_node_pt(0,0)->pin(i);
      mesh_pt()->boundary_node_pt(1,0)->pin(i);
     }

    //Assign the local equation numbers
    std::cout << assign_eqn_numbers() << " Equation numbers assigned "
              << std::endl;
   }

  /// Clean-up memory
  ~DGProblem()
   {
    delete Slope_limiter_pt;
   }


 /// Slope limit the solution, if required
 //void actions_after_explicit_timestep()
 // {dynamic_cast<DGMesh*>(this->mesh_pt())->limit_slopes(Slope_limiter_pt);}

}; //End of problem definition


//----------------------------MAIN FUNCTION-------------------------------

int main()
 {
  //One hunder elements
  unsigned n_element = 100;

  DGProblem<DGSpectralEulerElement<1,2> > problem(n_element);
  std::cout << "Sod problem\n";
  //Solve the sod problem
  problem.parameter_study(true);

  std::cout << "Lax problem\n";
  //Solve the lax problem
  problem.parameter_study(false);

  return 1;
 }
