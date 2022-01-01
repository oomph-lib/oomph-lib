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
//Demonstration code to solve a one-dimensional Poisson equation
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
 //Set a one dimensional constant wind in the x-direction
 void constant_wind(const oomph::Vector<double> &x,
                    oomph::Vector<double> &wind)
 {
  wind.resize(1);
  wind[0] = 1.0;
 }

  void initial_condition(const double &time,const Vector<double> &x,
                        Vector<double> &u)
 {
  u[0] = sin(8.0*atan(1.0)*x[0]);
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
  double X_min = 0.0;
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
    nod_pt->x(0) = e*el_length;
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
      nod_pt->x(0) = (e + s_fraction[0])*el_length;

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
         //Otherwie make periodic
         else
          {
           face_element_pt = 
            dynamic_cast<ELEMENT*>(Element_pt[n_element-1])->face_element_pt(1);
           //dynamic_cast<ELEMENT*>(Element_pt[e])->face_element_pt(0);
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
            dynamic_cast<ELEMENT*>(Element_pt[0])->face_element_pt(0);
           //dynamic_cast<ELEMENT*>(Element_pt[e])->face_element_pt(1);
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

class AdvectionProblem : public Problem 
{
  public:

  /// Problem constructor: 
 AdvectionProblem()
  {
   this->set_explicit_time_stepper_pt(new RungeKutta<4>);
   this->add_time_stepper_pt(new BDF<4>);
   this->disable_info_in_newton_solve();
   //the problem is linear
   Problem_is_nonlinear = false;
  }
 
 /// Destructor: clean up memory
 ~AdvectionProblem()
  {
   delete Problem::mesh_pt();
   delete this->time_stepper_pt();
   delete this->explicit_time_stepper_pt();
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
     dynamic_cast<ScalarAdvectionEquations<1>*>(
      Problem::mesh_pt()->element_pt(e))
      ->compute_error(std::cout,
                      Global::initial_condition,t,local_error,
                      local_norm);  
     error[0] += local_error[0];
    }
  }
 
 //Setup the initial conditions
  void set_initial_conditions(const double &dt)
   {
    this->initialise_dt(dt);
    const unsigned n_element = Problem::mesh_pt()->nelement();
    //Storage for the x-coordinate
    Vector<double> x(1);
    //Storage for the initial condition
    Vector<double> initial_u(1);

    //Cache the constant wind
    Vector<double> wind(1);
    Global::constant_wind(x,wind);

    //Set to a given function
    for(unsigned e=0;e<n_element;e++)
     {
      FiniteElement* const elem_pt = mesh_pt()->finite_element_pt(e);
      const unsigned n_node = elem_pt->nnode();
      for(unsigned n=0;n<n_node;n++)
       {  
        //Cache the node
        Node* nod_pt = elem_pt->node_pt(n);
        //Set the x-coordinate
        x[0] = nod_pt->x(0);
        //Get the initial condition
        Global::initial_condition(0.0,x,initial_u);
        nod_pt->set_value(0,initial_u[0]);
        //Set the previous values
        const unsigned n_prev = nod_pt->time_stepper_pt()->nprev_values();
        for(unsigned i=0;i<n_prev;i++)
         {
          //Move the x-coordinate backwards
          x[0] += wind[0]*dt;
          //Now get the new value
          Global::initial_condition(0.0,x,initial_u);
          nod_pt->set_value(i+1,0,initial_u[0]);
         }
       }
     }
   }

 /// Assign the equations number associated with the problem
 /// This virtual so that it can be overloaded in the discontinuous problem
 /// to setup the coupling degrees of freedom for implicit timestepping
 /// when required.
 virtual void assign_equation_numbers(const bool &explicit_timestepper)
  {
   //The default is simply to assign the equation numbers
   std::cout << assign_eqn_numbers() << " Equation numbers assigned "
             << std::endl;
  }

 void parameter_study(std::ostream &trace, const bool &explicit_timestepper,
                      const bool &disc)
  {
   //Defer assignment of the equation numbers until the parameter study
   //so that the inclusion of coupling terms between faces can be
   //determined from the type of timestep required.
   this->assign_equation_numbers(explicit_timestepper);
   this->enable_mass_matrix_reuse();
   this->enable_jacobian_reuse();
   double dt = 0.001;
   this->set_initial_conditions(dt);
   
   Vector<double> error(1,0.0);
    
   char filename[100];
   ofstream outfile;
   unsigned count=1;
   for(unsigned i=0;i<1000;i++)
     {
      if(explicit_timestepper)
       {
        explicit_timestep(dt);
       }
      else
       {
        unsteady_newton_solve(dt);
       }

      if(count==250)
       {
        if(disc)
         {
          sprintf(filename,"disc_%li_time%g.dat",mesh_pt()->nelement(),
                  this->time());
         }
        else
         {
          sprintf(filename,"cont_%li_time%g.dat",mesh_pt()->nelement(),
                  this->time());
         }
        outfile.open(filename);
        mesh_pt()->output(outfile);
        outfile.close();
        count=0;
       }
      
      compute_error(this->time(),error);
      trace << mesh_pt()->nelement() << " " << 
       this->time() << " " << sqrt(error[0]) << std::endl;
      ++count;
     }
    trace << "\n\n";
   }
};


//==========================================================================
/// Define the Problem which creates the mesh, applies the
/// boundary conditions, and assigns equation numbers.
//==========================================================================
template<class ELEMENT>
class DGProblem : public AdvectionProblem
 {
  public:

  /// Assign the equations number associated with the problem
  /// AFTER setting up any potential coupling between the faces.
  virtual void assign_equation_numbers(const bool &explicit_timestepper)
   {
    //Make the formulation discontinuous if we have an explicit timestepper
    if(explicit_timestepper)
     {this->enable_discontinuous_formulation();}

    //Setup the coupling between the faces
    dynamic_cast<OneDimMesh<ELEMENT>*>(this->mesh_pt())
     ->setup_face_neighbour_info(!explicit_timestepper);

    //The default is simply to assign the equation numbers
    std::cout << assign_eqn_numbers() << " Equation numbers assigned "
              << std::endl;
   }


  /// Problem constructor: 
  DGProblem(const unsigned &n_element) 
   {
    //Create a OneDimMesh Mesh object and set it to be the problem's mesh.
    //The element type, TwoNodePoissonElement, is passed  as a template 
    //parameter to the mesh. The argument to the constructor indicates
    //the number of elements in the mesh.
    Problem::mesh_pt() = new OneDimMesh<ELEMENT>(n_element,
                                                 this->time_stepper_pt());
    
    //Set the wind function
    for(unsigned e=0;e<n_element;e++)
     {
      dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
       wind_fct_pt() = &Global::constant_wind;
     }
   }

}; //End of problem definition


//-----------------------PROBLEM CLASS-------------------------------------

//==========================================================================
/// Define the Problem which creates the mesh, applies the
/// boundary conditions, and assigns equation numbers.
//==========================================================================
template<class ELEMENT>
class ContProblem : public AdvectionProblem
 {
  public:

  /// Problem constructor: 
  ContProblem(const unsigned &n_element) : AdvectionProblem()
   {
    //Create a OneDimMesh Mesh object and set it to be the problem's mesh.
    //The element type, TwoNodePoissonElement, is passed  as a template 
    //parameter to the mesh. The argument to the constructor indicates
    //the number of elements in the mesh.
    Problem::mesh_pt() = new OneDMesh<ELEMENT>(n_element,1.0,
                                               time_stepper_pt());

    //Make the boundaries periodic
    // This has been modified from the original (commented out below)
    // so that it passes the self tests
    mesh_pt()->boundary_node_pt(1,0)->
     make_periodic(mesh_pt()->boundary_node_pt(0,0));
//    mesh_pt()->boundary_node_pt(0,0)->
//     make_periodic(mesh_pt()->boundary_node_pt(1,0));

    //Set the wind function
    for(unsigned e=0;e<n_element;e++)
     {
      dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
       wind_fct_pt() = &Global::constant_wind;
     }
   }

}; //End of problem definition



//----------------------------MAIN FUNCTION-------------------------------

int main()
 {
  ofstream trace("trace_disc.dat");
  
  unsigned n_element = 2;

  for(unsigned i=0;i<3;i++)
   {
    DGProblem<DGSpectralScalarAdvectionElement<1,7> > problem(n_element);
    problem.parameter_study(trace,true,true);
    n_element *= 2;
   }
  trace.close();
  
  trace.open("trace_cont.dat");
  n_element = 2;
  for(unsigned i=0;i<3;i++)
   {
    ContProblem<QSpectralScalarAdvectionElement<1,7> > problem(n_element);
    problem.parameter_study(trace,true,false);
    n_element *= 2;
   }

  trace.close();

  return 1;
 }
