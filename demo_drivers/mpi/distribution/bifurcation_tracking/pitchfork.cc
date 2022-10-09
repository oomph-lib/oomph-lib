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
// Driver code that detects the PitchforkBifurcation in the Berman similarity
// solution for porous channel flow.

// Standard includes
#include <iostream>
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <cstdio>

//My own includes
#include "generic.h"

using namespace std;

using namespace oomph;

//Define a Global Reynolds number
namespace Global_Physical_Variables
{
 double Re=0.0;
}

//=====================================================================
//Now let's define a specific 1Dmesh class for the simple ODE problem
//=====================================================================
template <class ELEMENT>
class Mesh1D : public virtual Mesh
{
 unsigned N1, N2;
public:
 //Constructor, which is where all the work takes place
 Mesh1D(const unsigned &N1, const unsigned &N2, 
        TimeStepper* time_stepper_pt=&Mesh::Default_TimeStepper);
 //Access functions for private data
 const unsigned &n1() const {return N1;}
 const unsigned &n2() const {return N2;}
 const unsigned &nx() const {return N1+N2;}

 /// Return pointer to FiniteElement
 FiniteElement* finite_element_pt(const unsigned& ielem) 
  { return dynamic_cast<FiniteElement*>(Element_pt[ielem]);}

};

//Define the mesh constructor
//Argument list:
// N1 : number of elements in the x direction, lower region (-1 -> 0)
// N2 : number of elements in the x direction, upper region (0 -> 1)
// Ntime : number of time data that need to be stored
template <class ELEMENT>
Mesh1D<ELEMENT>::Mesh1D(const unsigned &n1,
                        const unsigned &n2, 
                        TimeStepper* time_stepper_pt)
{

 //Set the internal values
 N1 = n1; N2 = n2; 
 
 //Resize the boundary node pointer to the number of boundaries
 set_nboundary(2);

 unsigned Nx = N1 + N2;

 //Allocate the store for the elements
 Element_pt.resize(Nx);
 //Allocate the memory for the first element
 Element_pt[0] = new ELEMENT;
 //Read out the number of linear points in the element
 unsigned Np = dynamic_cast<ELEMENT*>(finite_element_pt(0))->nnode_1d();

 //Can now allocate the store for the nodes 
 Node_pt.resize(1 + (Np-1)*Nx);

 //Set up geometrical data
 unsigned long node_count=0;
 double x1init = -1.0, x2init = 0.0;
 //Set the values of the increments
 double x1step = 1.0/((Np-1)*N1), x2step = 1.0/((Np-1)*N2);

 //Now assign the topology
 //Boundaries are numbered 0 (bottom) 1 (top)

 //FIRST ELEMENT (bottom)
 //Allocate memory for the node
 Node_pt[node_count] = finite_element_pt(0)->
  construct_boundary_node(0,time_stepper_pt);
 //Set the position of the node
 Node_pt[node_count]->x(0) = x1init;
 //Push the node back onto boundary
 add_boundary_node(0,Node_pt[node_count]);
 //Increment the node number
 node_count++;

 //Loop over the other nodes in the first element
 for(unsigned l1=1;l1<Np;l1++)
  {
   //Allocate memory for the nodes
   Node_pt[node_count] = 
    finite_element_pt(0)->construct_node(l1,time_stepper_pt);
   //Set the position of the node
   Node_pt[node_count]->x(0) = x1init + x1step*l1;
   //Increment the node number
   node_count++;
  }
 //END OF FIRST ELEMENT

 //CENTRAL ELEMENTS
 //Loop over all other elements apart from the last element
 for(unsigned j=1;j<(Nx-1);j++)
  {
   //Allocate memory for new element
   Element_pt[j] = new ELEMENT;
   //First node is same as last node of previous element
   finite_element_pt(j)->node_pt(0) = finite_element_pt(j-1)->node_pt((Np-1));
   //New nodes 
   for(unsigned l1=1;l1<Np;l1++)
    {
     //Allocate memory for the nodes
     Node_pt[node_count] = 
      finite_element_pt(j)->construct_node(l1,time_stepper_pt);
     
     //Set the position of the node
     if(j<N1)
      {
       Node_pt[node_count]->x(0) = x1init + x1step*(l1+(Np-1)*j);
      }
     else
      {
       Node_pt[node_count]->x(0) = x2init + x2step*(l1+(Np-1)*(j-N1));
      }
     //Increment the node number
     node_count++;
    }
  }
 //END OF CENTRAL ELEMENTS

 //FINAL ELEMENT (top)
 //Allocate memory for element
 Element_pt[Nx-1] = new ELEMENT;
 //First nodes is same as last node of previous element
 finite_element_pt(Nx-1)->node_pt(0) = finite_element_pt(Nx-2)->node_pt(Np-1);
 
 //Assign the new middle nodes
 for(unsigned l1=1;l1<(Np-1);l1++)
  {
   //Allocate memory for node
   Node_pt[node_count] = finite_element_pt(Nx-1)->construct_node(l1,time_stepper_pt);
   //Set the position of the node
   Node_pt[node_count]->x(0) = x2init + x2step*(l1+(Np-1)*(N2-1));
   //Increment the node number
   node_count++;
  }

 //New final node
 //Allocate memory for the node
 Node_pt[node_count] = 
  finite_element_pt(Nx-1)->construct_boundary_node(Np-1,time_stepper_pt); 
 //Set the position of the node
 Node_pt[node_count]->x(0) = x2init + x2step*((Np-1)*N2);
 //Push the node back onto boundary
 add_boundary_node(1,Node_pt[node_count]);
 //Increment the node number
 node_count++;
 //END OF FINAL ELEMENT IN MESH
}


//======================================================================
/// A class for elements that integrate the Berman equations in a mixed
/// formulation. The single fourth-order equation is written as a coupled
/// system of two second-order equations. The element is implemented as
/// an EigenElement, so that the temporal stability problem can be solved.
//======================================================================
class SSPorousChannelEquations : public virtual FiniteElement
{
private:

 //Pointer to Global Reynolds number
 double *Re_pt;

protected:

 /// Function to compute the shape functions and derivatives w.r.t.
 /// global coords at local coordinates.
 virtual double dshape_and_dtest_eulerian_at_knot(const unsigned &ipt,
 Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx)=0;

public:

 /// Empty Constructor, set boolean flags to be false
 SSPorousChannelEquations() {}

 /// Access functions for the Reynolds number
 const double &re() const {return *Re_pt;}

 /// Access function for a pointer to the Reynolds number
 double* &re_pt() {return Re_pt;}

 /// Return the i-th ODE variable stored at the n-th local node
 virtual double f(const unsigned &b, const unsigned &i)const =0;

 /// Return the i-th ODE variable stored at the n-th local node at time t.
 virtual double f(const unsigned &t, const unsigned &n,
                  const unsigned &i) const=0;

 /// // i-th component of df/dt at local node n. 
 /// Uses suitably interpolated value for hanging nodes.
 double df_dt(const unsigned &n, const unsigned &i) const
  {
   // Get the data's timestepper
   TimeStepper* time_stepper_pt=node_pt(n)->time_stepper_pt();

   // Number of timsteps (past & present)
   unsigned n_time = time_stepper_pt->ntstorage();
   
   double dfdt=0.0;
   
   //Loop over the timesteps
   if (time_stepper_pt->type()!="Static")
    {
     for(unsigned t=0;t<n_time;t++)
      {
       dfdt+=time_stepper_pt->weight(1,t)*f(t,n,i);
      }
    }
   
   return dfdt;
  }

 /// This function returns just the residuals
 inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Create a dummy matrix
   DenseMatrix<double> dummy(1);
   //Call the generic residuals function with flag set to 0
   add_generic_residual_contribution(residuals,dummy,0);
  }

 /// This function returns the residuals and the jacobian
 inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                          DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   add_generic_residual_contribution(residuals,jacobian,1);
  }


 //----------------------------------------------------------------------
 /// This function returns the residuals for the ODE; J
 /// FLAG=1(or 0): do (or don't) compute the Jacobian as well. 
 //----------------------------------------------------------------------
 void add_generic_residual_contribution(Vector<double> &residuals, 
                                        DenseMatrix<double> &jacobian, 
                                        unsigned flag) 
  {
   //Find out how many nodes there are
   unsigned n_node = nnode();
           
   //Set up memory for the shape and test functions
   Shape psif(n_node), testf(n_node);
   DShape dpsifdx(n_node,1), dtestfdx(n_node,1);

   //Determine the number of integration points
   unsigned n_intpt = integral_pt()->nweight();

   //Get Reynolds number from the element.
   double Re = re();

   //Integers to store the local equation and unknown values
   int local_eqn=0, local_unknown=0;

   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Get the integral weight
     double w = integral_pt()->weight(ipt);

     //Find the shape and test functions and return the Jacobian
     //of the mapping.
     double J =
      dshape_and_dtest_eulerian_at_knot(ipt,psif,dpsifdx,testf,dtestfdx);

     //Premultiply the weights and the Jacobian
     double W = w*J;
     
     //Need to find position to feed into equations
     Vector<double> interpolated_x(1,0.0);
     //Assign memory for the ODE variables
     Vector<double> interpolated_f(2,0.0);
     Vector<double> interpolated_dfdx(2,0.0);
     double df1dt=0.0;


     //Calculate velocities and derivatives
     for(unsigned l=0;l<n_node;l++) 
      {
       //Set the value of df1dt
       df1dt += df_dt(l,1)*psif[l];           
       
       //Set value of x
       interpolated_x[0] += nodal_position(l,0)*psif[l];
       //Loop over the ODE variables
       for(unsigned i=0;i<2;i++)
        {
         interpolated_f[i] += f(l,i)*psif[l];
         //and their derivatives
         interpolated_dfdx[i] += f(l,i)*dpsifdx(l,0);
        }
      }

     //CALCULATE THE RESIDUALS --- THESE ARE THE GOVERNING EQUATIONS
     //FOR BERMAN FLOW!
 
     //Loop over the test functions
     for(unsigned l=0;l<n_node;l++)
      {
       //FIRST EQUATION
       local_eqn = nodal_local_eqn(l,0);
       //If not a boundary condition
       if(local_eqn >= 0)
        {
         //Time term
         residuals[local_eqn] -= Re*df1dt*testf[l]*W;
         //Viscous term
         residuals[local_eqn] -= 
          interpolated_dfdx[1]*dtestfdx(l,0)*W;
         //Inertial term
         residuals[local_eqn] += 
          Re*(interpolated_f[0]*interpolated_dfdx[1] 
              /*Comment out the next term for Iain's equations*/
              - interpolated_f[1]*interpolated_dfdx[0])*testf[l]*W;

         //If we are calculating the Jacobian
         if(flag)
          {
           //Loop over the nodes again
           for(unsigned l2=0;l2<n_node;l2++)
            {
             local_unknown = nodal_local_eqn(l2,0);
             //FIRST DOF: If at non-zero degree of freedom add in entry
             if(local_unknown >= 0)
              {
               jacobian(local_eqn,local_unknown)
                += Re*(psif[l2]*interpolated_dfdx[1] 
                       /*Comment out the next term for Iain's equations*/
                       - interpolated_f[1]*dpsifdx(l2,0))*testf[l]*W;
              }

             //SECOND DOF: If at non-zero degree of freedom add in entry
             local_unknown = nodal_local_eqn(l2,1);
             if(local_unknown >= 0)
              {
               jacobian(local_eqn,local_unknown)
                -= Re*node_pt(l2)->time_stepper_pt()->weight(1,0)*psif[l2]
                *testf[l]*W;
               jacobian(local_eqn,local_unknown)
                -= dpsifdx(l2,0)*dtestfdx(l,0)*W;
               jacobian(local_eqn,local_unknown)
                += Re*(interpolated_f[0]*dpsifdx(l2,0) 
                       /*Comment out the next term for Iains' equations*/
                       - psif[l2]*interpolated_dfdx[0])*testf[l]*W;
              }
            }
          } //End of Jacobian calculation
        } //End of first equation

       //SECOND EQUATION
       local_eqn = nodal_local_eqn(l,1);
       if(local_eqn >= 0)
        {
         residuals[local_eqn] += 
          (interpolated_f[1]*testf[l] + 
           interpolated_dfdx[0]*dtestfdx(l,0))*W;

         //If calculating the jacobian
         if(flag)
          {
           //Loop over the nodes again
           for(unsigned l2=0;l2<n_node;l2++)
            {
             local_unknown = nodal_local_eqn(l2,0);
             //FIRST DOF: If at non-zero degree of freedom, add in entry
             if(local_unknown >= 0)
              {
               jacobian(local_eqn,local_unknown) +=
                dpsifdx(l2,0)*dtestfdx(l,0)*W;
              }
             
             //SECOND DOF: If at non-zero degree of freedom, add in entry
             local_unknown = nodal_local_eqn(l2,1);
             if(local_unknown >= 0)
              {
               jacobian(local_eqn,local_unknown) +=
                psif[l2]*testf[l]*W;
              }
            }
          } //End of Jacobian calculation
        } //End of second equation
                
      } //End of loop over shape functions

    } //End of loop over integration points
  }

/// FE interpolated values of the arguments
 double interpolated_f(const Vector<double> &s, const unsigned &i)
  {
   //Find number of nodes
   unsigned n_node = nnode();
   //Local shape function
   Shape psi(n_node);
   //Find values of the shape function
   shape(s,psi);
   
   //Initialise value of f
   double interpolated_f = 0.0;
   //Loop over the nodes and add to the sum
   for(unsigned l=0;l<n_node;l++)
    {
     interpolated_f += f(l,i)*psi[l];
    }

   return(interpolated_f);
  }

 /// Overload the output function
 void output(ostream &outfile) {FiniteElement::output(outfile);}
 
/// Output function: x,y,[z],u,v,[w],p in tecplot format
 void output(ostream &outfile, const unsigned &Np)
  {FiniteElement::output(outfile,Np);}
 
}; 


//A Standard ODE Element for testing the mixed formulation
//==========================================================================
class SSPorousChannelElement : public virtual QElement<1,3>, 
 public virtual SSPorousChannelEquations
{
  private:
 //Static array of ints to hold number of variables at nodes
 static const unsigned Initial_Nvalue[];

 protected:

 //Function to do the whole derivatives and Jacobian thing
 double dshape_and_dtest_eulerian_at_knot(const unsigned &ipt, 
 Shape &psi, DShape &dpsidx, Shape &test, DShape &dtestdx) 
  {
   //Call the geometrical shape functions and derivatives
   double J = QElement<1,3>::dshape_eulerian_at_knot(ipt,psi,dpsidx);
 
  //Loop over the test functions and derivatives and set them equal
   //to the shape functions
   for(unsigned l=0;l<3;l++)
    {
     test[l] = psi[l];
     dtestdx(l,0) = dpsidx(l,0);
    }
   //Return the Jacobian
   return J;
  }
 
public:
 

 /// Constructor, there are no internal values 
 SSPorousChannelElement() : QElement<1,3>(), SSPorousChannelEquations() { }

 /// Required number of values at node n
 inline unsigned required_nvalue(const unsigned &n) const 
  {return Initial_Nvalue[n];}
 
 //Access function for the additionally stored ODEs
 double f(const unsigned &l, const unsigned &i) const
  {return *node_pt(l)->value_pt(i);}
 
 //Access function for the ODE variables at previous times
 double f(const unsigned &t, const unsigned &l, const unsigned &i)
 const {return *node_pt(l)->value_pt(t,i);}

 //Access function for the projected derivative function
 double df0dx(const unsigned &l) const
  {return *node_pt(l)->value_pt(2);}

 /// Overload the output function
 void output(ostream &outfile) {FiniteElement::output(outfile);}

 /// Output function: x,y[,z],u,v[,w],p at Nplot^DIM points
 void output(ostream &outfile, const unsigned &Nplot)
  {
   //Set output Vector
   Vector<double> s(1);
   //Set length of domain 
   double L=10.0;

   //Dummy loop over x
   for(unsigned e=0;e<20;e++)
    {
   //Tecplot head info
   outfile << "ZONE I=" << Nplot << ", J=" << Nplot << std::endl;
    
   //Loop over nodes in the element
   for(unsigned l2=0;l2<Nplot;l2++)
    {
     s[0] = -1.0 + l2*2.0/(Nplot-1);
   //Loop over the dummy nodes in the element
   for(unsigned l=0;l<Nplot;l++)
    {
     double x = e*L/20.0 + l*(L/20)/(Nplot-1);
     //Ouput x, y, u and v and f and f2
     outfile << x << " " << interpolated_x(s,0) << " ";
     outfile << -1.0*interpolated_f(s,0) << " ";
     for(unsigned i=0;i<2;i++) outfile << interpolated_f(s,i) << " ";
     outfile << std::endl;
    }
    }
   //Issue blank line
   outfile << std::endl;
    }
  }

};

//Non-inline functions for ODE_Element
const unsigned SSPorousChannelElement::Initial_Nvalue[3]={2,2,2};

//--------------------------------------------------------------------------
//Specific class for uniform transpiration problem
template<class ELEMENT>
class UniformTranspiration : public Problem
{
private:
 
 //Store the number of elements in the sections of the mesh
 unsigned NX1, NX2;

public:
 
//Constructor
 UniformTranspiration(const unsigned &N1, const unsigned &N2);
   
 //Access function for the mesh
 Mesh1D<ELEMENT>* mesh_pt() 
  {return dynamic_cast<Mesh1D<ELEMENT>*>(Problem::mesh_pt());}

 //Update function does nothing
 void actions_after_newton_solve() {}

 //Update before solve does nothing
 void actions_before_newton_solve() {}
};

//Constructor
template<class ELEMENT>
UniformTranspiration<ELEMENT>::UniformTranspiration
(const unsigned &N1, const unsigned &N2) : NX1(N1), NX2(N2)
{
 Max_residuals = 100.0;

 //Now create the mesh
 Problem::mesh_pt() = new Mesh1D<ELEMENT>(N1,N2); 

 //Complete the build of the elements
 using namespace Global_Physical_Variables;

 //Loop over all the "Normal" Elements and set the pointers
 unsigned long Nelement = mesh_pt()->nelement();
 for(unsigned long e=0;e<Nelement;e++)
  {
   //Cast to an ODE element
   ELEMENT *temp_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //Set the Reynolds number
   temp_pt->re_pt() = &Re;
  }

 //Top and bottom are dirichlet in first variable, but free in second
 mesh_pt()->boundary_node_pt(0,0)->pin(0);
 mesh_pt()->boundary_node_pt(0,0)->set_value(0,1.0);
 mesh_pt()->boundary_node_pt(1,0)->pin(0);
 mesh_pt()->boundary_node_pt(1,0)->set_value(0,-1.0);
 
 //Attach the boundary conditions to the mesh
 cout << assign_eqn_numbers() << std::endl; 
}

//Main driver loop
int main(int argc, char* argv[])
{
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif

 //Load the namespace
 using namespace Global_Physical_Variables;

 //Set the number of elements in each section of the mesh.
 unsigned Nx1=50, Nx2=50;
 //Set up the problem
 UniformTranspiration<SSPorousChannelElement> problem(Nx1,Nx2);

 //Open output file
 std::ostringstream trace_filename;
 trace_filename << "trace" << problem.communicator_pt()->my_rank()
                << ".dat";
 ofstream trace(trace_filename.str().c_str());

 //Track the node in the middle
 //The element Nx1 contains this node
 SSPorousChannelElement *Test_pt = 
  dynamic_cast<SSPorousChannelElement *>(problem.mesh_pt()->element_pt(Nx1));
  //Set the position of the node within the element
 Vector<double> s(1); s[0] = 1.0;

#ifdef OOMPH_HAS_MPI
 //Set up a dummy partition
 unsigned n_element = problem.mesh_pt()->nelement();
 Vector<unsigned> element_partition(n_element);
 for(unsigned e=0;e<n_element/2;e++) {element_partition[e]=0;}
 for(unsigned e=n_element/2;e<n_element;e++) {element_partition[e]=1;}

 //DocInfo mesh_doc_info;
 //bool report_stats=true;
 //mesh_doc_info.set_directory("RESLT_MESH");
 problem.distribute(element_partition);//,mesh_doc_info,report_stats);
 //problem.check_halo_schemes(mesh_doc_info);
 //problem.distribute();
#endif
 
 //Step up in Reynolds number
 for(unsigned i=0;i<6;i++)
  {
   //Increase the Reynolds number by 0.5 each time
   Re = i*0.5;
   problem.steady_newton_solve();
   //Output to the trace file
   trace << Re << " " 
        << Test_pt->interpolated_f(s,1) << " " 
        << Test_pt->interpolated_f(s,0) << std::endl;
  }


 //Temp to get the distribution
/* DoubleVector res;
 problem.get_residuals(res);

 Vector<unsigned> desired_global_eqns;
 for(unsigned i=0;i<400;i++)
  {
   desired_global_eqns.push_back(i);
  }

 DoubleVectorHaloScheme test(res.distribution_pt(),desired_global_eqns);

 DoubleVectorWithHaloEntries junk(res);

 junk.build_halo_scheme(&test);

 //Now let's check the synchronisation
 const unsigned n_row_local = junk.nrow_local();
 for(unsigned n=0;n<n_row_local;n++)
  {
   junk[n] = 100.0;
  }
 junk.global_value(200) = 100.0;
 junk.gather();

 for(unsigned n=0;n<problem.ndof();n++)
  {
   oomph_info << n << " " << junk.global_value(n) << "\n";
   }

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif
exit(1);*/
 


 //Specify the symmetry the hard way
 //Need to consider the possibility of distribution, so
 //only work with local dofs
 unsigned n_dof_local = problem.dof_distribution_pt()->nrow_local();
 Vector<double> backup(n_dof_local);
  for(unsigned n=0;n<n_dof_local;n++)
   {
    backup[n] = problem.dof(n);
   }
  
  //Now sort out the problem
  unsigned n_node=problem.mesh_pt()->nnode();
  for(unsigned n=0;n<n_node;n++)
   {
    Node* nod_pt = problem.mesh_pt()->node_pt(n);
    double x = nod_pt->x(0);
    if(!nod_pt->is_pinned(1))
     {
      nod_pt->set_value(1,x*x);
     }
    if(!nod_pt->is_pinned(0))
     {
      nod_pt->set_value(0,0.0);
     }
   }
  
  //The symmetry vector needs to have the DOF distribution

  //LinearAlgebraDistribution dist(problem.communicator_pt(),n_dof,false);
  DoubleVector symm(problem.dof_distribution_pt());
  for(unsigned n=0;n<n_dof_local;n++)
   {
    symm[n] = problem.dof(n);
    problem.dof(n) = backup[n];
   }


 //Let's try to find the pitchfork it
 problem.activate_pitchfork_tracking(&Re,symm,false);

 problem.steady_newton_solve();

 std::cout << "Pitchfork detected at " << Re  << std::endl;
 if(problem.communicator_pt()->my_rank()==0)
  {
   std::cout << "The slack parameter is " 
             << problem.dof(2*n_dof_local+1) << std::endl;
  }

 //Output to the trace file
 trace << Re << " " 
       << Test_pt->interpolated_f(s,1) << " " 
       << Test_pt->interpolated_f(s,0) << std::endl;
 

 // problem.deactivate_bifurcation_tracking();
 // Re += 0.01;
 // problem.steady_newton_solve();

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif
 }







