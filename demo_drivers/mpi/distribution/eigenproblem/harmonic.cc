//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
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
//Driver to solve the harmonic equation with homogeneous Dirichlet boundary
//conditions.

// Generic oomph-lib routines
#include "generic.h"

// Include the mesh
#include "meshes/one_d_mesh.h"

using namespace std;

using namespace oomph;

//===================================================================
/// Function-type-object to perform comparison of complex data types
/// Needed to sort the complex eigenvalues into order based on the
/// size of the real part.
//==================================================================
template <class T>
class ComplexLess
{
public:

 /// Comparison. Are the values identical or not?
 bool operator()(const complex<T> &x, const complex<T> &y) const
  {
   return x.real() < y.real();
  }
};


//=================================================================
/// A class for all elements that solve the simple one-dimensional
/// eigenvalue problem
/// \f[ 
/// \frac{\partial^2 u}{\partial x_i^2}  + \lambda u = 0
/// \f] 
/// These elements are very closely related to the Poisson
/// elements and could inherit from them. They are here developed
/// from scratch for pedagogical purposes.
/// This class  contains the generic maths. Shape functions, geometric
/// mapping etc. must get implemented in derived class.
//================================================================
class HarmonicEquations : public virtual FiniteElement
{

public:
 /// Empty Constructor
 HarmonicEquations() {}
 
 /// Access function: Eigenfunction value at local node n
 /// Note that solving the eigenproblem does not assign values
 /// to this storage space. It is used for output purposes only.
 virtual inline double u(const unsigned& n) const 
  {return nodal_value(n,0);}

 /// Output the eigenfunction with default number of plot points
 void output(ostream &outfile) 
  {
   unsigned nplot=5;
   output(outfile,nplot);
  }

 /// Output FE representation of soln: x,y,u or x,y,z,u at 
 /// Nplot  plot points
 void output(ostream &outfile, const unsigned &nplot)
  {
   //Vector of local coordinates
   Vector<double> s(1);

   // Tecplot header info
   outfile << tecplot_zone_string(nplot);
   
   // Loop over plot points
   unsigned num_plot_points=nplot_points(nplot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     // Get local coordinates of plot point
     get_s_plot(iplot,nplot,s);
     //Output the coordinate and the eigenfunction
     outfile << interpolated_x(s,0) << " " << interpolated_u(s) << std::endl;   
    }
   // Write tecplot footer (e.g. FE connectivity lists)
   write_tecplot_zone_footer(outfile,nplot);
  }

 /// Assemble the contributions to the jacobian and mass
 /// matrices
 void fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double> &residuals,
  DenseMatrix<double> &jacobian, DenseMatrix<double> &mass_matrix)
  {
   //Find out how many nodes there are
   unsigned n_node = nnode();
   
   //Set up memory for the shape functions and their derivatives
   Shape psi(n_node);
   DShape dpsidx(n_node,1);

   //Set the number of integration points
   unsigned n_intpt = integral_pt()->nweight();
   
   //Integers to store the local equation and unknown numbers
   int local_eqn=0, local_unknown=0;
   
   //Loop over the integration points
   for(unsigned ipt=0;ipt<n_intpt;ipt++)
    {
     //Get the integral weight
     double w = integral_pt()->weight(ipt);
     
     //Call the derivatives of the shape and test functions
     double J = dshape_eulerian_at_knot(ipt,psi,dpsidx);
     
     //Premultiply the weights and the Jacobian
     double W = w*J;

     //Assemble the contributions to the mass matrix
     //Loop over the test functions
     for(unsigned l=0;l<n_node;l++)
      {
       //Get the local equation number
       local_eqn = u_local_eqn(l);
       /*IF it's not a boundary condition*/
       if(local_eqn >= 0)
        {
         //Loop over the shape functions
         for(unsigned l2=0;l2<n_node;l2++)
          { 
           local_unknown = u_local_eqn(l2);
           //If at a non-zero degree of freedom add in the entry
           if(local_unknown >= 0)
            {
             jacobian(local_eqn,local_unknown) += dpsidx(l,0)*dpsidx(l2,0)*W;
             mass_matrix(local_eqn, local_unknown) += psi(l)*psi(l2)*W;
            }
          }
        }
      }
    }
  }
 
 /// Return FE representation of function value u(s) at local coordinate s
 inline double interpolated_u(const Vector<double> &s) const
  {
   unsigned n_node = nnode();

   //Local shape functions
   Shape psi(n_node);

   //Find values of basis function
   this->shape(s,psi);

   //Initialise value of u
   double interpolated_u = 0.0;

   //Loop over the local nodes and sum
   for(unsigned l=0;l<n_node;l++)  {interpolated_u+=u(l)*psi[l];}

   //Return the interpolated value of the eigenfunction
   return(interpolated_u);
  }

protected:

 /// Shape/test functions and derivs w.r.t. to global coords at 
 /// local coord. s; return  Jacobian of mapping
 virtual double dshape_eulerian(const Vector<double> &s, 
                                Shape &psi, 
                                DShape &dpsidx) const=0;

 /// Shape/test functions and derivs w.r.t. to global coords at 
 /// integration point ipt; return  Jacobian of mapping
 virtual double dshape_eulerian_at_knot(const unsigned &ipt, 
                                        Shape &psi, 
                                        DShape &dpsidx) const=0;
 
 /// Access function that returns the local equation number
 /// of the unknown in the problem. Default is to assume that it is the
 /// first (only) value stored at the node.
 virtual inline int u_local_eqn(const unsigned &n)
  {return nodal_local_eqn(n,0);}
 
    private:
 
};



//======================================================================
/// QHarmonicElement<NNODE_1D> elements are 1D  Elements with 
/// NNODE_1D nodal points that are used to solve the Harmonic eigenvalue
/// Problem described by HarmonicEquations.
//======================================================================
template <unsigned NNODE_1D>
class QHarmonicElement : public virtual QElement<1,NNODE_1D>, 
                         public HarmonicEquations
{
 
  public:

 /// Constructor: Call constructors for QElement and 
 /// Poisson equations
 QHarmonicElement() : QElement<1,NNODE_1D>(), HarmonicEquations()
  {}



 ///  Required  # of `values' (pinned or dofs) 
 /// at node n
 inline unsigned required_nvalue(const unsigned &n) const 
  {
   return Initial_Nvalue[n];
  }

 /// Output function:  
 ///  x,y,u   or    x,y,z,u
 void output(ostream &outfile)
  {
   HarmonicEquations::output(outfile);
  }



 ///  Output function:  
 ///   x,y,u   or    x,y,z,u at Nplot^DIM plot points
 void output(ostream &outfile, const unsigned &Nplot)
  {
   HarmonicEquations::output(outfile,Nplot);
  }


protected:

/// Shape, test functions & derivs. w.r.t. to global coords. Return Jacobian.
 inline double dshape_eulerian(const Vector<double> &s, 
                               Shape &psi, 
                               DShape &dpsidx) const;
 

 /// Shape, test functions & derivs. w.r.t. to global coords. at
 /// integration point ipt. Return Jacobian.
 inline double dshape_eulerian_at_knot(const unsigned& ipt,
                                       Shape &psi, 
                                       DShape &dpsidx) const;
private:
 
 /// Static array of ints to hold number of variables
 /// at nodes: Initial_Nvalue[n]
 static const unsigned Initial_Nvalue[];
 

};




//Inline functions:


//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned NNODE_1D>
double QHarmonicElement<NNODE_1D>::
dshape_eulerian(const Vector<double> &s,
                Shape &psi, 
                DShape &dpsidx) const
{
 //Call the geometrical shape functions and derivatives  
 double J = QElement<1,NNODE_1D>::dshape_eulerian(s,psi,dpsidx);

 //Return the jacobian
 return J;
}


//======================================================================
/// Define the shape functions and test functions and derivatives
/// w.r.t. global coordinates and return Jacobian of mapping.
///
/// Galerkin: Test functions = shape functions
//======================================================================
template<unsigned NNODE_1D>
double QHarmonicElement<NNODE_1D>::dshape_eulerian_at_knot(
 const unsigned &ipt,
 Shape &psi, 
 DShape &dpsidx) const
{
 //Call the geometrical shape functions and derivatives  
 double J = QElement<1,NNODE_1D>::dshape_eulerian_at_knot(ipt,psi,dpsidx);

 //Return the jacobian
 return J;
}


//======================================================================
// Set the data for the number of Variables at each node
//======================================================================
template<>
const unsigned QHarmonicElement<4>::Initial_Nvalue[4]={1,1,1,1};
template<>
const unsigned QHarmonicElement<3>::Initial_Nvalue[3]={1,1,1};
template<>
const unsigned QHarmonicElement<2>::Initial_Nvalue[2]={1,1};


//==start_of_problem_class============================================
/// 1D Harmonic problem in unit interval.
//====================================================================
template<class ELEMENT,class EIGEN_SOLVER> 
class HarmonicProblem : public Problem
{
public:

 /// Constructor: Pass number of elements and pointer to source function
 HarmonicProblem(const unsigned& n_element);

 /// Destructor (empty)
 ~HarmonicProblem(){delete this->mesh_pt(); delete this->eigen_solver_pt();}

 /// Update the problem specs before solve: (Re)set boundary conditions
 void actions_before_newton_solve() {}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}

 /// Solve the problem
 void solve(const unsigned &label);

 /// Doc the solution, pass the number of the case considered,
 /// so that output files can be distinguished.
 void doc_solution(const unsigned& label);

}; // end of problem class



//=====start_of_constructor===============================================
/// Constructor for 1D Harmonic problem in unit interval.
/// Discretise the 1D domain with n_element elements of type ELEMENT.
/// Specify function pointer to source function. 
//========================================================================
template<class ELEMENT,class EIGEN_SOLVER>
HarmonicProblem<ELEMENT,EIGEN_SOLVER>::HarmonicProblem(
 const unsigned& n_element)
{ 
 this->eigen_solver_pt() = new EIGEN_SOLVER;

 //Set domain length 
 double L=1.0;

 // Build mesh and store pointer in Problem
 Problem::mesh_pt() = new OneDMesh<ELEMENT>(n_element,L);
 
 // Set the boundary conditions for this problem: By default, all nodal
 // values are free -- we only need to pin the ones that have 
 // Dirichlet conditions. 

 // Pin the single nodal value at the single node on mesh 
 // boundary 0 (= the left domain boundary at x=0)
 mesh_pt()->boundary_node_pt(0,0)->pin(0);
 
 // Pin the single nodal value at the single node on mesh 
 // boundary 1 (= the right domain boundary at x=1)
 mesh_pt()->boundary_node_pt(1,0)->pin(0);
 
 // Setup equation numbering scheme
 assign_eqn_numbers();

} // end of constructor



//===start_of_doc=========================================================
/// Doc the solution in tecplot format. Label files with label.
//========================================================================
template<class ELEMENT,class EIGEN_SOLVER>
void HarmonicProblem<ELEMENT,EIGEN_SOLVER>::doc_solution(const unsigned& label)
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution with specified number of plot points per element
 sprintf(filename,"soln%i_on_proc%i.dat",label,
         this->communicator_pt()->my_rank());
 some_file.open(filename);
 mesh_pt()->output(some_file,npts);
 some_file.close();

} // end of doc

//===================================================================
/// Solve the eigenproblem 
//===================================================================
template<class ELEMENT,class EIGEN_SOLVER>
void HarmonicProblem<ELEMENT,EIGEN_SOLVER>::
solve(const unsigned& label)
{ 
 //Set external storage for the eigenvalues
 Vector<complex<double> > eigenvalue;
 //Set external storage for the eigenvectors
 Vector<DoubleVector> eigenvector_real;
 Vector<DoubleVector> eigenvector_imag;
 //Desired number eigenvalues
 unsigned n_eval=4;
 
 //Solve the eigenproblem
 this->solve_eigenproblem(n_eval,eigenvalue,eigenvector_real,eigenvector_imag);


 //We now need to sort the output based on the size of the real part
 //of the eigenvalues.
 //This is because the solver does not necessarily sort the eigenvalues
 Vector<complex<double> > sorted_eigenvalue = eigenvalue;
 sort(sorted_eigenvalue.begin(),sorted_eigenvalue.end(),
      ComplexLess<double>());

 //Read out the second smallest eigenvalue
 complex<double> temp_evalue = sorted_eigenvalue[1];
 unsigned second_smallest_index=0;
 //Loop over the unsorted eigenvalues and find the entry that corresponds
 //to our second smallest eigenvalue.
 for(unsigned i=0;i<eigenvalue.size();i++)
  {
   //Note that equality tests for doubles are bad, but it was just
   //sorted data, so should be fine
   if(eigenvalue[i] == temp_evalue) {second_smallest_index=i; break;}
  }
 
 //Normalise the eigenvector 
 {
  //Get the dimension of the eigenvector
  unsigned dim = eigenvector_real[second_smallest_index].nrow_local();
  double length=0.0;
  //Loop over all the entries
  for(unsigned i=0;i<dim;i++)
   {
    //Add the contribution to the length
    length += std::pow(eigenvector_real[second_smallest_index][i],2.0);
   }
//Add contributions from all processors
#ifdef OOMPH_HAS_MPI
  double length2 = length;
  if (eigenvector_real[second_smallest_index].distributed() && 
      eigenvector_real[second_smallest_index].distribution_pt()
      ->communicator_pt()->nproc() > 1)
   {
     MPI_Allreduce(&length,&length2,1,MPI_DOUBLE,MPI_SUM,
                   eigenvector_real[second_smallest_index].distribution_pt()
                   ->communicator_pt()->mpi_comm());
   }
  length = length2;
#endif
  //Now take the magnitude
  length = sqrt(length);
  
  //Fix the sign
  int sign=1;
  if(this->communicator_pt()->my_rank()==0)
   {
    if(eigenvector_real[second_smallest_index][0] < 0) {sign = -1;}
   }

#ifdef OOMPH_HAS_MPI
  //Broadcast to all other processes
  MPI_Bcast(&sign,1,MPI_INT,0,this->communicator_pt()->mpi_comm());
#endif

  //Finally normalise
  for(unsigned i=0;i<dim;i++)
   {
    eigenvector_real[second_smallest_index][i] /= length*sign;
   }
 }

 //Now assign the second eigenvector to the dofs of the problem
 this->assign_eigenvector_to_dofs(eigenvector_real[second_smallest_index]);
 //Output solution for this case (label output files with "1")
 this->doc_solution(label);

 char filename[100];
 sprintf(filename,"eigenvalues%i_on_proc%i.dat",label,
         this->communicator_pt()->my_rank());
 
 //Open an output file for the sorted eigenvalues
 ofstream evalues(filename);
 for(unsigned i=0;i<n_eval;i++)
  {
   //Print to screen
   oomph_info << sorted_eigenvalue[i].real() << " " 
              << sorted_eigenvalue[i].imag() << std::endl;
   //Send to file
   evalues << sorted_eigenvalue[i].real() << " " 
           << sorted_eigenvalue[i].imag() << std::endl;
  }
 
 evalues.close();
}
 

/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//======start_of_main==================================================
/// Driver for 1D Poisson problem
//=====================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
 MPI_Helpers::init(argc,argv);
#endif

#ifdef OOMPH_HAS_TRILINOS
 // Set up the problem: 
 unsigned n_element=100; //Number of elements

//Solve with Anasazi (non distributed)
//Note that the linear algebra will still be distributed
 clock_t t_start1 = clock();
 {
  HarmonicProblem<QHarmonicElement<3>,ANASAZI> problem(n_element);
  
  problem.solve(1);
 }
 clock_t t_end1 = clock();


 //Solve with Anasazi (distributed)
 clock_t t_start2 = clock();
 {

  HarmonicProblem<QHarmonicElement<3>,ANASAZI> problem(n_element);
  //In the validation case there should be two processors
  if(problem.communicator_pt()->nproc()==2)
   {
    //Setup a dummy partition
    unsigned n_element = problem.mesh_pt()->nelement();
    Vector<unsigned> element_partition(n_element);
    for(unsigned e=0;e<n_element/2;e++) {element_partition[e]=0;}
    for(unsigned e=n_element/2;e<n_element;e++) {element_partition[e]=1;}
    problem.distribute(element_partition);
   }
  else
   {
    problem.distribute();
   }
  problem.solve(2);
 }
 clock_t t_end2 = clock();

 oomph_info << "ANASAZI TIME (non-distributed): " 
            << (double)(t_end1 - t_start1)/CLOCKS_PER_SEC
            << std::endl;


 oomph_info << "ANASAZI TIME (distributed): " << 
  (double)(t_end2 - t_start2)/CLOCKS_PER_SEC
            << std::endl;
#endif

#ifdef OOMPH_HAS_MPI
 MPI_Helpers::finalize();
#endif

} // end of main










