//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2021 Matthias Heil and Andrew Hazel
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
// Demonstrate use of geometric object as GeneralisedElement

 
// Generic oomph-lib headers
#include "generic.h"

// Circle as generalised element:
#include "circle_as_generalised_element.h"

using namespace std;

using namespace oomph;

/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////


//======start_of_problem==============================================
/// Problem to demonstrate the use of a GeomObject as a 
/// GeneralisedElement: A geometric object (a Circle) is "upgraded"
/// to a GeneralisedElement. The position of the Circle is
/// determined by a balance of forces, assuming that the
/// Circle is mounted on an elastic spring of specified
/// stiffness and loaded by a vertical "load". 
//====================================================================
class GeomObjectAsGeneralisedElementProblem : public Problem
{

public:

 /// Constructor
 GeomObjectAsGeneralisedElementProblem();

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve(){}
  
 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}

 /// Doc the solution
 void doc_solution();

 /// Return value of the "load" on the elastically supported ring
 double& load()
  {
   return *Load_pt->value_pt(0);
  }

 /// Access to DocInfo object
 DocInfo& doc_info() {return Doc_info;}

private:

 /// Trace file
 ofstream Trace_file;

 /// Pointer to data item that stores the "load" on the ring
 Data* Load_pt;

 /// Doc info object
 DocInfo Doc_info;

};





//=============================start_of_problem_constructor===============
/// Constructor
//========================================================================
GeomObjectAsGeneralisedElementProblem::GeomObjectAsGeneralisedElementProblem()
{ 
 
 // Set coordinates and radius for the circle
 double x_c=0.5;
 double y_c=0.0;
 double R=1.0;

 // Build GeomObject that's been upgraded to a GeneralisedElement
 // GeneralisedElement* 
 ElasticallySupportedRingElement* geom_object_element_pt = 
  new ElasticallySupportedRingElement(x_c,y_c,R);

 // Set the stiffness of the elastic support
 geom_object_element_pt->k_stiff()=0.3;

 // Build mesh
 mesh_pt()=new Mesh;

 // So far, the mesh is completely empty. Let's add the 
 // one (and only!) GeneralisedElement to it:
 mesh_pt()->add_element_pt(geom_object_element_pt);

 // Create the load (a Data object with a single value)
 Load_pt=new Data(1);
   
 // The load is prescribed so its one-and-only value is pinned
 Load_pt->pin(0);

 // Set the pointer to the Data object that specifies the 
 // load on the ring
 geom_object_element_pt->set_load_pt(Load_pt);
 
 // Setup equation numbering scheme.
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

 // Set output directory
 Doc_info.set_directory("RESLT"); 
  
 // Open trace file
 char filename[100];
 sprintf(filename,"%s/trace.dat",Doc_info.directory().c_str());
 Trace_file.open(filename);
 Trace_file << "VARIABLES=\"load\",\"y<sub>circle</sub>\"" << std::endl;

} // end of constructor




//===========================start_of_doc_solution========================
/// Doc the solution in tecplot format.
//========================================================================
void GeomObjectAsGeneralisedElementProblem::doc_solution()
{ 

 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts=100;

 // Lagrangian coordinate and position vector (both as vectors)
 Vector<double> zeta(1);
 Vector<double> r(2);
 
 // Output solution 
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 for (unsigned i=0;i<npts;i++)
  {
   zeta[0]=2.0*MathematicalConstants::Pi*double(i)/double(npts-1);
   static_cast<ElasticallySupportedRingElement*>(mesh_pt()->element_pt(0))->
    position(zeta,r);  
   some_file << r[0] << " " << r[1] << std::endl;
  }
 some_file.close();


 // Write "load" and vertical position of the ring's centre
 Trace_file 
  << static_cast<ElasticallySupportedRingElement*>(
   mesh_pt()->element_pt(0))->load()
  << " "
  << static_cast<ElasticallySupportedRingElement*>(
   mesh_pt()->element_pt(0))->y_c()
  << " "
  << std::endl;

} // end of doc_solution


 





/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////





//===============start_of_driver==========================================
/// Driver 
//========================================================================
int main()
{

 // Set up the problem 
 GeomObjectAsGeneralisedElementProblem problem;
 
 // Initial value for the load 
 problem.load()=-0.3;
 
 // Loop for different loads
 //-------------------------
 
 // Number of steps
 unsigned nstep=5;
 
 // Increment in load
 double dp=0.6/double(nstep-1);
 
 for (unsigned istep=0;istep<nstep;istep++)
  {    
   // Solve/doc
   problem.newton_solve();
   problem.doc_solution();
   
   //Increment counter for solutions 
   problem.doc_info().number()++;    
   
   // Change load on ring
   problem.load()+=dp;
  } 
 
} // end of driver


