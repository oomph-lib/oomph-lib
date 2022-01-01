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
// Create an output file that can be read with geomview
// Change the infile name and the output name to use it 


#include<iostream>
#include<fstream>
#include<vector>
#include<math.h>
using namespace std;

int main(int argc, char* argv[])
 {

 // Convert argument to strings that specify the input file name
 string mesh_file_name(argv[1]);
  
// Read the output mesh file to find informations about the nodes
// and elements of the mesh 

ifstream infile(mesh_file_name.c_str(), ios_base::in);
unsigned n_node;
infile>>n_node;
vector<double> x(n_node);
vector<double> y(n_node);
vector<int> vertinfo(n_node);
 for(unsigned i=0;i<n_node;i++)
  { 
   infile>>x[i];
   infile>>y[i];
   infile>>vertinfo[i];
  }
unsigned n_vx;
infile>>n_vx;
vector<int> nodecode(n_vx);
vector<int> icurv(n_vx);
vector<double> ucurv(n_vx);
 for(unsigned i=0;i<n_vx;i++)
  {
   infile>>nodecode[i];
   infile>>icurv[i];
   infile>>ucurv[i];
  }
unsigned n_local_node;
infile>>n_local_node;
unsigned n_element;
infile>>n_element;
unsigned b=n_local_node*n_element;
vector<int> v(b);
vector<int> edgeinfo(b);
unsigned k=0;
for(unsigned i=0;i<n_element;i++)
 {
  for(unsigned j=0;j<n_local_node;j++)
   {
    infile>>v[k];
    k++;
   }
 }
unsigned l=0;
for(unsigned i=0;i<n_element;i++)
 {
  for(unsigned j=0;j<n_local_node;j++)
   {
    infile>>edgeinfo[l];
    l++;
   }
 }

infile.close();

// Create a file of type ".quad" to visualize the mesh with geomview


 unsigned nn=0;
 char result[20];
 sprintf(result,"%s","mesh.quad");
 ofstream outfile(result,ios_base::out);
 outfile<<"QUAD"<<'\n';
 for(unsigned i=0;i<n_element;i++)
   {
    for(unsigned j=0;j<n_local_node;j++) 
     {
      outfile<<x[v[nn]-1]<<" "<<y[v[nn]-1]<<" 0 ";
      nn++; 
     }
    outfile<<'\n';
   } 
 outfile.close(); 



}  //end of main
