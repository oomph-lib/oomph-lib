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
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>

using namespace std;



 /// Helper function -- ignores specified number of ints from input stream
 void ignore_ints(std::ifstream& fig_file, const unsigned& nints=1)
  {
   int dummy;
   for (unsigned i=0;i<nints;i++)
    {
     fig_file >> dummy;
    }
  }
 
 /// Helper function -- ignores specified number of floats from input stream
 void ignore_floats(std::ifstream& fig_file, const unsigned& nfloats=1)
  {
   float dummy;
   for (unsigned i=0;i<nfloats;i++)
    {
     fig_file >> dummy;
    }
  }

/// Helper function -- ignores specified number of floats from stringstream
 void ignore_ints(std::istream& fig_file, const unsigned& nints=1)
  {
   int dummy;
   for (unsigned i=0;i<nints;i++)
    {
     fig_file >> dummy;
    }
  }
 
 /// Helper function -- ignores specified number of floats from stringstream
 void ignore_floats(std::istream& fig_file, const unsigned& nfloats=1)
  {
   float dummy;
   for (unsigned i=0;i<nfloats;i++)
    {
     fig_file >> dummy;
    }
  }



/// /////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////



//====================================================================
/// fig2poly: Converts xfig drawings into input files for
/// Jonathan Shewczuk's triangle code.
//====================================================================
int main(int argc, char* argv[])
{
        
 // Check number of command line arguments: Need exactly one.
 if (argc!=2)
  {
   cout << "Wrong number of command line arguments. " << std::endl;
   cout << "Must specify the single file name that " << std::endl;
   cout << "specifies the xfig file" << std::endl;
   exit(1);
  }

 // Open fig file
 ifstream xfig_file(argv[1]);

 // Storage for the coordinates of the points on the various
 // boundaries: boundary_point_coords[i_boundary][i_point][i_coordinate]
 vector<vector<vector<double> > > boundary_point_coords;

 // Storage for the coordinates of the points in the 
 // holes: hole_coords[i_hole][i_coordinate]
 vector<vector<double> > hole_coords;
 
 // Temporary storage for input string
 std::string input_string;

 // *.fig file must conform to "Fig Format 3.2"
 // See http://en.wikipedia.org/wiki/Xfig [xfig is no
 // longer actively maintained (listed as dormant)
 // which is a shame. The user manual used to be at
 // http : // www . xfig . org / userman
 // (minus the spaces obviously -- added to stop the broken
 // link checker from complaining...)]

 // Read initial comment lines 
 //---------------------------

 // First line contains the version of xfig format
 getline(xfig_file,input_string,'\n');
 std::string version_string=input_string.substr(0,8);
 if (version_string!="#FIG 3.2")
  {
   std::cout << "*.fig file is in the wrong format!\n"
             << "It must conform to Fig Format 3.2\n"
             << "but the version string is: " << version_string 
             << std::endl;
   exit(1);
  }


 // Ignore a few other irrelevant header lines
 for (unsigned i=0;i<7;i++)
  {
   // Read entire line
   getline(xfig_file,input_string,'\n');
  }

 // The next two ints specify the resolution and some flag to do with the
 // coordinate system.
 int scaling_factor,dummy;
 xfig_file >> scaling_factor;
 xfig_file >> dummy;



 // Keep looping over objects and brutally write all lines into stringstream
 //-------------------------------------------------------------------------
 std::ostringstream buffer;
 while (!xfig_file.eof())
  {
   // Read the line
   getline(xfig_file,input_string,'\n');
   // Copy into buffer but do not end a "\n" at the end -- this will allow
   // us to detect the end of the data (the original fig file contains
   // a newline on the last line and this confuses the read loop)
   buffer << " " << input_string;
  }

 // Remove trailing space
 std::string aux_string(buffer.str(),0,buffer.str().length()-1);

  // Convert buffer into input string stream
 std::istringstream buffer2(aux_string);

 // Keep looping over objects
 //--------------------------
 while (!buffer2.eof())
  {

   // Now read the object identifier (first entry in object line)
   int object_id;
   buffer2 >> object_id;


   // Circle/ellipse object id 2
   //---------------------------
   if (object_id==1)
    {
     cout 
      << "I've read a circle -- it will be used as a hole identifier." 
      << std::endl;

     // Ignore what's not needed
     ignore_ints(buffer2,8);
     ignore_floats(buffer2,1);
     ignore_ints(buffer2,1);
     ignore_floats(buffer2,1);

     
     // Here's a new point in a hole
     vector<double> new_hole(2);
   
     // Now read in the x/y coordinate of centre
     int ix,iy;
     buffer2 >> ix;
     buffer2 >> iy;

     // Convert to scaled float
     new_hole[0]=double(ix)/double(scaling_factor);
     new_hole[1]=double(iy)/double(scaling_factor);

     // Stick into big container
     hole_coords.push_back(new_hole);
   
     // Ignore the remaining 6 ints
     ignore_ints(buffer2,6);

    }
   // Polyline has object id 2
   //-------------------------
   else if (object_id==2)
    {
     cout << "I've read a polyline -- it defines a boundary."
          << std::endl;
     
     // Next line: Type of polyline
     int sub_type;
     buffer2 >> sub_type;

     // Can only handle actual polylines (rather than polygons or boxes)
     if (sub_type==1)
      {
       // Ignore what's not needed
       ignore_ints(buffer2,7);
       ignore_floats(buffer2,1);
       ignore_ints(buffer2,5);

       // Finally: Number of points on line:
       int npoints;
       buffer2 >> npoints;

       // Here's another boundary
       vector<vector<double> > new_boundary(npoints);
              
       // Now read in the x/y coordinates
       int ix,iy;
       for (unsigned i=0;i<unsigned(npoints);i++)
        {
         buffer2 >> ix;
         buffer2 >> iy;
         new_boundary[i].resize(2);
         new_boundary[i][0]=double(ix)/double(scaling_factor);
         new_boundary[i][1]=double(iy)/double(scaling_factor);
        }

       // Stick into big container
       boundary_point_coords.push_back(new_boundary);
      }
     else
      {
       std::cout
        << "Can't handle this sub-type of polygon: " << sub_type << std::endl
        << "Can only do open polylines -- no closed boxes, etc. " << std::endl;
       exit(1);
      }

    }
   // Other object IDs can't be handled
   //----------------------------------
   else
    {
     std::cout
      << "Can't handle this object id: " << object_id 
      << "\nMake sure your figure only contains: \n"
      << "- polylines (which define boundaries) \n"
      << "- circles/ellipses (which define hole points) \n" 
      << "Also make sure that you've broken up any\n"
      << "compound objects. "  << std::endl;
     exit(1);
    }
   
  } // end of loop over objects



 // Get ready for output
 //---------------------
 
 // Get total number of points on all boundaries
 unsigned total_npoints=0;
 unsigned nbound=boundary_point_coords.size();
 for (unsigned b=0;b<nbound;b++)
  {
   unsigned npoints=boundary_point_coords[b].size();
   total_npoints+=npoints;
  }

 // Output points in pslg format
 //-----------------------------
 char filename[100];
 sprintf(filename,"%s.poly",argv[1]);
 ofstream poly_file(filename);
  {
   poly_file 
    << total_npoints 
    << " 2 0 0 # of pts, 2D, no attributes or boundary markers for points" 
    << std::endl;
  }
  

  unsigned count=1;
  for (unsigned b=0;b<nbound;b++)
   {
    unsigned npoints=boundary_point_coords[b].size();
    for (unsigned i=0;i<npoints;i++)
     {
      // Global node number, x, y, [additional comment that helps
      // to identify which boundary the node is located on -- mainly
      // for debugging
      poly_file << count << " " << boundary_point_coords[b][i][0]
                << " " << boundary_point_coords[b][i][1] 
                << " #  [located on boundary " << b << "]" << std::endl;
            
      // Increment counter
      count++;
     }
   }
  
  // End of point block
  poly_file << "# END_OF_NODE_BLOCK" << std::endl;

  // Output boundaries in poly format
  //---------------------------------
  
  
  // Number of edges, number of boundaries 
  poly_file << total_npoints << " 1 " //<< nbound 
            << " #   [number of segments (i.e. boundary edges),"
            << " flag for using boundary markers here]" 
            << std::endl;

  
  // Set segment bounded flag:
  
  // Initialise counters
  unsigned edge_count=1;
  unsigned npoints_in_all_previous_boundaries=0;
  
  // Loop over boundaries
  for (unsigned b=0;b<nbound;b++)
   {
    unsigned npoints=boundary_point_coords[b].size();
    // Up to the last point, each point on the boundary connects
    // to the next one
    for (unsigned i=0;i<npoints-1;i++)
     {
      // Edge number, first node index, second node index (in global 
      // node numbering scheme above), boundary id
      poly_file << edge_count << " " 
                << npoints_in_all_previous_boundaries+i+1 << " " 
                << npoints_in_all_previous_boundaries+i+2 << " " 
                << " " << b+1 << std::endl;

      // Increment counter
      edge_count++;
     }

    // Final point connects back to the first one
    // Segment number, first node index, second node index (in global 
    // node numbering scheme above), boundary id 
     poly_file << edge_count << " " 
               << npoints_in_all_previous_boundaries+npoints << " " 
               << npoints_in_all_previous_boundaries+1 << " " 
               << " " << b+1 << std::endl;
     
     // Bump up counter
     edge_count++;
     
     // Increment counter for points located in all previous boundaries
     npoints_in_all_previous_boundaries+=npoints;
   }
  
  // End of segment block
  poly_file << "# END_OF_SEGMENT_BLOCK" << std::endl;

  
  // Output "points in hole" in pslg format
  //---------------------------------------
  
  // Vector of points that are located in the holes in the domain
  // [NOTE: Historically (and in Gemma's thesis) this used to be 
  // called, somewhat confusingly, Hole_ID_pt]
  vector<vector<double> > Point_in_hole;
  
  // Number of holes
  unsigned nhole=hole_coords.size();
  Point_in_hole.resize(nhole);
  unsigned hole_count=1;
  poly_file << nhole << " #  [number of holes] " << std::endl;
  for (unsigned i=0;i<nhole;i++)
  {
   poly_file << hole_count << " " 
             << hole_coords[i][0] << " " 
             << hole_coords[i][1] << std::endl;
   
   //Add the point to the vector of hole IDs
   Point_in_hole[i].resize(2);
   Point_in_hole[i][0]= hole_coords[i][0];
   Point_in_hole[i][1]= hole_coords[i][1];
   
   // Increment counter
   hole_count++;
  }
  
  
  // End of hole block
  poly_file << "# END_OF_HOLE_BLOCK" << std::endl;

 poly_file.close();

   
} 

