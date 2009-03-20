#include<iostream>
#include<string>
#include<vector>
#include<set>
#include<list>
#include<cassert>
#include<fstream>
#include<sstream>
#include <stdlib.h>

using namespace std;

int main(int argc, char* argv[])
{



 // Geometry definition: At some point read this from an input file
 //================================================================

 double x_centre_inflow=0.0;
 double y_centre_inflow=0.0;
 double z_inflow_threshold=-0.9;
 double z_inflow=-1.0;

 double x_branch_divider=0.0;
 double z_outflow_threshold=1.9;
 double z_outflow=2.0;

 double x_centre_left_outflow=0.5;
 double y_centre_left_outflow=0.0;

 double x_centre_right_outflow=-0.5;
 double y_centre_right_outflow=0.0;


 // Flags for debugging of faces: Build them up slowly
 //===================================================
 bool use_inflow_faces=true;
 bool use_left_outflow_faces=true;
 bool use_right_outflow_faces=true;


 string input_string;
 char filename[100];
  

 // Get name of input file
 std::vector<string> file_name(2);
 if (argc==1)
  {
   file_name[0]="potential_boundary_outer_fluid";
   file_name[1]="potential_boundary_outer_solid";
  }
 else if (argc==3)
  {
   file_name[0]=argv[1];
   file_name[1]=argv[2];
  }
 else
  {
   std::cout << "Provide two (or no) input arguments, indicating stem of\n "
             << "tecplot boundary file" << std::endl;
   assert(false);
  }
 

 // Provide storage
 std::vector<unsigned> nnod(2);

 // Prepare storage for ordered nodes on in and outflow
 std::vector<std::list<std::pair<double,unsigned> > > 
  inflow_below_centre(2);
 std::vector<std::list<std::pair<double,unsigned> > > 
  inflow_above_centre(2);
 std::vector<std::list<std::pair<double,unsigned> > > 
  left_outflow_below_centre(2);
 std::vector<std::list<std::pair<double,unsigned> > >
  left_outflow_above_centre(2);
 std::vector<std::list<std::pair<double,unsigned> > > 
  right_outflow_below_centre(2);
 std::vector<std::list<std::pair<double,unsigned> > > 
  right_outflow_above_centre(2);
    
 // Nodal coordinates
 std::vector<std::vector<std::vector<double> > > x(2);

 // Counter for facets
 unsigned facet_count=0;
 unsigned fluid_facet_count=0;
 unsigned solid_facet_count=0;

 // Offset for node numbers once we've gone around
 unsigned node_offset=0;
  
 // Temporary output to string stream because we don't
 // know the number of factets yet.
 std::vector<std::ostringstream*> poly_file_stream_pt(2);
 
 // Stream for extra fluid/solid faces
 std::ostringstream* extra_fluid_poly_file_stream_pt= new ostringstream;
 std::ostringstream* extra_solid_poly_file_stream_pt= new ostringstream;
 

 // Loop over inner (fluid) and outer (solid) boundaries
 for (unsigned in_out=0;in_out<2;in_out++)
  {

   // Open temporary output stream
   poly_file_stream_pt[in_out]=new ostringstream;

   // Open input file
   sprintf(filename,"%s.dat",file_name[in_out].c_str());
   ifstream* input_file_pt=new ifstream(filename);
   
   
   // Check if it's been opened succesfully
   if (input_file_pt!=0)
    {
     cout << "Have opened input file: " << filename << std::endl;
    }  
   else
    {
     cout << "ERROR while trying to open input file: " << filename 
          << std::endl;
     assert(false);
    }
   
   
   // Process input file
   //===================
   
   
   // Skip first seven lines
   //-----------------------
   for (unsigned i=0;i<7;i++)
    {
     // Ignore rest of line
     input_file_pt->ignore(80,'\n');
    }  
   
   
   // Get Number of nodes
   //--------------------
   getline(*input_file_pt,input_string,'=');
   getline(*input_file_pt,input_string,',');
   nnod[in_out]=atoi(input_string.c_str());
   std::cout << "Number of nodes: " << nnod[in_out] << std::endl;
   
   
   // Get Number of elements
   //-----------------------
   getline(*input_file_pt,input_string,'=');
   getline(*input_file_pt,input_string,',');
   unsigned nelement=atoi(input_string.c_str());
   std::cout << "Number of elements: " << nelement << std::endl;
   
   // Ignore rest of line
   input_file_pt->ignore(80,'\n');
   
   
   
   // Skip next two lines
   //--------------------
   for (unsigned i=0;i<2;i++)
    {
     // Ignore rest of line
     input_file_pt->ignore(80,'\n');
    }  
   
   
   // Read in nodes
   //--------------
   x[in_out].resize(nnod[in_out]);
   for (unsigned j=0;j<nnod[in_out];j++)
    {
     x[in_out][j].resize(3);
     
     // Ignore leading blank
     getline(*input_file_pt,input_string,' ');
     
     // Read coordinates
     for (unsigned i=0;i<3;i++)
      {
       getline(*input_file_pt,input_string,' ');
       x[in_out][j][i]=atof(input_string.c_str());
      }
     
     // Ignore rest of line
     input_file_pt->ignore(80,'\n');
     
     
     // Inflow
     //-------
     double x_centre=x_centre_inflow;
     double y_centre=y_centre_inflow;
     if (x[in_out][j][2]<z_inflow_threshold)
      {
       // Below centre
       if (x[in_out][j][1]<y_centre)
        {
         inflow_below_centre[in_out].push_back(make_pair(x[in_out][j][0],j));
        }
       // Above centre
       else
        {
         inflow_above_centre[in_out].push_back(make_pair(x[in_out][j][0],j));
       }
  
      // Adjust to make sure it's exactly in outflow plane
      x[in_out][j][2]=z_inflow;
     }
  
   
    // Left outflow 
    //-------------
    x_centre=x_centre_left_outflow;
    y_centre=y_centre_left_outflow;
    if ((x[in_out][j][2]>z_outflow_threshold)&&
        (x[in_out][j][0]>x_branch_divider))
     {
      // Below centre
      if (x[in_out][j][1]<y_centre)
       {
        left_outflow_below_centre[in_out].push_back(
         make_pair(x[in_out][j][0],j));
       }
      // Above centre
      else
       {
        left_outflow_above_centre[in_out].push_back(
         make_pair(x[in_out][j][0],j));
       }
  
      // Adjust to make sure it's exactly in outflow plane
      x[in_out][j][2]=z_outflow;
     }

    // Right outflow
    //--------------
    x_centre=x_centre_right_outflow;
    y_centre=y_centre_right_outflow;
    if ((x[in_out][j][2]>z_outflow_threshold)&&
        (x[in_out][j][0]<x_branch_divider))
     {

      // Below centre
      if (x[in_out][j][1]<y_centre)
       {
        right_outflow_below_centre[in_out].push_back(
         make_pair(x[in_out][j][0],j));
       }
      // Above centre
      else
       {
        right_outflow_above_centre[in_out].push_back(
         make_pair(x[in_out][j][0],j));
       }
  
      // Adjust to make sure it's exactly in outflow plane
      x[in_out][j][2]=z_outflow;
     }
       
  }


 // Sort nodes on in/outflow boundaries
 inflow_below_centre[in_out].sort();
 inflow_above_centre[in_out].sort();
 inflow_above_centre[in_out].reverse();

 left_outflow_below_centre[in_out].sort();
 left_outflow_above_centre[in_out].sort();
 left_outflow_above_centre[in_out].reverse();
 
 right_outflow_below_centre[in_out].sort();
 right_outflow_above_centre[in_out].sort();
 right_outflow_above_centre[in_out].reverse();
 

 // Read in and doc connectivity
 //-----------------------------

 // Offset for "normal" faces -- in/outflow come first
 unsigned in_out_offset=4;
 
 std::set<unsigned> quad_node;
 for (unsigned j=0;j<nelement;j++)
  {
   // Ignore leading blank
   getline(*input_file_pt,input_string,' ');

   vector<unsigned> tmp(4);

   // Read node IDs
   for (unsigned i=0;i<3;i++)
    {
     getline(*input_file_pt,input_string,' ');
     unsigned k=atoi(input_string.c_str());
     tmp[i]=k;
     quad_node.insert(k);
    }
   
   getline(*input_file_pt,input_string,'\n');
   unsigned k=atoi(input_string.c_str());
   tmp[3]=k;
   quad_node.insert(k);


   
   bool use_bounding_faces=true;
   if (use_bounding_faces)
    {
     unsigned n_distinct=quad_node.size();
     switch (n_distinct)
      {

      case 4:
       (*poly_file_stream_pt[in_out]) 
        << "1 0 " << facet_count+in_out_offset 
        << " # one polygon, no hole, boundary marker\n3 ";
       (*poly_file_stream_pt[in_out]) << tmp[0]+node_offset << " " 
                                      << tmp[1]+node_offset << " " 
                                      << tmp[2]+node_offset << " " 
                                      << std::endl;
       facet_count++;
       
       (*poly_file_stream_pt[in_out]) 
        << "1 0 " 
        << facet_count+in_out_offset 
        << " # one polygon, no hole, boundary marker\n3 ";
       (*poly_file_stream_pt[in_out]) << tmp[0]+node_offset << " " 
                                      << tmp[2]+node_offset << " " 
                                      << tmp[3]+node_offset << " " 
                                      << std::endl;
       facet_count++;
       
       
       break;
       
      case 3:
       
       (*poly_file_stream_pt[in_out]) 
        << "1 0 " << facet_count+in_out_offset 
        << " # one polygon, no hole, boundary marker\n3 ";
       for (std::set<unsigned>::iterator it=quad_node.begin();
            it!=quad_node.end();it++)
        {
         (*poly_file_stream_pt[in_out]) << (*it)+node_offset << " ";
        }
       (*poly_file_stream_pt[in_out]) << std::endl;
       facet_count++;
       
       break;
       
       
      default:
       
       cout << "VERY ODD: n_distinct=" << n_distinct << " " << j << std::endl;
       assert(false);
       
      }
    }
   // Reset
   quad_node.clear();
   
  }
 
 // Keep track of number of facets and bump up node offset
 if (in_out==0)
  {
   fluid_facet_count=facet_count;
   node_offset=nnod[0];
  }
 else
  {
   solid_facet_count=facet_count;
  }
 
 input_file_pt->close();
 


  } // end over fluid/solid boundaries
 
 
 
 
 // Write poly file for fluid
 //==========================
 {
  unsigned in_out=0;
  ofstream poly_file("fsi_potential_bifurcation_fluid.poly");
  poly_file << "#Node list" << std::endl;
  poly_file << nnod[0] << " 3 0 1 # number of nodes; 3D "
            << " no attributes but boundary markers " << std::endl;
  for (unsigned j=0;j<nnod[in_out];j++)
   {
    poly_file << j+1 << " "; 
    for (unsigned i=0;i<3;i++)
     {
      poly_file << x[in_out][j][i] << " ";
     }
    poly_file << " 0 " << std::endl;
   }


  
  // Inflow faces
  //-------------
  if (use_inflow_faces)
   {
    // Inflow
    (*extra_fluid_poly_file_stream_pt) 
     << "1 0 1 # one polygon, no hole, boundary marker\n";
    unsigned n=inflow_below_centre[in_out].size()+
     inflow_above_centre[in_out].size();
    (*extra_fluid_poly_file_stream_pt) << n << " ";
    for (std::list<std::pair<double,unsigned> >::iterator 
          it=inflow_below_centre[in_out].begin();
         it!=inflow_below_centre[in_out].end();it++)
     {
      (*extra_fluid_poly_file_stream_pt) << (*it).second+1 << " ";
     }
    for (std::list<std::pair<double,unsigned> >::iterator 
          it=inflow_above_centre[in_out].begin();
         it!=inflow_above_centre[in_out].end();it++)
     {
      (*extra_fluid_poly_file_stream_pt) << (*it).second+1 << " ";
     }
    (*extra_fluid_poly_file_stream_pt) << std::endl;
    fluid_facet_count++;
   }
  
  
  // Left outflow
  //-------------
  if (use_left_outflow_faces)
   {
    // Left outflow
    (*extra_fluid_poly_file_stream_pt) 
     << "1 0 2 # one polygon, no hole, boundary marker\n";
    unsigned n=left_outflow_below_centre[in_out].size()+
     left_outflow_above_centre[in_out].size();
    (*extra_fluid_poly_file_stream_pt) << n << " ";
    for (std::list<std::pair<double,unsigned> >::iterator 
          it=left_outflow_below_centre[in_out].begin();
         it!=left_outflow_below_centre[in_out].end();it++)
     {
      (*extra_fluid_poly_file_stream_pt) << (*it).second+1 << " ";
     }
    for (std::list<std::pair<double,unsigned> >::iterator 
          it=left_outflow_above_centre[in_out].begin();
         it!=left_outflow_above_centre[in_out].end();it++)
     {
      (*extra_fluid_poly_file_stream_pt) << (*it).second+1 << " ";
     }
    
    (*extra_fluid_poly_file_stream_pt) << std::endl;
    fluid_facet_count++;
   }
  
  // Right outflow
  //--------------
  if (use_right_outflow_faces)
   {
    
    // Right outflow
    (*extra_fluid_poly_file_stream_pt) 
     << "1 0 3 # one polygon, no hole, boundary marker\n";
    unsigned n=right_outflow_below_centre[in_out].size()+
     right_outflow_above_centre[in_out].size();
    (*extra_fluid_poly_file_stream_pt) << n << " ";
    for (std::list<std::pair<double,unsigned> >::iterator 
          it=right_outflow_below_centre[in_out].begin();
         it!=right_outflow_below_centre[in_out].end();it++)
     {
      (*extra_fluid_poly_file_stream_pt) << (*it).second+1 << " ";
     }
    for (std::list<std::pair<double,unsigned> >::iterator 
          it=right_outflow_above_centre[in_out].begin();
         it!=right_outflow_above_centre[in_out].end();it++)
     {
      (*extra_fluid_poly_file_stream_pt) << (*it).second+1 << " ";
     }
    
    (*extra_fluid_poly_file_stream_pt) << std::endl;
    fluid_facet_count++;
   }
  
  
  // Now process the facet list
  poly_file << "# Facet list" << std::endl;
  poly_file << fluid_facet_count 
            << " 1 # Number of facets; boundary markers" 
            << std::endl;
  
  // Fill in accumulated output from "normal facets"
  poly_file << (*poly_file_stream_pt[in_out]).str();
  poly_file << (*extra_fluid_poly_file_stream_pt).str();
  
  // Finish it off
  poly_file << "# Hole list \n0 \n#Attributes list\n0" << std::endl;
  poly_file.close();
  
  std::cout << "Done -- total number of in fluid faces: " 
            << fluid_facet_count << std::endl;
  
 }  



 
 // Write poly file for solid
 //==========================
 {
  ofstream poly_file("fsi_potential_bifurcation_solid.poly");
  poly_file << "#Node list" << std::endl;
  poly_file << nnod[0]+nnod[1] << " 3 0 1 # number of nodes; 3D "
            << " no attributes but boundary markers " << std::endl;
  unsigned node_count=0;
  for (unsigned io=0;io<2;io++)
   {
    for (unsigned j=0;j<nnod[io];j++)
     {
      poly_file << node_count+1 << " "; 
      for (unsigned i=0;i<3;i++)
       {
        poly_file << x[io][j][i] << " ";
       }
      poly_file << " 0 " << std::endl;
      node_count++;
     }
   }

  
  // Inflow faces
  //-------------
  if (use_inflow_faces)
   {
    // Inflow
    (*extra_solid_poly_file_stream_pt) 
     << "2 1 1 # two polygons, one hole, boundary marker\n";

    // Two polygons: Inflow in fluid and solid
    {
     unsigned offset=0;
     for (unsigned i=0;i<2;i++)
      {
       unsigned n=inflow_below_centre[i].size()+
        inflow_above_centre[i].size();
       (*extra_solid_poly_file_stream_pt) << n << " ";
       for (std::list<std::pair<double,unsigned> >::iterator 
             it=inflow_below_centre[i].begin();
            it!=inflow_below_centre[i].end();it++)
        {
         (*extra_solid_poly_file_stream_pt) << (*it).second+1+offset << " ";
        }
       for (std::list<std::pair<double,unsigned> >::iterator 
             it=inflow_above_centre[i].begin();
            it!=inflow_above_centre[i].end();it++)
        {
         (*extra_solid_poly_file_stream_pt) << (*it).second+1+offset << " ";
        }
       (*extra_solid_poly_file_stream_pt) << std::endl;
       offset+=nnod[i];
      }
    }    
    // Here's the one and only hole
    (*extra_solid_poly_file_stream_pt) << "1 " 
                                   <<  x_centre_inflow  << " "
                                   <<  y_centre_inflow  << " "
                                   <<  z_inflow  << " "
                                   << std::endl;
    solid_facet_count++;
   }
  
  
  // Left outflow
  //-------------
  if (use_left_outflow_faces)
   {
    // Left outflow
    (*extra_solid_poly_file_stream_pt) 
     << "2 1 2 # two polygons, one hole, boundary marker\n";
    unsigned offset=0;
    for (unsigned i=0;i<2;i++)
     {
      
      unsigned n=left_outflow_below_centre[i].size()+
       left_outflow_above_centre[i].size();
      (*extra_solid_poly_file_stream_pt) << n << " ";
      for (std::list<std::pair<double,unsigned> >::iterator 
            it=left_outflow_below_centre[i].begin();
           it!=left_outflow_below_centre[i].end();it++)
       {
        (*extra_solid_poly_file_stream_pt) << (*it).second+1+offset << " ";
       }
      for (std::list<std::pair<double,unsigned> >::iterator 
            it=left_outflow_above_centre[i].begin();
           it!=left_outflow_above_centre[i].end();it++)
       {
        (*extra_solid_poly_file_stream_pt) << (*it).second+1+offset << " ";
       }
      
      (*extra_solid_poly_file_stream_pt) << std::endl;
       offset+=nnod[i];
     }

    // Here's the one and only hole
    (*extra_solid_poly_file_stream_pt) << "1 " 
                                       <<  x_centre_left_outflow  << " "
                                       <<  y_centre_left_outflow  << " "
                                       <<  z_outflow  << " "
                                       << std::endl;
    solid_facet_count++;
   }
  


  // Right outflow
  //--------------
  if (use_right_outflow_faces)
   {
    // Left outflow
    (*extra_solid_poly_file_stream_pt) 
     << "2 1 3 # two polygons, one hole, boundary marker\n";
    unsigned offset=0;
    for (unsigned i=0;i<2;i++)
     {
      
      unsigned n=right_outflow_below_centre[i].size()+
       right_outflow_above_centre[i].size();
      (*extra_solid_poly_file_stream_pt) << n << " ";
      for (std::list<std::pair<double,unsigned> >::iterator 
            it=right_outflow_below_centre[i].begin();
           it!=right_outflow_below_centre[i].end();it++)
       {
        (*extra_solid_poly_file_stream_pt) << (*it).second+1+offset << " ";
       }
      for (std::list<std::pair<double,unsigned> >::iterator 
            it=right_outflow_above_centre[i].begin();
           it!=right_outflow_above_centre[i].end();it++)
       {
        (*extra_solid_poly_file_stream_pt) << (*it).second+1+offset << " ";
       }
      
      (*extra_solid_poly_file_stream_pt) << std::endl;
       offset+=nnod[i];
     }

    // Here's the one and only hole
    (*extra_solid_poly_file_stream_pt) << "1 " 
                                       <<  x_centre_right_outflow  << " "
                                       <<  y_centre_right_outflow  << " "
                                       <<  z_outflow  << " "
                                       << std::endl;
    solid_facet_count++;
   }
  
















//   // Right outflow
//   //--------------
//   if (use_right_outflow_faces)
//    {
    
//     // Right outflow
//     (*extra_solid_poly_file_stream_pt) 
//      << "1 0 3 # one polygon, no hole, boundary marker\n";
//     unsigned n=right_outflow_below_centre[in_out].size()+
//      right_outflow_above_centre[in_out].size();
//     (*extra_solid_poly_file_stream_pt) << n << " ";
//     for (std::list<std::pair<double,unsigned> >::iterator 
//           it=right_outflow_below_centre[in_out].begin();
//          it!=right_outflow_below_centre[in_out].end();it++)
//      {
//       (*extra_solid_poly_file_stream_pt) << (*it).second+1 << " ";
//      }
//     for (std::list<std::pair<double,unsigned> >::iterator 
//           it=right_outflow_above_centre[in_out].begin();
//          it!=right_outflow_above_centre[in_out].end();it++)
//      {
//       (*extra_solid_poly_file_stream_pt) << (*it).second+1 << " ";
//      }
    
//     (*extra_solid_poly_file_stream_pt) << std::endl;
//     solid_facet_count++;
//    }
  
  
  // Now process the facet list
  poly_file << "# Facet list" << std::endl;
  poly_file << solid_facet_count 
            << " 1 # Number of facets; boundary markers" 
            << std::endl;
  
  // Fill in accumulated output from "normal facets": Fluid and solid!
  poly_file << (*poly_file_stream_pt[0]).str();
  poly_file << (*poly_file_stream_pt[1]).str();
  poly_file << (*extra_solid_poly_file_stream_pt).str();

  // One hole
  poly_file << "# Hole list \n1 " << std::endl;
  poly_file << "1 " 
            <<  x_centre_inflow  << " "
            <<  y_centre_inflow  << " "
            <<  z_inflow+0.1  << " " // hierher!
            << std::endl;  
  poly_file << "\n#Attributes list\n0" << std::endl;
  poly_file.close();
  
  std::cout << "Done -- total number of in solid faces: " 
            << solid_facet_count << std::endl;
  
 }  


};

