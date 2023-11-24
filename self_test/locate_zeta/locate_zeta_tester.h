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

// Check locate zeta machinery
void check_locate_zeta(Mesh* mesh_pt)
{
 
 // Spatial dimension
 unsigned n_dim=mesh_pt->finite_element_pt(0)->dim();


 //Check elemental locate zeta:
 {
  FiniteElement* el_pt=mesh_pt->finite_element_pt(0);
  Vector<double> cog(n_dim);
  double max_radius=0.0;
  el_pt->get_centre_of_gravity_and_max_radius_in_terms_of_zeta(cog,max_radius);
  std::ofstream outfile;
  outfile.open("RESLT/single_element.dat");
  el_pt->output(outfile,el_pt->nnode_1d());
  outfile.close();
  if (n_dim==2)
   {
    outfile.open("RESLT/enclosing_circle_for_single_element.dat");
    unsigned nplot=100;
    outfile << "ZONE I=" << nplot << std::endl;
    for (unsigned i=0;i<nplot;i++)
     {
      double phi=2.0*MathematicalConstants::Pi*double(i)/double(nplot-1);
      outfile << cog[0]+max_radius*cos(phi) << " " 
              << cog[1]+max_radius*sin(phi) << std::endl;
     }
    outfile.close();
   }
  else if (n_dim==3)
   {
    outfile.open("RESLT/enclosing_sphere_for_single_element.dat");
    unsigned nplot=100;
    outfile << "ZONE I=" << nplot << ", J=" << nplot << std::endl;
    for (unsigned i=0;i<nplot;i++)
     {
      double phi=2.0*MathematicalConstants::Pi*double(i)/double(nplot-1);
      for (unsigned j=0;j<nplot;j++)
       {
        double theta=MathematicalConstants::Pi*(-0.5+double(j)/double(nplot-1));
        outfile << cog[0]+max_radius*cos(phi)*cos(theta) << " " 
                << cog[1]+max_radius*sin(phi)*cos(theta) << " "
                << cog[2]+max_radius*sin(theta) << " "
                << std::endl;
       }
     }
    outfile.close();
   }
  else
   {
    oomph_info << "wtf?" << std::endl;
    abort();
   } 
  
  // Check elemental locate zeta
  Vector<double> zeta(cog);
  GeomObject* geom_object_pt=0;
  Vector<double> s(n_dim);
  bool use_coordinate_as_initial_guess=false;
  el_pt->locate_zeta(zeta,
                     geom_object_pt, 
                     s,
                     use_coordinate_as_initial_guess);

  if (geom_object_pt==0)
   {
    oomph_info << "elemental locate zeta failed\n";
    abort();
   }
  else
   {
    oomph_info << "elemental locate zeta succeeded\n";

   }
 }


 // Create Mesh as geometric object 
 RefineableBinArrayParameters ref_bin_array_parameters(mesh_pt);
 NonRefineableBinArrayParameters non_ref_bin_array_parameters(mesh_pt);
#ifdef OOMPH_HAS_CGAL
 CGALSamplePointContainerParameters sample_point_container_parameters(mesh_pt);
#endif
 

 // Loop over the various sample point containers
 unsigned n_sample_point_container_type=2;
#ifdef OOMPH_HAS_CGAL
 n_sample_point_container_type++;
#endif
 for (unsigned ref_flag=0;ref_flag<n_sample_point_container_type;ref_flag++)
  {
   
   SamplePointContainerParameters* sample_point_container_params_pt=0;
   if (ref_flag==0)
    {
     sample_point_container_params_pt=&non_ref_bin_array_parameters;
     oomph_info << "\n\n\nDOING NON-REFINEABLE BIN ARRAY" << std::endl;
     oomph_info <<       "==============================" << std::endl;
    }
   else if (ref_flag==1)
    {
     sample_point_container_params_pt=&ref_bin_array_parameters;
     oomph_info << "\n\n\nDOING REFINEABLE BIN ARRAY" << std::endl;
     oomph_info <<       "==========================" << std::endl;
    }
#ifdef OOMPH_HAS_CGAL
   else if (ref_flag==2)
    {
     sample_point_container_params_pt=&sample_point_container_parameters;
     oomph_info << "\n\n\nDOING CGAL BASED CONTAINER" << std::endl;
     oomph_info <<       "==========================" << std::endl;
    }
#endif
   else 
    {
     oomph_info << "wrong container\n";
     abort();
    }
   
   MeshAsGeomObject* mesh_geom_obj_pt=new MeshAsGeomObject(sample_point_container_params_pt);
   

// hierher kill
/*    if (ref_flag==0) */
/*     { */
/*      sample_point_container_params_pt=&non_ref_bin_array_parameters; */
/*      oomph_info << "\n\n\nDOING NON-REFINEABLE BIN ARRAY" << std::endl; */
/*      oomph_info <<       "==============================" << std::endl; */
/*     } */
/*    else if (ref_flag==1) */
/*     { */
/*      sample_point_container_params_pt=&ref_bin_array_parameters; */
/*      oomph_info << "\n\n\nDOING REFINEABLE BIN ARRAY" << std::endl; */
/*      oomph_info <<       "==========================" << std::endl; */
/*     } */
/* #ifdef OOMPH_HAS_CGAL */
/*    else if (ref_flag==2) */
/*     { */
/*      sample_point_container_params_pt=&sample_point_container_parameters; */
/*      oomph_info << "\n\n\nDOING CGAL BASED CONTAINER" << std::endl; */
/*      oomph_info <<       "==========================" << std::endl; */
/*     } */
/* #endif */
/*    else  */
/*     { */
/*      oomph_info << "wrong container\n"; */
/*      abort(); */
/*     } */


   // Output bin structure
   {

    if (ref_flag==0)
     {
      std::string filename="RESLT/non_ref_bin.dat";
      ofstream outfile;
      outfile.open(filename.c_str());
      dynamic_cast<BinArray*>(mesh_geom_obj_pt->sample_point_container_pt())->
       output_bin_vertices(outfile);
      outfile.close();
     }
    else if (ref_flag==1)
     {
      std::string filename="RESLT/ref_bin.dat";
      ofstream outfile;
      outfile.open(filename.c_str());
      dynamic_cast<BinArray*>(mesh_geom_obj_pt->sample_point_container_pt())->
       output_bin_vertices(outfile);
      outfile.close();
     }
   }
   

   // CHECK LOTS OF POINTS INSIDE THE MESH
   //-------------------------------------
   {
    oomph_info << "\n\n\nLots of zetas:" << std::endl;
    oomph_info << "--------------\n\n" << std::endl;


    // Do lots of locate zeta calls, check the outcome and doc how
    // long it took. 
    unsigned count=0;
    Vector<double> s(n_dim); 
    Vector<double> zeta(n_dim);
    GeomObject* geom_obj_pt=0;
    Vector<double> s_from_locate_zeta(n_dim);
    Vector<double> zeta_test(n_dim);
    double max_error=0.0;
    double t_locate_zeta=0.0;
    unsigned nelem=mesh_pt->nelement();
    for (unsigned e=0;e<nelem;e++)
     {
      FiniteElement* el_pt=mesh_pt->finite_element_pt(e);
      Integral* int_pt=el_pt->integral_pt();
      unsigned nint=int_pt->nweight();
      for (unsigned j=0;j<nint;j++)
       { 
        // Find integration point in element
        Vector<double> s(n_dim);
        for (unsigned i=0;i<n_dim;i++)
         {
          s[i]=int_pt->knot(j,i);
         }
        el_pt->interpolated_zeta(s,zeta);
      
        // Now reverse mapping
        double start_t = TimingHelpers::timer();
        mesh_geom_obj_pt->locate_zeta(zeta,geom_obj_pt,s_from_locate_zeta);
        double end_t = TimingHelpers::timer();
        t_locate_zeta+=end_t-start_t;
        count++;

        if (geom_obj_pt==0)
         {
          oomph_info << "Search for zeta = ( " 
                     << zeta[0] << " "
                     << zeta[0] << " "
                     << " failed "
                     << std::endl;
          abort();
         }

        // Test
        geom_obj_pt->interpolated_zeta(s_from_locate_zeta,zeta_test);
      
        double error=0.0;
        for (unsigned i=0;i<n_dim;i++)
         {
          error+=pow(zeta[i]-zeta_test[i],2);
         }
        error=sqrt(error);
        if (error>max_error) max_error=error;
       }
     }
    oomph_info << "done/tested " << count << " located zetas. Max. error = "
               << max_error << std::endl
               << "Total time for locate zeta: " << t_locate_zeta
               << " sec " << std::endl;

   }

   // NOW DO ONE POINT BUT MAKE LOCATE ZETA FAIL AND DOC SPIRAL
   //----------------------------------------------------------
   {
    oomph_info << "\n\n\nEnforced failure for one zeta:" << std::endl;
    oomph_info <<       "------------------------------\n\n" << std::endl;

    BinArray::Always_fail_elemental_locate_zeta=true;
    std::stringstream filename;
    filename << "RESLT/sample_points_" << n_dim << "d_";
    if (ref_flag==0)
     {
      filename << "nonref_bin";
     }
    else if (ref_flag==1)
     {
      filename << "ref_bin";
     }
    else if (ref_flag==2)
     {
      filename << "cgal";
     }
    else
     {
      oomph_info << "wrong container\n";
      abort();
     }
    filename << ".dat";
    BinArray::Visited_sample_points_file.open(filename.str().c_str());
    Vector<double> zeta(n_dim,0.5);
    GeomObject* geom_obj_pt=0;
    Vector<double> s_from_locate_zeta(n_dim);
    Vector<double> zeta_test(n_dim,0.5);
    double t_locate_zeta=0.0;
    double start_t = TimingHelpers::timer();
    mesh_geom_obj_pt->locate_zeta(zeta,geom_obj_pt,s_from_locate_zeta);
    double end_t = TimingHelpers::timer();
    t_locate_zeta+=end_t-start_t;

    // Test
    oomph_info << "It took " << t_locate_zeta << " secs to ";
    if (geom_obj_pt!=0)
     {
      geom_obj_pt->interpolated_zeta(s_from_locate_zeta,zeta_test);
      oomph_info << "find: " 
                 << zeta_test[0] << " "  
                 << zeta_test[1] << " "  
                 << std::endl;
      
      unsigned nvisit=mesh_geom_obj_pt->sample_point_container_pt()->
       total_number_of_sample_points_visited_during_locate_zeta_from_top_level();
      unsigned ntotal=mesh_geom_obj_pt->sample_point_container_pt()->
       total_number_of_sample_points_computed_recursively();
      
      oomph_info << "Visited " << nvisit << " sample points " 
                 << " out of a total of " << ntotal << std::endl;
     }
    else
     {
      oomph_info << "not find: " 
                 << zeta_test[0] << " "  
                 << zeta_test[1] << " "  
                 << std::endl;

      unsigned nvisit=mesh_geom_obj_pt->sample_point_container_pt()->
       total_number_of_sample_points_visited_during_locate_zeta_from_top_level();
      unsigned ntotal=mesh_geom_obj_pt->sample_point_container_pt()->
       total_number_of_sample_points_computed_recursively();

      oomph_info << "Visited " << nvisit << " sample points " 
                 << " out of a total of " << ntotal << std::endl;

      if (nvisit!=ntotal)
       {
        oomph_info << "ERROR\n";
        abort();
       }

     }
    
    // Reset and close
    BinArray::Always_fail_elemental_locate_zeta=false;
    BinArray::Visited_sample_points_file.close();
   }
   

   // NOW DO ONE POINT BUT MAKE LOCATE ZETA FAIL BECAUSE OF
   //------------------------------------------------------
   // LIMITED RADIUS; DOC SPIRAL
   // --------------------------
   {
    oomph_info 
     << "\n\n\nEnforced failure for one zeta because of limited radius:"
     << std::endl;
    oomph_info 
     << "--------------------------------------------------------\n\n" 
     << std::endl;

    BinArray::Always_fail_elemental_locate_zeta=true;
    std::stringstream filename;
    filename << "RESLT/sample_points_limited_radius" << n_dim << "d_";
    if (ref_flag==0)
     {
      filename << "nonref_bin";
     }
    else if (ref_flag==1)
     {
      filename << "ref_bin";
     }
    else if (ref_flag==2)
     {
      filename << "cgal";
     }
    else
     {
      oomph_info << "wrong container\n";
      abort();
     }
    filename << ".dat";
    BinArray::Visited_sample_points_file.open(filename.str().c_str());
    Vector<double> zeta(n_dim);
    GeomObject* geom_obj_pt=0;
    Vector<double> s_from_locate_zeta(n_dim);
    Vector<double> zeta_test(n_dim);
    Vector<std::pair<double, double> > min_and_max_coordinates=
    mesh_geom_obj_pt->sample_point_container_pt()->
     min_and_max_coordinates();
    double min_dim=DBL_MAX;
    for (unsigned i=0;i<n_dim;i++)
     {
      zeta[i]=0.5*(min_and_max_coordinates[i].first+
                   min_and_max_coordinates[i].second);
      double dist=fabs(min_and_max_coordinates[i].first-
                       min_and_max_coordinates[i].second);
      if (dist<min_dim) min_dim=dist;
     }
    mesh_geom_obj_pt->sample_point_container_pt()->max_search_radius()=
     0.1*min_dim;

    oomph_info << "Search radius: " 
               << mesh_geom_obj_pt->sample_point_container_pt()->
     max_search_radius() << std::endl;

    double t_locate_zeta=0.0;
    double start_t = TimingHelpers::timer();
    mesh_geom_obj_pt->locate_zeta(zeta,geom_obj_pt,s_from_locate_zeta);
    double end_t = TimingHelpers::timer();
    t_locate_zeta+=end_t-start_t;

    // Test
    oomph_info << "It took " << t_locate_zeta << " secs to ";
    if (geom_obj_pt!=0)
     {
      geom_obj_pt->interpolated_zeta(s_from_locate_zeta,zeta_test);
      oomph_info << "find: " 
                 << zeta_test[0] << " "  
                 << zeta_test[1] << " "  
                 << std::endl;
      
      unsigned nvisit=mesh_geom_obj_pt->sample_point_container_pt()->
       total_number_of_sample_points_visited_during_locate_zeta_from_top_level();
      unsigned ntotal=mesh_geom_obj_pt->sample_point_container_pt()->
       total_number_of_sample_points_computed_recursively();
      
      oomph_info << "Visited " << nvisit << " sample points " 
                 << " out of a total of " << ntotal << std::endl;
     }
    else
     {
      oomph_info << "not find: " 
                 << zeta_test[0] << " "  
                 << zeta_test[1] << " "  
                 << std::endl;

      unsigned nvisit=mesh_geom_obj_pt->sample_point_container_pt()->
       total_number_of_sample_points_visited_during_locate_zeta_from_top_level();
      unsigned ntotal=mesh_geom_obj_pt->sample_point_container_pt()->
       total_number_of_sample_points_computed_recursively();

      oomph_info << "Visited " << nvisit << " sample points " 
                 << " out of a total of " << ntotal << std::endl;

      /* // Here we should only visit fewer than the lot... */
      /* if (nvisit>=ntotal) */
      /*  { */
      /*   oomph_info << "ERROR\n"; */
      /*   abort(); */
      /*  } */

     }
    
    // Reset and close
    BinArray::Always_fail_elemental_locate_zeta=false;
    BinArray::Visited_sample_points_file.close();
   }

   // Reset
   mesh_geom_obj_pt->sample_point_container_pt()->max_search_radius()=
    DBL_MAX;


   // NOW DO ONE POINT BUT MAKE LOCATE ZETA FAIL AND DOC INCREMENTAL SPIRAL
   //---------------------------------------------------------------------
   {
    oomph_info << "\n\n\nEnforced failure for one zeta (spiraling):" << std::endl;
    oomph_info <<       "------------------------------------------\n\n" << std::endl;

    Vector<double> zeta(n_dim,0.5);
    GeomObject* geom_obj_pt=0;

    // Get number of sample points
    unsigned nsample_points_in_bin_array=mesh_geom_obj_pt->
     sample_point_container_pt()->total_number_of_sample_points_computed_recursively();
    
    // Maximum spiral level within the cartesian bin structure
    // (to be assigned below)
    unsigned n_max_level=0; 

    
    // Initialise:
    double t_locate_zeta=0.0;

    // Refineable
    if (mesh_geom_obj_pt->sample_point_container_version()==
        UseRefineableBinArray)
     {         
      RefineableBinArray* sample_point_container_pt=
       dynamic_cast<RefineableBinArray*>(mesh_geom_obj_pt->
                                         sample_point_container_pt());
      
      sample_point_container_pt->
       last_sample_point_to_actually_lookup_during_locate_zeta() =
       sample_point_container_pt->
       initial_last_sample_point_to_actually_lookup_during_locate_zeta();
      sample_point_container_pt->
       first_sample_point_to_actually_lookup_during_locate_zeta() = 0;      
     }
    // Non-refineable
    else if (mesh_geom_obj_pt->sample_point_container_version()==
             UseNonRefineableBinArray)
     {
      NonRefineableBinArray* sample_point_container_pt=
       dynamic_cast<NonRefineableBinArray*>(mesh_geom_obj_pt->
                                            sample_point_container_pt());
      
      // Initialise spiral levels
      sample_point_container_pt->current_min_spiral_level()=0;
      sample_point_container_pt->current_max_spiral_level()=
       sample_point_container_pt->n_spiral_chunk()-1; 
      
      // Find maximum spiral level within the cartesian bin structure
      n_max_level=sample_point_container_pt->max_bin_dimension();
      
      // Limit it 
      if (sample_point_container_pt->current_max_spiral_level()>n_max_level)
       {
        sample_point_container_pt->current_max_spiral_level()=n_max_level-1;
       }
     }
#ifdef OOMPH_HAS_CGAL
    // CGAL
    else if (mesh_geom_obj_pt->sample_point_container_version()==
        UseCGALSamplePointContainer)
     {         
      CGALSamplePointContainer* sample_point_container_pt=
       dynamic_cast<CGALSamplePointContainer*>(mesh_geom_obj_pt->
                                         sample_point_container_pt());

      sample_point_container_pt->
       last_sample_point_to_actually_lookup_during_locate_zeta() =
       sample_point_container_pt->
       initial_last_sample_point_to_actually_lookup_during_locate_zeta();
      sample_point_container_pt->
       first_sample_point_to_actually_lookup_during_locate_zeta() = 0;      
     }
#endif
    else
     {
      oomph_info << "wrong sample point container\n";
      abort();
     }

    // Loop over "spirals/levels" away from the current position
    unsigned i_level = 0;
    bool has_not_reached_max_level_of_search=true;
    while ((geom_obj_pt==0)&&(has_not_reached_max_level_of_search))
     {
      oomph_info << "\n\nSearch level: " << i_level << std::endl;
      
      // Refineable
      if (mesh_geom_obj_pt->sample_point_container_version()==
          UseRefineableBinArray)
       { 
        
        RefineableBinArray* sample_point_container_pt=
         dynamic_cast<RefineableBinArray*>(mesh_geom_obj_pt->
                                           sample_point_container_pt());

        oomph_info 
         << "First/last sample point: "
         << sample_point_container_pt->
         first_sample_point_to_actually_lookup_during_locate_zeta()  << " " 
         << sample_point_container_pt->
         last_sample_point_to_actually_lookup_during_locate_zeta() 
         << std::endl;
        
       }
      // Non-refineable
      else if (mesh_geom_obj_pt->sample_point_container_version()==
               UseNonRefineableBinArray)
       {
        NonRefineableBinArray* sample_point_container_pt=
         dynamic_cast<NonRefineableBinArray*>(mesh_geom_obj_pt->
                                              sample_point_container_pt());
        
        oomph_info << "First/last spiral level to be visited: "
                   << sample_point_container_pt->current_min_spiral_level() << " " 
                   << sample_point_container_pt->current_max_spiral_level() << " " 
                   << std::endl;
       }
#ifdef OOMPH_HAS_CGAL
      // CGAL
      else if (mesh_geom_obj_pt->sample_point_container_version()==
               UseCGALSamplePointContainer)
       {
        CGALSamplePointContainer* sample_point_container_pt=
         dynamic_cast<CGALSamplePointContainer*>(mesh_geom_obj_pt->
                                           sample_point_container_pt());
        oomph_info 
         << "First/last sample point: "
         << sample_point_container_pt->
         first_sample_point_to_actually_lookup_during_locate_zeta()  << " " 
         << sample_point_container_pt->
         last_sample_point_to_actually_lookup_during_locate_zeta() 
         << std::endl;
       }
#endif
      else
       {
        oomph_info << "wrong sample point container\n";
        abort();
       }

      BinArray::Always_fail_elemental_locate_zeta=true;
      std::stringstream filename;
      filename << "RESLT/sample_points_" << n_dim << "d";
      if (ref_flag==0)
       {
        filename << "_nonref_bin";
       }
      else if (ref_flag==1)
       {
        filename << "_ref_bin";
       }
      else if (ref_flag==2)
       {
        filename << "cgal";
       }
      else
       {
        oomph_info << "wrong container\n";
        abort();
       }
      filename << "_spiral_level"+StringConversion::to_string(i_level)+".dat";
      BinArray::Visited_sample_points_file.open(filename.str().c_str());      
      Vector<double> s_from_locate_zeta(n_dim);
      Vector<double> zeta_test(n_dim,0.5);
      double start_t = TimingHelpers::timer();
      mesh_geom_obj_pt->locate_zeta(zeta,geom_obj_pt,s_from_locate_zeta);
      double end_t = TimingHelpers::timer();
      t_locate_zeta+=end_t-start_t;
      
      
      // Reset and close
      BinArray::Always_fail_elemental_locate_zeta=false;
      BinArray::Visited_sample_points_file.close();
      


      // Bail
      if (geom_obj_pt!=0) 
       {
        oomph_info << "Found zeta in " << i_level << " search levels.\n";
        break;
       }
      
      // Bump

      // Refineable
      if (mesh_geom_obj_pt->sample_point_container_version()==
          UseRefineableBinArray)
       {        
        
        RefineableBinArray* sample_point_container_pt=
         dynamic_cast<RefineableBinArray*>(mesh_geom_obj_pt->
                                           sample_point_container_pt());
        sample_point_container_pt->
         first_sample_point_to_actually_lookup_during_locate_zeta() =
         sample_point_container_pt->
         last_sample_point_to_actually_lookup_during_locate_zeta();
        sample_point_container_pt->
         last_sample_point_to_actually_lookup_during_locate_zeta() *=
         sample_point_container_pt->
         multiplier_for_max_sample_point_to_actually_lookup_during_locate_zeta();
       }
      // Non refineable
      else if (mesh_geom_obj_pt->sample_point_container_version()==
               UseNonRefineableBinArray)
       {
        NonRefineableBinArray* sample_point_container_pt=
         dynamic_cast<NonRefineableBinArray*>(mesh_geom_obj_pt->
                                              sample_point_container_pt());
        
        sample_point_container_pt->current_min_spiral_level()+=sample_point_container_pt->n_spiral_chunk();
        sample_point_container_pt->current_max_spiral_level()+=sample_point_container_pt->n_spiral_chunk();
       }
#ifdef OOMPH_HAS_CGAL
      // CGAL 
      else if (mesh_geom_obj_pt->sample_point_container_version()==
               UseCGALSamplePointContainer)
       {        
        
        CGALSamplePointContainer* sample_point_container_pt=
         dynamic_cast<CGALSamplePointContainer*>(mesh_geom_obj_pt->
                                           sample_point_container_pt());
        sample_point_container_pt->
         first_sample_point_to_actually_lookup_during_locate_zeta() =
         sample_point_container_pt->
         last_sample_point_to_actually_lookup_during_locate_zeta();
        sample_point_container_pt->
         last_sample_point_to_actually_lookup_during_locate_zeta() *=
         sample_point_container_pt->
         multiplier_for_max_sample_point_to_actually_lookup_during_locate_zeta();
       }
#endif
      else
       {
        oomph_info << "wrong sample point container\n";
        abort();
       }


      // Stop it because we're out of the bin array?

      // Refineable
      if (mesh_geom_obj_pt->sample_point_container_version()==
          UseRefineableBinArray)
       {
        RefineableBinArray* sample_point_container_pt=
         dynamic_cast<RefineableBinArray*>(mesh_geom_obj_pt->
                                           sample_point_container_pt());
        if (sample_point_container_pt->
            first_sample_point_to_actually_lookup_during_locate_zeta()
            <= nsample_points_in_bin_array)
         {
          has_not_reached_max_level_of_search=true;
         }
        else
         {
          has_not_reached_max_level_of_search=false;
         }
       }
      // Non-refineable
      else if (mesh_geom_obj_pt->sample_point_container_version()==
               UseNonRefineableBinArray)
       {
        NonRefineableBinArray* sample_point_container_pt=
         dynamic_cast<NonRefineableBinArray*>(mesh_geom_obj_pt->
                                              sample_point_container_pt());
        
        if (sample_point_container_pt->current_max_spiral_level() < n_max_level)
         {
          has_not_reached_max_level_of_search=true;
         }
        else
         {
          has_not_reached_max_level_of_search=false;
         }
        
       }
#ifdef OOMPH_HAS_CGAL
      // CGAL
      else if (mesh_geom_obj_pt->sample_point_container_version()==
               UseCGALSamplePointContainer)
       {
        CGALSamplePointContainer* sample_point_container_pt=
         dynamic_cast<CGALSamplePointContainer*>(mesh_geom_obj_pt->
                                           sample_point_container_pt());
        if (sample_point_container_pt->
            first_sample_point_to_actually_lookup_during_locate_zeta()
            <= nsample_points_in_bin_array)
         {
          has_not_reached_max_level_of_search=true;
         }
        else
         {
          has_not_reached_max_level_of_search=false;
         }
       }
#endif

      // Bump
      i_level++;
     }
    
    
    // Test
    oomph_info << "\n\n\nIt took " << t_locate_zeta << " secs to ";
    if (geom_obj_pt!=0)
     {
      oomph_info << "find: " 
                 << zeta[0] << " "  
                 << zeta[1] << " "  
                 << std::endl;
     }
    else
     {
      oomph_info << "not find: " 
                 << zeta[0] << " "  
                 << zeta[1] << " "  
                 << std::endl;
     }

    unsigned nvisit=mesh_geom_obj_pt->sample_point_container_pt()->
     total_number_of_sample_points_visited_during_locate_zeta_from_top_level();
    unsigned ntotal=mesh_geom_obj_pt->sample_point_container_pt()->
     total_number_of_sample_points_computed_recursively();
    
    oomph_info << "Visited " << nvisit << " sample points " 
               << " out of a total of " << ntotal << std::endl;
    
    if (nvisit!=ntotal)
     {
      oomph_info << "ERROR\n";
      abort();
     }
   
  
   }
 
 
   // Clean up 
   delete mesh_geom_obj_pt;
   mesh_geom_obj_pt=0;
 
   // If we got to here without aborting; we're done
   
    if (ref_flag==0)
     {
      std::string filename="RESLT/success_non_ref_bin.dat";
      ofstream outfile;
      outfile.open(filename.c_str());
      outfile << "1" << std::endl;
      outfile.close();
     }
    else if (ref_flag==1)
     {
      std::string filename="RESLT/success_ref_bin.dat";
      ofstream outfile;
      outfile.open(filename.c_str());
      outfile << "1" << std::endl;
      outfile.close();
     }
    else if (ref_flag==2)
     {
      std::string filename="RESLT/success_cgal.dat";
      ofstream outfile;
      outfile.open(filename.c_str());
      outfile << "1" << std::endl;
      outfile.close();
     }
   else 
    {
     oomph_info << "wrong container\n";
     abort();
    }

  }


#ifndef OOMPH_HAS_CGAL
 // Dummy output
 { 
  std::string filename="RESLT/success_cgal.dat";
  ofstream outfile;
  outfile.open(filename.c_str());
  outfile << "1" << std::endl;
  outfile.close();
 }
#endif

}

   
