
  // hierher
  namespace VolumeConstraintNewtonSolver
  {
   
   /// Centre of gravity
   Vector<double> R_c;
   
   /// Original nodal positions
   Vector<Vector<double> > X_orig;
   
   /// Nodes
   Vector<Node*> Node_pt;
   
   /// Volume constraint mesh
   Mesh* Volume_constraint_mesh_pt=0;
   
   /// Target volume
   double Target_volume=0;
   
   /// Residual function for volume constraint
   void volume_constraint_residual(const Vector<double>& params,
                                   const Vector<double>& unknowns,
                                   Vector<double>& residuals)
   {
    // Push all nodes out to their new positions
    double h=unknowns[0];
    unsigned nnod=Node_pt.size();
    for (unsigned j=0;j<nnod;j++)
     {
      
      // Get unit vector in displacement direction
      Vector<double> unit_displacement(2,0.0);
      unit_displacement[0]=X_orig[j][0]-R_c[0];
      unit_displacement[1]=X_orig[j][1]-R_c[1];
      double length=sqrt(unit_displacement[0]*unit_displacement[0]+
                         unit_displacement[1]*unit_displacement[1]);
      unit_displacement[0]/=length;
      unit_displacement[1]/=length;
      
      // Shift it
      Node_pt[j]->x(0)=X_orig[j][0]+h*unit_displacement[0];
      Node_pt[j]->x(1)=X_orig[j][1]+h*unit_displacement[1];
      
     }
   
    // Get total volume enclosed by face elements 
    double vol=0.0;
    unsigned n_element=Volume_constraint_mesh_pt->nelement();
    for(unsigned e=0;e<n_element;e++) 
     {
      TemplateFreeVolumeConstraintBoundingElementBase* el_pt=
       dynamic_cast<TemplateFreeVolumeConstraintBoundingElementBase*>(
        Volume_constraint_mesh_pt->element_pt(e));
      
      if (el_pt!=0)
       {
        vol+=el_pt->contribution_to_enclosed_volume();
       }
      
      // One and only residual
      residuals[0]=vol-Target_volume;
      
     }
   }
   
   
   /// Update nodal positions at interface to achieve 
   /// desired volume
   template<class BULK_ELEMENT>
   void update_nodal_positions(Mesh* bulk_mesh_pt,
                               Vector<unsigned> bounding_boundary_id,
                               const double& target_volume)
   {


    // hierher dummy
    VolumeConstraintElement* Vol_constraint_el_pt= 
     new VolumeConstraintElement(&Problem_Parameter::Volume);

    // Clean up
    if (Volume_constraint_mesh_pt!=0)
     {
      delete Volume_constraint_mesh_pt;
     }
    
    // Make new mesh
    Volume_constraint_mesh_pt=new Mesh;

    //Loop over the bounding boundaries
    unsigned nb=bounding_boundary_id.size();
    for(unsigned i=0;i<nb;i++)
     {
      // What boundary?
      unsigned b=bounding_boundary_id[i];
      
      // How many bulk fluid elements are adjacent to boundary b?
      unsigned n_element = bulk_mesh_pt->nboundary_element(b);
      
      // Loop over the bulk fluid elements adjacent to boundary b?
      for(unsigned e=0;e<n_element;e++)
       {
        // Get pointer to the bulk fluid element that is 
        // adjacent to boundary b
        BULK_ELEMENT* bulk_elem_pt = dynamic_cast<BULK_ELEMENT*>(
         bulk_mesh_pt->boundary_element_pt(b,e));
        
        //Find the index of the face of element e along boundary b
        int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);
        
        // Create new element // hierher move into base class!
        LineVolumeConstraintBoundingSolidElement<BULK_ELEMENT>* el_pt =
         new LineVolumeConstraintBoundingSolidElement<BULK_ELEMENT>(
          bulk_elem_pt,face_index,Vol_constraint_el_pt);   
        
        // Add it to the mesh
        Volume_constraint_mesh_pt->add_element_pt(el_pt);     
       } 
     }
   
    // Store target volume
    Target_volume=target_volume;
    
    // Get centre of gravity and old nodal positions
    unsigned count=0;
    R_c.resize(2);
    R_c[0]=0.0;
    R_c[1]=0.0;
    std::map<Node*,bool> node_done;
    unsigned n_element=Volume_constraint_mesh_pt->nelement();
    for(unsigned e=0;e<n_element;e++) 
     {
      TemplateFreeVolumeConstraintBoundingElementBase* el_pt=
       dynamic_cast<TemplateFreeVolumeConstraintBoundingElementBase*>(
        Volume_constraint_mesh_pt->element_pt(e));         
      if (el_pt!=0)
       {
        unsigned nnod=el_pt->nnode();
        for (unsigned j=0;j<nnod;j++)
         {
          Node* nod_pt=el_pt->node_pt(j);
          if (!node_done[nod_pt])
           {
            R_c[0]+=nod_pt->x(0);
            R_c[1]+=nod_pt->x(1);
            count++;
            node_done[nod_pt]=true;
            
            // Store original position
            Vector<double> x_orig(2);
            x_orig[0]=nod_pt->x(0);
            x_orig[1]=nod_pt->x(1);
            X_orig.push_back(x_orig);
            Node_pt.push_back(nod_pt);
           }
         }
       }
     }
    R_c[0]/=double(count);
    R_c[1]/=double(count);


    
    BlackBoxFDNewtonSolver::Doc_Progress=true;

    // Call Newton solver
    Vector<double> params;
    Vector<double> unknowns(1,0.0);
    BlackBoxFDNewtonSolver::black_box_fd_newton_solve(
     volume_constraint_residual,params,unknowns);
    
    // Note: Final call to residual when asserting 
    // convergence of Newton method will already have
    // assigned the correct nodal positions so we just leave
    // them.
    
    // Clean up 
    Node_pt.clear();
    X_orig.clear();
    delete Volume_constraint_mesh_pt;
    Volume_constraint_mesh_pt=0;
   }
   
   
  } // End of namespace


  // Update function acting on mesh
  void mesh_update_for_volume_conservation(Mesh* mesh_pt)
  {
   ofstream some_file;
   char filename[100];

   sprintf(filename,"before.dat");
   some_file.open(filename);
   mesh_pt->output(some_file,5);   
   some_file.close();

   // Update nodal positions at interface to achieve 
   // desired volume
   Vector<unsigned> bounding_boundary_id;
   bounding_boundary_id.push_back(4);
   bounding_boundary_id.push_back(5);
   Problem_Parameter::VolumeConstraintNewtonSolver::
    update_nodal_positions<ProjectableTaylorHoodElement<MyTaylorHoodElement> >(
     mesh_pt,bounding_boundary_id,-Volume);



   sprintf(filename,"after.dat");
   some_file.open(filename);
   mesh_pt->output(some_file,5);   
   some_file.close();

  }
