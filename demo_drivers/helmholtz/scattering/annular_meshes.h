//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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

//Oomph-lib includes
#include "generic.h"

//The base meshes
#include "meshes/rectangular_quadmesh.h"
#include "meshes/backward_step_mesh.h"

namespace oomph
{


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//===================================================================
/// 2D annular mesh with a unit circle in the middle and a layer
/// of thickness h surrounding it.
//==================================================================
template<class ELEMENT>
class TwoDAnnularMesh : public virtual RectangularQuadMesh<ELEMENT>
{

public:


 ///Constructor
 TwoDAnnularMesh(const bool& periodic, const double& azimuthal_fraction,
                 const unsigned &ntheta, const unsigned &nr, 
                 const double& a, const double &h, TimeStepper* time_stepper_pt
                 = &Mesh::Default_TimeStepper) :
  RectangularQuadMesh<ELEMENT>(ntheta,nr,1.0,1.0,
                               periodic,time_stepper_pt)
  {
   // Wrap mesh into annular shape
   wrap_into_annular_shape(a,h,azimuthal_fraction);
  }


  private:
  
  /// Wrap mesh into annular shape
  void wrap_into_annular_shape(const double& a, const double& h,
                               const double& azimuthal_fraction)
  {
   //Create the hole
   Ellipse ellipse(a,a);
   
   //Set all the initial positions of the nodes
   Vector<double> xi(1);
   Vector<double> base(2);
   Vector<double> N(2);
   const unsigned n_node = this->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     // Pointer to node
     Node* nod_pt = this->node_pt(n);

     // Get the angle of the node -- rotate such that jump in angle
     // appears at periodic boundary. Shrink domain slightly
     // to keep angle unique
     xi[0] = (1.0-1.0e-10)*(-azimuthal_fraction*nod_pt->x(0))*
      2.0*MathematicalConstants::Pi+
      MathematicalConstants::Pi-
      (1.0-azimuthal_fraction)*2.0*MathematicalConstants::Pi;
   
     //Get the node's fraction in the radial direction
     double w = nod_pt->x(1);
     
     //Get the position on the ellipse base
     ellipse.position(xi,base);
     
     //Get the unit normal, if it were a circle , by normalising the base
     double norm = sqrt(base[0]*base[0] + base[1]*base[1]);
     N[0] = base[0]/norm; 
     N[1] = base[1]/norm;
     
     //Set circular film from the ellipse
     nod_pt->x(0) = base[0] + w*(h+a-norm)*N[0];
     nod_pt->x(1) = base[1] + w*(h+a-norm)*N[1];

     // Set boundary coordinates
     Vector<double> xi_bound(1);

     // Polar angle for boundary coordinate on boundary 0
     if(nod_pt->is_on_boundary(0))
      {       
       xi_bound[0]=atan2(nod_pt->x(1),nod_pt->x(0));
       nod_pt->set_coordinates_on_boundary(0,xi_bound);
      }


     // Radius for boundary coordinate on boundary 1
     if(nod_pt->is_on_boundary(1))
      {      
       xi_bound[0]=sqrt(pow(nod_pt->x(0),2)+pow(nod_pt->x(1),2));
       nod_pt->set_coordinates_on_boundary(1,xi_bound);
      }

     // Polar angle for boundary coordinate on boundary 2
     if(nod_pt->is_on_boundary(2))
      {
       xi_bound[0]=atan2(nod_pt->x(1),nod_pt->x(0));
       nod_pt->set_coordinates_on_boundary(2,xi_bound);
      }


     // Radius for boundary coordinate on boundary 3
     if(nod_pt->is_on_boundary(3))
      {      
       xi_bound[0]=sqrt(pow(nod_pt->x(0),2)+pow(nod_pt->x(1),2));
       nod_pt->set_coordinates_on_boundary(3,xi_bound);
      }
    }

   this->Boundary_coordinate_exists[0]=true;
   this->Boundary_coordinate_exists[1]=true;
   this->Boundary_coordinate_exists[2]=true;
   this->Boundary_coordinate_exists[3]=true;
  }
 
};


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//===================================================================
/// Refineable 2D annular mesh with a unit circle in the middle and a layer
/// of thickness h surrounding it.
//==================================================================
template<class ELEMENT>
 class RefineableTwoDAnnularMesh : public virtual TwoDAnnularMesh<ELEMENT>,
  public virtual RefineableQuadMesh<ELEMENT>
{

public:


 ///Constructor
 RefineableTwoDAnnularMesh(const bool& periodic, 
                           const double& azimuthal_fraction,
                           const unsigned &ntheta, 
                           const unsigned &nr, 
                           const double& a, 
                           const double &h, 
                           TimeStepper* time_stepper_pt
                           = &Mesh::Default_TimeStepper) :
 RectangularQuadMesh<ELEMENT>(ntheta,nr,1.0,1.0,
                              periodic,time_stepper_pt),
  TwoDAnnularMesh<ELEMENT>(periodic,azimuthal_fraction,ntheta, 
                           nr,a,h,time_stepper_pt)  
   {
    // Nodal positions etc. were created in constructor for
    // RectangularMesh<...>. Only need to setup quadtree forest
    this->setup_quadtree_forest();
    
    //Setup the periodic neighbour information of the TreeRoots
    //Cast to specific elements
    if (periodic)
     {
      Vector<TreeRoot*> left_root_pt(nr);
      Vector<TreeRoot*> right_root_pt(nr);
      for(unsigned i=0;i<nr;i++) 
       {
        left_root_pt[i] = 
         dynamic_cast<ELEMENT*>(this->element_pt(i*ntheta))->
         tree_pt()->root_pt();
        
        right_root_pt[i] = 
         dynamic_cast<ELEMENT*>(this->element_pt((i+1)*ntheta-1))->
         tree_pt()->root_pt();
       }
      
      using namespace QuadTreeNames;
      //Set the neighbour and periodicity
      for(unsigned i=0;i<nr;i++) 
       {
        left_root_pt[i]->neighbour_pt(W) = right_root_pt[i];
        left_root_pt[i]->set_neighbour_periodic(W); 
        
        right_root_pt[i]->neighbour_pt(E) = left_root_pt[i];
        right_root_pt[i]->set_neighbour_periodic(E); 
       }
     }
   }
 
};
 


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//===================================================================
/// 2D annular mesh with a unit circle in the middle and a layer
/// of thickness h surrounding it. hierher
//==================================================================
template<class ELEMENT>
class TwoDAnnularMesh2 : public RectangularQuadMesh<ELEMENT>,
                         public RefineableQuadMesh<ELEMENT>
{
 
public:
 
 ///Constructor
 TwoDAnnularMesh2(const unsigned &ntheta, const unsigned &nr, 
                  const double& a, const double &h, 
                  const double& h_nose, 
                  const double& azimuthal_fraction_of_nose,
                  const unsigned ntheta_nose, 
                  const unsigned nr_nose, 
                  TimeStepper* time_stepper_pt
                  = &Mesh::Default_TimeStepper) :
  RectangularQuadMesh<ELEMENT>(ntheta,nr,1.0,1.0,
                               true,time_stepper_pt)
  {
   //Create the hole
   Ellipse ellipse(a-h_nose,a-h_nose);
   
   // Radial fraction of changeover between nose and remaining layer
   double w_changeover=double(nr_nose)/double(nr);
   double w_theta_changeover=double(ntheta_nose)/double(ntheta);

   //Set all the initial positions of the nodes
   Vector<double> xi(1);
   Vector<double> base(2);
   Vector<double> N(2);
   const unsigned n_node = this->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     // Pointer to node
     Node* nod_pt = this->node_pt(n);
     
     // Redistribute nodes in azimuthal direction so that
     // ntheta_nose fit into nose
     double scale;
     if (nod_pt->x(0)<=w_theta_changeover)
      {
       scale=azimuthal_fraction_of_nose*nod_pt->x(0)/w_theta_changeover;
      }
     else
      {
       scale=azimuthal_fraction_of_nose+
        (1.0-azimuthal_fraction_of_nose)*
        (nod_pt->x(0)-w_theta_changeover)/(1.0-w_theta_changeover);
      }
     nod_pt->x(0)=scale;
     
     //Get the node's fraction in the radial direction
     double w = nod_pt->x(1);
     
     // Rescale
     double radial_push;
     if (w<=double(nr_nose)/(double(nr)))
      {
       radial_push=h_nose*w/w_changeover;
      }
     else
      {
       radial_push=h_nose+h*(w-w_changeover)/(1.0-w_changeover);
      }
          
     // Get the angle of the node -- rotate such that jump in angle
     // appears at periodic boundary. Shrink domain slightly
     // to keep angle unique
     xi[0]=(1.0-1.0e-10)*(-nod_pt->x(0))*
      2.0*MathematicalConstants::Pi+
      MathematicalConstants::Pi;
                          
     
     //Get the position on the ellipse base
     ellipse.position(xi,base);
     
     //Get the unit normal, if it were a circle , by normalising the base
     double norm = sqrt(base[0]*base[0] + base[1]*base[1]);
     N[0] = base[0]/norm; 
     N[1] = base[1]/norm;
     
     //Set circular film from the ellipse
     nod_pt->x(0) = base[0] + radial_push*N[0];
     nod_pt->x(1) = base[1] + radial_push*N[1];
    }
   

   //Setup the periodic neighbour information of the TreeRoots
   //Cast to specific elements
   Vector<ELEMENT*> left_el_pt(nr);
   Vector<ELEMENT*> right_el_pt(nr);
   for(unsigned i=0;i<nr;i++) 
    {
     left_el_pt[i] = dynamic_cast<ELEMENT*>(this->element_pt(i*ntheta));
     right_el_pt[i] =dynamic_cast<ELEMENT*>(this->element_pt((i+1)*ntheta-1));
    }
   

   // By default nobody's claiming any nodes
   std::map<Node*,bool> keep;
     
   // Get elements outside lower right block
   Vector<FiniteElement*> el_retained_pt;
   Vector<FiniteElement*> el_killed_pt;
   for (unsigned i=0;i<ntheta;i++)
    {     
     for (unsigned j=0;j<nr;j++)
      {
       FiniteElement* el_pt=this->finite_element_pt(i+ntheta*j);
       if ((i>(ntheta_nose-1))&&(j<nr_nose))
        {
         el_killed_pt.push_back(el_pt);
        }
       else
        {
         el_retained_pt.push_back(el_pt);
         unsigned nnod_el=el_pt->nnode();
         for (unsigned jj=0;jj<nnod_el;jj++)
          {
           keep[el_pt->node_pt(jj)]=true;
          }
        }
      }
    }
     
         
   // By default nobody's claiming the nodes; also store old
   // boundary ids
   std::map<Node*,std::set<unsigned> > boundaries; 
   unsigned nnod=this->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     std::set<unsigned>* aux_pt=0;
     this->node_pt(j)->get_boundaries_pt(aux_pt);
     if (aux_pt!=0)
      {
       boundaries[this->node_pt(j)]=(*aux_pt);
      }
    }
     
   // Remove information about boundary nodes
   this->remove_boundary_nodes();
     
   // Reset number of boundaries
   this->set_nboundary(7);
     
   // Kill superfluous nodes
   Vector<Node*> node_backup_pt(this->Node_pt);
   this->Node_pt.clear();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt=node_backup_pt[j];
     if (keep[nod_pt])
      {
       this->Node_pt.push_back(nod_pt);

       // Polar angle for boundary coordinate
       Vector<double> xi(1);
       xi[0]=atan2(nod_pt->x(1),nod_pt->x(0));

       // Check old boundaries of node
       std::set<unsigned> aux=boundaries[nod_pt];
       for (std::set<unsigned>::iterator it=boundaries[nod_pt].begin();
            it!=boundaries[nod_pt].end();it++)
        {
         unsigned b=(*it);
           
         // Node on old boundary 0 stays there
         if (b==0)
          {
           this->add_boundary_node(b,nod_pt);
           nod_pt->set_coordinates_on_boundary(b,xi);
          }
           
         // Node on old boundary 1 moves to boundary 5
         if (b==1)
          {
           this->add_boundary_node(5,nod_pt);
          }

         // Nodes on boundary 2 go to boundary 4 (outer boundary)
         if (b==2)
          {
           this->add_boundary_node(4,nod_pt);
           nod_pt->set_coordinates_on_boundary(4,xi);
          }

         // Some (but not all) nodes on old boundary 3 stay there
         // we'll deal with them separately below but we also
         // add all of them to boundary 6
         if (b==3)
          {
           this->add_boundary_node(6,nod_pt);
          }
           
        }
      }
     else
      {
       delete node_backup_pt[j];
      }
    }

 
   // Add nodes to new boundary 1
   unsigned i=ntheta_nose-1;
   for (unsigned j=0;j<nr_nose;j++)
    {
     FiniteElement* el_pt=this->finite_element_pt(i+ntheta*j);     
     unsigned nnod_1d=el_pt->nnode_1d();
     for (unsigned jj=0;jj<nnod_1d;jj++)
      {
       unsigned jnod=(nnod_1d-1)+jj*nnod_1d;
       Node* nod_pt=el_pt->node_pt(jnod);
       if (!(nod_pt->is_on_boundary()))
        {
         this->convert_to_boundary_node(nod_pt);
        }
       this->add_boundary_node(1,nod_pt);
       
       // Use radius as boundary coordinate
       Vector<double> xi(1);
       xi[0]=sqrt(pow(nod_pt->x(0),2)+pow(nod_pt->x(1),2));
       nod_pt->set_coordinates_on_boundary(1,xi);
      }
    }
   
   
   // Add nodes to new boundary 2
   unsigned j=nr_nose;
   for (unsigned i=ntheta_nose;i<ntheta;i++)
    {
     FiniteElement* el_pt=this->finite_element_pt(i+ntheta*j);     
     unsigned nnod_1d=el_pt->nnode_1d();
     for (unsigned jj=0;jj<nnod_1d;jj++)
      {
       unsigned jnod=jj;
       Node* nod_pt=el_pt->node_pt(jnod);
       if (!(nod_pt->is_on_boundary()))
        {         
         this->convert_to_boundary_node(nod_pt);
        }
       this->add_boundary_node(2,nod_pt);
       
       // Polar angle for boundary coordinate
       Vector<double> xi(1);
       xi[0]=atan2(nod_pt->x(1),nod_pt->x(0));
       nod_pt->set_coordinates_on_boundary(2,xi);
      }
    }
   
   // Add nodes to new boundary 3
   i=0;
   for (unsigned j=0;j<nr_nose;j++)
    {
     FiniteElement* el_pt=this->finite_element_pt(i+ntheta*j);     
     unsigned nnod_1d=el_pt->nnode_1d();
     for (unsigned jj=0;jj<nnod_1d;jj++)
      {
       unsigned jnod=jj*nnod_1d;
       Node* nod_pt=el_pt->node_pt(jnod);
       if (!(nod_pt->is_on_boundary()))
        {         
         this->convert_to_boundary_node(nod_pt);
        }
       this->add_boundary_node(3,nod_pt);
       if (nod_pt->is_on_boundary(6))
        {
         // Keep last one on two boundaries
         if ((j==(nr_nose-1))&&(jj==(nnod_1d-1)))
          {
           cout << "leaving node " 
                << nod_pt->x(0) << " " 
                << nod_pt->x(1) << "\n"; 
          }
         else
          {
           cout << "removing node " 
                << nod_pt->x(0) << " " 
                << nod_pt->x(1) << "\n"; 
           this->remove_boundary_node(6,nod_pt);
          }
        }

       // Use radius as boundary coordinate
       Vector<double> xi(1);
       xi[0]=sqrt(pow(nod_pt->x(0),2)+pow(nod_pt->x(1),2));
       nod_pt->set_coordinates_on_boundary(3,xi);
      }
    }
   
  
   // Kill redundant elements   
   this->Element_pt.clear();
   unsigned n_retained=el_retained_pt.size();
   for (unsigned e=0;e<n_retained;e++)
    {
     this->Element_pt.push_back(el_retained_pt[e]);
    }
   unsigned n_killed=el_killed_pt.size();
   for (unsigned e=0;e<n_killed;e++)
    {
     delete el_killed_pt[e];
    }
      
   // Setup quadtree forest again
   this->setup_quadtree_forest();
   
   using namespace QuadTreeNames;
   
   ofstream left;
   ofstream right;
   left.open("left.dat");
   right.open("right.dat");

   cout << "nr_nose nr " << nr_nose << " " << nr << std::endl;

   //Set the neighbour and periodicity
   for(unsigned i=nr_nose;i<nr;i++) 
    {
     left_el_pt[i]->output(left);

     left_el_pt[i]->tree_pt()->root_pt()->neighbour_pt(W) = 
      right_el_pt[i]->tree_pt()->root_pt();
     left_el_pt[i]->tree_pt()->root_pt()->set_neighbour_periodic(W); 
     
     right_el_pt[i]->output(right);
     right_el_pt[i]->tree_pt()->root_pt()->neighbour_pt(E) = 
      left_el_pt[i]->tree_pt()->root_pt();
     right_el_pt[i]->tree_pt()->root_pt()->set_neighbour_periodic(E); 
    }
   left.close();
   right.close();
   
   // Re-setup lookup scheme that establishes which elements are located
   // on the mesh boundaries
   this->setup_boundary_element_info();
   
   
   // Set flag to ensure boundary coordinates get interpolated
   // during mesh refinement
   this->Boundary_coordinate_exists[0]=true;
   this->Boundary_coordinate_exists[1]=true;
   this->Boundary_coordinate_exists[2]=true;
   this->Boundary_coordinate_exists[3]=true;
   this->Boundary_coordinate_exists[4]=true;
   
  }
 
 

 
};



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//===================================================================
/// "Wrapped around" backward facing step mesh, creating an annular 
/// mesh with a "nose" sticking into the inner hold.
//==================================================================
template<class ELEMENT>
class WrappedAroundBackwardStepQuadMesh : 
 public BackwardStepQuadMesh<ELEMENT>,
 public RefineableQuadMesh<ELEMENT>
{
 
public:
 
 
 ///Constructor
 WrappedAroundBackwardStepQuadMesh(
  const bool& periodic, const double& azimuthal_fraction,
  const unsigned &ntheta, const unsigned &nr, 
  const double& a, const double &h, TimeStepper* time_stepper_pt
  = &Mesh::Default_TimeStepper) :
  RectangularQuadMesh<ELEMENT>(ntheta,nr,8.0*atan(1.0),1.0,
                               periodic,time_stepper_pt),
  BackwardStepQuadMesh<ELEMENT>(ntheta,nr,1,1, // last two coutout
                                8.0*atan(1.0),1.0,time_stepper_pt)
  {
   // Nodal positions etc. were created in constructor for
   // RectangularMesh<...>. Only need to setup quadtree forest
   this->setup_quadtree_forest();

   //Setup the periodic neighbour information of the TreeRoots
   //Cast to specific elements
   Vector<TreeRoot*> left_root_pt(nr);
   Vector<TreeRoot*> right_root_pt(nr);
   for(unsigned i=0;i<nr;i++) 
    {
     left_root_pt[i] = 
      dynamic_cast<ELEMENT*>(this->element_pt(i*ntheta))->tree_pt()->root_pt();

     right_root_pt[i] = 
      dynamic_cast<ELEMENT*>(this->element_pt((i+1)*ntheta-1))->
      tree_pt()->root_pt();
    }
   
   using namespace QuadTreeNames;
   
   //Set the neighbour and periodicity
   for(unsigned i=0;i<nr;i++) 
    {
     left_root_pt[i]->neighbour_pt(W) = right_root_pt[i];
     left_root_pt[i]->set_neighbour_periodic(W); 
   
     right_root_pt[i]->neighbour_pt(E) = left_root_pt[i];
     right_root_pt[i]->set_neighbour_periodic(E); 
    }
   
   
   //Create the hole
   Ellipse ellipse(a,a);
   
   //Set all the initial positions of the nodes
   Vector<double> xi(1);
   Vector<double> base(2);
   Vector<double> N(2);
   const unsigned n_node = this->nnode();
   for(unsigned n=0;n<n_node;n++)
    {
     // Pointer to node
     Node* nod_pt = this->node_pt(n);

     //Get the angle of the node
     xi[0] = azimuthal_fraction*(8.0*atan(1.0) - nod_pt->x(0));

     //Get the node's fraction in the radial direction
     double w = nod_pt->x(1);

     //Get the position on the ellipse base
     ellipse.position(xi,base);

     //Get the unit normal, if it were a circle , by normalising the base
     double norm = sqrt(base[0]*base[0] + base[1]*base[1]);
     N[0] = base[0]/norm; 
     N[1] = base[1]/norm;
     
     //Set circular film from the ellipse
     nod_pt->x(0) = base[0] + w*(h+a-norm)*N[0];
     nod_pt->x(1) = base[1] + w*(h+a-norm)*N[1];

     //Set the boundary coordinate
     if(nod_pt->is_on_boundary(0))
      {
       nod_pt->set_coordinates_on_boundary(0,xi);
      }
     if(nod_pt->is_on_boundary(2))
      {
       nod_pt->set_coordinates_on_boundary(2,xi);
      }
    }
   this->Boundary_coordinate_exists[0]=true;
   this->Boundary_coordinate_exists[2]=true;
  }
 
};

}
