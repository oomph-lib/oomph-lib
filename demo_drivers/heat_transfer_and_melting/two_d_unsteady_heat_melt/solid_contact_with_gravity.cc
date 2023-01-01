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
// Driver for 2D contact problem -- switching from displacement to
// weight control
#include <fenv.h> 

//Generic routines
#include "generic.h"

// The solid elements
#include "solid.h"

// Mesh
#include "meshes/triangle_mesh.h"
#include "meshes/rectangular_quadmesh.h"


// Contact stuff
#include "contact_elements.h"

using namespace std;

using namespace oomph;

using namespace MathematicalConstants;





/// //////////////////////////////////////////////////////////////////// 
/// ////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////// 


//======================================================================
/// Penetrator that keeps circle in contact with a control node 
/// on target surface (made of solid contact face elements) -- centre
/// of the circular penetrator is located at
///  
///     {\bf r}_c = {\bf R}_p + R {\bf e}_alpha
///
/// where {\bf R}_p is the position of the control point, R the radius
/// of the circular penetrator, and {\bf e}_alpha is a unit vector
/// inclined at an angle \alpha against the vertical. 
/// Penetration can be driven in two ways.
/// (1) We impose the vertical position of the control point (by 
///     pseudo-hijacking the Lagrange-multiplier (representing the contact 
///     pressure) stored at the controlled node. This means that rather than
///     determining the contact pressure from the no-penetration constraint,
///     (which we know to be satisfied by construction) we determine it
///     from the condition that {\bf R}_p \cdot {\bf e}_y = Y_c which 
///     is prescribed. We also impose the angle \alpha (stored as an
///     internal Data value in the element) by solving it via the
///     equation \alpha-\alpha_{prescribed} = 0.
/// (2) We impose the weight (i.e. the vertical reaction force from the 
///     contact elements) by using the equation 
///
///       \int p_c {\bf n} \cdot {\bf e}_y ds - W = 0
///   
///     as the equation for the pseudo-hijacked contact pressure and similarly,
///     use the horizontal force balance
/// 
///       \int p_c {\bf n} \cdot {\bf e}_x ds - H = 0
///   
///     to determine the rotation angle. Here, W and H are prescribed
///     and the integral is computed from the contact elements that
///     potentially impact on the penetrator.
//======================================================================
class CircularPenetratorElement : public virtual GeneralisedElement,
                                  public virtual Penetrator
 {

   public:

  /// Constructor: Pass pointer to control node whose 
  /// index_of_contact_pressure-th value represents the Lagrange multiplier
  /// (the discrete contact pressure) that has been traded for the
  /// vertical displacement/weight constraint. 
  CircularPenetratorElement(SolidNode* control_node_pt, 
                            const unsigned& index_of_contact_pressure,
                            double* r_pt) 
   {
    // Create internal data, representing the angle of rotation about
    // contact point. Determined either directly via insisting
    // that the difference between this value and the target is zero
    // or by insisting that the horizontal force is zero.
    add_internal_data(new Data(1));
    internal_data_pt(0)->set_value(0,0.0);

    // Store pointer to radius
    Radius_pt=r_pt;
    
    // Control node 
    Control_node_pt=control_node_pt;

    // Where is the tradedd contact pressure stored?
    Index_of_contact_pressure=index_of_contact_pressure;

    //Pointer to target weight (null if vertical displacement of control
    //node is imposed)
    Target_weight_pt=0;

    // Pointer to target horizontal force (null if rotation angle angle 
    // about control node is imposed)
    Target_horizontal_force_pt=0;

    // Pointer to  target vertical displacement of control node (null if 
    // weight is imposed)
    Target_yc_pt=0;

    // Pointer to target rotation angle about control node (null 
    // if horizontal force is imposed)
    Target_rotation_angle_pt=0;

    // Pointer to mesh of contact elements that contribute to force
    Contact_element_mesh_pt=0;
   }
  

  /// Virtual destructor
  virtual ~CircularPenetratorElement(){};

  /// Vector of pairs identifying values (via a pair of pointer to 
  /// Data object and index within it) that correspond to the Data values 
  /// that are determined by the horizontal/vertical/... equilibrium equations.
  Vector<std::pair<Data*,unsigned> > equilibrium_data()
   {
    // We're in 2D
    Vector<std::pair<Data*,unsigned> > thingy(2);

    
    // Horizontal equilibrium determines the rotation angle
    // which is stored as the zero-th internal data
    if (Target_horizontal_force_pt==0)
     {
      thingy[0]=std::make_pair(static_cast<Data*>(0),0);
     }
    else
     {
      thingy[0]=std::make_pair(internal_data_pt(0),0);
     }


    // Vertical equilibrium determines the discrete contact pressure
    // (Lagrange multiplier) at control node
    if (Target_weight_pt==0)
     {
      thingy[1]=std::make_pair(static_cast<Data*>(0),0);
     }
    else
     {
      thingy[1]=std::make_pair(
       external_data_pt(External_data_index_of_traded_contact_pressure),
       Index_of_contact_pressure);
     }

    return thingy;
   }
  
  /// Angle of rotation around contact point
  double angle() const
   {
    return internal_data_pt(0)->value(0);
   }


  /// Set angle of rotation around contact point
  void set_angle(const double& angle)
   {
    internal_data_pt(0)->set_value(0,angle);
   }


  /// Access to pointer to mesh of contact elements that contribute to 
  /// force on penetrator
  Mesh* contact_element_mesh_pt() const
   {
    return Contact_element_mesh_pt;
   }

  /// Set pointer to mesh of contact elements and setup
  /// external Data, i.e. Data that affects the residuals in this
  /// element. Also set the node pointed to by Control_node_pt
  /// as external Data for the elements in the contact mesh
  /// (unless they contain this node already).
  void set_contact_element_mesh_pt(Mesh* contact_element_mesh_pt)
  {
   Contact_element_mesh_pt=contact_element_mesh_pt;
   flush_external_data();
   
   // Store Data associated with control node: It contains the traded
   // Lagrange multiplier (contact pressure) 
   External_data_index_of_traded_contact_pressure=
    add_external_data(Control_node_pt);
   
   // Store its position data
   add_external_data(Control_node_pt->variable_position_pt());
   
   // Loop over all the elements in the contact mesh
   // If they don't contain the contact node already, its position
   // is external data because it affects the penetration.
   unsigned nel=Contact_element_mesh_pt->nelement();
   for (unsigned e=0;e<nel;e++)
    {
     bool el_contains_control_node=false;
     FiniteElement* el_pt=Contact_element_mesh_pt->finite_element_pt(e);
     unsigned nnod=el_pt->nnode();
     for (unsigned j=0;j<nnod;j++)
      {
       SolidNode* nod_pt=dynamic_cast<SolidNode*>(el_pt->node_pt(j));
       if (nod_pt==Control_node_pt)
        {
         el_contains_control_node=true;
        }
      }
     
     // Position of control node affects position of penetrator and
     // therefore is external data for all contact elements (apart from
     // any that contain the control node as one of their own)
     if (!el_contains_control_node)
      {
       el_pt->add_external_data(Control_node_pt->variable_position_pt());
      }
     
     // Rotation angle angle affects position of penetrator and therefore
     // affects penetration at all contact elements
     el_pt->add_external_data(internal_data_pt(0));
    }
  }


  

  /// Set target horizontal and vertical force to be in equilibrium
  void set_equilibrium_target_forces()
   {
#ifdef PARANOID
    if (Target_horizontal_force_pt==0)
     {
      std::stringstream junk;
      junk << "Target_horizontal_force_pt==0\n"
           << "Please set it by call to impose_weight(...)\n";
      throw OomphLibError(
       junk.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
     } 
    if (Target_weight_pt==0)
     {
      std::stringstream junk;
      junk << "Target_weight_pt==0.\n"
           << "Please set it by call to impose_horizontal_force(...)\n";
      throw OomphLibError(
       junk.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
     }
#endif
    Vector<double> force(resulting_force());
    (*Target_horizontal_force_pt)=-force[0];
    (*Target_weight_pt)=-force[1];
   }

  /// Target weight (returns zero if not imposed)
  double target_weight()
   {
    if (Target_weight_pt==0)
     {
      return 0.0;
     }
    else
     {
      return *Target_weight_pt;
     }
   }

  /// Target horizontal force (returns zero if not imposed)
  double target_horizontal_force()
   {
    if (Target_horizontal_force_pt==0)
     {
      return 0.0;
     }
    else
     {
      return *Target_horizontal_force_pt;
     }
   }

  /// Target vertical position of control point (returns zero if not imposed)
  double target_yc()
   {
    if (Target_yc_pt==0)
     {
      return 0.0;
     }
    else
     {
      return *Target_yc_pt;
     }
   }

  /// Target rotation angle about contact point (returns zero if not imposed)
  double target_rotation_angle()
   {
    if (Target_rotation_angle_pt==0)
     {
      return 0.0;
     }
    else
     {
      return *Target_rotation_angle_pt;
     }
   }

  /// Is vertical positon of control node imposed? If false then weight imposed.
  bool yc_is_imposed()
   {
    return (Target_yc_pt!=0);
   }

  /// Impose weight (rather than imposed displacement). Target
  /// weight specified via pointer.
  void impose_weight(double* target_weight_pt)
   {
    Target_weight_pt=target_weight_pt;
    Target_yc_pt=0;
   }

  /// Impose vertical position of control node (rather than weight).
  /// Target vertical position of control node specified via pointer.
  void impose_yc(double* target_yc_pt)
   {
    Target_weight_pt=0;
    Target_yc_pt=target_yc_pt;
   }


  /// Is angle of rotation about control node imposed? If false then 
  /// horizontal force is imposed.
  bool rotation_angle_is_imposed()
   {
    return (Target_rotation_angle_pt!=0);
   }


  /// Impose horizontal force (rather than rotation about contact node). 
  /// Target force specified via pointer.
  void impose_horizontal_force(double* target_horizontal_force_pt)
   {
    Target_horizontal_force_pt=target_horizontal_force_pt;
    Target_rotation_angle_pt=0;
   }

  /// Impose rotation about contact node (rather than horizontal force)
  /// Target angle specified via pointer.
  void impose_rotation_angle(double* target_rotation_angle_pt)
   {
    Target_horizontal_force_pt=0;
    Target_rotation_angle_pt=target_rotation_angle_pt;
   }

  /// Fill in contribution to residuals
  void fill_in_contribution_to_residuals(Vector<double> &residuals) 
  {
   // Get resulting force from all associated PseudoContactElements
   // onto the elastic body
   Vector<double> force(resulting_force());

   // Equation for Lagrange multiplier (contact pressure) at controlled
   // node
   int local_eqn=external_local_eqn(
    External_data_index_of_traded_contact_pressure,Index_of_contact_pressure); 
   if (local_eqn>=0)
    {
     // Resulting force from all associated PseudoContactElements
     // onto the elastic body is equal and opposite to force on penetrator
     if (Target_weight_pt!=0)
      {
       residuals[local_eqn]+=(*Target_weight_pt);
      }
     // Impose vertical position of control node
     else
      {
#ifdef PARANOID
       if (Target_yc_pt!=0)
        {
#endif
         residuals[local_eqn]+=Control_node_pt->x(1)-(*Target_yc_pt);
#ifdef PARANOID
        }
       else
        {
         std::stringstream junk;
         junk << "Target_yc_pt=0\n"
              << "Set with impose_yc(...)\n";
         throw OomphLibError(
          junk.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
        }
#endif
      }
    }



   // Equation for rotation angle
   local_eqn=internal_local_eqn(0,0); 
   if (local_eqn>=0)
    {
     // Resulting force from all associated PseudoContactElements
     // onto the elastic body is equal and opposite to force on penetrator
     if (Target_horizontal_force_pt!=0)
      {
       residuals[local_eqn]+=(*Target_horizontal_force_pt);
      }
     // Set rotation angle 
     else
      {
#ifdef PARANOID
       if (Target_rotation_angle_pt!=0)
        {
#endif
         residuals[local_eqn]+=internal_data_pt(0)->value(0)-
          (*Target_rotation_angle_pt);
#ifdef PARANOID
        }
       else
        {
         std::stringstream junk;
         junk << "Target_rotation_angle_pt=0\n"
              << "Set with impose_rotation_angle(...)\n";
         throw OomphLibError(
          junk.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
        }
#endif
      }
    }
  }
    
  /// Get centre of penetrator
  Vector<double> centre() const
   {
    Vector<double> rc(2);
    rc[0]=centre(0); 
    rc[1]=centre(1); 
    return rc;
   }

  /// Get centre of penetrator 
  double centre(const unsigned& i) const
   {
    switch (i)
     {
     case 0:
      return Control_node_pt->x(0)+(*Radius_pt)*sin(angle());
      break;
      
     case 1:
      return Control_node_pt->x(1)+(*Radius_pt)*cos(angle());
      break;
      
     default:
      std::stringstream junk;
      junk << "Wrong index: " << i 
           << "\nCan only handle 0 or 1 (it's 2D!)\n";
      throw OomphLibError(
       junk.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
     }
   }
 
  /// Get penetration for given point x.
  void penetration(const Vector<double>& x, const Vector<double>& n, 
                   double& d, bool& intersection) const
  {
   // Vector from potential contact point to centre of penetrator
   Vector<double> l(2);
   l[0]=centre(0)-x[0];
   l[1]=centre(1)-x[1];
   
   // Distance from potential contact point to centre of penetrator
   double ll=sqrt(l[0]*l[0]+l[1]*l[1]);

   // Projection of vector from potential contact point to centre of penetrator
   // onto outer unit normal on potential contact point
   double project=n[0]*l[0]+n[1]*l[1];
   double project_squared=project*project;

   // Final term in square root
   double b_squared=ll*ll-(*Radius_pt)*(*Radius_pt);

   // Is square root negative? In this case we have no intersection
   // and we return penetration as -DBL_MAX
   if (project_squared<b_squared)
    {
     d= -DBL_MAX;
     intersection = false;
    }
   else
    {
     double sqr=sqrt(project_squared-b_squared);
     d= -std::min(project-sqr,project+sqr);
     intersection = true;
    }
  }


  /// Output coordinates of penetrator at nplot plot points
  void output(std::ostream &outfile, const unsigned& nplot) const
   {
    for (unsigned j=0;j<nplot;j++)
     {
      double phi=2.0*MathematicalConstants::Pi*double(j)/double(nplot-1);
      outfile << centre(0)+(*Radius_pt)*cos(phi) << " " 
              << centre(1)+(*Radius_pt)*sin(phi)
              << std::endl;
     }
   }
    
  /// Resulting force from all associated ContactElements
  Vector<double> resulting_force() const
  {
   Vector<double> contact_force(2,0.0);
   Vector<double> el_contact_force(2);
   unsigned nel=Contact_element_mesh_pt->nelement();
   for (unsigned e=0;e<nel;e++)
    {
     dynamic_cast<TemplateFreeContactElementBase*>(
      Contact_element_mesh_pt->element_pt(e))->
      resulting_contact_force(el_contact_force);
     for (unsigned i=0;i<2;i++)
      {
       contact_force[i]+=el_contact_force[i];
      }
    }
   return contact_force;
  }

  /// Radius of penetrator
  double radius() const {return *Radius_pt;}
 
 private:

  /// Pointer to radius of penetrator
  double* Radius_pt;

  /// Control node
  SolidNode* Control_node_pt;

  /// Index at which contact pressure (Lagr mult) is stored in nodal
  /// data associated with control node
  unsigned Index_of_contact_pressure;

  /// Index of external data that contains the the contact
  /// pressure in its Index_of_contact_pressure-th value
  unsigned External_data_index_of_traded_contact_pressure;
  
  /// Pointer to target weight (null if vertical displacement of control
  /// node is imposed)
  double* Target_weight_pt;

  /// Pointer to target horizontal force (null if rotation angle angle 
  /// about control node is imposed)
  double* Target_horizontal_force_pt;

  /// Pointer to  target vertical displacement of control node (null if 
  /// weight is imposed)
  double* Target_yc_pt;

  /// Pointer to target rotation angle about control node (null 
  /// if horizontal force is imposed)
  double* Target_rotation_angle_pt;

  /// Mesh of contact elements that contribute to weight/horizontal force
  Mesh* Contact_element_mesh_pt;

 };


/// //////////////////////////////////////////////////////////////////// 
/// ////////////////////////////////////////////////////////////////////
/// //////////////////////////////////////////////////////////////////// 

//======start_of_ProblemParameters=====================
/// Namespace for problem parameters
//=====================================================
namespace ProblemParameters
{
 
 /// Impose position of centre (i.e. a stand-alone penetrator with
 /// prescribed position) or indirectly via control node?
 bool Impose_position_of_centre=true;
 
 /// Non-dim density for solid
 double Lambda_sq=0.0;
 
 /// Poisson's ratio for solid
 double Nu=0.3;

 /// Pointer to constitutive law
 ConstitutiveLaw* Constitutive_law_pt=0;

 /// Radius of penetrator
 double Radius=0.5;

 /// Radius of elastic body
 double Radius_of_elastic_body=2.0;

 /// Penetrator
 Penetrator* Penetrator_pt=0;

 /// NOTE: WE IMPOSE EITHER THESE ...

 /// Weight of penetrator
 double Weight=0.0;

 /// Horizontal force of penetrator
 double Horizontal_force=0.0;

 /// ... OR THESE...

 /// Target vertical position of control node
 double Y_c=0.0;

 /// Target rotation angle about control node
 double Rotation_angle=0.0;

 /// ...OR THIS
 
 /// Position of centre of penetrator
 Vector<double> Centre;

 /// Initial/max element area
 double El_area=0.02;

 /// Factor for element length on contact boundary
 double Element_length_factor=1.0;


 /// Body force magnitude
 double Body_force_amplitude=0.0;
 
 // Sharpness of body force
 double Body_force_alpha=1.0e4;

 /// The body force function
 void body_force(const double &time,
                 const Vector<double> &x,
                 Vector<double> &result)
 {
 result[0] = 0.0;
 result[1] = -Body_force_amplitude*(1.0-x[0])*x[0]*
  exp(-Body_force_alpha*(x[0]-0.5)*(x[0]-0.5));
 }
} // end of ProblemParameters


/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////


//=====start_of_problem_class=========================================
/// Problem class
//====================================================================
template<class ELEMENT>
class ContactProblem : public Problem
{

public:

 /// Constructor
 ContactProblem();
 
 /// Destructor (empty)
 ~ContactProblem(){}
 
 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}
 
 /// Update the problem specs before solve (empty)
 void actions_before_newton_solve(){}
  

 /// Actions before adapt: wipe contact elements
 void actions_before_adapt() 
  {
   // Backup x coordinate of old control node
   Xc_old=Control_node_pt->x(0);

#ifdef STRUCTURED_MESH

   // Make backup of surface mesh
   Backed_up_surface_contact_mesh_pt=
    new BackupMeshForProjection<QElement<1,3> >(Surface_contact_mesh_pt,
                                                Contact_boundary_id);
#else

   Bulk_mesh_pt->boundary_polyline_pt(Contact_boundary_id)->set_maximum_length(
    Maximum_element_length_on_contact_boundary*
    ProblemParameters::Element_length_factor);

   // Make backup of surface mesh
   Backed_up_surface_contact_mesh_pt=
    new BackupMeshForProjection<TElement<1,3> >(Surface_contact_mesh_pt,
                                                Contact_boundary_id);
#endif

   // // Output contact elements
   // ofstream some_file;
   // char filename[100];
   // sprintf(filename,"contact_before.dat");
   // some_file.open(filename);
   // unsigned nel=Surface_contact_mesh_pt->nelement();
   // for (unsigned e=0;e<nel;e++)
   //  {
   //   dynamic_cast<NonlinearSurfaceContactElement<ELEMENT>* >(
   //    Surface_contact_mesh_pt->element_pt(e))->output(some_file);
   //  }
   // some_file.close();
   
   // // Kill the  elements and wipe surface mesh
   delete_contact_elements();
   
   // Wipe the mesh
   Penetrator_mesh_pt->flush_element_and_node_storage();

   // Rebuild the Problem's global mesh from its various sub-meshes
   rebuild_global_mesh();
  }
 
 /// Actions after adapt: 
 /// Setup the problem again -- remember that the mesh has been
 /// completely rebuilt and its element's don't have any
 /// pointers to source fcts etc. yet
 void actions_after_adapt()
  {
   // Create contact elements
   create_contact_elements();
   
   // Now project from backup of original contact mesh to new one
   Backed_up_surface_contact_mesh_pt->project_onto_new_mesh(
    Surface_contact_mesh_pt);

   // Rebuild elements
   complete_problem_setup();
   
   // Rebuild the Problem's global mesh from its various sub-meshes
   rebuild_global_mesh();
   
   // Kill backed up mesh
   delete Backed_up_surface_contact_mesh_pt;
   Backed_up_surface_contact_mesh_pt=0;

   // // Output contact elements
   // ofstream some_file;
   // char filename[100];
   // sprintf(filename,"contact_after.dat");
   // some_file.open(filename);
   // unsigned nel=Surface_contact_mesh_pt->nelement();
   // for (unsigned e=0;e<nel;e++)
   //  {
   //   dynamic_cast<NonlinearSurfaceContactElement<ELEMENT>* >(
   //    Surface_contact_mesh_pt->element_pt(e))->output(some_file);
   //  }
   // some_file.close();

  }

 /// Switch to displ control
 void switch_to_displ_control()
  {
   dynamic_cast<CircularPenetratorElement*>(
    ProblemParameters::Penetrator_pt)->impose_yc(&ProblemParameters::Y_c);
   dynamic_cast<CircularPenetratorElement*>(
    ProblemParameters::Penetrator_pt)->impose_rotation_angle(
    &ProblemParameters::Rotation_angle);
   ProblemParameters::Y_c=Control_node_pt->x(1);
  }


 /// Switch to force control
 void switch_to_force_control()
  {
   dynamic_cast<CircularPenetratorElement*>(
    ProblemParameters::Penetrator_pt)->impose_weight(
     &ProblemParameters::Weight);
   dynamic_cast<CircularPenetratorElement*>(
    ProblemParameters::Penetrator_pt)->impose_horizontal_force(
     &ProblemParameters::Horizontal_force); 
   dynamic_cast<CircularPenetratorElement*>(
    ProblemParameters::Penetrator_pt)->set_equilibrium_target_forces();

   // Re-set contact mesh -- we now need to treat the positions
   // and penalty pressures of all nodes as external data of the
   // penetrator element!
   dynamic_cast<CircularPenetratorElement*>(
    ProblemParameters::Penetrator_pt)->set_contact_element_mesh_pt(
    Surface_contact_mesh_pt);

   // Reset penetrator because its equilibrium data (which is external
   // data for contact elements) has changed
   unsigned n_element=Surface_contact_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement 
     NonlinearSurfaceContactElement<ELEMENT> *el_pt = 
      dynamic_cast<NonlinearSurfaceContactElement<ELEMENT>*>(
       Surface_contact_mesh_pt->element_pt(e));
     
     // Set pointer to penetrator
     el_pt->set_penetrator_pt(ProblemParameters::Penetrator_pt);
    }
   
   cout <<"New number of equations: " << assign_eqn_numbers() << std::endl; 

  }

 /// Doc the solution
 void doc_solution();
 

private:

 /// Create contact elements
 void create_contact_elements()
  {
   // How many bulk elements are adjacent to boundary b?
   unsigned b=Contact_boundary_id; 
   unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
   
   // Loop over the bulk elements adjacent to boundary b?
   for(unsigned e=0;e<n_element;e++)
    {
     // Get pointer to the bulk element that is adjacent to boundary b
     ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->boundary_element_pt(b,e));
     
     //What is the face index of element e along boundary b
     int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
     
     // Build the corresponding contact element
     NonlinearSurfaceContactElement<ELEMENT>* contact_element_pt = new 
      NonlinearSurfaceContactElement<ELEMENT>(bulk_elem_pt,face_index,
                                                Contact_id);
     
     //Add the contact element to the surface mesh
     Surface_contact_mesh_pt->add_element_pt(contact_element_pt);

    } //end of loop over bulk elements adjacent to boundary b    
  }



 /// Delete contact elements
 void delete_contact_elements()
  {
   // How many surface elements are in the surface mesh
   unsigned n_element = Surface_contact_mesh_pt->nelement();
   
   // Loop over the surface elements
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill surface element
     delete Surface_contact_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Surface_contact_mesh_pt->flush_element_and_node_storage();
  }


 /// Helper function to (re-)set boundary condition
 /// and complete the build of  all elements
 void complete_problem_setup()
  {

   // Set (pseudo-)solid mechanics properties for all elements
   //---------------------------------------------------------
   unsigned n_element = Bulk_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     //Cast to a solid element
     ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
     
     // Set the constitutive law
     el_pt->constitutive_law_pt() =
      ProblemParameters::Constitutive_law_pt;
     
     // Set density to zero
     el_pt->lambda_sq_pt()=&ProblemParameters::Lambda_sq;
     
     // Disable inertia
     el_pt->disable_inertia();

     // Set body force
     el_pt->body_force_fct_pt() = &ProblemParameters::body_force;

    }

   // Apply boundary conditions for solid
   //------------------------------------

   // Bottom: completely pinned
   unsigned b=Bottom_boundary_id;
   unsigned nnod=Bulk_mesh_pt->nboundary_node(b);
   for (unsigned j=0;j<nnod;j++)
    {
     SolidNode* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,j);
     nod_pt->pin_position(0);
     nod_pt->pin_position(1);
    }

   // Sides: Symmetry bcs
   b=Left_boundary_id;
   nnod=Bulk_mesh_pt->nboundary_node(b);
   for (unsigned j=0;j<nnod;j++)
    {
     SolidNode* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,j);
     nod_pt->pin_position(0);
    }
   b=Right_boundary_id;
   nnod=Bulk_mesh_pt->nboundary_node(b);
   for (unsigned j=0;j<nnod;j++)
    {
     SolidNode* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,j);
     nod_pt->pin_position(0);
    }


   // Enforce contact at most central or most heavily loaded node
   //------------------------------------------------------------
   {
    
    // Update angle?
    bool update_angle=false;
    double phi_old=0.0;
    CircularPenetratorElement* pen_el_pt=
     dynamic_cast<CircularPenetratorElement*>(
      ProblemParameters::Penetrator_pt);
    if (pen_el_pt!=0)
     {
      update_angle=true;
      phi_old=pen_el_pt->angle();
     }
    
    double x_c=0.5;
    Control_node_pt=0;
    SolidNode* most_central_node_pt=0;
    SolidNode* most_loaded_node_pt=0;
    double dist_min=DBL_MAX;
    double load_max=0.0;
    unsigned b=Contact_boundary_id;
    unsigned nnod=Bulk_mesh_pt->nboundary_node(b);
    for (unsigned j=0;j<nnod;j++)
     {
      SolidNode* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,j);
      
      // Find closest node
      double dist=std::fabs(nod_pt->x(0)-x_c);
      if (dist<dist_min) 
       {
        dist_min=dist;
        most_central_node_pt=nod_pt;
       }

      // Find most loaded node
      BoundaryNodeBase *bnod_pt=dynamic_cast<BoundaryNodeBase*>(nod_pt);
      unsigned index_of_contact_pressure=
       bnod_pt->index_of_first_value_assigned_by_face_element(Contact_id);
      if (nod_pt->value(index_of_contact_pressure)>load_max)
       {
        load_max=nod_pt->value(index_of_contact_pressure);
        most_loaded_node_pt=nod_pt;
       }
     }

    // No load
    if (most_loaded_node_pt==0)
     {
      oomph_info << "Choosing most central node as control node\n";
      Control_node_pt=most_central_node_pt;
     }
    else
     {
      oomph_info << "Choosing most loaded node as control node\n";
      Control_node_pt=most_loaded_node_pt;
     }
    oomph_info << "Control node located at: "
               << Control_node_pt->x(0) << " " 
               << Control_node_pt->x(1) << " "
               << std::endl;


    // Update angle
    if (update_angle)
     {
      double phi_new=asin((Xc_old-Control_node_pt->x(0))/
                          ProblemParameters::Radius+sin(phi_old));
      dynamic_cast<CircularPenetratorElement*>(
       ProblemParameters::Penetrator_pt)->set_angle(phi_new);
      oomph_info << "Old/new angle: " << phi_old << " "
                 << phi_new << std::endl;
     }

    
    // Set target vertical position of control node to its current
    // position. This is fine as initial assignment; it's then over-written
    // before the next solve (if the displacement is controlled). 
    // If/when this function is called during adaptation it's fine too
    // because it doesn't impose another displacement increment -- it simply
    // maintains the position that was obtained on the previous solve.
    ProblemParameters::Y_c=Control_node_pt->x(1);
    
    // Index of nodal value at control node that stores the traded
    // contact pressure
    unsigned index_of_traded_contact_pressure=UINT_MAX;

    // Impose position of centre when penetrator position directly imposed
    //--------------------------------------------------------------------
    if (ProblemParameters::Impose_position_of_centre)
     {
      // Delete old one
      delete ProblemParameters::Penetrator_pt;
      
      // Make new one
      ProblemParameters::Penetrator_pt =
       new CircularPenetrator(&ProblemParameters::Centre,
                              ProblemParameters::Radius);
     }
    // Compute penetrator position as part of the solution, either by
    // --------------------------------------------------------------
    // prescribing nodal position or weight
    // ------------------------------------
    else
     {
      // Loop over face elements to identify the ones that contain 
      // the control node
      unsigned nel=Surface_contact_mesh_pt->nelement();
      bool found=false;
      for (unsigned e=0;e<nel;e++)
       {
        NonlinearSurfaceContactElement<ELEMENT>* el_pt=
         dynamic_cast<NonlinearSurfaceContactElement<ELEMENT>*>(
          Surface_contact_mesh_pt->element_pt(e));
        unsigned nnod=el_pt->nnode();
        for (unsigned j=0;j<nnod;j++)
         {
          SolidNode* solid_nod_pt=dynamic_cast<SolidNode*>(el_pt->node_pt(j));
          if (solid_nod_pt==Control_node_pt)
           {
            // Got it!
            found=true;
            
            // Find index at which contact pressure/Lagrange multiplier
            // is stored
            BoundaryNodeBase *bnod_pt = 
             dynamic_cast<BoundaryNodeBase*>(solid_nod_pt); 
            
            // Get the index of the first nodal value associated with
            // this FaceElement
            unsigned new_index_of_traded_contact_pressure=
             bnod_pt->index_of_first_value_assigned_by_face_element(Contact_id);
            
            // Copy across
#ifdef PARANOID
            if (index_of_traded_contact_pressure!=UINT_MAX)
             {
              if (new_index_of_traded_contact_pressure!=
                  index_of_traded_contact_pressure)
               {
                std::stringstream junk;
                junk << "Inconsistency in identification of index of traded"
                     << "contact pressure: " 
                     << new_index_of_traded_contact_pressure
                     << " != " << index_of_traded_contact_pressure;
                throw OomphLibError(
                 junk.str(),
                 OOMPH_CURRENT_FUNCTION,
                 OOMPH_EXCEPTION_LOCATION);
               }
             }
#endif
            index_of_traded_contact_pressure=
             new_index_of_traded_contact_pressure;
           }
         }
       }
      if (!found) 
       {
        std::stringstream junk;
        junk << "Control node not found!";
        throw OomphLibError(
         junk.str(),
         OOMPH_CURRENT_FUNCTION,
         OOMPH_EXCEPTION_LOCATION);
       }
      
      // Back up old penetrator (if it existed)
      CircularPenetratorElement* old_penetrator_pt=
       dynamic_cast<CircularPenetratorElement*>(
        ProblemParameters::Penetrator_pt);
      
      // Make new one
      ProblemParameters::Penetrator_pt =
       new CircularPenetratorElement(Control_node_pt,
                                     index_of_traded_contact_pressure,
                                     &ProblemParameters::Radius);
      
      // Displacement or weight imposed?
      bool impose_displ=true;
      if (old_penetrator_pt!=0)
       {
        if (!old_penetrator_pt->yc_is_imposed())
         {
          impose_displ=false;
         }
       }
      if (impose_displ)
       {
        dynamic_cast<CircularPenetratorElement*>(
         ProblemParameters::Penetrator_pt)->impose_yc(&ProblemParameters::Y_c);
       }
      else
       {
        dynamic_cast<CircularPenetratorElement*>(
         ProblemParameters::Penetrator_pt)->impose_weight(
          &ProblemParameters::Weight);
       }
      
      
      // Angle or horizontal force imposed?
      bool impose_angle=true;
      if (old_penetrator_pt!=0)
       {
        if (!old_penetrator_pt->rotation_angle_is_imposed())
         {
          impose_angle=false;
         }
       }
      if (impose_angle)
       {
        dynamic_cast<CircularPenetratorElement*>(
         ProblemParameters::Penetrator_pt)->impose_rotation_angle(
          &ProblemParameters::Rotation_angle);
       }
      else
       {
        dynamic_cast<CircularPenetratorElement*>(
         ProblemParameters::Penetrator_pt)->impose_horizontal_force(
          &ProblemParameters::Horizontal_force);
        dynamic_cast<CircularPenetratorElement*>(
         ProblemParameters::Penetrator_pt)->set_angle(
          old_penetrator_pt->angle());
       }
      
      // Now kill the old one!
      delete old_penetrator_pt;
      
      // Add to mesh
      Penetrator_mesh_pt->add_element_pt(
       dynamic_cast<CircularPenetratorElement*>(
        ProblemParameters::Penetrator_pt));
     }
    
    // Pass contact elements to penetrator element and declare their
    // positions and Lagrange multiplier (contact pressure) values
    // to be external data.
    CircularPenetratorElement* el_pt=dynamic_cast<CircularPenetratorElement*>(
     ProblemParameters::Penetrator_pt);
    if (el_pt!=0)
     {
      el_pt->set_contact_element_mesh_pt(Surface_contact_mesh_pt);
     }
    
   } // end of penetrator position computed as part of the solution
     // (directly or indirectly)


   // Loop over the contact elements to pass pointer to penetrator
   //-------------------------------------------------------------
   n_element=Surface_contact_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement 
     NonlinearSurfaceContactElement<ELEMENT> *el_pt = 
      dynamic_cast<NonlinearSurfaceContactElement<ELEMENT>*>(
       Surface_contact_mesh_pt->element_pt(e));
     
     // Set pointer to penetrator
     el_pt->set_penetrator_pt(ProblemParameters::Penetrator_pt);
    }
   
  }

#ifdef STRUCTURED_MESH

 /// Pointer to bulk mesh
 ElasticRefineableRectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;

#else

 /// Pointer to bulk mesh
 RefineableSolidTriangleMesh<ELEMENT>* Bulk_mesh_pt;

#endif
 
 /// Pointer to the "surface" mesh
 Mesh* Surface_contact_mesh_pt;
 
 /// Penetrator mesh
 Mesh* Penetrator_mesh_pt;

 /// ID of contact boundary
 unsigned Contact_boundary_id;

 /// ID of top left boundary
 unsigned Left_top_boundary_id;

 /// ID of top right boundary
 unsigned Right_top_boundary_id;

 /// ID of bottom boundary
 unsigned Bottom_boundary_id;

 /// ID of left boundary
 unsigned Left_boundary_id;

 /// ID of right boundary
 unsigned Right_boundary_id;

 /// Trace file
 ofstream Trace_file;

 // Setup labels for output
 DocInfo Doc_info;

#ifdef STRUCTURED_MESH

 /// Backup of Surface_contact_mesh_pt so the Lagrange multipliers
 /// can be projected across
 BackupMeshForProjection<QElement<1,3> >* Backed_up_surface_contact_mesh_pt;

#else

 /// Backup of Surface_contact_mesh_pt so the Lagrange multipliers
 /// can be projected across
 BackupMeshForProjection<TElement<1,3> >* Backed_up_surface_contact_mesh_pt;

#endif

 /// Pointer to control node where Lagrange multiplier (contact pressure)
 /// is "pseudo-hijacked" to impose either displacement or weight constraint.
 SolidNode* Control_node_pt;

 /// x coordinate of old control node
 double Xc_old;

 /// ID of additional nodal values created by contact elements to store
 /// contact pressure/Lagr. mult.
 unsigned Contact_id;

 /// x coordinate of lower left corner 
 double X_ll;

 /// x coordinate of upper right corner 
 double X_ur;
 
 /// y coordinate of lower left corner 
 double Y_ll;

 /// y coordinate of upper right corner 
 double Y_ur;

 /// Contact boundary in its poly line representation
 TriangleMeshPolyLine* Contact_boundary_pt;

 /// Max. element length on contact boundary
 double Maximum_element_length_on_contact_boundary;

}; // end of problem class


//========start_of_constructor============================================
/// Constructor for contact problem in square domain
//========================================================================
template<class ELEMENT>
ContactProblem<ELEMENT>::ContactProblem()
{ 

 // Allow for crap initial guess or convergence problems...
 Problem::Max_newton_iterations=50;

 // Initialise
 Control_node_pt=0;

 // Initialise x coordinate of old control node (will be overwritten 
 // when needed)
 Xc_old=0.0;

 // ID of additional nodal values created by contact elements to store
 // contact pressure/Lagr. mult.
 Contact_id=1;

 // Initialise
 Backed_up_surface_contact_mesh_pt=0;

 // Output directory
 Doc_info.set_directory("RESLT");
 
 // Output number
 Doc_info.number()=0;

 // Open trace file
 Trace_file.open("RESLT/trace.dat");
 
 // Set boundaries of domain
 X_ll=0.0;;
 X_ur=1.0;
 Y_ll=0.0;
 Y_ur=1.0;


 // Parameters for deflection of upper boundary to concave circular arc
 double radius=ProblemParameters::Radius_of_elastic_body;
 double alpha=asin(0.5*(X_ur+X_ll)/radius);
 double y_c=Y_ur;
 if (ProblemParameters::Radius_of_elastic_body<0.0)
  {
   alpha=0.0;
  }
 else
  {
   y_c+=radius*cos(alpha);
  }

#ifdef STRUCTURED_MESH

 // Create the mesh

 // # of elements in x-direction
 unsigned n_x=11;

 // # of elements in y-direction
 unsigned n_y=10;

 // Domain length in x-direction
 double l_x=1.0;

 // Domain length in y-direction
 double l_y=1.0;

 //Now create the mesh 
 Bulk_mesh_pt = new ElasticRefineableRectangularQuadMesh<ELEMENT>
  (n_x,n_y,l_x,l_y);

 // Contact boundary
 Contact_boundary_id=2;

 // ID of top left boundary -- not used
 Left_top_boundary_id=UINT_MAX;

 // ID of top right boundary -- not used
 Right_top_boundary_id=UINT_MAX;

 // ID of bottom boundary
 Bottom_boundary_id=0;

 // ID of left boundary
 Left_boundary_id=3;

 // ID of right boundary
 Right_boundary_id=1;

 // Move nodes on contact boundary on circle
  if (ProblemParameters::Radius_of_elastic_body>0.0)
   {
    unsigned nnod=Bulk_mesh_pt->nnode();
    for (unsigned j=0;j<nnod;j++)
     {
      SolidNode* nod_pt=Bulk_mesh_pt->node_pt(j);
      double x=nod_pt->x(0);
      double y_fract=nod_pt->x(1)/Y_ur;
      double y_b=y_c-sqrt(radius*radius-(x-0.5*(X_ll+X_ur))*
                          (x-0.5*(X_ll+X_ur)));
      nod_pt->x(1)=y_b*y_fract;
     }
    Bulk_mesh_pt->set_lagrangian_nodal_coordinates();
   }

#else

 // Pointer to the closed curve that defines the outer boundary
 TriangleMeshClosedCurve* closed_curve_pt=0;

 // Build outer boundary as Polygon
  
 // The boundary is bounded by five distinct boundaries, each
 // represented by its own polyline
 Vector<TriangleMeshCurveSection*> boundary_polyline_pt(6);
 
 // Vertex coordinates on boundary
 Vector<Vector<double> > bound_coords(2);
 
 // Left boundary
 bound_coords[0].resize(2);
 bound_coords[0][0]=X_ll;
 bound_coords[0][1]=Y_ur;

 bound_coords[1].resize(2);
 bound_coords[1][0]=X_ll;
 bound_coords[1][1]=Y_ll;

 // Build the boundary polyline
 Left_boundary_id=0;
 boundary_polyline_pt[0]=new TriangleMeshPolyLine(bound_coords,
                                                  Left_boundary_id);

 // Bottom boundary
 bound_coords[0].resize(2);
 bound_coords[0][0]=X_ll;
 bound_coords[0][1]=Y_ll;

 bound_coords[1].resize(2);
 bound_coords[1][0]=X_ur;
 bound_coords[1][1]=Y_ll;

 // Build the boundary polyline
 Bottom_boundary_id=1;
 boundary_polyline_pt[1]=new TriangleMeshPolyLine(bound_coords,
                                                  Bottom_boundary_id);

 // Right boundary
 bound_coords[0].resize(2);
 bound_coords[0][0]=X_ur;
 bound_coords[0][1]=Y_ll;

 bound_coords[1].resize(2);
 bound_coords[1][0]=X_ur;
 bound_coords[1][1]=Y_ur;
 
 // Build the boundary polyline
 Right_boundary_id=2;
 boundary_polyline_pt[2]=new TriangleMeshPolyLine(bound_coords,
                                                  Right_boundary_id);


 // Right top boundary
 unsigned npt_right=10;
 double x_contact_end_right=0.7;
 double x_contact_end_left=0.3;
 Vector<Vector<double> > right_top_bound_coords(npt_right);
 right_top_bound_coords[0].resize(2);
 right_top_bound_coords[0][0]=X_ur;
 right_top_bound_coords[0][1]=Y_ur;
 for (unsigned j=1;j<npt_right;j++)
  {  
   right_top_bound_coords[j].resize(2);
   double x=X_ur-(X_ur-x_contact_end_right)*double(j)/double(npt_right-1);
   double y=y_c; 
   if (ProblemParameters::Radius_of_elastic_body>0.0)
    {
     y-=sqrt(radius*radius-(x-0.5*(X_ll+X_ur))*(x-0.5*(X_ll+X_ur)));
    }
   right_top_bound_coords[j][0]=x;
   right_top_bound_coords[j][1]=y;
  }

 
 // Build boundary poly line
 Right_top_boundary_id=3;
 TriangleMeshPolyLine* right_top_boundary_pt=
  new TriangleMeshPolyLine(right_top_bound_coords,
                           Right_top_boundary_id);
 boundary_polyline_pt[3]=right_top_boundary_pt;
  

 // Contact boundary
 unsigned npt_contact=50; 
 Vector<Vector<double> > contact_bound_coords(npt_contact);
 contact_bound_coords[0].resize(2);
 contact_bound_coords[0][0]=right_top_bound_coords[npt_right-1][0];
 contact_bound_coords[0][1]=right_top_bound_coords[npt_right-1][1];
 for (unsigned j=1;j<npt_contact;j++)
  {  
   contact_bound_coords[j].resize(2);
   double x=x_contact_end_right-
    (x_contact_end_right-x_contact_end_left)*double(j)/double(npt_contact-1);
   double y=y_c;
   if (ProblemParameters::Radius_of_elastic_body>0.0)
    {
     y-=sqrt(radius*radius-(x-0.5*(X_ll+X_ur))*(x-0.5*(X_ll+X_ur)));
    }
   contact_bound_coords[j][0]=x;
   contact_bound_coords[j][1]=y;
  }
 
 // Build boundary poly line
 Contact_boundary_id=4;
 Contact_boundary_pt=
  new TriangleMeshPolyLine(contact_bound_coords,
                           Contact_boundary_id);
 boundary_polyline_pt[4]=Contact_boundary_pt;
  
 // Keep elements near boundary nice and small
 Maximum_element_length_on_contact_boundary=(X_ur-X_ll)/double(npt_contact);
 Contact_boundary_pt->set_maximum_length(
  Maximum_element_length_on_contact_boundary);



 // Left top boundary
 unsigned npt_left=15; 
 Vector<Vector<double> > top_left_bound_coords(npt_left);
 top_left_bound_coords[0].resize(2);
 top_left_bound_coords[0][0]=contact_bound_coords[npt_contact-1][0];
 top_left_bound_coords[0][1]=contact_bound_coords[npt_contact-1][1];
 for (unsigned j=1;j<npt_left-1;j++)
  {  
   top_left_bound_coords[j].resize(2);
   double x=x_contact_end_left-(x_contact_end_left-X_ll)*
    double(j)/double(npt_left-1);
   double y=y_c; 
   if (ProblemParameters::Radius_of_elastic_body>0.0)
    {
     y-=sqrt(radius*radius-(x-0.5*(X_ll+X_ur))*(x-0.5*(X_ll+X_ur)));
    }
   top_left_bound_coords[j][0]=x;
   top_left_bound_coords[j][1]=y;
  }
 top_left_bound_coords[npt_left-1].resize(2);
 top_left_bound_coords[npt_left-1][0]=X_ll;
 top_left_bound_coords[npt_left-1][1]=Y_ur;
 
 // Build boundary poly line
 Left_top_boundary_id=5;
 TriangleMeshPolyLine* top_left_boundary_pt=
  new TriangleMeshPolyLine(top_left_bound_coords,
                           Left_top_boundary_id);
 boundary_polyline_pt[5]=top_left_boundary_pt;
  
 // Create the triangle mesh polygon for outer boundary
 //----------------------------------------------------
 TriangleMeshPolygon *outer_polygon =
  new TriangleMeshPolygon(boundary_polyline_pt);
  
 // Set the pointer
 closed_curve_pt = outer_polygon;
 
 // Now build the mesh
 //===================

 // Use the TriangleMeshParameters object for helping on the manage of the
 // TriangleMesh parameters
 TriangleMeshParameters triangle_mesh_parameters(closed_curve_pt);

 // Specify the maximum area element
 double uniform_element_area=ProblemParameters::El_area;
 triangle_mesh_parameters.element_area() = uniform_element_area;
 
 // Create the mesh
 Bulk_mesh_pt=
  new RefineableSolidTriangleMesh<ELEMENT>(triangle_mesh_parameters);
  

 // Move nodes on contact boundary on circle
 //-----------------------------------------
 if (ProblemParameters::Radius_of_elastic_body>0.0)
  {
   unsigned b=Contact_boundary_id;
   unsigned nnod=Bulk_mesh_pt->nboundary_node(b);
  for (unsigned j=0;j<nnod;j++)
   {
    SolidNode* nod_pt=Bulk_mesh_pt->boundary_node_pt(b,j);
    double x=nod_pt->x(0);
    double y=y_c-sqrt(radius*radius-(x-0.5*(X_ll+X_ur))*(x-0.5*(X_ll+X_ur)));
    nod_pt->x(1)=y;
   }
  Bulk_mesh_pt->set_lagrangian_nodal_coordinates();
  }

#endif


#ifdef STRUCTURED_MESH
 
 // Setup fake error estimator: 
 unsigned n_extra=1;
 Vector<unsigned> elements_to_refine;
 if (n_x%2==0)
  {
   elements_to_refine.push_back(n_x*n_y-1-n_x/2+1);
   elements_to_refine.push_back(n_x*n_y-1-n_x/2);
   for (unsigned i=0;i<n_extra;i++)
    {
     elements_to_refine.push_back(n_x*n_y-1-n_x/2+1+(i+1));
     elements_to_refine.push_back(n_x*n_y-1-n_x/2+ -(i+1));
    }
  }
 else
  {
   elements_to_refine.push_back(n_x*n_y-1-(n_x+1)/2+1);
   for (unsigned i=0;i<n_extra;i++)
    {
     elements_to_refine.push_back(n_x*n_y-1-(n_x+1)/2+1+(i+1));
     elements_to_refine.push_back(n_x*n_y-1-(n_x+1)/2+1-(i+1));
    }
  }
 unsigned central_node_number=4;
 bool use_lagrangian_coordinates=true;
 Bulk_mesh_pt->spatial_error_estimator_pt()=
  new DummyErrorEstimator(Bulk_mesh_pt,elements_to_refine,
                          central_node_number,
                          use_lagrangian_coordinates);
 
#else

 //hierher
 // Disable the use of an iterative solver for the projection
 // stage during mesh adaptation
 Bulk_mesh_pt->disable_iterative_solver_for_projection();

 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 Bulk_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;

 // Set element size limits
 Bulk_mesh_pt->max_element_size()=ProblemParameters::El_area;
 Bulk_mesh_pt->min_element_size()=0.1*ProblemParameters::El_area;
  
#endif

 // Create the surface mesh as an empty mesh
 Surface_contact_mesh_pt=new Mesh;
 
 // Build 'em 
 create_contact_elements();

 // Create mesh for penetrator element
 Penetrator_mesh_pt=new Mesh;

 // Set boundary condition and complete the build of all elements
 complete_problem_setup();

 // Initial assigment
 ProblemParameters::Y_c=Control_node_pt->x(1);

 // Add the sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_contact_mesh_pt); 
 add_sub_mesh(Penetrator_mesh_pt);

 // Combine all submeshes into a single global Mesh
 build_global_mesh();

 // Do equation numbering
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

} // end of constructor



//=======start_of_doc_solution============================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void ContactProblem<ELEMENT>::doc_solution()
{ 

 oomph_info << "Outputting for step: " << Doc_info.number() << std::endl;
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5;
 
 // Output solution 
 sprintf(filename,"%s/soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();
 
 // Output solution coarsely (only element vertices for easier
 // mesh visualisation)
 sprintf(filename,"%s/coarse_soln%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,2);
 some_file.close();
 
 
 // Output contact elements
 sprintf(filename,"%s/contact%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned nel=Surface_contact_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   dynamic_cast<NonlinearSurfaceContactElement<ELEMENT>* >(
    Surface_contact_mesh_pt->element_pt(e))->output(some_file,20);
  }
 some_file.close();



 // Output integration points of contact elements
 sprintf(filename,"%s/contact_integration_points%i.dat",
         Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 Vector<double> s(1);
 Vector<double> x(2);
 nel=Surface_contact_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   NonlinearSurfaceContactElement<ELEMENT>* el_pt=
    dynamic_cast<NonlinearSurfaceContactElement<ELEMENT>* >(
     Surface_contact_mesh_pt->element_pt(e));
   unsigned nint=el_pt->integral_pt()->nweight();
   for (unsigned j=0;j<nint;j++)
    {
     s[0]=el_pt->integral_pt()->knot(j,0);
     oomph_info << "s_int[ " << j << " " << s[0] << std::endl;
     el_pt->interpolated_x(s,x);
     some_file << x[0] << " " << x[1] << std::endl;
    }
   oomph_info << std::endl;
  }
 some_file.close();

 
 // Output penetrator
 sprintf(filename,"%s/penetrator%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned nplot=500;
 ProblemParameters::Penetrator_pt->output(some_file,nplot);
 some_file.close();
  
 // Output contact elements and assemble total resulting force
 Vector<double> total_contact_force(2,0.0);
 Vector<double> contact_force(2,0.0);
 sprintf(filename,"%s/contact%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 nel=Surface_contact_mesh_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   NonlinearSurfaceContactElement<ELEMENT>* el_pt=
    dynamic_cast<NonlinearSurfaceContactElement<ELEMENT>*>(
     Surface_contact_mesh_pt->element_pt(e));
   el_pt->output(some_file,3);
   el_pt->resulting_contact_force(contact_force);
   total_contact_force[0]+=contact_force[0];
   total_contact_force[1]+=contact_force[1];
  }
 some_file.close();
 
 // Get half-width of Hertz contact region
 double b_hertz=sqrt(4.0*(1.0-ProblemParameters::Nu*ProblemParameters::Nu)/
                     (MathematicalConstants::Pi*
                      (-1.0/ProblemParameters::Radius_of_elastic_body+
                       1.0/ProblemParameters::Radius))*
                     (-total_contact_force[1]));
 oomph_info << "b_hertz " << b_hertz <<  std::endl;
 
 double p_max_hertz = 0.0;
 if(b_hertz!=0.0)
  {
   p_max_hertz=2.0*total_contact_force[1]/
   (MathematicalConstants::Pi*b_hertz);
  }

 // Output Hertzian pressure contact distribution
 sprintf(filename,"%s/hertz%i.dat",Doc_info.directory().c_str(),
         Doc_info.number());
 some_file.open(filename);
 unsigned n=500;
 CircularPenetratorElement* pen_el_pt=
  dynamic_cast<CircularPenetratorElement*>(ProblemParameters::Penetrator_pt);
 Vector<double> centre(2);
 if (pen_el_pt!=0)
  {
   Vector<double> my_centre(pen_el_pt->centre());
   centre[0]=my_centre[0];
   centre[1]=my_centre[1];
  }
 else
  {
   centre[0]=ProblemParameters::Centre[0];
   centre[1]=ProblemParameters::Centre[1];
  }
 double x_c=centre[0];
 double width=2.0*b_hertz;
 for (unsigned j=0;j<n;j++)
  {
   double x=x_c-0.5*width+width*double(j)/double(n-1);
   double p=0.0;
   if (abs((x-x_c))<b_hertz)
    {
     p=p_max_hertz*sqrt(1.0-pow((x-x_c)/b_hertz,2));
    }
   some_file << x << " 0.0 " << p << std::endl;
  }
 some_file.close();


 double target_weight=0.0;
 if (pen_el_pt!=0)
  {
   target_weight=pen_el_pt->target_weight();
  }
 // double angle=0.0;
 // if (pen_el_pt!=0)
 //  {
 //   angle=pen_el_pt->angle();
 //  }
 

 // Write trace file
 Trace_file << total_contact_force[0] << " " 
            << total_contact_force[1] << " " 
            << centre[0] << " "
            << centre[1] << " "
            << target_weight << " " 
  //ALH/hierher Suppressed so that validation tests pass
  //For small changes in loading a different contact point can be selected
  //which means that different angles are computed
  // << angle << " " 
            << std::endl;
 
 //Increment counter for solutions 
 Doc_info.number()++;

} // end of doc_solution



/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////
/// /////////////////////////////////////////////////////////////////////



//=======start_of_main====================================================
/// Driver code
//========================================================================
int main(int argc, char* argv[])
{

 FiniteElement::Accept_negative_jacobian=true;

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);

 // Define possible command line arguments and parse the ones that
 // were actually specified
 
 // Suppress adaptation
 CommandLineArgs::specify_command_line_flag("--no_adapt");
    
 // Initial element size
 CommandLineArgs::specify_command_line_flag("--el_area",
                                            &ProblemParameters::El_area);
    
 // Suppress adaptation
 CommandLineArgs::specify_command_line_flag("--validate");
    
 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 // Create generalised Hookean constitutive equations
 ProblemParameters::Constitutive_law_pt = 
  new GeneralisedHookean(&ProblemParameters::Nu);
      
 // Define centre of penetrator
 ProblemParameters::Centre.resize(2);
 ProblemParameters::Centre[0]=0.5;
 ProblemParameters::Centre[1]=1.0+ProblemParameters::Radius;

#ifdef STRUCTURED_MESH

 // Build problem
 ContactProblem<RefineableQPVDElement<2,3> > problem;

#else

 // Build problem
 ContactProblem<ProjectablePVDElement<TPVDElement<2,3> > > problem;

#endif


 // Update position of centre: Y_c has been determined as the 
 // vertical position of a (control) node on the contact boundary.
 // We place the penetrator vertically above this.
 ProblemParameters::Centre[1]=ProblemParameters::Y_c
  +ProblemParameters::Radius;

  //Output initial condition
 problem.doc_solution();
 
 unsigned max_adapt=1;
 if (CommandLineArgs::command_line_flag_has_been_set("--no_adapt"))
  {
   max_adapt=0;
  }


 // Move position of centre directly
 //----------------------------------
 double dyc=0.00024;
 unsigned nstep=1;
 for (unsigned i=0;i<nstep;i++)
  {
   ProblemParameters::Centre[1]-=dyc;
   oomph_info << "Re-solving imposed circle pos for yc=" 
              << ProblemParameters::Centre[1]
              << std::endl;

   // Solve
   problem.newton_solve(max_adapt);
   
   //Output solution
   problem.doc_solution();
  }

 // Now increase resolution
 //------------------------ 
 nstep=1;
 for (unsigned i=0;i<nstep;i++)
 {

  oomph_info << "Re-solving imposed circle pos for element length factor: "
             << ProblemParameters::Element_length_factor
             << std::endl;

  // Re-solve
  problem.newton_solve(max_adapt);
  
  //Output solution
  problem.doc_solution();

  // Refine 
  ProblemParameters::Element_length_factor/=2.0;
 }

 // Switch to node control
 ProblemParameters::Impose_position_of_centre=false;

 // Use displacement control initially
 //-----------------------------------
 problem.adapt();
 problem.switch_to_displ_control();
 dyc=0.0003; // 0.00006; 
 nstep=1;
 for (unsigned i=0;i<nstep;i++)
  {
   ProblemParameters::Y_c-=dyc;
   oomph_info << "Re-solving for yc=" 
              << ProblemParameters::Y_c
              << std::endl;

   // Solve
   problem.newton_solve(max_adapt);
   
   //Output solution
   problem.doc_solution();
  }
 

 // Switch to force control
 problem.switch_to_force_control();
 
 // Balance it in the horizontal direction
 //---------------------------------------
 ProblemParameters::Horizontal_force=0.0;
 
 oomph_info << "RE-solving for weight=" 
            << dynamic_cast<CircularPenetratorElement*>(
             ProblemParameters::Penetrator_pt)->target_weight() 
            << " and horizontal force: "
            << dynamic_cast<CircularPenetratorElement*>(
             ProblemParameters::Penetrator_pt)->target_horizontal_force() 
            << std::endl;
 
 // Re-solve
 problem.newton_solve(max_adapt);
 
 //Output solution
 problem.doc_solution();
 

 // Now increase weight
 //-------------------- 
 double dweight=0.0003; // 0.0001; 
 nstep=1;
 for (unsigned i=0;i<nstep;i++)
  {
   oomph_info << "Re-solving for weight=" 
              << dynamic_cast<CircularPenetratorElement*>(
               ProblemParameters::Penetrator_pt)->target_weight() 
              << " and horizontal force: "
              << dynamic_cast<CircularPenetratorElement*>(
               ProblemParameters::Penetrator_pt)->target_horizontal_force() 
              << std::endl;
   
   // Re-solve
   problem.newton_solve(max_adapt);
   
   //Output solution
   problem.doc_solution();
   
   // Increase weight
   ProblemParameters::Weight+=dweight;
  }
 
 
// Now increase body force
//------------------------
 nstep=4;
 for (unsigned i=0;i<nstep;i++)
  {
  // Bump up
  ProblemParameters::Body_force_amplitude+=1.0;

  oomph_info << "Re-solving for weight=" 
              << dynamic_cast<CircularPenetratorElement*>(
               ProblemParameters::Penetrator_pt)->target_weight() 
              << " and horizontal force: "
              << dynamic_cast<CircularPenetratorElement*>(
               ProblemParameters::Penetrator_pt)->target_horizontal_force() 
              << " and body force ampl: "
              << ProblemParameters::Body_force_amplitude
              << std::endl;
  
  // Re-solve
  problem.newton_solve(max_adapt);
  
  //Output solution
  problem.doc_solution();
  
 }

 exit(0);

 // Now increase resolution
 //-------------------- 
 nstep=3;
 for (unsigned i=0;i<nstep;i++)
 {

  oomph_info << "Re-solving for weight=" 
             << dynamic_cast<CircularPenetratorElement*>(
              ProblemParameters::Penetrator_pt)->target_weight() 
             << " and horizontal force: "
             << dynamic_cast<CircularPenetratorElement*>(
              ProblemParameters::Penetrator_pt)->target_horizontal_force() 
             << " and body force ampl: "
             << ProblemParameters::Body_force_amplitude
             << " and element length factor: "
             << ProblemParameters::Element_length_factor
             << std::endl;

  // Re-solve
  problem.newton_solve(max_adapt);
  
  //Output solution
  problem.doc_solution();

  // Refine 
  ProblemParameters::Element_length_factor/=2.0;
 }


}
