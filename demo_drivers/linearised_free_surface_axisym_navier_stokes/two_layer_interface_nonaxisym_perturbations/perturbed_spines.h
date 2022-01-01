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
// Header file for perturbed spine nodes, elements and meshes

#ifndef OOMPH_PERTURBED_SPINES_HEADER
#define OOMPH_PERTURBED_SPINES_HEADER

// oomph-lib includes
#include "generic.h"

namespace oomph
{
 
 // Forward declaration
 class PerturbedSpineNode;

//=======================================================================
/// PerturbedSpines are used for algebraic node update operations
/// in free-surface fluid problems in which we wish to model a
/// perturbation to the `base' position of the free-surface. Perturbed
/// spines work in an analogous manner to their `standard' counterparts.
/// In addition, each PerturbedSpine stores a pointer to its
/// corresponding `base' Spine.
//=======================================================================
class PerturbedSpine
{
  
  public:
  
 /// Default constructor: Create the PerturbedSpine and initialise
 /// its two "heights" to zero. Pass a pointer to the corresponding (real)
 /// spine in the base problem.
 PerturbedSpine(Spine*& base_spine_pt) : Base_spine_pt(base_spine_pt)
  {
   Geom_data_pt.resize(1);
   
   // Create Data for the two "heights". By default they are free
   Geom_data_pt[0] = new Data(2);
  }
  
 /// Default constructor (unsteady version): Create the
 /// PerturbedSpine and initialise its two "heights" to zero. Pass a
 /// pointer to the corresponding (real) spine in the base problem.
 PerturbedSpine(TimeStepper* const &time_stepper_pt,
                Spine*& base_spine_pt) : Base_spine_pt(base_spine_pt)
  {
   Geom_data_pt.resize(1);
   
   // Create Data for the two "heights". By default they are free
   Geom_data_pt[0] = new Data(time_stepper_pt,2);
  }
    
 /// Constructor: Create the PerturbedSpine and initialise its
 /// heights to the specified value. Pass a pointer to the corresponding
 /// (real) spine in the base problem.
 PerturbedSpine(Spine*& base_spine_pt,
                const double& height) : Base_spine_pt(base_spine_pt)
  {
   Geom_data_pt.resize(1);
   
   // Create Data for the two "heights". By default they are free
   Geom_data_pt[0] = new Data(2);
   
   // Set value
   for(unsigned i=0;i<2;i++) { Geom_data_pt[0]->set_value(i,height); }
  }
     
 /// Destructor: Wipe Data object that stores the
 /// PerturbedSpine height. All other objects (geometric Data and
 /// geometric objects) were created outside the PerturbedSpine
 /// and must be deleted there.
 ~PerturbedSpine() { delete Geom_data_pt[0]; }

 /// Access function to pointer to base spine
 Spine* &base_spine_pt() { return Base_spine_pt; }

 /// Access function to i-th component of spine "height"
 double& height(const unsigned& i)
  {
   return *(Geom_data_pt[0]->value_pt(i));
  }

 /// Access function to t-th history value of the i-th component of
 /// the perturbed spine "height"
 double& height(const unsigned& t, const unsigned& i)
  {
   return *(Geom_data_pt[0]->value_pt(t,i));
  }

 /// Access function to Data object that stores the spine "heights"
 Data*& height_pt() { return Geom_data_pt[0]; }

 /// Access function to Data object that stores the spine "heights"
 /// (const version)
 Data* height_pt() const { return Geom_data_pt[0]; }

 /// Access function to SpineNode at bottom of spine
 PerturbedSpineNode*& node_at_bottom_of_spine_pt()
  {
   return Node_at_bottom_of_spine_pt;
  }

 /// Access function to SpineNode at bottom of spine
 /// (const version)
 PerturbedSpineNode* node_at_bottom_of_spine_pt() const
  {
   return Node_at_bottom_of_spine_pt;
  }

 /// Access function to SpineNode at top of spine
 PerturbedSpineNode*& node_at_top_of_spine_pt()
  {
   return Node_at_top_of_spine_pt;
  }

 /// Access function to SpineNode at top of spine
 /// (const version)
 PerturbedSpineNode* node_at_top_of_spine_pt() const
  {
   return Node_at_top_of_spine_pt;
  }

 /// Number of geometric Data that is involved in the
 /// node update operations for this PerturbedSpine
 unsigned ngeom_data() { return Geom_data_pt.size(); }

 /// Set vector of (pointers to) geometric Data that is
 /// involved in the node update operations for this PerturbedSpine.
 /// Wipes any previously existing geometric Data.
 void set_geom_data_pt(const Vector<Data*>& geom_data_pt)
  {
   unsigned n_geom_data=geom_data_pt.size();
   Geom_data_pt.resize(n_geom_data+1);
   for (unsigned i=1;i<n_geom_data;i++)
    {
     Geom_data_pt[i+1]=geom_data_pt[i];
    }
  }

 /// Add (pointer to) geometric Data that is
 /// involved in the node update operations for this PerturbedSpine
 void add_geom_data_pt(Data* geom_data_pt)
  {
   Geom_data_pt.push_back(geom_data_pt);
  }

 /// Return i-th geometric Data that is involved in the
 /// node update operations for this PerturbedSpine
 Data*& geom_data_pt(const unsigned& i){return Geom_data_pt[i];}

 /// Return i-th geometric Data that is involved in the
 /// node update operations for this PerturbedSpine. Const version
 Data* geom_data_pt(const unsigned& i) const {return Geom_data_pt[i];}
 
 /// Return the vector of geometric data
 Vector<Data*> &vector_geom_data_pt() {return Geom_data_pt;}

 /// Number of geometric objects that is involved in the
 /// node update operations for this PerturbedSpine
 unsigned ngeom_object(){return Geom_object_pt.size();}

 /// Set vector of (pointers to) geometric objects that is
 /// involved in the node update operations for this PerturbedSpine
 void set_geom_object_pt(const Vector<GeomObject*>& geom_object_pt)
  {
   unsigned n_geom_object=geom_object_pt.size();
   Geom_object_pt.resize(n_geom_object);
   for (unsigned i=0;i<n_geom_object;i++)
    {
     Geom_object_pt[i]=geom_object_pt[i];
    }
  }

 /// Add (pointer to) geometric object that is
 /// involved in the node update  operations for this PerturbedSpine
 void add_geom_object_pt(GeomObject* geom_object_pt)
  {
   Geom_object_pt.push_back(geom_object_pt);
  }

 /// Return i-th geometric object that is involved in the
 /// node update operations for this PerturbedSpine
 GeomObject*& geom_object_pt(const unsigned& i) { return Geom_object_pt[i]; }

 /// Return i-th geometric object that is involved in the
 /// node update operations for this PerturbedSpine. Const version
 GeomObject* geom_object_pt(const unsigned& i) const
  {
   return Geom_object_pt[i];
  }

 /// Return vector of all geometric objects that affect this spine
 Vector<GeomObject*> &vector_geom_object_pt() { return Geom_object_pt; }

 /// Number of geometric parameters that are involved in the
 /// node update operations for this PerturbedSpine
 unsigned ngeom_parameter() { return Geom_parameter.size(); }

 /// Set vector of geometric parameters that are
 /// involved in the node update operations for this PerturbedSpine.
 /// Wipes any previously existing geometric parameters
 void set_geom_parameter(const Vector<double>& geom_parameter)
  { Geom_parameter = geom_parameter; }

 /// Add geometric parameter 
 /// involved in the node update operations for this PerturbedSpine
 void add_geom_parameter(const double &geom_parameter)
  { Geom_parameter.push_back(geom_parameter); }

 /// Return i-th geometric parameter that is involved in the
 /// node update operations for this PerturbedSpine
 double& geom_parameter(const unsigned& i)
  { return Geom_parameter[i]; }

 /// Return i-th geometric parameter that is involved in the
 /// node update operations for this PerturbedSpine. Const version
 const double &geom_parameter(const unsigned& i) const 
  { return Geom_parameter[i]; }


  private:

 /// Pointer to corresponding (real) spine in the base problem
 Spine* Base_spine_pt;

 /// Vector that stores the pointers to additional geometric Data
 Vector<Data*> Geom_data_pt;

 /// Vector that stores the pointers to geometric objects that is
 /// involved in the node update operation
 Vector<GeomObject*> Geom_object_pt;

 /// Vector that stores doubles that are used in the geometric
 /// updates
 Vector<double> Geom_parameter;

 /// Store a pointer to the nodes at either end of the spine so that a
 /// unit vector in the direction of the spine can be calculated
 PerturbedSpineNode* Node_at_bottom_of_spine_pt;
 PerturbedSpineNode* Node_at_top_of_spine_pt;

}; // End of PerturbedSpine class definition



// Forward declaration
class PerturbedSpineMesh;

//=======================================================================
/// Class for nodes that `live' on perturbed spines.
//=======================================================================
class PerturbedSpineNode : public Node
{
  private:

 /// Private internal data pointer to a perturbed spine
 PerturbedSpine* PerturbedSpine_pt;

 /// Private double that represents the fixed fraction along the spine
 double Fraction;

 /// Pointer to PerturbedSpineMesh that this node is a part of.
 /// (The mesh implements the node update function(s))
 PerturbedSpineMesh* PerturbedSpine_mesh_pt;

 /// ID of node update function (within specific mesh -- useful if there
 /// are multiple node update functions, e.g. in two-layer problems.
 unsigned Node_update_fct_id;

 
  public:

 /// Steady Constructor, initialise pointers to zero
 PerturbedSpineNode(const unsigned &n_dim,
                    const unsigned &n_position_type,
                    const unsigned &initial_nvalue) :
  Node(n_dim,n_position_type,initial_nvalue),
  PerturbedSpine_pt(0), Fraction(0),
  PerturbedSpine_mesh_pt(0), Node_update_fct_id(0) {}
  
 /// Unsteady Constructor, initialise pointers to zero
 PerturbedSpineNode(TimeStepper* const &time_stepper_pt,
                    const unsigned &n_dim,
                    const unsigned &n_position_type, 
                    const unsigned &initial_nvalue) :
  Node(time_stepper_pt,n_dim,n_position_type,initial_nvalue),
  PerturbedSpine_pt(0),
  Fraction(0), PerturbedSpine_mesh_pt(0), Node_update_fct_id(0) {}
   
 /// Access function to perturbed spine
 PerturbedSpine* &perturbed_spine_pt() { return PerturbedSpine_pt; }
   
 /// Set reference to fraction along spine
 double &fraction() { return Fraction; }
   
 /// Access function to ID of node update function (within specific mesh)
 unsigned& node_update_fct_id() { return Node_update_fct_id; }
   
 /// Access function to Pointer to PerturbedSpineMesh that this   
 /// node is a part of and which implements the node update function(s)
 PerturbedSpineMesh*& spine_mesh_pt() { return PerturbedSpine_mesh_pt; }

 /// Access function to i-th component of spine "height"
 double &h(const unsigned& i) { return PerturbedSpine_pt->height(i); }
   
 /// Overload the node update function, call 
 /// the update function in the Node's PerturbedSpineMesh
 void node_update(const bool& update_all_time_levels_for_new_node=false);
   
 /// Return the number of geometric data, zero if no spine.
 unsigned ngeom_data() const
 {
  if(PerturbedSpine_pt) { return PerturbedSpine_pt->ngeom_data(); }
  else { return 0; }
 }
 
 /// Return the number of geometric objects, zero if no spine.
 unsigned ngeom_object() const
 {
  if(PerturbedSpine_pt) { return PerturbedSpine_pt->ngeom_object(); }
  else { return 0; }
 }

 /// Return the vector of all geometric data
 Data** all_geom_data_pt() { return &(PerturbedSpine_pt->geom_data_pt(0)); }
   
 /// Return the vector of all geometric objects
 GeomObject** all_geom_object_pt()
  {
   return &(PerturbedSpine_pt->geom_object_pt(0));
  }
   
}; // End of PerturbedSpineNode class definition



/// ///////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////



//=======================================================================
/// A policy class that serves only to establish the interface for 
/// assigning the perturbed spine equation numbers
//=======================================================================
class PerturbedSpineFiniteElement
{
  public:

 /// Empty constructor
 PerturbedSpineFiniteElement() {}

 /// Emtpty virtual destructor
 virtual ~PerturbedSpineFiniteElement(){}

};



//=======================================================================
/// The PerturbedSpineElement<ELEMENT> class takes an existing
/// element as a template parameter and adds the necessary additional
/// functionality to allow the element to be update using the method of
/// spines. A vector of pointers to spines and storage for the local
/// equation numbers associated with the spines are added to the element. 
//=======================================================================
template<class ELEMENT>
class PerturbedSpineElement : 
public ElementWithSpecificMovingNodes<ELEMENT,PerturbedSpineNode>,
 public PerturbedSpineFiniteElement
{
  private:
 
 /// Array to hold the index of the geometric data associated with
 /// the spine height of the spine that affects the n-th node
 unsigned *PerturbedSpine_geometric_index;

 /// Complete the setup of additional dependencies. Overloads
 /// empty virtual function in GeneralisedElement to determine the "geometric 
 /// Data", i.e. the Data that affects the element's shape.
 /// This function is called (for all elements) at the very beginning of the
 /// equation numbering procedure to ensure that all dependencies
 /// are accounted for.
 void complete_setup_of_dependencies();

 
  public:

 /// Constructor, call the constructor of the base element
 PerturbedSpineElement() : 
  ElementWithSpecificMovingNodes<ELEMENT,PerturbedSpineNode>(),  
  PerturbedSpineFiniteElement(),
  PerturbedSpine_geometric_index(0) {}
 
 /// Constructor used for spine face elements
 PerturbedSpineElement(FiniteElement* const &element_pt,
                       const int &face_index) : 
  ElementWithSpecificMovingNodes<ELEMENT,PerturbedSpineNode>(element_pt,
                                                             face_index),
  PerturbedSpineFiniteElement(),
  PerturbedSpine_geometric_index(0) {}

  /// Destructor, clean up the storage allocated to the local equation numbers
 ~PerturbedSpineElement()
  {
   if(PerturbedSpine_geometric_index) 
    {
     delete[] PerturbedSpine_geometric_index;
    }
  }

 /// Return the local equation number corresponding to the i-th
 /// component of the "height" of the perturbed spine at the n-th node
 inline int spine_local_eqn(const unsigned &n, const unsigned &i)
  {
#ifdef RANGE_CHECKING
   const unsigned n_node = this->nnode();
   if(n >= n_node)
    {
     std::ostringstream error_message;
     error_message << "Range Error:  Node number " << n
                   << " is not in the range (0,"
                   << n_node-1 << ")";
     throw OomphLibError(error_message.str(),
                         "PerturbedSpineElement::spine_local_eqn()",
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif

#ifdef PARANOID
   // The perturbed spine "height" only has two components
   if(i>=2)
    {
     std::ostringstream error_message;
     error_message << "Perturbed spines only have two height components:"
                   << "\nComponent i=" << i
                   << " is not in the range (0,2)";
     throw OomphLibError(error_message.str(),
                         "PerturbedSpineElement::spine_local_eqn()",
                         OOMPH_EXCEPTION_LOCATION);
    }
   // If there is no spine then you can't get the local equation
   if(PerturbedSpine_geometric_index[n]==this->ngeom_data())
    {
     std::ostringstream error_stream;
     error_stream 
      << "PerturbedSpineNode " << n << " does not have a PerturbedSpine\n"
      << "attached, so you can't get its local equation number. Check that\n"
      << "the Mesh is correctly associating PerturbedSpines with is Nodes\n";
     throw OomphLibError(error_stream.str(),
                         "PerturbedSpineElement<ELEMENT>::spine_local_eqn()",
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif
   return this->geometric_data_local_eqn(PerturbedSpine_geometric_index[n],i);
  }


 /// Return the local equation number associated with the i-th `height' of
 /// the n-th local perturbed spine
 virtual int get_local_eqn_number_corresponding_to_geometric_dofs(
  const unsigned &n, const unsigned &i)
 {
  return this->spine_local_eqn(n,i);
 }

}; // End of PerturbedSpineElement<ELEMENT> class definition



//=======================================================================
/// Explicit definition of the face geometry for spine elements:
/// The same as the face geometry of the underlying element
//=======================================================================
template<class ELEMENT>
class FaceGeometry<PerturbedSpineElement<ELEMENT> >:
public virtual FaceGeometry<ELEMENT>  
{
  public:

 /// Constructor
 FaceGeometry() : FaceGeometry<ELEMENT>() {}
};

//=======================================================================
/// Explicit definition of the face geometry for spine elements:
/// The same as the face geometry of the underlying element
//=======================================================================
template<class ELEMENT>
class FaceGeometry<FaceGeometry<PerturbedSpineElement<ELEMENT> > >:
public virtual FaceGeometry<FaceGeometry<ELEMENT> > 
{
  public:

 /// Constructor
 FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT> >() {}
};

//=======================================================================
/// Explicit definition of the face geometry for spine elements:
/// The same as the face geometry of the underlying element
//=======================================================================
template<class ELEMENT>
class FaceGeometry<PerturbedSpineElement<FaceGeometry<ELEMENT> > >:
public virtual FaceGeometry<FaceGeometry<ELEMENT> > 
{
  public:

 /// Constructor
 FaceGeometry() : FaceGeometry<FaceGeometry<ELEMENT> >() {}
};



/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////



//========================================================================
/// General PerturbedSpineMesh class.
///
/// Derived from Mesh with virtual so that perturbed
/// spine meshes can be derived from general meshes, without
/// multiple copies of Mesh objects.
//========================================================================
class PerturbedSpineMesh : public virtual Mesh
{
  protected:
 
 /// A PerturbedSpine mesh contains a Vector of pointers to perturbed spines
 Vector<PerturbedSpine*> PerturbedSpine_pt;

  public:
 
 /// Destructor to clean up the memory allocated to the perturbed spines
 virtual ~PerturbedSpineMesh()
  {
   // Set the range of PerturbedSpine_pt
   unsigned long PerturbedSpine_pt_range = PerturbedSpine_pt.size();

   // Loop over the entries in reverse and free memory
   for(unsigned long i=PerturbedSpine_pt_range;i>0;i--)
    {
     delete PerturbedSpine_pt[i-1]; PerturbedSpine_pt[i-1] = 0;
    }
  }

 /// Return the i-th perturbed spine in the mesh
 PerturbedSpine* const & perturbed_spine_pt(const unsigned long &i) const
  {
   return PerturbedSpine_pt[i];
  }

 /// Return the number of perturbed spines in the mesh
 unsigned long nspine() const { return PerturbedSpine_pt.size(); }

 /// Add a perturbed spine to the mesh
 void add_perturbed_spine_pt(PerturbedSpine* const &spine_pt)
 {
  PerturbedSpine_pt.push_back(spine_pt);
 }
 
 /// Return a pointer to the n-th global PerturbedSpineNode 
 PerturbedSpineNode* node_pt(const unsigned long &n) 
  {
#ifdef PARANOID
   if(!dynamic_cast<PerturbedSpineNode*>(Node_pt[n]))
    {
     std::ostringstream error_message;
     error_message << "Node " << n << "is a "
                   << typeid(Node_pt[n]).name() 
                   << ", not a PerturbedSpineNode" << std::endl;
     
     throw OomphLibError(error_message.str(),
                         "PerturbedSpineMesh::node_pt()",
                         OOMPH_EXCEPTION_LOCATION);
    } 
#endif

   // Return a cast to the pointer to the node
   return (dynamic_cast<PerturbedSpineNode*>(Node_pt[n]));
  }

 /// Return the n-th local PerturbedSpineNode in element e. This is
 /// required to cast the nodes in a spine mesh to be PerturbedSpineNodes
 /// and therefore allow access to the extra PerturbedSpineNode data.
 PerturbedSpineNode* element_node_pt(const unsigned long &e,
                                     const unsigned &n)
  {
#ifdef PARANOID
   // Try to cast to FiniteElement
   FiniteElement* el_pt=dynamic_cast<FiniteElement*>(Element_pt[e]);
   if (el_pt==0)
    {
     throw OomphLibError(
      "Can't execute element_node_pt(...) for non FiniteElements",
      "PerturbedSpineMesh::element_node_pt()",
      OOMPH_EXCEPTION_LOCATION);
    }
   if(!dynamic_cast<PerturbedSpineNode*>(el_pt->node_pt(n)))
    {
     std::ostringstream error_message;
     error_message << "Node " << n << "is a "
                   << typeid(el_pt->node_pt(n)).name() 
                   << ", not a PerturbedSpineNode" << std::endl;
     
     throw OomphLibError(error_message.str(),
                         "PerturbedSpineMesh::node_pt()",
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif
   // Return a cast to the node pointer
   return(dynamic_cast<PerturbedSpineNode*>(
           dynamic_cast<FiniteElement*>(Element_pt[e])->node_pt(n)));
  }
 
 /// Assign equation numbers for perturbed spines
 unsigned long assign_global_eqn_numbers(Vector<double *> &Dof_pt);
 
 /// Update function to update all nodes of mesh
 /// [Doesn't make sense to use this mesh with SolidElements anyway,
 /// so we buffer the case if update_all_solid_nodes is set to 
 /// true.]
 void node_update(const bool& update_all_solid_nodes=false);

 /// Update function for given spine node -- this must be implemented
 /// by all specific PerturbedSpineMeshes.
 virtual void perturbed_spine_node_update(PerturbedSpineNode* spine_node_pt)=0;

 /// Overload the dump function so that the spine data is dumped
 void dump(std::ofstream &dump_file,const bool &use_old_ordering=true) const;

 /// Overload the read function so that the spine data is read 
 /// from the restart file
 void read(std::ifstream &restart_file);

}; // End of PerturbedSpineMesh class definition




} // End of oomph namespace

#endif
