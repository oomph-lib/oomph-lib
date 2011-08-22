#ifndef OOMPH_QUARTER_PIPE_MESH_TEMPLATE_CC
#define OOMPH_QUARTER_PIPE_MESH_TEMPLATE_CC

#include "quarter_pipe_mesh.template.h"



namespace oomph
{

 //====================================================================
 /// Constructor: Pass number of elements in various directions,
 /// the inner and outer radius and the length of the tube
 //====================================================================
 template<class ELEMENT>
 QuarterPipeMesh<ELEMENT>::QuarterPipeMesh(const unsigned &ntheta,
                                           const unsigned &nr,
                                           const unsigned &nz,  
                                           const double &rmin,
                                           const double &rmax,
                                           const double &length,
                                           TimeStepper* time_stepper_pt) :
  SimpleCubicMesh<ELEMENT>(ntheta,nr,nz,1.0,1.0,length,time_stepper_pt)
 {
  //Variables declaration
  Ntheta=ntheta;
  Nr=nr;
  Nz=nz; 
  Rmin=rmin;
  Rmax=rmax;
  Length=length;
  
  //Build macro element-based domain
  Domain_pt = new QuarterPipeDomain(ntheta,nr,nz,rmin,rmax,length);

  // Loop over all elements 
  unsigned nel=this->nelement();
  for (unsigned e=0;e<nel;e++)
  {
   // Try to cast to FiniteElement
   FiniteElement* el_pt = dynamic_cast<FiniteElement*>(this->element_pt(e));

   // Set macro element pointer
   el_pt->set_macro_elem_pt(Domain_pt->macro_element_pt(e));
  }

  // Update node coordinates with macroelement coordinates,
  // updating solid coordinates too.
  this->node_update(true);
 }
}  


#endif
