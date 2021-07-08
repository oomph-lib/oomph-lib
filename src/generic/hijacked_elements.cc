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
//Non inline functions for Hijacked elements
#include "hijacked_elements.h"

namespace oomph
{

//======================================================================
/// Static default value of the multiplier for the original residuals
/// The contribution to the total residuals is the product of the
/// multiplier and the original value of the residuals
//======================================================================
double HijackedElementBase::Default_residual_multiplier=0.0;

//====================================================================
/// Destructor that cleans up any memory allocated by the class
//===================================================================
HijackedElementBase::~HijackedElementBase() 
{
 //If the hijacked equation number storage has been alloacted
 //clear it
 if(Hijacked_global_eqn_number_pt!=0) 
  {delete Hijacked_global_eqn_number_pt;}
 
 //If the hijacked equation number storage has been alloacted
 //clear it
 if(Hijacked_local_eqn_number_pt!=0) 
  {delete Hijacked_local_eqn_number_pt;}
}

//======================================================================
/// Mark the global equation, addressed by global_eqn_pt, 
/// as hijacked by this element.
//======================================================================
void HijackedElementBase::hijack_global_eqn(long* const &global_eqn_pt)
{
 //If the storage has not been allocated, allocate it
 if(Hijacked_global_eqn_number_pt==0) 
  {
   Hijacked_global_eqn_number_pt = new std::set<long*>;
  }
 
 //Now insert the value, note that this prevents multiple inclusions,
 //which is neater, but possibly inefficient.
 Hijacked_global_eqn_number_pt->insert(global_eqn_pt);
}
 
//=====================================================================
/// The global equation, addressed by global_eqn_pt,
/// is no longer hijacked by this element.
//====================================================================
void HijackedElementBase::unhijack_global_eqn(long* const &global_eqn_pt)
{
 //Check that the storage has been allocated
 if(Hijacked_global_eqn_number_pt!=0)
  {
   Hijacked_global_eqn_number_pt->erase(global_eqn_pt);
  }
}
 
}
