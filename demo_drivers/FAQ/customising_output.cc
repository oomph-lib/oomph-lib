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
//Demo code to document customisation of output

// oomph-lib includes
#include "generic.h"
#include "poisson.h"


// The wrapper class for the element has to be included into 
// the oomph-lib namespace
namespace oomph
{

//======customised_poisson=============================================
/// Customised Poisson element -- simply overloads the output function.
/// All other functionality is retained.
//=====================================================================
template<unsigned DIM, unsigned NNODE_1D>
class CustomisedQPoissonElement : public virtual QPoissonElement<DIM,NNODE_1D>
{

public:

 /// Empty constructor
 CustomisedQPoissonElement(){};

 /// Empty virtual constructor
 ~CustomisedQPoissonElement(){};

 /// Overload output function
 void output(std::ostream& output_file)
  {
   output_file << "Hello world" << std::endl;
  }

};

} //end extension of oomph-lib namespace



//======start_of_main==================================================
/// Driver
//=====================================================================
int main()
{

 using namespace oomph;

 // Build the templated object:
 CustomisedQPoissonElement<2,2> element;

 // Call the element's (customised) output function and dump to screen
 element.output(std::cout);

} // end of main









