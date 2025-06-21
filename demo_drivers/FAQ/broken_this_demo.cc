//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2025 Matthias Heil and Andrew Hazel
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
//Demo code to document problem with missing "this->"
#include<iostream>



//==========Templated_base_class=======================================
/// Some Non-templated base class
//=====================================================================
template<unsigned TEMPLATE_PARAMETER>
class TemplatedBaseClass
{
public:

 /// Empty constructor
 TemplatedBaseClass(){};

 /// Empty virtual constructor
 ~TemplatedBaseClass(){};

 /// Some member function
 void say_hello_world()
  {
   std::cout << "Hello world from base class " << std::endl;
  }

};



//======templated_derived_class=======================================
/// Some templated derived class
//=====================================================================
template<unsigned TEMPLATE_PARAMETER>
class SomeDerivedClass : public virtual TemplatedBaseClass<TEMPLATE_PARAMETER>
{

public:

 // Empty constructor
 SomeDerivedClass(){};

 // Virtual empty constructor
 virtual ~SomeDerivedClass(){};

 /// Some member function
 void output_template_parameter()
  {
   std::cout << "My template parameter is: " 
             << TEMPLATE_PARAMETER << std::endl;

   // Now call the function in the base class

#ifdef USE_BROKEN_VERSION

   // This is illegal according to the C++ standard 
   say_hello_world();

#else 

   // This is stupid but in line with the C++ standard
   this->say_hello_world();
 
#endif

  }

};


//======start_of_main==================================================
/// Driver
//=====================================================================
int main()
{

 // Build the templated object:
 SomeDerivedClass<2> object;

 // Get it to output its template parameter and say hello:
 object.output_template_parameter();

} // end of main

