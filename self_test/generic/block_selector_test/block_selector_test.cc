//LIC//====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC//           Version 0.85. June 9, 2008.
//LIC//
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
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

// Oomph-lib includes
#include "generic.h"

using namespace oomph;

//===start_of_main======================================================
/// Driver code: Testing BlockSelector class
//======================================================================
int main(int argc, char* argv[])
{
  // Create the output file.
  std::ostringstream out_stream;
  out_stream << "OUTPUT";
  std::ofstream out_file;
  out_file.open(out_stream.str().c_str());

  // Testing the default constructor.
  // This sets the Block_pt to 0
  // then calls BlockSelector::build(0,0,false)
  //
  // Which sets Row_index = 0
  //            Column_index = 0
  //            Wanted = false
  //
  // We check this.
  BlockSelector block_selector_default;
  out_file << "Test 1: default constructor:\n";
  out_file << block_selector_default.block_pt() << "\n";
  out_file << block_selector_default.row_index() << "\n";
  out_file << block_selector_default.column_index() << "\n";
  out_file << block_selector_default.wanted() << "\n";
  out_file << std::endl;

  //////////////////////////////////////////////////////////////////////////

  // Testing the second constructor (with parameters)
  BlockSelector block_selector_without_replacement(1,2,true);
  out_file << "Test 2: constructor with param, no replacement:\n";
  out_file << block_selector_without_replacement.block_pt() << "\n";
  out_file << block_selector_without_replacement.row_index() << "\n";
  out_file << block_selector_without_replacement.column_index() << "\n";
  out_file << block_selector_without_replacement.wanted() << "\n";
  out_file << std::endl;
 
  //////////////////////////////////////////////////////////////////////////

  // Testing second constructor given a replacement matrix.
  CRDoubleMatrix testmat;
  BlockSelector block_selector_with_replacement(3,4,true,&testmat);
  out_file << "Test 3: constructor with param, with replacement:\n";
  if(block_selector_with_replacement.block_pt() == &testmat)
  {
    // This worked.
    out_file << "1" << std::endl;
  }
  else
  {
    // This did not work.
    out_file << "0" << std::endl;
  }
  out_file << block_selector_with_replacement.row_index() << "\n";
  out_file << block_selector_with_replacement.column_index() << "\n";
  out_file << block_selector_with_replacement.wanted() << "\n";
  out_file << std::endl;

  //////////////////////////////////////////////////////////////////////////

  // Testing  BlockSelector::want_block() from block_selector_default.
  // This is currently set to false, we print this out to make sure.
  out_file << "Test 4: Test Block_selector::want_block()\n";
  out_file << "block_selector_default.wanted(): "
           << block_selector_default.wanted() << std::endl;
  // Now we want the block.
  block_selector_default.want_block();
  out_file << "block_selector_default.wanted(): "
           << block_selector_default.wanted() << std::endl;
  out_file << std::endl;

  // Set it back so the destructor does not complain when it goes out of 
  // scope.
  block_selector_default.do_not_want_block();

  //////////////////////////////////////////////////////////////////////////
  
  // Now testing the do_not_want_block functionality.
  // The do_not_want_block() gives a warning if a replacement block is set.
  //   it then nulls the block pointer.
  //   then sets the boolean Wanted to false.
  // If no replacement is set, i.e. block_pt() is 0, then it simply set the
  // boolean Wanted to false without complaining. We test both these actions.
  
  // No replacement set:
  out_file << "Test 5.1: Calling do_not_want_block() "
           << "with no replacement set:" << std::endl;
  out_file << "block_selector_without_replacement.wanted(): " 
           << block_selector_without_replacement.wanted() << std::endl;
  block_selector_without_replacement.do_not_want_block();
  out_file << "block_selector_without_replacement.wanted(): " 
           << block_selector_without_replacement.wanted() << std::endl;

  out_file << "Test 5.2: Calling do_not_want_block() "
           << "with replacement block set:" << std::endl;
  out_file << "block_selector_with_replacement.wanted(): "
           << block_selector_with_replacement.wanted() << std::endl;
  block_selector_with_replacement.do_not_want_block();
  out_file << "block_selector_with_replacement.wanted(): "
           << block_selector_with_replacement.wanted() << std::endl;
  out_file << std::endl;
  
  
  //////////////////////////////////////////////////////////////////////////

  // Test 6
  out_file << "Test 6: Test BlockSelector::set_row_index(...)\n";
  // This is currently set to 1.
  out_file << "block_selector_without_replacement.row_index(): "
           << block_selector_without_replacement.row_index() << std::endl;
  // Now change it to 5
  block_selector_without_replacement.set_row_index(5);
  // Output
  out_file << "block_selector_without_replacement.row_index(): "
           << block_selector_without_replacement.row_index() << std::endl;
  out_file << std::endl;

  //////////////////////////////////////////////////////////////////////////

  // Test 7
  out_file << "Test 7: Test BlockSelector::set_column_index(...)\n";
  // This is currently set to 2.
  out_file << "block_selector_without_replacement.column_index(): "
           << block_selector_without_replacement.column_index() << std::endl;
  // Now change it to 6
  block_selector_without_replacement.set_column_index(6);
  // Output
  out_file << "block_selector_without_replacement.column_index(): "
           << block_selector_without_replacement.column_index() << std::endl;
  out_file << std::endl;


  out_file.close();

  // Cleanup: Null the block_pt from block_selector_with_replacement.
  // (Otherwise it will complain when the destructor is called)
  block_selector_with_replacement.null_block_pt();

  return(EXIT_SUCCESS);
} // end_of_main
