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

// Oomph-lib includes
#include "generic.h"

using namespace oomph;

//===start_of_main======================================================
/// Driver code: Testing BlockSelector class
//======================================================================
int main(int argc, char* argv[])
{
  // Store command line arguments
  CommandLineArgs::setup(argc,argv);

  CommandLineArgs::specify_command_line_flag(
      "--test_default_constructor");
  CommandLineArgs::specify_command_line_flag(
      "--test_constructor_with_param");
  CommandLineArgs::specify_command_line_flag(
      "--test_select_block_function");
  CommandLineArgs::specify_command_line_flag(
      "--test_want_block_function");
  CommandLineArgs::specify_command_line_flag(
      "--test_do_not_want_block_function");
  CommandLineArgs::specify_command_line_flag(
      "--test_do_not_want_block_function_replace");
  CommandLineArgs::specify_command_line_flag(
      "--test_set_row_index_function");
  CommandLineArgs::specify_command_line_flag(
      "--test_set_column_index_function");


  // Parse the above flags.
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();




  // Testing the default constructor.
  // This sets the Replacement_block_pt to 0
  // then calls BlockSelector::build(0,0,false)
  //
  // Which sets Row_index = 0
  //            Column_index = 0
  //            Wanted = false
  //
  // We check this.
  if(CommandLineArgs::command_line_flag_has_been_set(
     "--test_default_constructor"))
  {
    // Create the output file.
    std::ostringstream out_stream_name;
    out_stream_name << "OUTFILE_default_constructor_test";
    std::ofstream out_file;
    out_file.open(out_stream_name.str().c_str());

    BlockSelector block_selector_default;
    out_file << block_selector_default << "\n";

    out_file.close();
  }
  // Testing the constructor with parameters but no replacement block.
  // We test both wanted and not wanted versions.
  else if(CommandLineArgs::command_line_flag_has_been_set(
          "--test_constructor_with_param"))
  {
    // Create the output file.
    std::ostringstream out_stream_name1;
    out_stream_name1 << "OUTFILE_constructor_with_param_wanted";
    std::ofstream out_file1;
    out_file1.open(out_stream_name1.str().c_str());
    
    BlockSelector block_selector_with_param_wanted(1,2,true);
    out_file1 << block_selector_with_param_wanted << "\n";
    
    out_file1.close();


    // Create the output file.
    std::ostringstream out_stream_name2;
    out_stream_name2 << "OUTFILE_constructor_with_param_not_wanted";
    std::ofstream out_file2;
    out_file2.open(out_stream_name2.str().c_str());
    
    BlockSelector block_selector_with_param_not_wanted(3,4,false);
    out_file2 << block_selector_with_param_not_wanted << "\n";
    
    out_file2.close();


    // Create the output file.
    std::ostringstream out_stream_name3;
    out_stream_name3 << "OUTFILE_constructor_with_param_replace";
    std::ofstream out_file3;
    out_file3.open(out_stream_name3.str().c_str());
    
    CRDoubleMatrix testmat;
    BlockSelector block_selector_with_param_replace(5,6,true,&testmat);
    out_file3 << block_selector_with_param_replace << "\n";
    block_selector_with_param_replace.null_replacement_block_pt();
    
    out_file3.close();
  }
  // Once the BlockSelector is constructed, you can use the select_block()
  // function to select the block you want. We test this.
  else if(CommandLineArgs::command_line_flag_has_been_set(
          "--test_select_block_function"))
  {
    // Create the output file.
    std::ostringstream out_stream_name1;
    out_stream_name1 << "OUTFILE_select_block_wanted";
    std::ofstream out_file1;
    out_file1.open(out_stream_name1.str().c_str());
    
    BlockSelector block_selector_wanted;
    block_selector_wanted.select_block(1,2,true);
    out_file1 << block_selector_wanted << "\n";
    out_file1.close();


    // Create the output file.
    std::ostringstream out_stream_name2;
    out_stream_name2 << "OUTFILE_select_block_not_wanted";
    std::ofstream out_file2;
    out_file2.open(out_stream_name2.str().c_str());
    
    BlockSelector block_selector_not_wanted;
    block_selector_not_wanted.select_block(3,4,false);
    out_file2 << block_selector_not_wanted << "\n";
    out_file2.close();


    // Create the output file.
    std::ostringstream out_stream_name3;
    out_stream_name3 << "OUTFILE_select_block_replace";
    std::ofstream out_file3;
    out_file3.open(out_stream_name3.str().c_str());
    
    CRDoubleMatrix testmat;
    BlockSelector block_selector_replace;
    block_selector_replace.select_block(5,6,true,&testmat);
    out_file3 << block_selector_replace << "\n";
    block_selector_replace.null_replacement_block_pt();
    out_file3.close();
  }
  // Testing  BlockSelector::want_block().
  else if(CommandLineArgs::command_line_flag_has_been_set(
          "--test_want_block_function"))
  {
    // Create the output file.
    std::ostringstream out_stream_name;
    out_stream_name << "OUTFILE_want_block";
    std::ofstream out_file;
    out_file.open(out_stream_name.str().c_str());
    
    BlockSelector block_selector(1,2,false);
    block_selector.want_block();
    out_file << block_selector << "\n";
    out_file.close();
  }
  // Testing  BlockSelector::do_not_want_block().
  // There are two tests:
  // If no replacement block is set, then calling do_not_want_block simply
  // sets Wanted to false.
  // If a replacement block is set, it nulls the pointer, gives a warning 
  // (if paranoid is turned on), then sets Wanted to false.
  else if(CommandLineArgs::command_line_flag_has_been_set(
          "--test_do_not_want_block_function"))
  {
    // Create the output file.
    std::ostringstream out_stream_name;
    out_stream_name << "OUTFILE_do_not_want_block";
    std::ofstream out_file;
    out_file.open(out_stream_name.str().c_str());
    
    BlockSelector block_selector(1,2,true);
    block_selector.do_not_want_block();
    out_file << block_selector << "\n";
    out_file.close(); 
  }
  else if(CommandLineArgs::command_line_flag_has_been_set(
          "--test_do_not_want_block_function_replace"))
  {
    // Create the output file.
    std::ostringstream out_stream_name;
    out_stream_name << "OUTFILE_do_not_want_block_replace";
    std::ofstream out_file;
    out_file.open(out_stream_name.str().c_str());
   
    CRDoubleMatrix testmat; 
    BlockSelector block_selector(1,2,true,&testmat);
    block_selector.do_not_want_block();
    out_file << block_selector << "\n";
    out_file.close(); 
  }
  // Testing  BlockSelector::set_row_index().
  else if(CommandLineArgs::command_line_flag_has_been_set(
          "--test_set_row_index_function"))
  {
    // Create the output file.
    std::ostringstream out_stream_name;
    out_stream_name << "OUTFILE_set_row_index";
    std::ofstream out_file;
    out_file.open(out_stream_name.str().c_str());
   
    BlockSelector block_selector;
    block_selector.set_row_index(42);
    out_file << block_selector << "\n";
    out_file.close(); 
  }
  // Testing  BlockSelector::set_row_index().
  else if(CommandLineArgs::command_line_flag_has_been_set(
          "--test_set_column_index_function"))
  {
    // Create the output file.
    std::ostringstream out_stream_name;
    out_stream_name << "OUTFILE_set_column_index";
    std::ofstream out_file;
    out_file.open(out_stream_name.str().c_str());
   
    BlockSelector block_selector;
    block_selector.set_column_index(42);
    out_file << block_selector << "\n";
    out_file.close(); 
  }
  else
  {
    std::ostringstream err_msg;
    err_msg << "No test is recognised.\n"
            << "Please set the appropriate command line flag." 
            << std::endl;

    throw OomphLibError(err_msg.str(),
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION); 
  }

  return(EXIT_SUCCESS);
} // end_of_main
