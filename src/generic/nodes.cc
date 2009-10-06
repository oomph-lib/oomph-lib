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
//Functions for the Node/Data/SolidNode classes

#include<algorithm>

//oomph-lib headers
#include "nodes.h"
#include "timesteppers.h"


namespace oomph
{

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
//Functions for the Data class
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
 
 //=================================================================
 /// \short Private function to check that the arguments are within
 /// the range of the stored data values and timesteps.
 //=================================================================
 void Data::range_check(const unsigned &t, const unsigned &i) const
 {
  //If either the value or the time history value are out of range
  if((i>= Nvalue) || (t >= ntstorage()))
   {
    std::ostringstream error_message;
    //Value range check
    if(i>= Nvalue)
     {
      error_message << "Range Error: Value " << i
                    << " is not in the range (0,"
                    << Nvalue-1 << ")";
     }
    //Time range check
    if(t >= ntstorage())
     {
      error_message << "Range Error: Time Value " << t
                    << " is not in the range (0,"
                    << ntstorage() - 1 << ")";
     }
    //Throw the error
    throw OomphLibError(error_message.str(),
                        "Data::range_check()",
                        OOMPH_EXCEPTION_LOCATION);
   }
 }


//====================================================================
/// Add the pointer data_pt to the internal storage used to keep track
/// of copies of the Data object.
//====================================================================
 void Data::add_copy(Data* const &data_pt)
 {
  //Find the current number of copies
  const unsigned n_copies = Ncopies;
  //Allocate storage for the pointers to the new number of copies
  Data** new_copy_of_data_pt = new Data*[n_copies+1];
  //Copy over the exisiting pointers
  for(unsigned i=0;i<n_copies;i++)
   {new_copy_of_data_pt[i] = Copy_of_data_pt[i];}
  //Add the new pointer to the end
  new_copy_of_data_pt[n_copies] = data_pt;

  //Delete the old storage 
  delete[] Copy_of_data_pt;
  //Allocate the new storage
  Copy_of_data_pt = new_copy_of_data_pt;
  //Increase the number of copies
  ++Ncopies;
 }

//=====================================================================
/// Remove the pointer data_pt from the internal storage used to keep
/// track of copies
//=====================================================================
 void Data::remove_copy(Data* const &data_pt)
 {
  //Find the current number of copies
  const unsigned n_copies = Ncopies; 
  //Index of the copy
  unsigned data_index = n_copies;
  //Check that the existing data is actually a copy
  for(unsigned i=0;i<n_copies;i++)
   {
    if(Copy_of_data_pt[i]==data_pt) {data_index = i; break;}
   }

  //If we have not found an index throw an error
  if(data_index==n_copies)
   {
    std::ostringstream error_stream;
    error_stream << "Data pointer " << data_pt 
                 << " is not stored as a copy of the data object " << this 
                 << std::endl;
    throw OomphLibError(error_stream.str(),
                        "Data::remove_copy()",
                        OOMPH_EXCEPTION_LOCATION);
   }

  //If we still here remove the data
  //Allocate storage for the pointers to the new number of copies
  Data** new_copy_of_data_pt = new Data*[n_copies-1];
  
  unsigned index=0;
  //Copy over the exisiting pointers
  for(unsigned i=0;i<n_copies;i++)
   {
    //If we are not at the copied index
    if(i!=data_index)
     {
      //Copy the data across
      new_copy_of_data_pt[index] = Copy_of_data_pt[i];
      //Increase the index
      ++index;
     }
   }

  //Delete the old storage 
  delete[] Copy_of_data_pt;
  //Allocate the new storage
  Copy_of_data_pt = new_copy_of_data_pt;
  //Set the new number of copies
  --Ncopies;
 }

//================================================================
/// \short Helper function that should be overloaded in classes
/// that contain copies of Data. The function must 
/// reset the internal pointers to the copied data. This is used
/// when resizing data to ensure that all the pointers remain valid.
/// The base Data class cannot be a copy, so throw an error 
 //==================================================================
 void Data::reset_copied_pointers() 
 {
  throw OomphLibError("Data can never be a copy",
                      "Data::reset_copied_pointers()",
                      OOMPH_EXCEPTION_LOCATION);
 }

//=======================================================================
/// \short Helper function that should be overloaded classes
/// that contain copies of data. The function must
/// unset (NULL out) the internal pointers to the copied data.
/// This is used when destructing data to ensure that all pointers remain
/// valid.
//======================================================================
 void Data::clear_copied_pointers() 
 {
  throw OomphLibError("Data can never be a copy",
                      "Data::clear_copied_pointers()",
                      OOMPH_EXCEPTION_LOCATION);
 }

//================================================================
/// Delete all the storage allocated by the Data object and 
/// set its pointers to NULL
//================================================================
 void Data::delete_value_storage()
 {
  //If we have nulled out the storage already return immediately
  if((Value==0) && (Eqn_number==0)) {return;}
  
  //Delete the double storage arrays at once (they were allocated at once)
  delete[] Value[0];
  //Delete the pointers to the arrays.
  delete[] Value; delete[] Eqn_number;
  //Null out the pointers
  Value = 0; Eqn_number = 0;
 }
 
//================================================================
/// Default (steady) timestepper for steady Data
//================================================================
TimeStepper* Data::Default_static_time_stepper_pt=new Steady<0>();

//================================================================
/// Static "Magic number" to indicate pinned values
//================================================================
 long Data::Is_pinned=-1;

//================================================================
/// \short Static "Magic number" to indicate values that haven't 
/// been classified as pinned or free
//================================================================
long Data::Is_unclassified=-10;

//================================================================
/// Static "Magic number" to indicate that the value is constrained,
/// usually because is it associated with non-conforming data,
/// otherwise known as hanging nodes
//================================================================
long Data::Is_constrained=-2;

//================================================================
/// Default constructor.
//================================================================
 Data::Data() : Value(0), Eqn_number(0), 
                Time_stepper_pt(Data::Default_static_time_stepper_pt),
                Copy_of_data_pt(0),
                Nvalue(0), Ncopies(0)
#ifdef OOMPH_HAS_MPI
                , Is_halo(false)
#endif

 {}

//================================================================
/// Default constructor for steady problems. Memory is assigned for a given 
/// number of values, which are assumed to be free (not pinned)
//================================================================
 Data::Data(const unsigned &initial_n_value) : 
  Value(0), Eqn_number(0),
  Time_stepper_pt(Data::Default_static_time_stepper_pt),
  Copy_of_data_pt(0),
  Nvalue(initial_n_value),
  Ncopies(0)
#ifdef OOMPH_HAS_MPI
  , Is_halo(false)
#endif
 {
  //Only bother to do something if there are values
  if(initial_n_value > 0)
   {
   //Allocate initial_n_value values in the value and equation number
   //storage schemes.
   Value = new double*[initial_n_value];
   Eqn_number = new long[initial_n_value]; 
  
  //Allocate contiguous arrays of doubles and longs to
  //hold the data values.
  double *values = new double[initial_n_value];

  //Set the pointers to the data values and equation numbers
  //and initialise the actual values.
  for(unsigned i=0;i<initial_n_value;i++) 
   {
    //Set the pointer from the address in the contiguous array
    Value[i] = &values[i];
    //Initialise the value to zero
    Value[i][0] = 0.0;
    //Initialise the equation number to Is_unclassified
    Eqn_number[i] = Is_unclassified;
   }
 }
 }

//================================================================
/// Constructor for unsteady problems. Memory is assigned for a given 
/// number of values; and the additional storage required by the Timestepper.
/// The values are assumed to be free (not pinned).
//================================================================
Data::Data(TimeStepper* const &time_stepper_pt, 
           const unsigned &initial_n_value,
           const bool &allocate_storage)
 : 
 Value(0), Eqn_number(0), Time_stepper_pt(time_stepper_pt),
 Copy_of_data_pt(0),
 Nvalue(initial_n_value), 
 Ncopies(0)
#ifdef OOMPH_HAS_MPI
 , Is_halo(false)
#endif
{
 //If we are in charge of allocating the storage,
 //and there are data to allocate, do so
 if((allocate_storage) && (initial_n_value > 0))
  {
   //Allocate storage for initial_n_value equation numbers
   Eqn_number = new long[initial_n_value];
   
   //Locally cache the number of time history values
   const unsigned n_tstorage = ntstorage();
   
   //There will be initial_n_value pointers each addressing and array
   //of n_tstorage doubles.
   Value = new double*[initial_n_value];
   
   //Allocate all the data values in one big array for data locality.
   double *values = new double[initial_n_value*n_tstorage];

   //Set the pointers to the data values and equation numbers
   for(unsigned i=0;i<initial_n_value;i++)
    {
     //Set the pointers to the start of the time history values
     //allocated for each value.
     Value[i] = &values[i*n_tstorage];
     //Initialise all values to zero
     for(unsigned t=0;t<n_tstorage;t++) {Value[i][t] = 0.0;}
     //Initialise the equation number to be unclassified.
     Eqn_number[i] = Is_unclassified;
    }
  }
}

//================================================================
///Virtual destructor, deallocates memory assigned for data
//================================================================
Data::~Data()
{
 //If we have any copies clear their pointers 
 for(unsigned i=0;i<Ncopies;i++) 
  {
   Copy_of_data_pt[i]->clear_copied_pointers();
  }

 //Now delete the storage allocated for pointers to the copies
 delete[] Copy_of_data_pt; Copy_of_data_pt=0;

 //Clean up the allocated storage
 delete_value_storage();
}


//==================================================================
/// Compute Vector of values (dofs or pinned) at this Data object
//==================================================================
void Data::value(Vector<double>& values) const
{
 //Loop over all the values and set the appropriate value
 const unsigned n_value = nvalue();
 for(unsigned i=0;i<n_value;i++) {values[i] = value(i);}
}

//==================================================================
/// Compute Vector of values (dofs or pinned) at this node
/// at time level t (t=0: present; t>0: previous)
//==================================================================
void Data::value(const unsigned& t, Vector<double>& values) const
{
 //Loop over all the values and set the value at time level t
 const unsigned n_value = nvalue();
 for(unsigned i=0;i<n_value;i++) {values[i] = value(t,i);}
}

//================================================================
/// Copy Data values from specified Data object
//================================================================
void Data::copy(Data* orig_data_pt)
{

 //Find the amount of data stored 
 const unsigned n_value = nvalue();
 const unsigned n_time = ntstorage();

 // Check # of values:
 const unsigned long n_value_orig=orig_data_pt->nvalue();
 if (n_value!=n_value_orig)
  {
   std::ostringstream error_stream;
   error_stream << "The number of values, " << n_value 
                << " is not the same of those in the original data "
                << n_value_orig << std::endl;
   
   throw OomphLibError(error_stream.str(),
                       "Data::copy()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 const unsigned long n_time_orig=orig_data_pt->ntstorage();
 if (n_time!=n_time_orig)
  {
   std::ostringstream error_stream;
   error_stream << "The number of time history values, " << n_time 
                << " is not the same of those in the original data "
                << n_time_orig << std::endl;
   
   throw OomphLibError(error_stream.str(),
                       "Data::copy()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Read data
 for(unsigned t=0;t<n_time;t++)
  {
   for(unsigned j=0;j<n_value;j++) 
    {
     set_value(t,j,orig_data_pt->value(t,j));
    }
  }
}


//================================================================
///Dump data object to a file 
//================================================================
void Data::dump(std::ostream& dump_file)
{
 //Find the amount of storage used
 const unsigned value_pt_range = nvalue();
 const unsigned time_steps_range = ntstorage();

 //Only write data if there is some stored
 if (value_pt_range*time_steps_range > 0)
  {
   dump_file << value_pt_range << " # number of data values" << std::endl;
   dump_file << time_steps_range << " # number of doubles for time history" 
           << std::endl;

   // Write data
   for(unsigned t=0;t<time_steps_range;t++)
    {
     for(unsigned j=0;j<value_pt_range;j++) 
      {
       dump_file << value(t,j) << std::endl ;
      }
    }
  }
}

//================================================================
///Read data object from file 
//================================================================
void Data::read(std::ifstream& restart_file)
{
 std::string input_string;

 //Find the amount of data stored 
 const unsigned value_pt_range = nvalue();
 const unsigned time_steps_range = ntstorage();

 //Only read in data if there is some storage available
 if (value_pt_range*time_steps_range > 0)
  {
   // Read line up to termination sign
   getline(restart_file,input_string,'#');
   // Ignore rest of line
   restart_file.ignore(80,'\n');
   // Check # of values:
   const unsigned long check_nvalues=atoi(input_string.c_str());
   if (check_nvalues!=value_pt_range)
    {
     std::ostringstream error_stream;
     error_stream 
      << "Number of values stored in dump file is not equal to the amount "
      << "of storage allocated in Data object "
      <<  check_nvalues << " " << value_pt_range << std::endl;
     
     throw OomphLibError(error_stream.str(),
                         "Data::read()",
                         OOMPH_EXCEPTION_LOCATION);
    }
   
   // Read line up to termination sign
   getline(restart_file,input_string,'#');
   // Ignore rest of line
   restart_file.ignore(80,'\n');
   // Check # of values:
   const unsigned long check_ntvalues=atoi(input_string.c_str());

   // Dynamic run restarted from steady run
   if (check_ntvalues<time_steps_range)
    {
     std::ostringstream warning_stream;
     warning_stream 
      << "Number of time history values in dump file is less "
      << "than the storage allocated in Data object: "
      <<  check_ntvalues << " " << time_steps_range << std::endl;
     warning_stream
      << "We're using steady data as initial data for unsteady \n"
      << "run. I'll fill in the remaining history values with zeroes. \n"
      << "If you don't like this \n"
      << "you'll have to overwrite this yourself with whatever is \n "
      << "appropriate for your timestepping scheme. \n";
     //Issue the warning
     OomphLibWarning(warning_stream.str(),
                     "Data::read()",
                     OOMPH_EXCEPTION_LOCATION);
     
     // Read data
     for(unsigned t=0;t<time_steps_range;t++)
      {
       for(unsigned j=0;j<value_pt_range;j++) 
        {
         if (t==0)
          {
           // Read line
           getline(restart_file,input_string);
           // Transform to double
           set_value(t,j,atof(input_string.c_str()));
          }
         else
          {
           set_value(t,j,0.0);
          }
        }
      }
    }
   // Static run restarted from unsteady run
   else if (check_ntvalues>time_steps_range)
    {
     std::ostringstream warning_stream;
     warning_stream 
      << "Warning: number of time history values in dump file is greater "
      << "than the storage allocated in Data object: " 
      <<  check_ntvalues << " " << time_steps_range << std::endl;
     warning_stream << "We're using the current values from an unsteady \n"
                    << "restart file to initialise a static run. \n";
     //Issue the warning
     OomphLibWarning(warning_stream.str(),
                     "Data::read()",
                     OOMPH_EXCEPTION_LOCATION);

     // Read data
     for(unsigned t=0;t<check_ntvalues;t++)
      {
       for(unsigned j=0;j<value_pt_range;j++) 
        {
         // Read line
         getline(restart_file,input_string);
         if (t==0)
          {
           // Transform to double
           set_value(t,j,atof(input_string.c_str()));
          }
        }
      }
    }
   // Proper dynamic restart
   else
    {
     // Read data
     for(unsigned t=0;t<time_steps_range;t++)
      {
       for(unsigned j=0;j<value_pt_range;j++) 
        {
         // Read line
         getline(restart_file,input_string);
         // Transform to double
         set_value(t,j,atof(input_string.c_str()));
        }
      }
    }
  }
}

//===================================================================
/// Return the total number of doubles stored per value to record
/// the time history of ecah value. The information is read from the
/// time stepper
//===================================================================
unsigned Data::ntstorage() const {return Time_stepper_pt->ntstorage();}

//================================================================
///Assign (global) equation number.
/// This function does NOT initialise the value because
/// if we're using things like node position as variables in the problem
/// they will have been set before the call to assign equation numbers
/// and setting it to zero will wipe it out :(.
///
/// Pass: 
/// - current number of global dofs global_number (which gets incremented)
/// - the Vector of pointers to global dofs (to which new dofs
///   get added)
//================================================================
void Data::assign_eqn_numbers(unsigned long &global_number, 
                              Vector<double *> &dof_pt)
{
 //Loop over the number of variables
 //Set temporary to hold range
 const unsigned eqn_number_range = Nvalue;
 for(unsigned i=0;i<eqn_number_range;i++)
  {
#ifdef OOMPH_HAS_MPI
   // Is the node a halo? If so, treat it as pinned for now
   // This will be overwritten with the actual equation number
   // during the synchronisation phase.
   if (is_halo())
    {
     eqn_number(i) = Is_pinned; 
    }
   else
#endif
    {
     //Boundary conditions test: if it's not a pinned or constrained variable,
     //The assign a new global equation number
     if((!is_pinned(i)) && (!is_constrained(i))) 
      {
       //Assign the equation number and increment global equation number
       Eqn_number[i] = global_number++;
       //Add pointer to global dof vector 
       dof_pt.push_back(value_pt(i));
      }
    }
  }
}

//================================================================
///  Self-test: Have all values been classified as pinned/unpinned?
/// Return 0 if OK. 
//================================================================
unsigned Data::self_test()
{
 //Initialise test flag
 bool passed=true;

 //Loop over all equation numbers
 const unsigned eqn_number_range = Nvalue;
 for(unsigned i=0;i<eqn_number_range;i++)
  {
   //If the equation number has not been assigned, issue an error
   if (Eqn_number[i]==Is_unclassified)
    {
     passed=false;
     oomph_info 
      << "\n ERROR: Failed Data::self_test() for i=" << i << std::endl;
     oomph_info 
      << "          (Value is not classified as pinned or free)" << std::endl;
    }
  }

 //Return verdict
 if (passed) {return 0;}
 else {return 1;}
}

//================================================================
/// Increase the number of data values stored, useful when adding
/// additional data at a node, almost always Lagrange multipliers.
/// Note if any of the unresized data is copied, then we assume all the 
/// resized data is copied from the same node as the unresized data.
//================================================================
void Data::resize(const unsigned &n_value)
{
 //Find current number of values
 const unsigned n_value_old = nvalue();
 //Set the desired number of values
 const unsigned n_value_new = n_value;

 //If the number of values hasn't changed, do nothing
 if(n_value_new==n_value_old) {return;}

 //Put in a little safely check here
#ifdef PARANOID
 if(n_value_new < n_value_old) 
  {
   std::ostringstream error_stream;
   error_stream 
    << "Warning : Data cannot be resized to a smaller value!" << std::endl;
   throw OomphLibError(error_stream.str(),
                       "Data::resize()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 //Find amount of additional time storage required
 //N.B. We can't change timesteppers in this process
 const unsigned t_storage = ntstorage();
 
 //Create new sets of pointers of the appropriate (new) size
 double **value_new_pt = new double*[n_value_new];
 long *eqn_number_new = new long[n_value_new];

 //Create new array of values that is contiguous in memory
 double *values = new double[n_value_new*t_storage];

 //Copy the old values over into the new storage scheme
 for(unsigned i=0;i<n_value_old;i++)
  {
   //Set pointer for the new values
   value_new_pt[i] = &values[i*t_storage];
   //Copy value
   for(unsigned t=0;t<t_storage;t++)
    {value_new_pt[i][t] = Value[i][t];}

   //Copy equation number
   eqn_number_new[i] = Eqn_number[i];
  }

 //Loop over the new entries, set pointers and initialise data
 for(unsigned i=n_value_old;i<n_value_new;i++)
  {
   //Set the pointer
   value_new_pt[i] = &values[i*t_storage];
   //Initialise the new data values to zero
   for(unsigned t=0;t<t_storage;t++) {value_new_pt[i][t] = 0.0;}

   //Initialise the equation number to Is_unclassified
   eqn_number_new[i] = Is_unclassified;
  }

 //Set the number of new values
 Nvalue = n_value_new;

 //Now delete the old storage and set the new pointers
 if (n_value_old!=0) delete[] Value[0];
 delete[] Value;
 Value = value_new_pt;
 delete[] Eqn_number;
 Eqn_number = eqn_number_new;

 //Now update pointers in any copies of this data
 for(unsigned i=0;i<Ncopies;i++)
  {
   Copy_of_data_pt[i]->reset_copied_pointers();
  }
}

//================================================================
/// Reset the pointers to the copied data
//===============================================================
void HijackedData::reset_copied_pointers()
{
 //Copy the pointer to the value. This will give the appropriate
 //"slice" of the array
 Value = &Copied_data_pt->Value[Copied_index];
 
 //Copy the pointer to the equation number
 Eqn_number = &Copied_data_pt->Eqn_number[Copied_index];
}


//===============================================================
/// Clear ther pointers to the copied data
//===============================================================
void HijackedData::clear_copied_pointers()
{
 Copied_data_pt = 0;
 Value = 0; Eqn_number = 0;
}

//================================================================
/// Constructor, creates a HijackedData object with a single value 
/// that is copied from another Data object. 
/// The ordering of the aguments is used to 
/// distinguish this case from that of copying all data values, except one
/// independent value.
//================================================================
HijackedData::
HijackedData(const unsigned &copied_index, Data* const &data_pt) : 
 Data(data_pt->time_stepper_pt(),1,false),
 Copied_data_pt(data_pt),
 Copied_index(copied_index)
{
 //Don't allow copying of a copy
 if(data_pt->is_a_copy(copied_index))
  {
   std::ostringstream error_stream;
   error_stream << "The data you are trying to hijack is already a copy"
                << std::endl;
   error_stream << "Please copy the original data" << std::endl;
   error_stream << "In a later version, I might do this for you,"
                << " but not today" << std::endl;

   throw OomphLibError(error_stream.str(),
                       "HijackedData::HijackedData()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 //Copy the pointer to the value. This will give the appropriate
 //"slice" of the array
 Value = &data_pt->Value[copied_index];
 //Copy the pointer to the equation number
 Eqn_number = &data_pt->Eqn_number[copied_index];
 //Inform the original data that it has been copied
 data_pt->add_copy(this);
}

//=================================================================
/// We do not allow Hijacked Data to be resized
//=================================================================
void HijackedData::resize(const unsigned &n_value)
{
 throw OomphLibError("HijackedData cannot be resized",
                     "HijackedData::resize()",
                     OOMPH_EXCEPTION_LOCATION);
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//Functions for the CopiedData class
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////


//================================================================
/// Reset the pointers to the copied data
//===============================================================
void CopiedData::reset_copied_pointers()
{
 //Copy the pointer to the value. This will give the appropriate
 //"slice" of the array
 Value = Copied_data_pt->Value;
 
 //Copy the pointer to the equation numbers
 Eqn_number = Copied_data_pt->Eqn_number;
}


//===============================================================
/// Clear ther pointers to the copied data
//===============================================================
void CopiedData::clear_copied_pointers()
{
 Copied_data_pt = 0;
 Value = 0; Eqn_number = 0;
}

//================================================================
/// Constructor, creates a CopiedData object with all values 
/// copied from another Data object. 
//================================================================
CopiedData::CopiedData(Data* const &data_pt) : 
 Data(data_pt->time_stepper_pt(),1,false),
 Copied_data_pt(data_pt)
{
 //Don't allow copying of a copy
 if(data_pt->is_a_copy())
  {
   std::ostringstream error_stream;
   error_stream << "The data you are trying to copy is already a copy"
                << std::endl;
   error_stream << "Please copy the original data" << std::endl;
   error_stream << "In a later version, I might do this for you,"
                << " but not today" << std::endl;

   throw OomphLibError(error_stream.str(),
                       "CopiedData::CopiedData()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 //Copy the pointer to the value. 
 Value = data_pt->Value;
 //Copy the pointer to the equation number
 Eqn_number = data_pt->Eqn_number;
 //Inform the original data that it has been copied
 data_pt->add_copy(this);
}

//=================================================================
/// We do not allow Copied Data to be resized
//=================================================================
void CopiedData::resize(const unsigned &n_value)
{
 throw OomphLibError("CopiedData cannot be resized",
                     "CopiedData::resize()",
                     OOMPH_EXCEPTION_LOCATION);
}



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//Functions for the HangInfo class
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//================================================================
/// Check that the argument is within the range of 
/// stored master nodes
//=================================================================
void HangInfo::range_check(const unsigned &i) const
{
 //If the argument is negative or greater than the number of stored
 //values die
 if(i >= Nmaster)
  {
   std::ostringstream error_message;
   error_message << "Range Error: the index " << i
                 << " is not in the range (0,"
                 << Nmaster-1 << ")";
   throw OomphLibError(error_message.str(),
                       "HangInfo::range_check()",
                       OOMPH_EXCEPTION_LOCATION);
  }
}


//=====================================================================
/// Set the pointer to the i-th master node and its weight
//=====================================================================
void HangInfo::set_master_node_pt(const unsigned &i, 
                                  Node* const &master_node_pt,
                                  const double &weight)
{
#ifdef RANGE_CHECKING
 range_check(i);
#endif
 Master_nodes_pt[i] = master_node_pt;
 Master_weights[i] = weight;
}

//====================================================================
/// Add (pointer to) master node and corresponding weight to 
/// the internally stored  (pointers to) master nodes and weights.
//====================================================================
 void HangInfo::add_master_node_pt(Node* const &master_node_pt,
                                   const double &weight) 
{
 //Find the present number of master nodes
 const unsigned n_master = Nmaster;
 //Make new data
 Node* *new_master_nodes_pt = new Node*[n_master+1];
 double *new_master_weights = new double[n_master+1];
 
 //Copy the old values over to the new data
 for(unsigned i=0;i<n_master;i++)
  {
   new_master_nodes_pt[i] = Master_nodes_pt[i];
   new_master_weights[i] = Master_weights[i];
  }
 //Add the new values at the end
 new_master_nodes_pt[n_master] = master_node_pt;
 new_master_weights[n_master] = weight;
   
 //Reset the pointers
 delete[] Master_nodes_pt; Master_nodes_pt = new_master_nodes_pt;
 delete[] Master_weights; Master_weights = new_master_weights;
 //Increase the number of master nodes
 ++Nmaster;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//Functions for the Node class
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
 
//=================================================================
/// \short Private function to check that the arguments are within
/// the range of the stored coordinates, position types and time history 
/// values.
//=================================================================
void Node::x_gen_range_check(const unsigned &t, const unsigned &k, 
                             const unsigned &i) const
{
 //Number of stored history values
 const unsigned position_ntstorage = Position_time_stepper_pt->ntstorage();
 //If any of the coordinates or time values are out of range
 if((i >= Ndim) || (k >= Nposition_type) || 
    (t >= position_ntstorage))
  {
   std::ostringstream error_message;
   //If it's the dimension
   if(i >= Ndim)
    {
     error_message << "Range Error: X coordinate " << i
                   << " is not in the range (0,"
                   << Ndim-1 << ")";
    }
   //If it's the position type
   if(k >= Nposition_type)
    {
     error_message << "Range Error: Position type " << k
                   << " is not in the range (0,"
                   << Nposition_type-1 << ")";
    }
   //If it's the time
   if(t >= position_ntstorage)
    {
     error_message << "Range Error: Position Time Value " << t
                   << " is not in the range (0,"
                   << position_ntstorage - 1 << ")";
    }
   //Throw the error
   throw OomphLibError(error_message.str(),
                       "Node::x_gen_range_check()",
                       OOMPH_EXCEPTION_LOCATION);
  }
}

//========================================================================
/// Static "Magic number" passed as independent_position when there is
/// no independent position in the periodic node. For example, in a periodic
/// mesh.
//=======================================================================
unsigned Node::No_independent_position=10;

//========================================================================
/// Default constructor.
//========================================================================
Node::Node() : Data(), 
               Position_time_stepper_pt(Data::Default_static_time_stepper_pt),
               Hanging_pt(0),
               Ndim(0), Nposition_type(0), Aux_node_update_fct_pt(0)
{
#ifdef LEAK_CHECK
 LeakCheckNames::Node_build+=1;
#endif
}

//========================================================================
/// Steady Constructor, allocates storage for initial_n_value values 
/// at a node of spatial dimension NDim. nposition_type: # of  coordinate 
/// types needed in the mapping between local and global coordinates 
/// (e.g. 1 for Lagrange-type elements; 2 for 1D Hermite elements; 4 for
/// 2D Hermite elements, etc).
//========================================================================
Node::Node(const unsigned &n_dim, 
           const unsigned &n_position_type,
           const unsigned &initial_n_value,
           const bool &allocate_x_position) : 
 Data(initial_n_value), 
 X_position(0),
 Position_time_stepper_pt(Data::Default_static_time_stepper_pt),
 Hanging_pt(0),
 Ndim(n_dim), 
 Nposition_type(n_position_type), Aux_node_update_fct_pt(0)
{
#ifdef LEAK_CHECK
 LeakCheckNames::Node_build+=1;
#endif

 //Determine the total amount of storage required for position variables
 const unsigned n_storage = n_dim*n_position_type;

 //If we are in charge of the x coordinates (non-solid node)
 //the allocate storage
 if(allocate_x_position)
  {
   //Allocate the pointers to each coordinate and coordinate type
   X_position = new double*[n_storage];

   //Create one big array of positions
   double *x_positions = new double[n_storage];

   //Set pointers from the contiguous array
   for(unsigned j=0;j<n_storage;j++) 
    {
     X_position[j] = &x_positions[j];
     //Initialise value to zero
     X_position[j][0] = 0.0;
    }
  }
 
}

//========================================================================
/// Unsteady Constructor for a node of spatial dimension n_dim. 
/// Allocates storage
/// for initial_n_value values with history values as required
/// by timestepper. n_position_type: # of coordinate 
/// types needed in the mapping between local and global coordinates  
/// (e.g. 1 for Lagrange-type elements; 2 for 1D Hermite elements; 4 for 
/// 2D Hermite elements)
//========================================================================
Node::Node(TimeStepper* const &time_stepper_pt, 
           const unsigned &n_dim, 
           const unsigned &n_position_type, 
           const unsigned &initial_n_value,
           const bool &allocate_x_position) 
 : Data(time_stepper_pt,initial_n_value),
   X_position(0),
   Position_time_stepper_pt(time_stepper_pt), 
   Hanging_pt(0),
   Ndim(n_dim), Nposition_type(n_position_type), Aux_node_update_fct_pt(0)
{
#ifdef LEAK_CHECK
 LeakCheckNames::Node_build+=1;
#endif

 //Determine the total amount of storage required for position variables
 const unsigned n_storage = n_dim*n_position_type;
 
 //If we are allocating the storage (non-solid node)
 if(allocate_x_position)
  {
   //Amount of storage required for history values
   const unsigned n_tstorage = Position_time_stepper_pt->ntstorage();
   
   //Allocate the pointers to each coordinate and coordinate type
   X_position = new double*[n_storage];

   //Allocate the positions in one big array
   double *x_positions = new double[n_storage*n_tstorage];

   //Set the pointers to the contiguous memory
   for(unsigned j=0;j<n_storage;j++) 
    {
     //Set the pointer from the bug array
     X_position[j] = &x_positions[j*n_tstorage];
     //Initialise all values to zero
     for(unsigned t=0;t<n_tstorage;t++) {X_position[j][t] = 0.0;}
    }
  }
}


//========================================================================
/// Destructor to clean up the memory allocated for nodal position
//========================================================================
Node::~Node()
{
#ifdef LEAK_CHECK
 LeakCheckNames::Node_build-=1;
#endif

 //Clean up memory allocated to hanging nodes
 if(Hanging_pt!=0)
  {
   //The number of hanging pointers is the number of values plus one
   const unsigned nhang = nvalue() + 1;
   for(unsigned ival=1;ival<nhang;ival++)
    {
     //If the ival-th HangInfo object is not the same as the geometrical 
     //one, delete it
     if(Hanging_pt[ival]!=Hanging_pt[0]) {delete Hanging_pt[ival];}
     //Always NULL out the HangInfo pointer
     Hanging_pt[ival]=0;
    }

   //Delete the Geometrical HangInfo pointer
   delete Hanging_pt[0];
   Hanging_pt[0]=0;
 
   //Delete the Hanging_pt
   delete[] Hanging_pt;
   Hanging_pt=0;
  }

 //Free the memory allocated
 
 //If we did not allocate then the memory must have been freed by the
 //destructor of the object that did the allocating and X_position MUST
 //have been set back to zero
 //Test this and if so, we're done
 if(X_position==0) {return;}

 //If we're still here we must free our own memory which was allocated
 //in one block
 delete[] X_position[0];

 //Now delete the pointer
 delete[] X_position; X_position=0;
}


//================================================================
///  Return the i-th component of nodal velocity: dx/dt
//================================================================
double Node::dx_dt(const unsigned &i) const
{
 //Number of timsteps (past & present)
 const unsigned n_time = Position_time_stepper_pt->ntstorage();

 double dxdt=0.0;
 
 //If the timestepper is not steady
 if (!Position_time_stepper_pt->is_steady())
  {
   //Loop over the additional storage and add the appropriate contributions
   for(unsigned t=0;t<n_time;t++)
    {
     dxdt+=Position_time_stepper_pt->weight(1,t)*x(t,i);
    }
  }
 
 return dxdt;
}

//================================================================
/// Return the i-th component of j-th derivative of nodal position: 
/// d^jx/dt^j.
//================================================================
double Node::dx_dt(const unsigned &j, const unsigned &i) const
{
 // Number of timsteps (past & present)
 const unsigned n_time = Position_time_stepper_pt->ntstorage();
 
 double dxdt=0.0;
 
 //If the timestepper is not steady
 if ((!Position_time_stepper_pt->is_steady()) || (j==0))
  {
   //Loop over the additional storage and add the appropriate contributions 
   for(unsigned t=0;t<n_time;t++)
    {
     dxdt+=Position_time_stepper_pt->weight(j,t)*x(t,i);
    }
  }
 
 return dxdt;
}

//================================================================
/// \short  i-th component of time derivative (velocity) of the 
/// generalised position, dx(k,i)/dt. `Type': k; Coordinate direction: i.
//================================================================
double Node::dx_gen_dt(const unsigned &k, const unsigned &i) const
{
 // Number of timsteps (past & present)
 const unsigned n_time = Position_time_stepper_pt->ntstorage();
 
 double dxdt=0.0;
 
 //If the timestepper is not steady
 if (!Position_time_stepper_pt->is_steady())
  {
   //Loop over the additional time storage and add the appropriate 
   //contributions
   for(unsigned t=0;t<n_time;t++)
    {
     dxdt+=Position_time_stepper_pt->weight(1,t)*x_gen(t,k,i);
    }
  }
 
 return dxdt;
}
 
//================================================================
/// \short  i-th component of j-th time derivative (velocity) of the 
/// generalised position, d^jx(k,i)/dt^j. `Type': k; Coordinate direction: i.
//================================================================
double Node::dx_gen_dt(const unsigned &j, const unsigned &k, 
                        const unsigned &i) const
{
 // Number of timsteps (past & present)
 const unsigned n_time = Position_time_stepper_pt->ntstorage();
 
 double dxdt=0.0;
 
 //If the timestepper is not steady
 if ((!Position_time_stepper_pt->is_steady()) || (j==0))
  {
   //Loop over the additional storage and add the appropriate contributions
   for(unsigned t=0;t<n_time;t++)
    {
     dxdt+=Position_time_stepper_pt->weight(j,t)*x_gen(t,k,i);
    }
  }
 
 return dxdt;
}



//================================================================
/// Copy all nodal data from specified Node object
//================================================================
void Node::copy(Node* orig_node_pt)
{

 // Number of positional values
 const unsigned npos_storage = Ndim*Nposition_type;

 // Check # of values:
 const unsigned long npos_storage_orig
  = orig_node_pt->ndim()*orig_node_pt->nposition_type();
 if (npos_storage!=npos_storage_orig)
  {
   std::ostringstream error_stream;
   error_stream << "The allocated positional storage " 
                << npos_storage << " is not the same as the original Node "
                << npos_storage_orig << std::endl;

   throw OomphLibError(error_stream.str(),
                       "Node::copy()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 // Number of time values (incl present)
 const unsigned n_time = Position_time_stepper_pt->ntstorage();

 // Check # of values:
 const unsigned long n_time_orig =
  orig_node_pt->position_time_stepper_pt()->ntstorage();
 if (n_time!=n_time_orig)
    {
     std::ostringstream error_stream;
     error_stream << "The number of positional time history values, " 
                  << n_time 
                  << " is not the same of those in the original node "
                  << n_time_orig << std::endl;
     
     throw OomphLibError(error_stream.str(),
                         "Node::copy()",
                         OOMPH_EXCEPTION_LOCATION);
    }

 // Copy fixed nodal positions
 for(unsigned t=0;t<n_time;t++)
  {
   for(unsigned j=0;j<npos_storage;j++) 
    {
     X_position[j][t] = orig_node_pt->X_position[j][t];
    }
  }

 //  Read associated data
 Data::copy(orig_node_pt);

}


//================================================================
///Dump nodal positions and associated data to file for restart
//================================================================
void Node::dump(std::ostream& dump_file)
{
 // Number of positional values
 const unsigned npos_storage = Ndim*Nposition_type;
 dump_file << npos_storage 
           << " # number of fixed position variables" << std::endl;

 const unsigned Time_steps_range = Position_time_stepper_pt->ntstorage();
 dump_file << Time_steps_range 
           << " # total number of doubles for time history (incl present)" 
           << std::endl;

 for(unsigned t=0;t<Time_steps_range;t++)
  {
   for(unsigned j=0;j<npos_storage;j++) 
    {
     dump_file << X_position[j][t] << std::endl;
    }
  }

 // Dump out data
 Data::dump(dump_file);
}

//================================================================
///Read nodal positions and associated data from file for restart
//================================================================
void Node::read(std::ifstream& restart_file)
{

 std::string input_string;

 // Number of positional values
 const unsigned npos_storage = Ndim*Nposition_type;

 // Read line up to termination sign
 getline(restart_file,input_string,'#');
 // Ignore rest of line
 restart_file.ignore(80,'\n');
 // Check # of values:
 const unsigned long check_npos_storage=atoi(input_string.c_str());
 if (check_npos_storage!=npos_storage)
    {
     std::ostringstream error_stream;
     error_stream << "The allocated positional storage " 
                  << npos_storage << 
      " is not the same as that in the input file"
                  << check_npos_storage << std::endl;
     
     throw OomphLibError(error_stream.str(),
                         "Node::read()",
                         OOMPH_EXCEPTION_LOCATION);
   
    }
 
 // Number of time values (incl present)
 const unsigned time_steps_range = Position_time_stepper_pt->ntstorage();

 // Read line up to termination sign
 getline(restart_file,input_string,'#');
 // Ignore rest of line
 restart_file.ignore(80,'\n');
 // Check # of values:
 const unsigned long check_time_steps_range=atoi(input_string.c_str());
 if (check_time_steps_range!=time_steps_range)
  {
   std::ostringstream error_stream;
   error_stream
    << "Number of positional history values in dump file is less "
    << "than the storage allocated in Node object: "
    <<  check_time_steps_range 
    << " " << time_steps_range << std::endl;

   throw OomphLibError(error_stream.str(),
                       "Node::read()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 // Read fixed nodal positions
 for(unsigned t=0;t<time_steps_range;t++)
  {
   for(unsigned j=0;j<npos_storage;j++) 
    {
     // Read line
     getline(restart_file,input_string);
     
     // Transform to double
     X_position[j][t] = atof(input_string.c_str());
    }
  }

 //  Read associated data
 Data::read(restart_file);
}

//=====================================================================
/// Set the hanging data for the i-th  value. 
/// If node is already hanging, simply overwrite the appropriate entry.
/// If the node isn't hanging (because it might not be hanging
/// geometrically), create the Vector of hanging pointers
/// and make the other entries point to the node's geometric
/// hanging data.
//=====================================================================
void Node::set_hanging_pt(HangInfo* const &hang_pt, const int &i)
{
 //The number of hanging values is the number of stored values plus
 //one (geometry)
 unsigned n_hang = nvalue() + 1;

 //Has the vector of pointers to the HangInfo objects already been created?
 //If not create it
 if(Hanging_pt==0) 
  {
   Hanging_pt = new HangInfo*[n_hang];
   //Initialise all entries to zero
   for(unsigned i=0;i<n_hang;i++) {Hanging_pt[i] = 0;}
  }
 
 //Geometric hanging data
 if(i==-1) 
  {
   //Setup boolean array to find which pointers match the geometric pointer
   std::vector<bool> Same_as_geometric(n_hang,true);
   
   //Mark up any values that DON'T use the geometric hanging scheme
   for(unsigned n=1;n<n_hang;n++)
    {if(Hanging_pt[n] != Hanging_pt[0]) {Same_as_geometric[n] = false;}}
   
   //Remove the old geometric HangInfo
   delete Hanging_pt[0];
   //Assign the new geometric hanging data 
   Hanging_pt[0] = hang_pt;
   
   //Constrain the geometric data (virtual function that is
   //overladed in solid nodes)
   constrain_positions();

   //Loop over the entries again and update all pointers that pointed to
   //the geometric data
   for(unsigned n=1;n<n_hang;n++)
    {
     if(Same_as_geometric[n]==true) 
      {
       Hanging_pt[n] = Hanging_pt[0];
       //In addition set the corresponding value to be constrained (hanging)
       constrain(n-1);
      }
    }

  }
 //Value data
 else
  {
   //If the data is different from geometric, delete it
   if(Hanging_pt[i+1] != Hanging_pt[0]) 
    {
     delete Hanging_pt[i+1];
     Hanging_pt[i+1] = 0;
    }

   //Overwrite hanging data for the required value
   //Do not need to delete previous value, because it is assigned outside
   //the Node class
   Hanging_pt[i+1]=hang_pt;

   //In addition set the value to be constrained (hanging)
   constrain(i);
  }
}

//=====================================================================
/// Make the node periodic by copying values from node_pt. 
/// Broken virtual (only implemented in BoundaryNodes)
//===================================================================== 
void Node::make_periodic(Node* const &node_pt)
{
 throw OomphLibError("Only BoundaryNodes can be made periodic",
                     "Node::make_periodic()",
                     OOMPH_EXCEPTION_LOCATION);
}

//=====================================================================
/// Make the nodes passed in periodic_nodes_pt periodic by copying values
/// across from this node. At present all the positions will be assumed
/// to be independent.
/// Broken virtual (only implemented in BoundaryNodes)
//===================================================================== 
void Node::make_periodic_nodes(const Vector<Node*> &periodic_nodes_pt)
{
 throw OomphLibError("Only BoundaryNodes can make periodic nodes",
                     "Node::make_periodic_nodes()",
                     OOMPH_EXCEPTION_LOCATION);
}


//====================================================================
/// Label node as non-hanging node by removing all hanging node data.
//====================================================================
void Node::set_nonhanging()
{
 if(Hanging_pt!=0)
  {
   //Kill any additional hanging data for values
   const unsigned nhang = nvalue() + 1;
   for(unsigned ival=1;ival<nhang;ival++)
    {
     // Only kill it if it's different from the geometric hanging node data
     if (Hanging_pt[ival]!= Hanging_pt[0]) {delete Hanging_pt[ival];}
     //Always zero the entry
     Hanging_pt[ival] = 0;

     //Unconstrain any values that were constrained only because they
     //were hanging
     unconstrain(ival-1);
    }
   
   //Unconstrain the positions (virtual function that is overloaded for
   //solid nodes)
   unconstrain_positions();
   
   //Kill the geometric hanging node data
   delete Hanging_pt[0];
   Hanging_pt[0]=0;

   //Kill the pointer to all hanging data
   delete[]  Hanging_pt;
   Hanging_pt=0;
  }
}



//=======================================================================
/// Interface for function to add the node to the mesh boundary b.
/// Broken here in order to report run-time errors. Must be overloaded
/// by all boundary nodes
//=======================================================================
void Node::add_to_boundary(const unsigned &b)
{
 throw OomphLibError("Cannot add non BoundaryNode<NODE> to boundary",
                     "Node::add_to_boundary()",
                     OOMPH_EXCEPTION_LOCATION);
}


//=======================================================================
/// Interface for function to remove the node from the mesh boundary b.
/// Broken here in order to report run-time erorrs. Must be overloaded
/// by all boundary nodes
//=======================================================================
void Node::remove_from_boundary(const unsigned &b)
{
 throw OomphLibError("Cannot remove non BoundaryNode<NODE> to boundary",
                     "Node::remove_from_boundary()",
                     OOMPH_EXCEPTION_LOCATION);
}


//=========================================================================
///  Interface to get the number of boundary coordinates on mesh boundary b. 
/// Broken here in order to provide run-time error reporting. Must 
/// be overloaded by all boundary nodes.
//=========================================================================
unsigned Node::ncoordinates_on_boundary(const unsigned &b)
{
 throw OomphLibError("Non-boundary Node cannot have boundary coordinates",
                     "Node::ncoordinates_on_boundary()",
                     OOMPH_EXCEPTION_LOCATION);
 // dummy return
 return 0;
}


//=========================================================================
/// Interface for function to get the k-th generalised boundary coordinate
/// of the node on boundary b. Broken here in order to 
/// provide run-time error reporting. Must be overloaded by all boundary
/// nodes.
//=========================================================================
void Node::get_coordinates_on_boundary(const unsigned &b, const unsigned& k,
                                       Vector<double> &boundary_zeta)
{
 throw OomphLibError("Non-boundary Node cannot have boundary coordinates",
                     "Node::get_coordinate_on_boundary()",
                     OOMPH_EXCEPTION_LOCATION);
}


//=========================================================================
/// Interface for function to set the k-th generalised boundary coordinate
///  of the node on boundary b. Broken here to provide 
/// run-time error reports. Must be overloaded by all boundary nodes.
//=========================================================================
void Node::set_coordinates_on_boundary(const unsigned &b, const unsigned& k,
                                       const Vector<double> &boundary_zeta)
{
 throw OomphLibError("Non-boundary Node cannot have boundary coordinates",
                     "Node::set_coordinate_on_boundary()",
                     OOMPH_EXCEPTION_LOCATION);
}


//=================================================================
/// Return i-th value (free or pinned) at this node
/// either directly or via hanging node representation.
//================================================================
double Node::value(const unsigned &i) const
{
 //If value is not hanging, just return the underlying value
 if(!is_hanging(i)) {return raw_value(i);}
 // Hanging node: Use hanging node representation
 else
  {
   // Initialise
   double sum=0.0;
   // Add contribution from master nodes
   const unsigned n_master = hanging_pt(i)->nmaster();
   for(unsigned m=0;m<n_master;m++)
    {
     //A master node cannot be hanging by definition.
     //so we get the raw value to avoid an unnecessary it
     sum += hanging_pt(i)->master_node_pt(m)->raw_value(i)*
      hanging_pt(i)->master_weight(m);
    }
   return sum;
  }
}

//=================================================================
/// Return i-th value (free or pinned) at this node at time level t
/// either directly or via hanging node representation.
//================================================================
double Node::value(const unsigned &t, const unsigned &i) const
{
 //If value is not hanging, just return the raw value
 if(!is_hanging(i)) {return raw_value(t,i);}
 // Hanging node: Use hanging node representation
 else
  {
   // Initialise
   double sum=0.0;

   // Add contribution from master nodes
   const unsigned n_master=hanging_pt(i)->nmaster();
   for(unsigned m=0;m<n_master;m++)
    {
     //Get the raw nodal values at each master to avoid un-necessary ifs
     sum += hanging_pt(i)->master_node_pt(m)->raw_value(t,i)*
      hanging_pt(i)->master_weight(m);
    }
   return sum;
  }
}

//==================================================================
/// Compute Vector of values (dofs or pinned) at this Data object
/// either directly or via hanging node representation.
//==================================================================
void Node::value(Vector<double>& values) const
{
 //Loop over all the values
 const unsigned n_value = nvalue();
 for(unsigned i=0;i<n_value;i++)
  {
   //Set the value, using the hanging node representation if necessary
   values[i] = value(i);
  }
}

//==================================================================
/// Compute Vector of values (dofs or pinned) at this node
/// at time level t (t=0: present; t>0: previous)
/// either directly or via hanging node representation.
//==================================================================
void Node::value(const unsigned& t, Vector<double>& values) const
{
 //Loop over all the values
 const unsigned n_value = nvalue();
 for(unsigned i=0;i<n_value;i++)
  {
   //Set the value at the time-level t, using the hanging node representation
   //if necessary
   values[i] = value(t,i);
  }
}


//============================================================
/// Compute Vector of nodal positions
/// either directly or via hanging node representation
//===========================================================
void Node::position(Vector<double>& pos) const
{
 //Assign all positions using hanging node representation where necessary
 const unsigned n_dim = ndim();
 for(unsigned i=0;i<n_dim;i++) {pos[i] = position(i);}
}

//===============================================================
/// Compute Vector of nodal position at timestep t
/// (t=0: current; t>0: previous timestep),
/// either directly or via hanging node representation.
//==============================================================
void Node::position(const unsigned &t, Vector<double>& pos) const
{
 //Assign all positions, using hanging node representation where necessary
 const unsigned n_dim = ndim();
 for(unsigned i=0;i<n_dim;i++) {pos[i] = position(t,i);}
}

//=======================================================================
/// Return i-th nodal coordinate
/// either directly or via hanging node representation.
//======================================================================
double Node::position(const unsigned &i) const
{
 double posn=0.0;

 // Non-hanging node: just return value
 if (!is_hanging()) {posn = x(i);}
 // Hanging node: Use hanging node representation
 else
  {
   // Initialise
   double position=0.0;
   
   // Add contribution from master nodes
   const unsigned n_master=hanging_pt()->nmaster();
   for (unsigned m=0;m<n_master;m++)
    {
     position += hanging_pt()->master_node_pt(m)->x(i)*
      hanging_pt()->master_weight(m);
    } 
   posn=position;
  }
 return posn;
}

//================================================================
/// Return i-th nodal coordinate at time level t
/// (t=0: current; t>0: previous time level),
/// either directly or via hanging node representation.
//================================================================
double Node::position(const unsigned &t, const unsigned &i) const
{
 double posn=0.0;
 
 // Non-hanging node: just return value
 if(!is_hanging()) {posn = x(t,i);}
 // Hanging node: Use hanging node representation
 else
  {
   // Initialise
   double position=0.0;
   
   // Add contribution from master nodes
   const unsigned n_master=hanging_pt()->nmaster();
   for (unsigned m=0;m<n_master;m++)
    {
     position+=
      hanging_pt()->master_node_pt(m)->x(t,i)*
      hanging_pt()->master_weight(m);
    }  
   posn=position;
  }
 
 return posn;
}

//=======================================================================
/// Return generalised nodal coordinate
/// either directly or via hanging node representation.
//======================================================================
double Node::position_gen(const unsigned &k, const unsigned &i) const
{
 double posn=0.0;

 // Non-hanging node: just return value
 if(!is_hanging()) {posn=x_gen(k,i);}
 // Hanging node: Use hanging node representation
 else
  {
   // Initialise
   double position=0.0;
   
   // Add contribution from master nodes
   const unsigned n_master=hanging_pt()->nmaster();
   for (unsigned m=0;m<n_master;m++)
    {
     position+=
      hanging_pt()->master_node_pt(m)->x_gen(k,i)*
      hanging_pt()->master_weight(m);
    } 
   posn=position;
  }
 return posn;
}

//================================================================
/// Return generalised nodal coordinate at time level t
/// (t=0: current; t>0: previous time level),
/// either directly or via hanging node representation.
//================================================================
double Node::position_gen(const unsigned &t, const unsigned &k,
                          const unsigned &i) const
{
 double posn=0.0;
 
 // Non-hanging node: just return value
 if(!is_hanging()) {posn=x_gen(t,k,i);}
 // Hanging node: Use hanging node representation
 else
  {
   // Initialise
   double position=0.0;
   
   // Add contribution from master nodes
   const unsigned n_master=hanging_pt()->nmaster();
   for (unsigned m=0;m<n_master;m++)
    {
     position+=
      hanging_pt()->master_node_pt(m)->x_gen(t,k,i)*
      hanging_pt()->master_weight(m);
    }  
   posn=position;
  }
 
 return posn;
}

//================================================================
///  Return the i-th component of nodal velocity: dx/dt
//// Use the hanging node representation if required.
//================================================================
double Node::dposition_dt(const unsigned &i) const
{
 //Number of timsteps (past & present)
 const unsigned n_time = Position_time_stepper_pt->ntstorage();

 double dxdt=0.0;
 
 //If the timestepper is not steady
 if (!Position_time_stepper_pt->is_steady())
  {
   //Loop over the additional storage and add the appropriate contributions
   for(unsigned t=0;t<n_time;t++)
    {
     dxdt+=Position_time_stepper_pt->weight(1,t)*position(t,i);
    }
  }
 
 return dxdt;
}

//================================================================
/// Return the i-th component of j-th derivative of nodal position: 
/// d^jx/dt^j. Use the hanging node representation.
//================================================================
double Node::dposition_dt(const unsigned &j, const unsigned &i) const
{
 // Number of timsteps (past & present)
 const unsigned n_time = Position_time_stepper_pt->ntstorage();
 
 double dxdt=0.0;
 
 //If the timestepper is not steady
 if ((!Position_time_stepper_pt->is_steady()) || (j==0))
  {
   //Loop over the additional storage and add the appropriate contributions 
   for(unsigned t=0;t<n_time;t++)
    {
     dxdt+=Position_time_stepper_pt->weight(j,t)*position(t,i);
    }
  }
 
 return dxdt;
}

//================================================================
/// \short  i-th component of time derivative (velocity) of the 
/// generalised position, dx(k,i)/dt. `Type': k; Coordinate direction: i.
/// Use the hanging node representation
//================================================================
double Node::dposition_gen_dt(const unsigned &k, const unsigned &i) const
{
 // Number of timsteps (past & present)
 const unsigned n_time = Position_time_stepper_pt->ntstorage();
 
 double dxdt=0.0;

 //If the timestepper is not steady
 if (!Position_time_stepper_pt->is_steady())
  {
   //Loop over the additional time storage and add the appropriate 
   //contributions
   for(unsigned t=0;t<n_time;t++)
    {
     dxdt+=Position_time_stepper_pt->weight(1,t)*position_gen(t,k,i);
    }
  }
 
 return dxdt;
}
 
//================================================================
/// \short  i-th component of j-th time derivative (velocity) of the 
/// generalised position, d^jx(k,i)/dt^j. `Type': k; Coordinate direction: i.
/// Use the hanging node representation.
//================================================================
double  Node::dposition_gen_dt(const unsigned &j, const unsigned &k, 
                               const unsigned &i) const
{
 // Number of timsteps (past & present)
 const unsigned n_time = Position_time_stepper_pt->ntstorage();
 
 double dxdt=0.0;
 
 //If the timestepper is not steady
 if ((!Position_time_stepper_pt->is_steady()) || (j==0))
  {
   //Loop over the additional storage and add the appropriate contributions
   for(unsigned t=0;t<n_time;t++)
    {
     dxdt+=Position_time_stepper_pt->weight(j,t)*position_gen(t,k,i);
    }
  }
 
 return dxdt;
}


//========================================================================
/// Output nodal coordinates
//========================================================================
void Node::output(std::ostream &outfile)
{
 //Loop over the number of dimensions of the node
 const unsigned ndim = this->ndim();
 for (unsigned i=0;i<ndim;i++) {outfile << x(i) << " ";}  
 outfile << std::endl;
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
//Functions for the BoundaryNodeBase class
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

//==================================================================
/// \short Helper function that is used to turn BoundaryNodes into
/// periodic boundary nodes by setting the data values of the nodes
/// in the vector periodic_copies_pt to be the same as those
/// in copied_node_pt. This function should be used when making doubly periodic
/// sets of nodes.
//==================================================================
void BoundaryNodeBase::make_nodes_periodic(
 Node* const &copied_node_pt,
 Vector<Node*> const &periodic_copies_pt)
{
 //Don't allow copying if the original or periodic nodes are already
 //periodic
 bool already_a_copy = false;
 already_a_copy |= copied_node_pt->is_a_copy();
 const unsigned n_periodic = periodic_copies_pt.size();
 for(unsigned n=0;n<n_periodic;n++)
  {
   already_a_copy |= periodic_copies_pt[n]->is_a_copy();
  }

 //If we have a copy bail
 if(already_a_copy)
  {
   std::ostringstream error_stream;
   error_stream << 
    "The nodes you are trying to make periodic are already periodic\n"
                << 
    "Or you are trying to make a copy of another already periodic node\n";
   error_stream << "Please copy the original data if you can\n";
   throw OomphLibError(error_stream.str(),
                       "BoundaryNodeBase::make_nodes_periodic()",
                       OOMPH_EXCEPTION_LOCATION);

  }

 //Now we simply delete and copy over for each node
 for(unsigned n=0;n<n_periodic;n++)
  {
   //Local cache of the node
   Node* const nod_pt = periodic_copies_pt[n];
   //Miss out the node itself if it's in the list
   if(nod_pt != copied_node_pt)
    {
     //Delete the storage allocated in the copy
     nod_pt->delete_value_storage();
     //Now set the Value and Equation number pointers to be the same
     nod_pt->Value = copied_node_pt->Value;
     nod_pt->Eqn_number = copied_node_pt->Eqn_number;
     
     //Set the copied node pointer in the copy
     BoundaryNodeBase* cast_nod_pt = dynamic_cast<BoundaryNodeBase*>(nod_pt);
     cast_nod_pt->Copied_node_pt = copied_node_pt;
     //Inform the node that it has been copied
     copied_node_pt->add_copy(nod_pt);
    }
  }

}


//====================================================================
/// Helper function that is used to turn BoundaryNodes into
/// peridic boundary nodes by setting the data values of
/// copy_of_node_pt to those of copied_node_pt. 
//=====================================================================
void BoundaryNodeBase::make_node_periodic(Node* const &node_pt, 
                                          Node* const &copied_node_pt)
{
 //Don't allow the copying of a copy or an orignal node
 if((node_pt->is_a_copy()) || (copied_node_pt->is_a_copy()))
  {
   std::ostringstream error_stream;
   error_stream << 
    "The node you are trying to make periodic is already periodic\n"
                << 
    "Or you are trying to make a copy of another already periodic node\n";
   error_stream << "Please copy the original data if you can\n";
   error_stream << 
    "If you wish to make doubly periodic nodes, e.g. in the corners of a cell\n";
   error_stream << "You should use the function \n\n"
                << "Node::make_periodic_nodes()\n";
   throw OomphLibError(error_stream.str(),
                       "BoundaryNodeBase::make_node_periodic()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 //Set the copied node pointer
 Copied_node_pt = copied_node_pt;

 //First copy the data values
 //Delete the storage allocated in the copy
 node_pt->delete_value_storage();
 //Now set the Value and Equation number pointers to be the same
 node_pt->Value = copied_node_pt->Value;
 node_pt->Eqn_number = copied_node_pt->Eqn_number;

 //Inform the node that it has been copied
 copied_node_pt->add_copy(node_pt);
}

//======================================================================
/// Destructor to clean up any memory that might have been allocated.
//=======================================================================
BoundaryNodeBase::~BoundaryNodeBase()
{
 //Delete the set of boundaries on which the Node lies
 delete Boundaries_pt;
 Boundaries_pt=0;

 //If the Boundary coordinates have been set then delete them
 if(Boundary_coordinates_pt != 0)
  {
   //Loop over the boundary coordinate entries and delete them
   for(std::map<unsigned,DenseMatrix<double> *>::iterator it 
        = Boundary_coordinates_pt->begin();
       it!=Boundary_coordinates_pt->end();++it)
    {
     //Delete the vectors that have been allocated for the storage
     //of the boundary coordinates
     delete it->second;
    }
   
   //Now delete the Boundary coordinates map itself
   delete Boundary_coordinates_pt;
   //Set the pointer to null to be on the safe side
   Boundary_coordinates_pt=0;
  }

 //Delete the map of face element's first value 
 delete Index_of_first_value_assigned_by_face_element_pt;
 Index_of_first_value_assigned_by_face_element_pt=0;
}


//=======================================================================
/// Add the node to the mesh boundary b
//=======================================================================
void BoundaryNodeBase::add_to_boundary(const unsigned &b)
{
 //If there is not storage then create storage and set the entry
 //to be the boundary b
 if(Boundaries_pt==0) {Boundaries_pt = new std::set<unsigned>;}
 
//  //If the boundary is already stored in the node issue a warning
//  if(find(Boundaries_pt->begin(),Boundaries_pt->end(),b) != 
//     Boundaries_pt->end())
//   {
// // MH: who cares?
// //    oomph_info << std::endl << "============================================"
// //         << std::endl << "Warning in Node::add_to_boundary()          "
// //         << std::endl << "Node is already marked as being on boundary " << b  
// //         << std::endl << "============================================"
// //         << std::endl << std::endl;
//   }
//  else
//   {

   Boundaries_pt->insert(b);

//   }
}


//=======================================================================
/// Remove the node from the mesh boundary b
//=======================================================================
void BoundaryNodeBase::remove_from_boundary(const unsigned &b)
{
#ifdef PARANOID
 if(is_on_boundary(b)==false)
  {
   std::ostringstream error_stream;
   error_stream << "Node is not on boundary " << b << std::endl;

   throw OomphLibError(error_stream.str(),
                       "Node::remove_from_boundary()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 //Remove the boundary from the set 
 Boundaries_pt->erase(b);

 //Need to delete the equivalent entry in the Boundary coordinate
 //map, if the storage has actually been allocated
 if(Boundary_coordinates_pt!=0)
  {
   //Delete the vector storage that has been allocated
   delete (*Boundary_coordinates_pt)[b];
   //Delete the entry in the map
   Boundary_coordinates_pt->erase(b);
  }

 //If all entries have been removed, delete the storage
 if(Boundaries_pt->size()==0) 
  {
   delete Boundaries_pt;
   Boundaries_pt=0;
  }
}

//========================================================================
/// Test whether the node lies on the mesh boundary b
//========================================================================
bool BoundaryNodeBase::is_on_boundary(const unsigned &b)
{
 //If the node lies on any boundary
 if(Boundaries_pt!=0)
  {
   if(find(Boundaries_pt->begin(),Boundaries_pt->end(),b) 
      != Boundaries_pt->end()) {return true;}
  }

 //If we haven't returned yet, then the node does not lie on the boundary
 return false;
}


//=========================================================================
/// Get the number of boundary coordinates on mesh boundary b
//=========================================================================
unsigned BoundaryNodeBase::ncoordinates_on_boundary(const unsigned &b)
{
 //Check that the node lies on a boundary
#ifdef PARANOID
 if(Boundaries_pt==0)
  {
   throw OomphLibError("Node does not lie on any boundary",
                       "Node::ncoordinates_on_boundary()",
                       OOMPH_EXCEPTION_LOCATION);
  }


 //Does the node lie on the mesh boundary b
 if(!is_on_boundary(b))
  {
   std::ostringstream error_stream;
   error_stream << "Node is not on boundary " << b << std::endl;

   throw OomphLibError(error_stream.str(),
                       "Node::get_coordinates_on_boundary()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 //Check that the boundary coordinates have been set
 if(Boundary_coordinates_pt == 0)
  {
   std::ostringstream error_stream;
   error_stream << "Boundary coordinates have not been set\n"
                << "[Note: In refineable problems, the boundary coordinates\n"
                << "       will only be interpolated to newly created nodes\n"
                << "       if Mesh::Boundary_coordinate_exists[...] has been\n"
                << "       set to true!]\n";
   throw OomphLibError(error_stream.str(),
                       "Node::ncoordinates_on_boundary()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
  
 //Find out how may coordinates there are from the map
 return (*Boundary_coordinates_pt)[b]->nrow();

}



//=========================================================================
/// Given the mesh boundary b, return the k-th generalised boundary 
/// coordinates of the node in the vector boundary_zeta
//=========================================================================
void BoundaryNodeBase::
get_coordinates_on_boundary(const unsigned &b, const unsigned& k,
                            Vector<double> &boundary_zeta)
{
 //Check that the node lies on a boundary
#ifdef PARANOID
 if(Boundaries_pt==0)
  {
   throw OomphLibError("Node does not lie on any boundary",
                       "Node::coordinate_on_boundary()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 //Does the node lie on the mesh boundary b
 if(!is_on_boundary(b))
  {
   std::ostringstream error_stream;
   error_stream << "Node is not on boundary " << b << std::endl;

   throw OomphLibError(error_stream.str(),
                       "Node::get_coordinates_on_boundary()",
                       OOMPH_EXCEPTION_LOCATION);
  }


#ifdef PARANOID
 //Check that the boundary coordinates have been set
 if(Boundary_coordinates_pt == 0)
  {
   std::ostringstream error_stream;
   error_stream << "Boundary coordinates have not been set\n"
                << "[Note: In refineable problems, the boundary coordinates\n"
                << "       will only be interpolated to newly created nodes\n"
                << "       if Mesh::Boundary_coordinate_exists[...] has been\n"
                << "       set to true!]\n";
   throw OomphLibError(error_stream.str(),
                       "Node::get_coordinates_on_boundary()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif
 
 
 //Find out how may coordinates there are from the map
 const unsigned nboundary_coord = (*Boundary_coordinates_pt)[b]->nrow();
#ifdef PARANOID
 if(nboundary_coord != boundary_zeta.size())
  {
   std::ostringstream error_stream;
   error_stream 
    << "Wrong number of coordinates in the vector boundary_zeta"
    << std::endl << "There are " << nboundary_coord 
    << " boundary coordinates"
    << std::endl << "But bounday_zeta() has size " << boundary_zeta.size()
    << std::endl;

   throw OomphLibError(error_stream.str(),
                       "Node::coordinate_on_boundary()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 //Loop over and assign the coordinates
 for(unsigned i=0;i<nboundary_coord;i++)
  {boundary_zeta[i] = (*(*Boundary_coordinates_pt)[b])(i,k);}
}


//=========================================================================
/// Given the mesh boundary b, set the k-th generalised boundary 
/// coordinates of the node from the vector boundary_zeta
//=========================================================================
void BoundaryNodeBase::
set_coordinates_on_boundary(const unsigned &b,  const unsigned& k,
                            const Vector<double> &boundary_zeta)
{
 //Check that the node lies on a boundary
#ifdef PARANOID
 if(Boundaries_pt==0)
  {
   throw OomphLibError("Node does not lie on any boundary",
                       "Node::set_coordinates_on_boundary()",
                       OOMPH_EXCEPTION_LOCATION);
  }
#endif

 //Does the node lie on the mesh boundary b
 if(!is_on_boundary(b))
  {
   std::ostringstream error_stream;
   error_stream << "Node is not on boundary " << b << std::endl;

   throw OomphLibError(error_stream.str(),
                       "Node::set_coordinates_on_boundary()",
                       OOMPH_EXCEPTION_LOCATION);
  }

 //If the storage has not been assigned, then assign it
 if(Boundary_coordinates_pt == 0)
  {
   Boundary_coordinates_pt = new std::map<unsigned,DenseMatrix<double> *>;
  }

  
 //Find the number of boundary coordinates
 const unsigned nboundary_coord = boundary_zeta.size();
 
 //Allocate the vector for the boundary coordinates, if we haven't already
 if((*Boundary_coordinates_pt)[b] == 0)
  {
   // Need k+1 columns initially
   (*Boundary_coordinates_pt)[b] = 
    new DenseMatrix<double>(nboundary_coord,k+1);
  }
 //Otherwise resize it, in case the number of boundary coordinates 
 //or the number of types has changed
 else
  {
   // Adjust number of boundary coordinates -- retain number of types
   unsigned ncol=(*Boundary_coordinates_pt)[b]->ncol();
   {
    (*Boundary_coordinates_pt)[b]->resize(nboundary_coord,ncol);
   }
   
   // Resize number of types if required
   if ((k+1)>(*Boundary_coordinates_pt)[b]->ncol())
    {
     (*Boundary_coordinates_pt)[b]->resize(nboundary_coord,k+1);
    }
  }
 
 //Loop over and assign the coordinates
 for(unsigned i=0;i<nboundary_coord;i++)
  {(*(*Boundary_coordinates_pt)[b])(i,k) = boundary_zeta[i];}
}




///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//Functions for the SolidNode class
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


//=================================================================
/// Private function to check that the argument is within the
/// range of stored Lagrangain coordinates and position types.
//=================================================================
void SolidNode::xi_gen_range_check(const unsigned &k, const unsigned &i) const
{
 //If either the coordinate or type are out of range
 if((i >= Nlagrangian) || (k >= Nlagrangian_type))
  {
   std::ostringstream error_message;
   //If it's the lagrangian coordinate
   if(i >= Nlagrangian)
    {
     error_message << "Range Error: Xi coordinate " << i
                   << " is not in the range (0,"
                   << Nlagrangian-1 << ")";
    }
   //If it's the position type
   if(k >= Nlagrangian_type)
    {
     error_message << "Range Error: Lagrangian type " << k
                   << " is not in the range (0,"
                   << Nlagrangian_type-1 << ")";
    }
   
   throw OomphLibError(error_message.str(),
                       "SolidNode::xi_gen_range_check()",
                       OOMPH_EXCEPTION_LOCATION);
  }
}
  


//========================================================================
///  Steady constructor. The node has NLgrangian Lagrangian 
/// coordinates of n_lagrangian_type types (1 for Lagrange elements, 
/// 2 for 1D Hermite etc.).
/// The Eulerian dimension of the Node is n_dim and we have n_position_type
/// (generalised) Eulerian coordinates. There are 
/// initial_n_value values stored at
/// this node and NAdditional_solid additional values associated with the 
/// solid equations are stored in a separate Data object at the node.
//========================================================================
SolidNode::SolidNode(const unsigned &n_lagrangian, 
                     const unsigned &n_lagrangian_type,
                     const unsigned &n_dim, 
                     const unsigned &n_position_type,
                     const unsigned &initial_n_value)
 : Node(n_dim,n_position_type,initial_n_value,false),
   Nlagrangian(n_lagrangian), Nlagrangian_type(n_lagrangian_type)
{
 //Calculate the total storage required for positions
 const unsigned n_storage = Ndim*Nposition_type;

 //Allocate a data object with exactly the coorect number of coordinates
 Variable_position_pt = new Data(n_storage);
 //Set X_position to point to the data's positions
 X_position = Variable_position_pt->Value;

 //Setup the lagrangian storage
 const unsigned n_lagrangian_storage = n_lagrangian*n_lagrangian_type;
 Xi_position = new double[n_lagrangian_storage];
 //Initialise lagrangian positions to zero
 for(unsigned j=0;j<n_lagrangian_storage;j++) {Xi_position[j] = 0.0;}
}

//========================================================================
/// Unsteady constructor.  
/// Allocates storage for initial_n_value nodal values with history values
/// as required by timestepper.
/// The node has NLgrangian Lagrangian coordinates of
/// n_lagrangian_type types (1 for Lagrange elements, 2 for 1D Hermite etc.)/
/// The Eulerian dimension of the Node is n_dim and we have n_position_type
/// generalised Eulerian coordinates. 
//========================================================================
SolidNode::SolidNode(TimeStepper* const &time_stepper_pt,
                     const unsigned &n_lagrangian,
                     const unsigned &n_lagrangian_type,
                     const unsigned &n_dim, 
                     const unsigned &n_position_type,
                     const unsigned &initial_n_value)
 : Node(time_stepper_pt,n_dim,n_position_type,initial_n_value,false),
   Nlagrangian(n_lagrangian), Nlagrangian_type(n_lagrangian_type)
{
 //Calculate the total storage required for positions
 const unsigned n_storage = n_dim*n_position_type;

 //Allocate a Data value to each element of the Vector
 Variable_position_pt = new Data(time_stepper_pt,n_storage);
 //Set the pointer
 X_position = Variable_position_pt->Value;

 //Setup the lagrangian storage
 const unsigned n_lagrangian_storage = n_lagrangian*n_lagrangian_type;
 Xi_position = new double[n_lagrangian_storage];
 //Initialise lagrangian positions to zero
 for(unsigned j=0;j<n_lagrangian_storage;j++) {Xi_position[j] = 0.0;}
}

//========================================================================
///Destructor to clean up the memory allocated for nodal positions and
///additional solid variables
//========================================================================
SolidNode::~SolidNode()
{
 //Null out X_position so that the Node destructor doesn't delete it
 X_position=0;
 //Delete the position data
 delete Variable_position_pt;  Variable_position_pt = 0;
 //Now clean up lagrangian position data
 delete[] Xi_position; Xi_position=0;
}


//================================================================
/// Copy nodal positions and associated data from specified
/// node object
//================================================================
void SolidNode::copy(SolidNode* orig_node_pt)
{
 // Eulerian positions are stored as Data, so copy the data values
 // from one data to another
 Variable_position_pt->copy(orig_node_pt->variable_position_pt());
 
 //Copy the Lagrangian coordinates
 const unsigned nlagrangian_storage = Nlagrangian*Nlagrangian_type;

 // Check # of values:
 const unsigned long nlagrangian_storage_orig
  = orig_node_pt->nlagrangian()*orig_node_pt->nlagrangian_type();
 if (nlagrangian_storage!=nlagrangian_storage_orig)
  {
   std::ostringstream error_stream;
   error_stream << "The allocated lagrangian storage " 
                << nlagrangian_storage 
                << " is not the same as the original Solid Node "
                << nlagrangian_storage_orig << std::endl;

   throw OomphLibError(error_stream.str(),
                       "Node::copy()",
                       OOMPH_EXCEPTION_LOCATION);
  }
 
 // Copy lagrangian positions
 for(unsigned j=0;j<nlagrangian_storage;j++) 
  {
   Xi_position[j] = orig_node_pt->Xi_position[j];
  }

 // Copy the associated data
 Data::copy(orig_node_pt);

}


//================================================================
///Dump nodal positions and associated data to file for restart
//================================================================
void SolidNode::dump(std::ostream& dump_file)
{
 //Dump out the Lagrangian coordinates
 // Number of lagrangian values
 const unsigned nlagrangian_storage = Nlagrangian*Nlagrangian_type;
 dump_file << nlagrangian_storage 
           << " # number of Lagrangian position variables" << std::endl;

 for(unsigned j=0;j<nlagrangian_storage;j++) 
  {
   dump_file << Xi_position[j] << std::endl;
  }

 // Dump out Eulerian positions and nodal data
 Node::dump(dump_file);
}

//================================================================
/// Read nodal positions and associated data to file for restart
//================================================================
void SolidNode::read(std::ifstream& restart_file)
{
 std::string input_string;
 
 // Number of lagrangian values
 const unsigned nlagrangian_storage = Nlagrangian*Nlagrangian_type;

  // Read line up to termination sign
  getline(restart_file,input_string,'#');
 // Ignore rest of line
 restart_file.ignore(80,'\n');
 // Check # of values:
 const unsigned long check_nlagrangian_storage=atoi(input_string.c_str());
 if(check_nlagrangian_storage!=nlagrangian_storage)
  {
   std::ostringstream error_stream;
   error_stream << "The allocated Lagrangian storage " 
                << nlagrangian_storage << 
    " is not the same as that in the input file"
                << check_nlagrangian_storage << std::endl;
     
   throw OomphLibError(error_stream.str(),
                       "SolidNode::read()",
                       OOMPH_EXCEPTION_LOCATION);
   
  }
 
 // Read Lagrangian positions
 for(unsigned j=0;j<nlagrangian_storage;j++) 
  {
   // Read line
   getline(restart_file,input_string);
   
   // Transform to double
   Xi_position[j] = atof(input_string.c_str());
  }
 
 // Read Eulerian positions and nodal data
 Node::read(restart_file);
}

//===================================================================
/// Set the variable position data from an external source. 
/// This is mainly used when setting periodic solid problems.
//==================================================================
void SolidNode::set_external_variable_position_pt(Data* const &data_pt)
{
 //Wipe the existing value
 delete Variable_position_pt;
 //Set the new value
 Variable_position_pt = new CopiedData(data_pt);
 //Set the new value of x
 X_position = Variable_position_pt->Value;
}


//=======================================================================
/// Return lagrangian coordinate
/// either directly or via hanging node representation.
//======================================================================
double SolidNode::lagrangian_position(const unsigned &i) const
{
 double posn=0.0;

 // Non-hanging node: just return value
 if(!is_hanging()) {posn=xi(i);}
 // Hanging node: Use hanging node representation
 else
  {
   // Initialise
   double position=0.0;
   
   // Add contribution from master nodes
   const unsigned nmaster=hanging_pt()->nmaster();
   for (unsigned imaster=0;imaster<nmaster;imaster++)
    {
     position+=
      static_cast<SolidNode*>(hanging_pt()->master_node_pt(imaster))->xi(i)*
      hanging_pt()->master_weight(imaster);
    } 
   posn=position;
  }
 return posn;
}

//=======================================================================
/// Return generalised lagrangian coordinate
/// either directly or via hanging node representation.
//======================================================================
double SolidNode::lagrangian_position_gen(const unsigned &k,
                                          const unsigned &i) const
{
 double posn=0.0;

 // Non-hanging node: just return value
 if(!is_hanging()) {posn=xi_gen(k,i);}
 // Hanging node: Use hanging node representation
 else
  {
   // Initialise
   double position=0.0;
   
   // Add contribution from master nodes
   const unsigned nmaster=hanging_pt()->nmaster();
   for (unsigned imaster=0;imaster<nmaster;imaster++)
    {
     position+=
      static_cast<SolidNode*>(hanging_pt()->master_node_pt(imaster))
      ->xi_gen(k,i)*
      hanging_pt()->master_weight(imaster);
    } 
   posn=position;
  }
 return posn;
}

//================================================================
/// Assign (global) equation number, for SolidNodes
//================================================================
void SolidNode::assign_eqn_numbers(unsigned long &global_number, 
                                   Vector<double*>& dof_pt)
{
 //Let's call position equations first
 Variable_position_pt->assign_eqn_numbers(global_number,dof_pt);
 //Then call standard Data assign_eqn_numbers 
 Data::assign_eqn_numbers(global_number,dof_pt);
} 

}
