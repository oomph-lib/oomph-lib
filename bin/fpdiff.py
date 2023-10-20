#!/usr/bin/env python3

import sys


def gettype(a): 
 """ Distinguish between a number and a string:
     
     Returns integer 1 if the argument is a number,
                     2 if the argument is a string.
 """
 import re
 type1 = re.compile("(^[+-]?[0-9]*[.]?[0-9]+$)|(^[+-]?[0-9]+[.]?[0-9]*$)|(^[+-]?[0-9]?[.]?[0-9]+[EeDd][+-][0-9]+$)")
 if type1.match(a):
  return 1
 else:
  return 2

def stuff(string,symbol,number): 
 """ Add number copies of symbol to  string

     Returns modified string.
 """
 for i in range(number):
  string += symbol
 string += " "
 return string 

def read_file(filename):
 """ Read file into a list of strings, uncluding direct reading of ".gz" files
 Different operations required for gzip files in Python 3 because they are
 read as binary (rather than text) files
 """
 import gzip
    
 #If the first file is a gzipped file, open it via the gzip module
 if(filename.find(".gz") != -1):
   F=gzip.open(filename)
   filedata = [l.decode() for l in F.readlines() ] 
   F.close()
 #Otherwise open as normal
 else:
   F=open(filename);
   filedata = F.readlines();
   F.close() 

 return filedata

def fpdiff_helper(filename1,filename2,relative_error,small,
                  outstream, details_stream):
 """ Calculate the floating-point difference between two data files.
     The idea is to use a looser tolerance than the UNIX diff command,
     so that if two entries have a relative error less than the argument
     relative_error, they are counted as the same.

     Note that the relative error is percentage!
     
     Information on pass/failure is written to outstream. Details on which
     lines failed are written to details_stream. Warning: if run on
     large files the details_stream may be overwhelmingly long.
     
     First return value: 0 if the two files are the same, 1 if they are
     different or 5 if the files cannot be opened.

     Second return value: the maximum relative error.

     Third return value: the largest entry that caused an error
     (i.e. what "small" would need to be set as for there to be no
     differences).
 """

 import math

 # Storage for worst case error sizes
 max_rel_diff = 0
 max_wrong_entry = 0 
 
 # Open the files (if run as a script then open failures are handled higher
 # up by catching the error, otherwise it is the parent program's job to
 # handle the error and so we shouldn't do anything weird here).
 tmpfile1 = read_file(filename1)
 tmpfile2 = read_file(filename2)
 # this line catches the error when an empty selftest reference data file is provided and will return a FAIL
 for file_length, filename in [(len(tmpfile1), filename1), (len(tmpfile2), filename2)]:
  if file_length == 0:
   details_stream.write(f"\nWarning: file {filename} is empty and the test will therefore be treated as a fail.")
   outstream.write("\n\n   [FAILED]\n")
   return 1

 #Find the number of lines in each file
 n1 = len(tmpfile1)
 n2 = len(tmpfile2)

 #If file1 has more lines than file2, keep order the same 
 if n1 >= n2:
  file1 = tmpfile1; file2 = tmpfile2; n = n2
 #Otherwise swap the order of the files
 else:
  file1 = tmpfile2; file2 = tmpfile1; n = n1

 #Counter for the number of errors
 nerr = 0
 #Counter for the number of lines
 count = -1
 #Counter for the number of lines with errors
 nline_error = 0

 #Loop over the lines in file1 (the file with the most lines!)
 for line1 in file1:
  #Increase the counter
  count += 1
  #If we've run over the end of the file2, issue a warning and end the loop
  if count >= n:
   details_stream.write("\nWarning: files have different numbers of lines")
   details_stream.write("\nResults are for first %d lines of both files\n" % count)
   nerr += 1
   break
  #Read the next line from file2
  line2 = file2[count]

  #If the lines are the same, we're done
  if(line1 == line2):
   continue
  #If not need to do more work
  else:
   #Split each line into its separate fields
   fields1 = line1.split(); fields2 = line2.split()
   #Find the number of fields in each line
   nfields1 = len(fields1); nfields2 = len(fields2)

   #If the number of fields is not the same, report it as an error
   if nfields1 != nfields2:
     details_stream.write("\n =====> line %d: different number of fields\n" \
     % (count+1))
     details_stream.write("%s fields: %s" % (nfields1, line1)) 
     details_stream.write("%s fields: %s" % (nfields2, line2))
     nerr += 1
     continue
    
   #Otherwise, we now compare field by field
   else:
    #Flag to indicate whether there has been a problem in the field
    problem = 0
    #Strings that will hold the output data
    outputline1 = ""; outputline2 = ""; outputline3 = ""

    #Loop over the fields
    for i in range(nfields1):
     #Start by loading the fields into the outputlines (plus whitespace)
     outputline1 += fields1[i] + " "; outputline3 += fields2[i] + " "

     #Find the lengths of the fields
     length1 = len(fields1[i]); length2 = len(fields2[i])

     #Pad the shortest field so the lengths are the same
     if length1 < length2:
      fieldlength = length2
      for j in range(length2-length1):
       outputline1 += " "
     else:
      fieldlength = length1
      for j in range(length1 - length2):
       outputline3 += " "
      
     #If the fields are identical, we are fine
     if fields1[i] == fields2[i]:
      #Put spaces into the error line
      outputline2 = stuff(outputline2," ",fieldlength)
     #Otherwise time for yet more analysis
     else:
      #Find the type (numeric or string) of each field
      type1 = gettype(fields1[i]); type2 = gettype(fields2[i])

      #If the data-types aren't the same issue an error
      if type1 != type2:
       problem = 1
       nerr += 1
       #Put the appropriate symbol into the error line
       outputline2 = stuff(outputline2,"*",fieldlength)
      #Otherwise more analysis
      #If the types are both strings then report the error
      elif type1 == 2:
        problem = 1
        nerr += 1
        #Put the appropriate symbol into the error line
        outputline2 = stuff(outputline2,"%",fieldlength)
      else:
       #Convert strings to floating point number
       x1 = float(fields1[i].lower().replace("d","e"))
       x2 = float(fields2[i].lower().replace("d","e"))

       #If both numbers are very small, that's fine
       if math.fabs(x1) <= small and math.fabs(x2) <= small:
        #Put spaces into the error line
        outputline2 = stuff(outputline2," ",fieldlength)
       else:
        #Find the relative difference based on the largest number
	#Note that this "minimises" the relative error (in some sense)
	#but means that I don't have to separately trap the cases
	#when x1, x2 are zero
        if math.fabs(x1) > math.fabs(x2) :
         diff = 100.0*(math.fabs(x1 - x2) / math.fabs(x1))
        else:
         diff = 100.0*(math.fabs(x1 - x2) / math.fabs(x2))

        #If the relative error is smaller than the tolerance, that's fine
        if diff <= relative_error:
         #Put spaces into the error line
         outputline2 = stuff(outputline2," ",fieldlength)
        #Otherwise issue an error
        else:
         problem = 1
         nerr += 1
         #Put the appropriate symbols into the error line 
         outputline2 = stuff(outputline2,"-",fieldlength)

         # Record any changes in the worst case values
         if diff > max_rel_diff:
            max_rel_diff = diff

         if math.fabs(x1) > max_wrong_entry:
            max_wrong_entry = math.fabs(x1)
         elif  math.fabs(x2) > max_wrong_entry:
            max_wrong_entry = math.fabs(x2)

   #If there has been any sort of error, print it
   if problem == 1:
    nline_error += 1
    details_stream.write("\n =====> line %d\n" % (count+1))
    details_stream.write("%s\n%s\n%s\n" % (outputline1, outputline2, outputline3))
 
 #Final print out, once loop over lines is complete
 if nerr > 0:
  outstream.write("\n In files %s %s" % (filename1, filename2))
  outstream.write("\n number of lines processed: %d" % count)
  outstream.write("\n number of lines containing errors: %d" % nline_error)
  outstream.write("\n number of errors: %d " % nerr)
  outstream.write("\n largest relative error: %g " % max_rel_diff)
  outstream.write("\n largest abs value of an entry which caused an error: %g "
                   % max_wrong_entry)
  outstream.write("\n========================================================")
  outstream.write("\n    Parameters used:")
  outstream.write("\n        threshold for numerical zero : %g" % small)
  outstream.write("\n        maximum rel. difference [percent] : %g" % relative_error)
  outstream.write("\n    Legend: ")
  outstream.write("\n        *******  means differences in data type (string vs number)")
  outstream.write("\n        -------  means real data exceeded the relative difference maximum") 
  outstream.write("\n        %%%%%%%  means that two strings are different")
  outstream.write("\n========================================================")
  outstream.write("\n\n   [FAILED]\n")
  # Return failure
  return 2, max_rel_diff, max_wrong_entry

 else:
  outstream.write("\n\n In files %s %s" % (filename1, filename2))
  outstream.write(\
  "\n   [OK] for fpdiff.py parameters: - max. rel. error = %g " % relative_error)
  outstream.write("%")
  outstream.write(\
  "\n                                  - numerical zero  = %g\n" % small)
  # Return success
  return 0, max_rel_diff, max_wrong_entry


def fpdiff(filename1, filename2, relative_error=0.1, small=1e-14,
           outstream=sys.stdout, details_stream=sys.stdout):
    """Wrapper for using fpdiff inside python. Has default args and returns a
    bool."""
    ok, max_rel_diff, max_wrong_entry = \
       fpdiff_helper(filename1, filename2, relative_error, small,
                         outstream, details_stream)

    return (ok == 0), max_rel_diff, max_wrong_entry


def run_as_script(argv):
 """Run fpdiff as a script (handles argument parsing, output as error codes
    and some helpful messages).
 """
 # Note that we shouldn't just put this code this under 'if __name__ ==
 # "__main__":' because variables created there are global. This resulted
 # in some bugs before.
    
 #Set the defaults
 maxreld = 1.0e-1 # max relative difference in percent
 small = 1.0e-14  # small number -- essentially round-off error 

 #Remove the program name from the front of the argument list
 argv.pop(0)
 
 #Let's find the number of command line arguments
 narg = len(argv)

 #If we're out of range, issue a usage message 
 if narg < 2 or narg > 4:
  sys.stdout.write("\n      *********   ERROR   **********\n")
  sys.stdout.write("\nMust specify 2, 3 or 4 keywords on the command line. ")
  sys.stdout.write("You have specified %d" % narg) 
  sys.stdout.write("\n   Proper usage:  ")
  sys.stdout.write("\n         fpdiff file1 file2 [max_rel_diff_percent] [small]\n")
  sys.stdout.write("\n      *********  PROGRAM TERMINATING   ***********")
  sys.stdout.write("\n   [FAILED] \n")
  sys.exit(4)
 
 #Read any optional arguments
 if narg >= 3:
  maxreld = float(argv[2])
  if narg == 4:
   small = float(argv[3])
 
 # Run the diff
 try:
    error_code, _, _ = fpdiff_helper(argv[0], argv[1], maxreld, small,
                                  sys.stdout, sys.stdout)

 # If there has been an IO error then fail with a useful message
 except IOError:

     # Get the exception that was raised
     _, err, _ = sys.exc_info()
     # We have to get exceptions manually using this function instead of the
     # usual `except IOError as err:` or `except IOError, err:` for
     # compatibility with both the ancient version of python on the wulfling
     # and python3+.

     # Write the message
     sys.stderr.write("\n   [FAILED] I/O error(%d): %s \"%s\"\n"
                      % (err.errno, err.strerror, err.filename))

     return 5

 return error_code




# What to do if this is run as a script, rather than loaded as a module
if __name__ == "__main__":

 # Run and return whether it succeeded or not
 sys.exit(run_as_script(sys.argv))
