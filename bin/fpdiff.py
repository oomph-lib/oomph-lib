#!/usr/bin/env python

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

def fpdiff_helper(filename1,filename2,relative_error,small):
 """ Calculate the floating-point difference between two data files.
     The idea is to use a looser tolerance than the UNIX diff command,
     so that if two entries have a relative error less than the argument
     relative_error, they are counted as the same.

     Note that the relative error is percentage!

     Returns 0 if the two files are the same, 1 if they are different or
     5 if the files cannot be opened.
 """

 import math

 
 #Open the files
 try:
   tmpfile1 = read_file(filename1)
 except IOError:
   #If there has been an IO error fail
   sys.stdout.write("\n   [FAILED] : Unable to open the input file %s \n\n"
      % filename1)
   return 5

 try:
   tmpfile2 = read_file(filename2)
 except IOError:
   #If there has been an IO error fail
   sys.stdout.write("\n   [FAILED] : Unable to open the input file %s \n\n"
      % filename2)
   return 5

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
   sys.stdout.write("\nWarning: files have different numbers of lines")
   sys.stdout.write("\nResults are for first %d lines of both files\n" % count)
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
     sys.stdout.write("\n =====> line %d: different number of fields\n" \
     % (count+1))
     sys.stdout.write("%s fields: %s" % (nfields1, line1)) 
     sys.stdout.write("%s fields: %s" % (nfields2, line2))
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
   
   #If there has been any sort of error, print it
   if problem == 1:
    nline_error += 1
    sys.stdout.write("\n =====> line %d\n" % (count+1))
    sys.stdout.write("%s\n%s\n%s\n" % (outputline1, outputline2, outputline3))
 
 #Final print out, once loop over lines is complete
 if nerr > 0:
  sys.stdout.write("\n In files %s %s" % (filename1, filename2))
  sys.stdout.write("\n number of lines processed: %d" % count)
  sys.stdout.write("\n number of lines containing errors: %d" % nline_error)
  sys.stdout.write("\n number of errors: %d " % nerr)
  sys.stdout.write("\n========================================================")
  sys.stdout.write("\n    Parameters used:")
  sys.stdout.write("\n        threshold for numerical zero : %g" % small)
  sys.stdout.write("\n        maximum rel. difference [percent] : %g" % relative_error)
  sys.stdout.write("\n    Legend: ")
  sys.stdout.write("\n        *******  means differences in data type (string vs number)")
  sys.stdout.write("\n        -------  means real data exceeded the relative difference maximum") 
  sys.stdout.write("\n        %%%%%%%  means that two strings are different")
  sys.stdout.write("\n========================================================")
  sys.stdout.write("\n\n   [FAILED]\n")
  # Return failure
  return 2

 else:
  sys.stdout.write("\n\n In files %s %s" % (filename1, filename2))
  sys.stdout.write(\
  "\n   [OK] for fpdiff.py parameters: - max. rel. error = %g " % relative_error)
  sys.stdout.write("%")
  sys.stdout.write(\
  "\n                                  - numerical zero  = %g\n" % small)
  # Return success
  return 0


def fpdiff(filename1, filename2, relative_error=0.1, small=1e-14):
    """Wrapper for using fpdiff inside python. Has default args and returns a
    bool."""
    return fpdiff_helper(filename1, filename2, relative_error, small) == 0


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
 error_code = fpdiff_helper(argv[0], argv[1], maxreld, small)

 return error_code




# What to do if this is run as a script, rather than loaded as a module
if __name__ == "__main__":

 # Run and return whether it succeeded or not
 sys.exit(run_as_script(sys.argv))
