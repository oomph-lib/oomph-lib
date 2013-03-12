#!/usr/bin/env python

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

def fpdiff(filename1,filename2,relative_error,small):
 """ Calculate the floating-point difference between two data files.
     The idea is to use a looser tolerance than the UNIX diff command,
     so that if two entries have a relative error less than the argument
     relative_error, they are counted as the same.

     Returns nothing
 """

 import math
 import gzip
 
 #Open the files

 #If the first file is a gzipped file, open it via the gzip module
 try:
  if(filename1.find(".gz") != -1):
   F1=gzip.open(filename1)
 #Otherwise open as normal
  else:
   F1=open(filename1);
 #If there has been an IO error fail
 except IOError:
   print 
   print "   [FAILED] : Unable to open the input file", filename1
   print

 #read contents of file1 into a list and close
 tmpfile1 = F1.readlines(); F1.close()

 #If the second file is a gzipped file, open it via the gzip module
 try:
  if(filename2.find(".gz") != -1):
   F2=gzip.open(filename2)
 #Otherwise open as normal
  else:
   F2=open(filename2);
 #If there has been an IO error, fail
 except IOError:
   print
   print "   [FAILED] : Unable to open the input file", filename2
   print

 #read contents of file2 into a list and close
 tmpfile2 = F2.readlines(); F2.close()

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
   print
   print "Warning: files have different numbers of lines"
   print "Results are for first", count, "lines of both files" 
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
     print "\n =====> line", count+1,": different number of fields"
     print nfields1, "fields:", line1; 
     print nfields2, "fields:", line2
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
        if diff <= maxreld:
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
    print "\n =====> line", count+1
    print outputline1, "\n", outputline2, "\n", outputline3	
 
 #Final print out, once loop over lines is complete
 if nerr > 0: 
  print
  print "number of lines processed: ", count
  print "number of lines containing errors: ", nline_error
  print "number of errors: ", nerr
  print "========================================================"
  print "    Parameters used:"
  print "        threshold for numerical zero : ", small
  print "        maximum rel. difference [percent] : ", maxreld
  print "    Legend: "
  print "        *******  means differences in data type (string vs number)"
  print "        -------  means real data exceeded the relative difference maximum" 
  print "        %%%%%%%  means that two strings are different"
  print "========================================================"
  print
  print "   [FAILED]"
  print
  # Return non-zero since it failed
  return 2
 else:
  print
  print "   [OK] for fpdiff.py parameters: - max. rel. error = ",maxreld,"%"
  print "                                  - numerical zero  = ",small
  print
  # Success: return 0
  return 0

#What to do if this is run as a script, rather than loaded as a module
if __name__ == "__main__":
 
 import sys

 #Set the defaults
 maxreld = 1.0e-1 # max relative difference in percent
 small = 1.0e-14  # small number -- essentially round-off error 

 #Remove the program name from the front of the argument list
 sys.argv.pop(0)
 
 #Let's find the number of command line arguments
 narg = len(sys.argv)

 #If we're out of range, issue a usage message 
 if narg < 2 or narg > 4:
  print "\n      *********   ERROR   **********\n"
  print "Must specify 2, 3 or 4 keywords on the command line. ",
  print "You have specified", narg, "\n" 
  print "   Proper usage:  "
  print "         fpdiff file1 file2 [max_rel_diff_percent] [small]\n"
  print "      *********  PROGRAM TERMINATING   ***********"
  print "   [FAILED] "
  assert 0
 
 #Read any optional arguments
 if narg >= 3:
  maxreld = float(sys.argv[2])
  if narg == 4:
   small = float(sys.argv[3])
 
 # Run the diff
 return_val = fpdiff(sys.argv[0],sys.argv[1],maxreld,small)

 # Return whether it succeeded or not
 sys.exit(return_val)
