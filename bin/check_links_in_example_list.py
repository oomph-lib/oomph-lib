#!/usr/bin/env python

#Little script to check for dead links
def search_dead_links(doc_root_directory):
 #Open the list of example codes
 filename=doc_root_directory + "/example_code_list/html/index.html"
 index_file=open(filename)
 entire_file = index_file.read()

 #Pull out all urls
 import re
 urls = re.findall(r'[Hh][Rr][Ee][Ff] *= *[\'"]?([^\'" >]+)',entire_file)

 import string
 #Now form the list of local urls 
 #Need to do some mangling to get the correct start, which
 #may be specific to UNIX installations, but only we are going to be
 #using this script so it should be OK
 local_urls = [] 
 for url in urls:
  #Replace items in the header with the same string as those in .txt
  new_url = string.replace(url,"../../../doc","../..")
  if new_url.startswith('..'):
   #Strip out anchors
   if "#" in new_url:
    local_urls.append(new_url[4:new_url.index("#")])
   else:
    local_urls.append(new_url[4:])

 #Now find all index.html files in the documentation tree
 find_command = "find " + doc_root_directory + " -name \'index.html\'"
 import subprocess
 proc = subprocess.Popen(find_command, shell=True, stdin=subprocess.PIPE,\
 stdout=subprocess.PIPE)
 proc.wait()
 local_files = proc.stdout.read()
 new_local_files = local_files.split()

 #Need to remove the documentation root and replace it by a .
 local_files = []
 for url in new_local_files:
  local_file = string.replace(url,doc_root_directory,"./")
  local_files.append(local_file)

 new_local_files = local_files

 #Finally print out ones that don't match
 print "Files not linked in example list file\n"

 for url in new_local_files:
  if url not in local_urls:
   print url


 # for url,file in map(None,local_urls,local_files):
 #  print url,file

if __name__ == "__main__":
 import sys
 #Pop the program name from the front of the lis
 sys.argv.pop(0)
 #Check the number of argumens
 narg = len(sys.argv)
 if narg != 1:
  print "\n Must specify path to doc directory on the command line"
  assert 0

 search_dead_links(sys.argv[0])

