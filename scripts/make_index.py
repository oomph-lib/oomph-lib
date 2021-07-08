#!/usr/bin/env python

from __future__ import print_function
#Script to generate an index page for oomph-lib from a very simple
#input file

# Global (bit bad) counter to give unique ids for each division for use
# when hiding sections of text
div_counter=0

#Function that is used to convert the entries of a list into a key
#and a link or cross-reference, provided that the second entry has
#the appropriate symbol % (link) ^ (cross-ref)
def get_key_and_link(entry):
 #Initialise strings to be empty
 key=''
 link=''
 synonym=''
 n = len(entry)
 #If there is only one entry, we only have a key
 if n==1:
  key = entry[0]
 #Otherwise we could have a link or cross reference as well
 elif n > 1:
  key = entry[0]
  #If the entry starts with a % it's a link, find the relative link name
  if entry[1][0]=='%':
   link = entry[1].lstrip('%')
  #If it starts with a ^ it's a cross-ref
  #Determine the relative hyperlink and
  #also provide the synonym to which we are referring in synonym
  elif entry[1][0]=='^':
   anchor = entry[1].lstrip('^').strip().replace(' ','')
   synonym = "see entry on " + entry[1].lstrip('^').strip().replace('.',':')
   link = "index/html/index.html#" + anchor[0].upper() + '.' + anchor 
 return (key,link,synonym)  

#Output the key and link or cross-reference in pretty format
def print_key(doc_root_directory,key,level,anchor):
 #Increment the global counter to give unique text block IDs
 global div_counter
 div_counter += 1
 style=""
 #Set the html styles for the different levels
 if level == 1:
  begin_tag = "<h2>"
  end_tag = "</h2>"
 elif level == 2:
  begin_tag = "<h3>"
  end_tag = "</h3>"
  #Uncomment to collapse anything below level 1 tags by default
  #  style = " style=\"display:none\""
 else: 
  begin_tag = "<li>"
  end_tag = "</li>"
  #Uncomment to collapse anything below level 1 tags by default
  #  style = " style=\"display:none\""

 #Create unique tags for html anchors based on the tags in the file
 #with spaces removed and '.' instead of '@'
 new_anchor = anchor.replace(' ','')
 anchor = new_anchor
 anchor_string = "<a class=\"anchor\" id=\"" \
 + anchor + "\"></a>"
 anchor_tag = "<a href=\"#" + anchor + "\">"

 #If we don't have a link then it's a submenu heading, which we allow to
 #collapse and uncollapse by the following html magic
 #The link is to its own anchor so it shouldn't move much, but has the
 #same look as real links (In the end I found this confusing, so have removed 
 #it)
 if key[1]=='':
  print(begin_tag, anchor_string) #, \
  #Uncomment for self anchor-tag 
  #"\htmlonly",\
  #anchor_tag, key[0], "</a>"
  print(key[0])
  #Uncomment for the collapsing section stuff
  #print "<a onclick=\"index_toggle(\'div_id",div_counter,"\',\'im_id",\
  #div_counter,"\')\">", \
  #"<img src=\"../collapse.png\" id=\"im_id",div_counter,"\">", "</a>",
  #print  "\endhtmlonly"
  print(end_tag, \
  "<div id=\"div_id",div_counter,"\"",style,">")

 #Otherwise just print the link and key
 else:
  link_target = doc_root_directory + "/" + key[1]
  if key[2]=='':
   output_key = key[0]
   print(begin_tag, \
   anchor_string, \
   "<a href=\"",link_target,"\">",output_key,"</a>",end_tag)
  else:
   output_key = key[2]
   print(begin_tag, \
   anchor_string, \
   key[0], "<br> &nbsp; &nbsp; &nbsp;  <a href=\"", link_target,"\">",output_key,"</a>",end_tag)

#The main nested output function. The idea is to loop through the list
#Output the keys and subkeys and keep going until we hit a link
#or cross-reference. The index must be sorted before first call of
#this function
def nested_output(doc_root_directory,index,level=0,anchor=""):
 level += 1
 #Return when we have no subindex
 if index == []: 
  return 
 if level>2 :
  print("<ul>")

 #The key is the first entry in the first list
 key = get_key_and_link(index[0])
 #Initialise subindex
 subindex=[]
 #Loop over entries in index
 for entry in index:
  #If the key has changed
  if entry[0] != key[0]:
   #Make unique anchor from the fields
   if level == 1:
    anchor_stem = key[0]
    anchor = key[0]
   else:
    anchor_stem = anchor
    anchor += "." + key[0]

   #Output the list we have so far
   print_key(doc_root_directory,key,level,anchor)
   nested_output(doc_root_directory,subindex,level,anchor)
   #Restore to the previous stem anchor
   anchor = anchor_stem

   #If there is no key the entry is not a link, so we end the text div
   if key[1]=='':
    print("</div>")
   #Get the new key
   key = get_key_and_link(entry)
   #There will only be a subindex entry if we have no link
   if key[1]=='' and len(entry) > 1:
    subindex = [entry[1:len(entry)]]
   else:
    subindex = []
  #Othewise the key is the same so add to the subindex list
  else:
   if len(entry) > 1:
    if entry[1][0]=='%':
     if entry[1].lstrip('%') != key[1]:
      print("Multiple links for the same key", key[1], "\n")
      assert 0
    elif entry[1][0]=='^':
     #Error trapping here will be a pain, but should be done
     x=1
    #Otherwise not a link so add to the subindex
    else:
     subindex.append(entry[1:len(entry)]) 

 #Output final entry
 #Make unique anchor from the fields
 if level == 1:
  anchor_stem = key[0]
  anchor = key[0]
 else:
  anchor_stem = anchor
  anchor += "." + key[0]
 print_key(doc_root_directory,key,level,anchor)
 nested_output(doc_root_directory,subindex,level,anchor)
 #Restore to the prevous stem
 anchor = anchor_stem

 if key[1]=='':
  print("</div>")
 if level>2 :
  print("</ul>")


def make_index(filename,doc_root_directory):

 main_index = []
 index_file=open(filename)
 entire_file = index_file.read().split("@end")

 #Loop over lines in file split by @end
 for line in entire_file:
  #Split into tags
  entry = line.split('@')
  # Remove trailing whitespace
  out = [field.strip() for field in entry]
  #Add stripped entries to the main list, if they have 
  #Non-null first entry
  if len(out[0]) > 1:
   main_index.append(out)

 #Sort the list
# main_index.sort()

 #Dump the output to a sorted file
 #sorted_index=open("sorted.idx",'w')
 #for entry in main_index:
 # n = len(entry)
 # for index,field in enumerate(entry):
 #  #If we have a link write a newline first
 #  if field[0]=='^' or field[0]=='%':
 #   sorted_index.write("\n")
 #  #Write the entry
 #  sorted_index.write(field)
 #  #Write the separator
 #  if index!=n-1:
 #   sorted_index.write(" @ ")
 #  else:
 #   sorted_index.write(" @end \n\n")
 #sorted_index.close()

 #Add first letters of each key as the first entry in the list
 for line in main_index:
  capital = line[0].upper()
  line.insert(0,capital[0])

 #Output Quick links to the Alphabet
 import string
 first_letters = "<h1>"
 for char in string.ascii_uppercase:
  first_letters += "<a href=" + doc_root_directory + \
  "/index/html/index.html#" + char + ">" + char + "</a>"
 first_letters += "</h1>"
 print("\htmlonly")
 print(first_letters)

 #Output the list
 nested_output(doc_root_directory,main_index)
 print("\endhtmlonly")
  

if __name__ == "__main__":
 import sys
 #Pop the program name from the front of the lis
 sys.argv.pop(0)
 #Check the number of argumens
 narg = len(sys.argv)
 if narg != 2:
  print("\n Must specify input filename and path to doc directory on the command line")
  assert 0

 make_index(sys.argv[0],sys.argv[1])

