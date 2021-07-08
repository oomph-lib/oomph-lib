#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=demo_code

# Name of sample makefile to built
noinst_SCRIPTS=makefile.sample

# Sources for executable
demo_code_SOURCES = demo_code.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
demo_code_LDADD = -L@libdir@ -lpoisson -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-------------------------------------------------------------------------


# Produce same executable (with different name) this time by
# using make outside the automake/autoconf enviroment
my_demo_code: makefile.sample demo_code.cc Makefile.am   
	make -f makefile.sample my_demo_code



#-------------------------------------------------------------------------


# Manually create makefile to demonstrate linking against oomph-lib
# from outside the automake/autoconf framework:
makefile.sample: 
	echo "###############################################################"> makefile.sample
	echo "# Automatically-generated sample makefile to illustrate how to  ">> makefile.sample
	echo "# link against oomph-lib from outside the automake/autoconf">> makefile.sample
	echo "# framework. Do not edit this -- make a copy first" >> makefile.sample
	echo "# " >> makefile.sample
	echo "# When customising this makefile, you should only have to change" >> makefile.sample
	echo "# " >> makefile.sample
	echo "# - the variable OOMPH-LIB_LIBS:" >> makefile.sample
	echo "#         Add any additional oomph-lib sub-libraries that " >> makefile.sample
	echo "#         you may wish to use in your code. " >> makefile.sample
	echo "# " >> makefile.sample
	echo "# - the specific dependencies for your driver code: " >> makefile.sample
	echo "#         Include any additional local dependencies such as " >> makefile.sample
	echo "#         your own header files etc." >> makefile.sample
	echo "# " >> makefile.sample
	echo "###############################################################" >> makefile.sample
	echo " "  >> makefile.sample
	echo " "  >> makefile.sample
	echo "# Installation-specific information -- don't change any of this! "  >> makefile.sample
	echo "#-------------------------------------------------------------- "  >> makefile.sample
	echo " "  >> makefile.sample
	echo "# Flags for C pre-processor " >> makefile.sample
	echo "AM_CPPFLAGS=-DHAVE_CONFIG_H -I. -I../.. " $(AM_CPPFLAGS) >> makefile.sample
	echo " "  >> makefile.sample
	echo "# Flags for C++ compiler " >> makefile.sample
	echo "CXXFLAGS=" $(CXXFLAGS) >> makefile.sample
	echo " "  >> makefile.sample
	echo "# Library include directory: This is where all the header files live" >> makefile.sample
	echo "OOMPH-LIB_INCLUDE_DIR="@includedir@  >> makefile.sample
	echo " "  >> makefile.sample
	echo "# Library directory: This is where all of oomph-lib's sub-libraries live"  >> makefile.sample
	echo "OOMPH-LIB_LIB_DIR="@libdir@ >> makefile.sample
	echo " "  >> makefile.sample
	echo "# These are the external (3rd party) libraries that are distributed" >> makefile.sample
	echo "# with oomph-lib and that we always link against"  >> makefile.sample
	echo "OOMPH-LIB_EXTERNAL_LIBS="$(EXTERNAL_LIBS)  >> makefile.sample
	echo " "  >> makefile.sample


	echo "# This specifies where libraries built from third party "  >> makefile.sample
	echo "# distributions can be found"  >> makefile.sample
	echo "EXTERNAL_DIST_LIBRARIES="$(AM_LDFLAGS)  >> makefile.sample
	echo " "  >> makefile.sample



	echo "# This is additional machine-specific linking information that " >> makefile.sample
	echo "# allows mixed-language compilation/linking"  >> makefile.sample
	echo "FLIBS="$(FLIBS)  >> makefile.sample
	echo " "  >> makefile.sample
	echo "# Flags required for the use of shared libraries "  >> makefile.sample
if HAVE_DARWIN
	echo "# Mac OSX with clang: Replace the above line with the following" >> makefile.sample	
	echo "SHARED_LIBRARY_FLAGS= -Wl,-bind_at_load" >> makefile.sample
else
	echo "# Linux: Replace the above line with the following" >> makefile.sample
	echo "SHARED_LIBRARY_FLAGS=-Wl,--rpath -Wl,"@libdir@  >> makefile.sample
endif
	echo " "  >> makefile.sample
	echo "# Mac OSX without clang(?): Replace the above line with the following" >> makefile.sample	
	echo "#SHARED_LIBRARY_FLAGS= --rpath="@libdir@  >> makefile.sample	
	echo " "  >> makefile.sample
	echo "# Problem-specific information -- edit this for your driver code "  >> makefile.sample
	echo " " >> makefile.sample
	echo "#---------------------------------------------------------------"  >> makefile.sample
	echo "# These are the specific oomph-lib sub-libraries that we have to link against"  >> makefile.sample
	echo "# for this driver code -- edit this according to your requirements"  >> makefile.sample
	echo "# but remember that the order of the libraries matters: List the"  >> makefile.sample
	echo "# the more specific ones before the more general ones!"  >> makefile.sample
	echo "OOMPH-LIB_LIBS=-lpoisson -lgeneric"  >> makefile.sample
	echo " "  >> makefile.sample
	echo "# Dependencies for this driver code and compile instructions: " >> makefile.sample
	echo "# Which local source (usually *.cc or *.h) files does the"  >> makefile.sample 
	echo "# driver code depend on?"  >> makefile.sample
	echo "my_demo_code.o: demo_code.cc"  >> makefile.sample
#
#
# NOTE \044 IS OCTAL ASCII CODE FOR $ TO STOP EXPANSION OF VARIABLES IN HERE
#
# ALSO NOTE THAT THERE MUST NOT BE ANY WHITESPACE AFTER CONTINUATION "\"
#
#

#old	echo -e "\t $(CXX) \044(AM_CPPFLAGS) -c demo_code.cc -o my_demo_code.o \\"   >> makefile.sample
	printf  "\t $(CXX) \044(AM_CPPFLAGS)  \044(CXXFLAGS) -c demo_code.cc -o my_demo_code.o \\"   >> makefile.sample
	echo ""   >> makefile.sample
#old	echo -e "\t       -I\044(OOMPH-LIB_INCLUDE_DIR)  " >> makefile.sample
	printf  "\t       -I\044(OOMPH-LIB_INCLUDE_DIR) \n " >> makefile.sample
	echo " "  >> makefile.sample
	echo "# Linking instructions: Just declare the target (i.e. the name of the executable)"  >> makefile.sample
	echo "# and the dependencies (i.e. the object files created above). The rest "  >> makefile.sample
	echo "# should not have to be changed. " >> makefile.sample
	echo "my_demo_code: my_demo_code.o"  >> makefile.sample
#old	echo -e "\t $(CXX) \044(SHARED_LIBRARY_FLAGS) \044< -o \044@ \\"   >> makefile.sample
	printf  "\t $(CXX) \044(SHARED_LIBRARY_FLAGS) \044< -o \044@ \\"   >> makefile.sample	
	echo ""   >> makefile.sample	
#old	echo -e "\t       -L\044(OOMPH-LIB_LIB_DIR) \044(EXTERNAL_DIST_LIBRARIES) \044(OOMPH-LIB_LIBS) \\" >> makefile.sample
	printf  "\t       -L\044(OOMPH-LIB_LIB_DIR) \044(EXTERNAL_DIST_LIBRARIES) \044(OOMPH-LIB_LIBS) \\" >> makefile.sample
	echo ""   >> makefile.sample	
#old	echo -e "\t        \044(OOMPH-LIB_EXTERNAL_LIBS) \044(FLIBS)"  >> makefile.sample
	printf  "\t        \044(OOMPH-LIB_EXTERNAL_LIBS) \044(FLIBS) \n"  >> makefile.sample


#-------------------------------------------------------------------------

#Additional cleanup commands
CLEANFILES=makefile.sample my_demo_code
