Hlib installation instructions
====================

Hlib is a library for hierarchical matrices.

Unfortunately it is not available under a completely open source license so
we cannot redistribute the code. However the developers do allow a limited
license for academic use if you contact them, details are on their
[website](http://www.hlib.org/).

Because Hlib is a C library some changes need to be made before it can be
linked directly into C++ code. The patches in this folder contain the
required changes, along with some additional patches which contain code for
use in micromagnetics simulations (written by Andreas Knittel and Matteo
Franchin at the University of Southampton).

To build oomph-lib with hlib follow these instructions:

1. Apply the patches. Extract hlib into a folder then cd into that folder
   and run the command

     patch -p1 < "$OOMPH_DIR/external_distributions/hlib/cpp_hlib.patch"
    
  where $OOMPH_DIR is your oomph-lib directory.
  
2. Build hlib using the command 

     ./configure && make && sudo make install
     
   it should install the library to "/usr/local/lib" and the headers to
   "/usr/local/include/HLib".
   
3. Build oomph-lib as ususal (using autogen) with the additional configure
   option:

    --with-hlib
    
   (Configure options in oomph-lib are usually specified by modifying an
   options file in "$OOMPH_DIR/config/configure_options/".)

