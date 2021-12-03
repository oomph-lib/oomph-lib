#!/usr/bin/python

# -----------------------------------------
# ScaLAPACK installer
# University of Tennessee Knoxville
# October 16, 2007
# ----------------------------------------

from utils import runShellCommand, ynask, reverse, writefile, killfiles, downloader, getURLName, fixpaths
import sys
import os
import getopt
import string

class Frame:
    """ This class takes care of the BLACS installation. """
    
    # set default values
    prefix = ''
    ccflags, fcflags, noopt = ['-O3','-O3','']
    mpibindir = ''               # directory containing mpicc and mpif77
    mpiincdir = ''               # directory containing mpi.h
    mpicc, mpif77 = ['', '']     # the mpicc and mpif77 commands
    make         = 'make'        # the "make" command
    plat         = 'LINUX'       
    mangling     = ''            # sets the mangling type for Fortran to C calls
    blaslib      = ''            # the blas library
    blacslib     = ''            # the blacs library
    blacsClib    = ''            # the blacs C library
    blacsF77lib  = ''            # the blacs F77 library
    lapacklib    = ''            # lapack library
    scalapacklib = ''            # scalapack library
    proxy        = 0
    downblas     = 0             # whether or not to download reference blas
    downblacs    = 1             # whether or not to download reference blacs
    downlapack   = 0             # whether or not to download reference lapack
    ldflags_c    = ''            # loader flags when main program is in C
    ldflags_f77  = ''            # loader flags when main program is in Fortran
    testing      = 1             # whether on not to compile and run LAPACK and ScaLAPACK test programs
    downcmd      = ''            # the command used to download stuff
    blasurl      = 'http://netlib.org/blas/blas.tgz'
    lapackurl    = 'http://netlib.org/lapack/lapack.tgz'
    blacsurl     = 'http://netlib.org/blacs/mpiblacs.tgz'
    scalapackurl = 'http://netlib.org/scalapack/scalapack-1.8.0.tgz'
    ranlib       = ''            # the ranlib command
    clean        = 0

    def __init__(self, argv):
        
        print('='*40)
        print('Setting up the framework\n')


        # parse input arguments
        self.parse_args(argv)
        
        if(Frame.clean == 1):
            self.cleanup()
            sys.exit()

        if(Frame.prefix != ''):
	        if(not os.path.isdir(os.path.join(os.getcwd(),Frame.prefix))):
	           print("Creating directory",Frame.prefix)
	           os.mkdir(os.path.join(os.getcwd(),Frame.prefix))
	        os.chdir(Frame.prefix)
	
        if(Frame.mpicc=='' or Frame.mpif77==''):
            #  if no MPI compilers are provided, look for them in mpibindir or in PATH
            if(not self.look_for_mpibinaries()):
                print("""
MPI C and Fortran77 compilers are needed  to complete
the installation. You can specify the using the --mpibindir
or the --mpicc and --mpif77 flags.""")
                sys.exit()
                
        if(Frame.mpiincdir==''):
            # the path to the mpi.h file is needed
            if(not self.look_for_mpih()):
                print('Please provide the location for mpi.h file using the --mpiincdir flag')
                sys.exit()

        # if testing has to be done, BLAS, BLACS andLAPACK are required
        if(Frame.testing == 1):
            if(Frame.blaslib == '' and Frame.downblas == 0):
                print("""
Please provide a working BLAS library. If a BLAS library
is not present on the system, the reference BLAS library it can be
automatically downloaded and installed by adding the --downblas flag.
Be aware that a reference BLAS library will be installed with the --downblas
flag so don't expect performance.
The BLAS library is not needed in the case where testing is disabled
by means of the --notesting flag. 
                """)
                sys.exit()
            if((Frame.blacslib == '' or Frame.blacsClib == '' or Frame.blacsF77lib == '') and Frame.downblacs == 0):
                print("""
Please provide a working BLACS library. If a BLACS library
is not present on the system, it can be automatically downloaded and
installed by adding the --downblacs flag. The BLAS library is not
needed in the case where testing is disabled by means of the
--notesting flag.
                """)
                sys.exit()
            if(Frame.lapacklib == '' and Frame.downlapack == 0):
                print("""
Please provide a working LAPACK library. If a LAPACK library
is not present on the system, it can be automatically downloaded and
installed by adding the --downlapack flag. The BLAS library is not
needed in the case where testing is disabled by means of the
--notesting flag.
                """)
                sys.exit()
        
          
        # CUSTOM CHECKS
        self.check_mpicc()
        self.check_mpif77()
        self.set_mangling()
        self.set_download()
        self.set_ranlib()
        self.detect_compilers()
        self.check_linking()
        
        return


    def usage(self):
          print('='*40)
          print("""
          ScaLAPACK configuration script
          The script will try to figure out some of the features of your system
          like the mpi compiler and the location of few other tools required for
          the installation of the ScaLAPACK package and the other software
          packages that it requires.
          It is strongly recommended that you help the script by providing info 
          through the flags listed below

          
          -h or --help         : prints this message
 
          --prefix             : path where to create the libraries, build and log of the installer

          --mpibindir=[DIR]    : the path to where the mpi
                                 binaries (mpicc and mpif77)
                                 are contained

          --mpicc=[CMD]        : the mpicc command.
                                 (default: will take the first available mpicc in the path)

          --mpif77=[CMD]       : the mpif77 command.
                                 (default: will take the first available mpif77 in the path)

          --mpiincdir=[DIR]    : the path to the directory containing mpi.h
                                 (default: will take the corresponding mpi.h to the mpicc/mpif77 found)

          --ccflags=[FLAGS]    : the flags for the C compiler
                                 (default -O3)

          --fcflags=[FLAGS]    : the flags for the F77 compiler
                                 (default -O3)

          --noopt=[FLAGS]      : compilation flags to be used
                                 on machine dependent subroutines
                                 in LAPACK and ScaLAPACK.
                                 See LAPACK and ScaLAPACK documentation.

          --ldflags_c=[flags]  : loader flags when main program is in C. Some compilers (e.g. PGI) require 
                                 different options when linking C main programs to Fortran subroutines
                                 and vice-versa

          --ldflags_f77=[flags]: loader flags when main program is in Fortran. Some compilers (e.g. PGI) require 
                                 different options when linking Fortran main programs to C subroutines
                                 and vice-versa

          --makecmd=[CMD]      : the make command
                                 (make by default)

          --blacslib=[LIB]     : the BLACS library
                                 (path should be absolute if --prefix is used)

          --blacsclib=[LIB]    : the BLACS C interface library
                                 (path should be absolute if --prefix is used)

          --blacsf77lib=[LIB]  : the BLACS F77 interface library
                                 (path should be absolute if --prefix is used)

          --blaslib=[LIB]      : a BLAS library
                                 (path should be absolute if --prefix is used)

          --lapacklib=[LIB]    : a LAPACK library
                                 (path should be absolute if --prefix is used)

          --downblas           : if you do not want to provide a BLAS
                                 we can download and install it for you

          --downblacs          : if you do not want to provide a BLACS
                                 we can download and install it for you

          --downlapack         : if you do not want to provide a LAPACK
                                 we can download and install it for you

          --notesting          : disables the ScaLAPACK testing. The
                                 BLAS, BLACS and LAPACK libraries are not
                                 required in this case.

          --clean              : cleans up the installer directory.

Note: If you want use a proxy for downloading, the http_proxy environment variable has to be set.
          """)
          print('='*40)
          sys.exit()
      


    def parse_args(self, argv):
        """ Parse input argument to get compilers etc. from command line. """

        try:
          opts, args = getopt.getopt(sys.argv[1:], "h", ["help","prefix=",
          "ccflags=", "fcflags=", "noopt=", "makecmd=", "mpibindir=", "mpiincdir=", "blacslib=",
          "blacsclib=", "blacsf77lib=", "blaslib=", "lapacklib=", "ldflags_c=", "ldflags_f77=","mpicc=","mpif77=",
          "downblas", "downblacs", "downlapack", "notesting","clean"])
        except getopt.error as msg:
          print(msg)
          print("for help use --help")
          sys.exit(2)
        # process options
        for o, a in opts:
          if o in ("-h", "--help"):
            self.usage()
            sys.exit(0)
          else:
            if(o == '--clean'):
                Frame.clean = 1
                return;
            elif(o == '--prefix'):
                Frame.prefix = fixpaths(a)
                print('Install directory is...', Frame.prefix)
            elif(o == '--ccflags'):
                Frame.ccflags = a
                print('C flags are ', a)
            elif(o=='--fcflags'):
                Frame.fcflags = a
                print('Fortran flags are ', a)
            elif(o=='--noopt'):
                Frame.noopt = a
                print('NOOPT flags are ', a)
            elif(o=='--makecmd'):
                Frame.make = a
            elif(o=='--mpibindir'):
                Frame.mpibindir = fixpaths(a)
                print('MPI bin dir is ', Frame.mpibindir)
            elif(o=='--mpiincdir'):
                Frame.mpiincdir = fixpaths(a)
                print('MPI include dir is ', Frame.mpiincdir)
            elif(o=='--mpicc'):
                Frame.mpicc = a
                print('mpicc is ', a)
            elif(o=='--mpif77'):
                Frame.mpif77 = a
                print('mpif77 is ', a)
            elif(o == '--blacslib'):
                Frame.blacslib = fixpaths(a)
                Frame.downblacs = 0
            elif(o == '--blacsclib'):
                Frame.blacsClib = fixpaths(a)
                Frame.downblacs = 0
            elif(o == '--blacsf77lib'):
                Frame.blacsF77lib = fixpaths(a)
                Frame.downblacs = 0
            elif (o == '--blaslib'):
                Frame.blaslib = fixpaths(a)
            elif (o == '--lapacklib'):
                Frame.lapacklib = fixpaths(a)
            elif (o == '--downblas'):
                Frame.downblas = 1
            elif (o == '--downblacs'):
                Frame.downblacs = 1
            elif (o == '--downlapack'):
                Frame.downlapack = 1
            elif (o == '--notesting'):
                Frame.testing = 0
            elif (o == '--ldflags_c'):
                Frame.ldflags_c = a
            elif (o == '--ldflags_f77'):
                Frame.ldflags_f77 = a

 
    # look for MPI
    def look_for_mpibinaries(self):
        """ looks for MPI compilers in mpibindir or in PATH """
        # This function is only able to find mpicc/mpif77. If the name of the 
        # compilers are different (like mpixlc or mpixlf) the user must explicitly provide them
        # with the related flags        
        if(Frame.mpibindir != ''):
            # look for mpicc and mpif77 in mpibindir  
            if(os.path.isfile(os.path.join(Frame.mpibindir,'mpicc'))):
                Frame.mpicc = os.path.join(Frame.mpibindir,'mpicc')
                print('mpicc :',Frame.mpicc)
            else:
                print('Could not find mpicc in ', Frame.mpibindir)
            if(os.path.isfile(os.path.join(Frame.mpibindir,'mpif77'))):
                Frame.mpif77 = os.path.join(Frame.mpibindir,'mpif77')
                print('mpif77 :',Frame.mpif77)
            else:
                print('Could not find mpif77 in ',Frame.mpibindir)

        # is mpicc and mpif77 haven't been found
        if(Frame.mpicc=='' and Frame.mpif77==''):

            path=str(os.getenv('PATH')).split(os.pathsep)
            print('Looking for MPI binaries...', end=' ')
            sys.stdout.flush()
            for i in path:
                mpicc  = os.path.join(i,'mpicc')
                mpif77 = os.path.join(i,'mpif77')
                if (Frame.mpicc == '' and os.path.isfile(mpicc)):
                    Frame.mpicc = mpicc
                if (Frame.mpif77 == '' and os.path.isfile(mpif77)):
                    Frame.mpif77 = mpif77
            if (Frame.mpicc!='' or Frame.mpif77!=''):
                print("mpicc and mpif77 found.\nmpicc : "+Frame.mpicc+"\nmpif77 :"+Frame.mpif77)
                return 1;
                    
            elif (Frame.mpicc=='' or Frame.mpif77==''):
                print("""
The information about the location of MPI commands is incomplete.
Please, either provide the path to the directory containing mpicc
and mpif77 with the --mpibindir flag or the full path to both 
commands with the --mpicc AND --mpif77 flags. In case none of these
flags are provided, the installer will look for mpicc and mpif77
in the PATH environment variable.
""")
                return 0;

        return 1;




    def look_for_mpih(self):
        """ looks for mpi.h close to mpibindir """
        import re
      
        print('Looking for mpi.h...', end=' ')
        for i in [Frame.mpibindir[0:Frame.mpibindir.find('bin')], Frame.mpicc[0:Frame.mpicc.find('bin')], Frame.mpif77[0:Frame.mpif77.find('bin')]]:
            sys.stdout.flush()
            tmp = i+'include'
            if(os.path.isfile(os.path.join(tmp,'mpi.h'))):
                Frame.mpiincdir = tmp
                print('found in '+i+'include')
                return 1;

        print('not found')

        return 0;


            
    def check_mpicc(self):
        """ checks if mpicc works """
        # simply generates a C program containing a couple of calls
        # to MPI routines and checks if the compilation and execution
        # are succesful
        print('Checking if mpicc works...', end=' ')
        sys.stdout.flush()
        # generate
        writefile('tmpc.c',"""
            #include \"mpi.h\"
            #include <stdio.h>
            int main(int argc, char **argv){
            int iam;
            MPI_Init( &argc, &argv );
            MPI_Comm_rank( MPI_COMM_WORLD, &iam );
            if(iam==0){fprintf(stdout, \"success\" );fflush(stdout);}
            MPI_Finalize();
            return 0;
            }\n""")

        # compile
        ccomm = Frame.mpicc+' -o tmpc '+os.path.join(os.getcwd(),'tmpc.c')
        (output, error, retz) = runShellCommand(ccomm)
      
        if(retz != 0):
            print('\n\nCOMMON: mpicc not working! aborting... retz = ',retz)
            print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()



        # run
        comm = 'date'
        print("before runShellcommand for date for mpicc")
        (output, error, retz) = runShellCommand(comm)
        print("back from runShellcommand for date for mpicc")
        if(retz != 0):
            print('\n\nCOMMON: mpicc not working! aborting...')
            print('error is:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()
        print("date for mpicc")

        # cleanup
        killfiles(['tmpc.c','tmpc'])
        print('yes')
        print("about to leave after mpicc")

        return 0;



    def check_mpif77(self):
        """ check if mpif77 works """
        # simply generates a F77 program containing a couple of calls
        # to MPI routines and checks if the compilation and execution
        # are succesful
        print('Checking if mpif77 works...', end=' ')
        sys.stdout.flush()
        # generate        
        writefile('tmpf.f',"""
      program ftest
      include 'mpif.h'
      integer Iam, i, ierr
      call mpi_init(ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, Iam, ierr)
      if (Iam .eq. 0) then
            print*,'success'
      end if
      call mpi_finalize(ierr)
      stop
      end\n""")

        # compile
        fcomm = Frame.mpif77+' -o tmpf '+'tmpf.f'
        (output, error, retz) = runShellCommand(fcomm)

        if(retz != 0):
            print('\n\nCOMMON: mpif77 not working! aborting...')
            print('error is:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()
        
        # run
        comm = 'date'
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nCOMMON: mpif77 not working! aborting...')
            print('error is:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()

        # cleanup        
        killfiles(['tmpf.f','tmpf','tmpf.o'])
        print('yes')

        return 0;
        

    def set_mangling(self):
        """ Sets the INTFACE variable in Bmake.inc """
        # This one generates a program equivalent to that in BLACS/INSTALL
        # that checks the mangling in FORTRAN function symbols
        sys.stdout.flush()
        writefile('tmpf.f',"""
      program intface
      external c_intface
      integer i
      call c_intface(i)
      stop
      end\n""")
        writefile('tmpc.c',"""
      #include <stdio.h>
      void c_intface_(int *i){fprintf(stdout, \"-DAdd_\");fflush(stdout);}
      void c_intface(int *i){fprintf(stdout, \"-DNoChange\");fflush(stdout);}
      void c_intface__(int *i){fprintf(stdout, \"-Df77IsF2C\");fflush(stdout);}
      void C_INTFACE(int *i){fprintf(stdout, \"-DUpCase\");fflush(stdout);}\n""")

        ccomm = Frame.mpicc+' '+Frame.ccflags+' -c tmpc.c -o tmpc.o'
        fcomm = Frame.mpif77+' '+Frame.fcflags+' tmpf.f tmpc.o -o xintface'

        (output, error, retz) = runShellCommand(ccomm)
        if(retz != 0):
            print('\n\nCOMMON: in set_mangling: cannot compile')
            print('error is:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()

        (output, error, retz) = runShellCommand(fcomm)
        if(retz != 0):
            print('\n\nCOMMON: in set_mangling: cannot compile')
            print('error is:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()

        comm = os.path.join(os.getcwd(),'xintface')
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nCOMMON: in set_mangling: cannot run xintface')
            print('error is:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()

        Frame.mangling = str(output)
        killfiles(['xintface', 'tmpf.f', 'tmpf.o', 'tmpc.c', 'tmpc.o'])
        
        return 1;


    def check_linking(self):
        """ Check if C main can be linked to F77 subroutine """

        # This one checks if the linking command works out of the box or
        # if any specific flag is required. For example if the linker if the
        # Intel FORTRAN compiler, then the "-nofor_main" is usually required.
        # This function only checks if linker works but does not automatically
        # detect the required flags
        print('Checking loader...', end=' ')
        sys.stdout.flush()
        writefile('tmpf.f',"""
      subroutine fsub()
      write(*,*)'success'
      stop
      end\n""")
        writefile('tmpc.c',"""
      #if defined Add_
      #define fsub fsub_
      #elif defined NoChange
      #define fsub fsub
      #elif defined f77IsF2C
      #define fsub fsub_
      #elif defined UpCase
      #define fsub FSUB
      #endif
      void main(){
      fsub();}\n""")

        ccomm = Frame.mpicc+' '+Frame.ccflags+' '+str(Frame.mangling)+' -c -o tmpc.o tmpc.c'
        fcomm = Frame.mpif77+' '+Frame.fcflags+' -c -o tmpf.o tmpf.f'
        lcomm = Frame.mpif77+' '+Frame.ldflags_c+' -o lnk tmpf.o tmpc.o'
        (output, error, retz) = runShellCommand(ccomm)
        if(retz != 0):
            print('\n\nCOMMON: in check_linking: cannot compile')
            print('command is: ',ccomm)
            print('error is:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()

        (output, error, retz) = runShellCommand(fcomm)
        if(retz != 0):
            print('\n\nCOMMON: in check_linking: cannot compile')
            print('command is: ',fcomm)
            print('error is:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()

        (output, error, retz) = runShellCommand(lcomm)
        if(retz != 0):
            print("""\n\nCOMMON: in check_linking: cannot link
            Cannot link a C main program to a Fortran77 subroutine
            Make sure that the appropriate flags are passed to the linker.""")
            print('command is: ',lcomm)
            print('error is:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()


        killfiles(['lnk', 'tmpf.f', 'tmpf.o', 'tmpc.c', 'tmpc.o'])
        
        print('works')
        return 1;


          
    def set_download(self):
        """ Figures out how to download files """
        print('Setting download command...')
        wget = 0
        urllib = 0
        # look for proxy setting in the environnement variables
        proxy_http=str(os.getenv('http_proxy'))
        # if proxy_http found in the environnement, we are going to use wget as download command
        if (proxy_http != 'None'):
            Frame.proxy=1
            urllib=0
        else:    
            try:
                # check if the urllib2 module is included in the
                # python installation. If yes files are downloaded with urllib2            
                print("Checking availability of urllib...", end=' ')
                import urllib.request, urllib.error, urllib.parse
                urllib=1
                print("available")
                print("Testing urllib...", end=' ')
                try:
#                     name = getURLName('http://www.netlib.org/lapack/index')
                    url = urllib.request.urlopen('http://www.netlib.org/lapack/index')
                    f = open(name,'w')
                    for line in url.readlines():
                        f.write(line)
                    url.close()
                    f.close()
                except:
                    print("not working")
                    urllib = -1
                else:
                    print("working")
                    Frame.downcmd = 'urllib2'
                    os.remove('index')
                    return
            except ImportError:
                print("not available")
                urllib=0

        if(urllib <= 0):
            # if urllib2 is not present checks if wget is present
            # in the PATH and if yes it sets the download command
            # to be wget
            print("Checking availablility of wget...", end=' ')
            path=str(os.getenv('PATH')).split(os.pathsep)
            for i in path:
                if (os.path.isfile(os.path.join(i,'wget'))):
                    print("available")
                    wget = 1
                    break

            if(wget == 1):
                # test wget
                if (proxy_http != 'None'):
                    print("Testing wget through the "+proxy_http+" proxy...", end=' ')
                else:
                    print("Testing wget...", end=' ')
                try:
                    comm = 'wget http://www.netlib.org/lapack/index'
                    (output, error, retz) = runShellCommand(comm)
                    if(retz != 0):
                        raise
                except:
                    print('not working.')
                    wget = -1
                else:
                    print('working')
                    Frame.downcmd='wget'
                    os.remove('index')
                    return
            else:
                # wget not available
                print("not available")
                wget=0


    def set_ranlib(self):
        """ Sets the ranlib command """
        # Some systems don't have the ranlib command (e.g. SGIs). 
        # In the case where ranlib is not present in the PATH,
        # echo is used instead of ranlib        
        print("Setting ranlib command...", end=' ')

        path=str(os.getenv('PATH')).split(os.pathsep)
        for i in path:
            if (os.path.isfile(os.path.join(i,'ranlib'))):
                Frame.ranlib=os.path.join(i,'ranlib')
                print(Frame.ranlib)
                return

        for i in path:
            if (os.path.isfile(os.path.join(i,'echo'))):
                Frame.ranlib=os.path.join(i,'echo')
                print(Frame.ranlib)
                return





    def detect_compilers(self):
        """ Tries to detect the compilers type """
        # By users experience it is known which compiler flags are required
        # in some cases. This function tries to detect which compilers are used
        # and sets the flags accordingly

        print('Detecting Fortran compiler...', end=' ')
        if (self.fc_is_intel()):
            # The Intel FORTRAN compiler requires -nofor_main flag
            # for the linking and the -mp flag to maintain the 
            # floating-point precision
            Frame.ldflags_c   += ' -nofor_main'
            Frame.ldflags_f77 += ''
            Frame.noopt += ' -mp'
            print('Intel')
        elif (self.fc_is_gnu()):
            print('GNU')
        elif (self.fc_is_pgi()):
            Frame.ldflags_c   += ' -Mnomain'
            Frame.ldflags_f77 += ''
        else:
            print('unknown')


        print('Detecting C compiler...', end=' ')
        if (self.cc_is_intel()):
            print('Intel')
        elif (self.cc_is_gnu()):
            print('GNU')
        elif (self.cc_is_pgi()):
            print('PGI')
        else:
            print('unknown')


        print('Selected C compiler flags: '+self.ccflags)
        print('Selected Fortran77 compiler flags: '+self.fcflags)
        print('Selected loader flags (C main): '+self.ldflags_c)
        print('Selected loader flags (F77 main): '+self.ldflags_f77)
        print('Selected NOOPT flags: '+self.noopt)

        return


    def fc_is_intel(self):
        comm = self.mpif77+' -V'
        (output, error, retz) = runShellCommand(comm)
        isifort = error.find('Intel(R) Fortran Compiler')
        if (isifort != -1):
            return 1

        return 0


    def fc_is_gnu(self):
        comm = self.mpif77+' --help'
        (output, error, retz) = runShellCommand(comm)
        isifort = error.find('gnu.org')
        if (isifort != -1):
            return 1

        return 0

    def fc_is_pgi(self):
        comm = self.mpif77+' -V'
        (output, error, retz) = runShellCommand(comm)
        isifort = error.find('pgif77')
        if (isifort != -1):
            return 1
        isifort = error.find('pgif90')
        if (isifort != -1):
            return 1
        isifort = error.find('pgif95')
        if (isifort != -1):
            return 1
        isifort = error.find('Portland')
        if (isifort != -1):
            return 1
        isifort = output.find('pgif77')
        if (isifort != -1):
            return 1
        isifort = output.find('pgif90')
        if (isifort != -1):
            return 1
        isifort = output.find('pgif95')
        if (isifort != -1):
            return 1
        isifort = output.find('Portland')
        if (isifort != -1):
            return 1

        return 0



    def cc_is_intel(self):
        comm = self.mpicc+' -V'
        (output, error, retz) = runShellCommand(comm)
        isifort = error.find('Intel(R) C Compiler')
        if (isifort != -1):
            return 1

        return 0


    def cc_is_gnu(self):
        comm = self.mpicc+' --help'
        (output, error, retz) = runShellCommand(comm)
        isifort = error.find('gnu.org')
        if (isifort != -1):
            return 1

        return 0
        

    def cc_is_pgi(self):
        comm = self.mpicc+' -V'
        (output, error, retz) = runShellCommand(comm)
        isifort = error.find('pgicc')
        if (isifort != -1):
            return 1
        isifort = error.find('Portland')
        if (isifort != -1):
            return 1
        isifort = output.find('pgicc')
        if (isifort != -1):
            return 1
        isifort = output.find('Portland')
        if (isifort != -1):
            return 1

        return 0


    def resume(self):
        
        cwd = os.getcwd()
        print("""
******************************************************
******************************************************

ScaLAPACK installation completed.


Your BLAS library is:\n"""+self.blaslib+"""
\nYour BLACS libraries are:\n"""+self.blacslib+'\n'+self.blacsClib+'\n'+self.blacsF77lib+"""
\nYour LAPACK library is:\n"""+self.lapacklib+"""\n
Your ScaLAPACK library is:\n"""+self.scalapacklib+"""

Log messages are in the\n"""+os.path.join(cwd,'log')+"""
directory.

The Scalapack testing programs are in:\n"""+os.path.join(cwd,'build/scalapack-1.8.0/TESTING')+"""


The\n"""+os.path.join(cwd,'build')+"""
directory contains the source code of the libraries 
that have been installed. It can be removed at this time.


******************************************************
******************************************************
""")



    def cleanup(self):
        """ Cleans up the installer directory """

        print("Cleaning up...", end=' ')
        sys.stdout.flush()
        
        cwd = os.getcwd()
        builddir = os.path.join(cwd,'build')
        libdir   = os.path.join(cwd,'lib')
        logdir   = os.path.join(cwd,'log')

        comm = 'rm -rf '+builddir+' '+libdir+' '+logdir
        (output, error, retz) = runShellCommand(comm)
        
        print("done.")
