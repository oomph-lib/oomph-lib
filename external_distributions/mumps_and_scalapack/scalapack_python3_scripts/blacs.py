#!/usr/bin/python

# -----------------------------------------
# ScaLAPACK installer
# University of Tennessee Knoxville
# October 16, 2007
# ----------------------------------------

from utils import writefile, runShellCommand, killfiles, downloader, getURLName
import sys
import os
from framework import Frame


class Blacs(Frame):
    """ This class takes care of the BLACS installation. """
    def __init__(self):
        print('\n','='*40)
        print('  BLACS installation/verification')
        print('='*40)

        if (self.downblacs):
            self.down_install_blacs()
        elif(self.blacslib != '' and self.blacsClib!='' and self.blacsF77lib !=''):
            self.check_blacs()
        else:
            print('need blacs')
            sys.exit()



    def check_blacs(self):
        """ Checks if BLACS libraries are good """

        # This function simply generates a FORTRAN program
        # that contains few calls to BLACS routine and then
        # checks if compilation, linking and execution are succesful
        print('Checking if provided BLACS works...', end=' ')
        sys.stdout.flush()
        writefile('tmpf.f',"""
      program ftest
      integer itmp, iam, nnodes
      call blacs_pinfo(iam, nnodes)
      if(nnodes.gt.0) then
         call blacs_get( 0, 0, itmp )
         call blacs_gridinit(itmp, 'c', 1, nnodes)
         call blacs_gridexit(itmp)
      end if
      stop
      end\n""")

        fcomm = self.mpif77+' '+self.ldflags_f77+' -o tmpf '+'tmpf.f '+self.blacsF77lib+' '+self.blacslib+' '+self.blacsClib
        (output, error, retz) = runShellCommand(fcomm)
        
        if(retz != 0):
            print('\n\nBLACS: provided BLACS cannot be used! aborting...')
            print('error is:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()
        
        comm = 'date'
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nBLACS: provided BLACS cannot be used! aborting...')
            print('error is:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()
        
        killfiles(['tmpf.f','tmpf'])
        print('yes')

        return 0;



    def down_install_blacs(self):
        """ Downloads and installs BLACS libraries """

        savecwd = os.getcwd()

        # creating the build and lib dirs if don't exist
        if(not os.path.isdir(os.path.join(os.getcwd(),'build'))):
            os.mkdir(os.path.join(os.getcwd(),'build'))

        if(not os.path.isdir(os.path.join(os.getcwd(),'lib'))):
            os.mkdir(os.path.join(os.getcwd(),'lib'))

        if(not os.path.isdir(os.path.join(os.getcwd(),'log'))):
            os.mkdir(os.path.join(os.getcwd(),'log'))

        # chdir into the build directory            
        os.chdir(os.path.join(os.getcwd(),'build'))
        # checks if mpiblacs.tgz is already present in the current directory
        # and if not it downloads it
        if(not os.path.isfile(os.path.join(os.getcwd(),getURLName(self.blacsurl)))):
            print('Downloading BLACS...', end=' ')
            downloader(self.blacsurl, self.downcmd)
            print('done')

        # unzip and untar            
        print('Unzip and untar BLACS...', end=' ')
        comm = 'gunzip -f mpiblacs.tgz'
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nBLACS: cannot unzip mpiblacs.tgz')
            print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()
        

        comm = 'tar xf mpiblacs.tar'
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nBLACS: cannot untar mpiblacs.tgz')
            print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()
        os.remove('mpiblacs.tar')
        print('done')
        
        # goes in the BLACS dir
        os.chdir(os.path.join(os.getcwd(),'BLACS'))

        # set TRANSCOMM
        self.set_transcomm()
        # write Bmake.inc file
        self.write_bmake()

        # compile
        print('Compiling BLACS...', end=' ')
        sys.stdout.flush()
        comm = self.make+' mpi'
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nBLACS: error building BLACS')
            print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()
        
        # write a log file containing the BLACS compilation output
        fulllog = os.path.join(savecwd,'log/blacslog')
        writefile(fulllog, output+error)
        print('done.')
        print(' (log is in ',fulllog,')')

        os.rename('LIB/blacs.a',os.path.join(savecwd,'lib/blacs.a'))
        os.rename('LIB/blacsC.a',os.path.join(savecwd,'lib/blacsC.a'))
        os.rename('LIB/blacsF77.a',os.path.join(savecwd,'lib/blacsF77.a'))

        # sets framework variables to point to the freshly installed BLACS libraries
        self.blacslib     = os.path.join(savecwd,'lib/blacs.a ')
        self.blacsClib    = os.path.join(savecwd,'lib/blacsC.a ')
        self.blacsF77lib  = os.path.join(savecwd,'lib/blacsF77.a ')
        Frame.blacslib    = os.path.join(savecwd,'lib/blacs.a ')
        Frame.blacsClib   = os.path.join(savecwd,'lib/blacsC.a ')
        Frame.blacsF77lib = os.path.join(savecwd,'lib/blacsF77.a ')
        os.chdir(savecwd)






    def set_transcomm(self):
        """ Sets the TRANSCOMM variable in Bmake.inc """
        # This one compiles few programs equivalent to those in BLACS/INSTALL
        # in order to set the TRANSCOMM value

        print('Setting TRANSCOMM...', end=' ') 
        sys.stdout.flush()

        writefile('tmpf.f',"""
      program tctst
      include 'mpif.h'
      integer i, ierr
      external Ccommcheck
      integer  Ccommcheck

      call mpi_init(ierr)

      i = Ccommcheck(MPI_COMM_WORLD, MPI_COMM_SELF)
      if (i .ne. 0) then
         print*,'-DCSameF77'
      end if

      call mpi_finalize(ierr)

      stop
      end\n""")

        writefile('tmpc.c',"""
#include <mpi.h>

int Ccommcheck(int F77World, int f77comm){
   int Np, Iam, i, OK=1;

   if (sizeof(int) != sizeof(MPI_Comm)) OK=0;
   else if ((MPI_Comm) F77World != MPI_COMM_WORLD) OK=0;
   else {
      MPI_Comm_rank(MPI_COMM_WORLD, &Iam);
      i = MPI_Comm_size((MPI_Comm) f77comm, &Np);
      if (i != MPI_SUCCESS) OK = 0;
      else if (Np != 1) OK = 0;
   }
   return(OK);
}
int CCOMMCHECK(int *F77World, int *f77comm){ return(Ccommcheck(*F77World, *f77comm));}
int ccommcheck_(int *F77World, int *f77comm){return(Ccommcheck(*F77World, *f77comm));}
int ccommcheck(int *F77World, int *f77comm){return(Ccommcheck(*F77World, *f77comm));}\n""")
        
        ccomm = self.mpicc+' '+self.ccflags+' -c tmpc.c -o tmpc.o'
        fcomm = self.mpif77+' '+self.fcflags+' tmpf.f tmpc.o -o xtc_CsameF77'

        (output, error, retz) = runShellCommand(ccomm)
        if(retz != 0):
            print('\n\nBLACS: Error in transcomm setting! cannot compile')
            print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()

        (output, error, retz) = runShellCommand(fcomm)
        if(retz != 0):
            print('\n\n1BLACS: Error in transcomm setting! cannot compile')
            print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()

        comm = os.path.join(os.getcwd(),'xtc_CsameF77')
        (output, error, retz) = runShellCommand(comm)
        #if(retz != 0):
        #    print '\n\nBLACS: Error in transcomm setting! cannot run xtc_CsameF77'
        #    print 'stderr:\n','*'*40,'\n',error,'\n','*'*40
        #    sys.exit()

        killfiles(['tmpf.f', 'tmpc.c', 'tmpc.o', 'xtc_CsameF77'])
        
        self.transcomm = output


        if(self.transcomm == ''):
            writefile('tmpc.c',"""
            #include <stdio.h>
            #include <mpi.h>
            main(){
               MPI_Comm ccomm;
               int fcomm;
               extern void *MPIR_ToPointer();
               extern int   MPIR_FromPointer();
               extern void *MPIR_RmPointer();
            
               if (sizeof(int) < sizeof(int*))   {
                  fcomm = MPIR_FromPointer(MPI_COMM_WORLD);
                  ccomm = (MPI_Comm) MPIR_ToPointer(fcomm);
                  if (ccomm == MPI_COMM_WORLD)
                     fprintf(stdout,\" -DUseMpich -DPOINTER_64_BITS=1\");
               }
               return 0;
            }\n""")
        
            ccomm = self.mpicc+' '+self.ccflags+' tmpc.c -o UseMpich'
            (output, error, retz) = runShellCommand(ccomm)
            if(retz != 0):
                print(self.transcomm)
                return 1;
            
            comm = os.path.join(os.getcwd(),'UseMpich')
            (output, error, retz) = runShellCommand(comm)
            if(retz != 0):
                print('\n\nBLACS: Error in transcomm setting! cannot run UseMpich')
                print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
                sys.exit()

            killfiles(['tmpc.c','UseMpich'])
            self.transcomm = output

# Testing MPI2 implementation
        if(self.transcomm == ''):
	            writefile('tmpc.c',"""
	            #include <stdio.h>
	            #include <mpi.h>
				int main(int argc, char** argv)
				{
				  int version, subversion;
				  MPI_Init(&argc, &argv);
				  MPI_Get_version(&version, &subversion);
				  if (version==2) fprintf(stdout,\" -DUseMpi2\");
				  MPI_Finalize();
				  return 0;
	            }\n""")

	            ccomm = self.mpicc+' '+self.ccflags+' tmpc.c -o UseMPI2'
	            (output, error, retz) = runShellCommand(ccomm)
	            if(retz != 0):
	                print(self.transcomm)
	                print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
	                sys.exit()

	            comm = os.path.join(os.getcwd(),'UseMPI2')
	            (output, error, retz) = runShellCommand(comm)
	            if(retz != 0):
	                print('\n\nBLACS: Error in transcomm setting! cannot run UseMPI2')
	                print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
	                sys.exit()

	            killfiles(['tmpc.c','UseMPI2'])
	            self.transcomm = output
	            
        print(self.transcomm)

            
        return 1;
        

    def write_bmake(self):
      """ Writes the Bmake.inc file for BLACS installation """
      print('Writing Bmake.inc...', end=' ')
      sys.stdout.flush()
      writefile('Bmake.inc',"""
SHELL = /bin/sh
BTOPdir = """+os.getcwd()+"""
COMMLIB = MPI
PLAT = 
BLACSdir    = $(BTOPdir)/LIB
BLACSDBGLVL = 0
BLACSFINIT  = $(BLACSdir)/blacsF77.a
BLACSCINIT  = $(BLACSdir)/blacsC.a
BLACSLIB    = $(BLACSdir)/blacs.a

MPIINCdir = """+self.mpiincdir+"""
MPILIB = 

BTLIBS = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT) $(MPILIB)
INSTdir = $(BTOPdir)/INSTALL/EXE
TESTdir = $(BTOPdir)/TESTING/EXE
FTESTexe = $(TESTdir)/xFbtest_$(COMMLIB)-$(PLAT)-$(BLACSDBGLVL)
CTESTexe = $(TESTdir)/xCbtest_$(COMMLIB)-$(PLAT)-$(BLACSDBGLVL)
SYSINC = -I$(MPIINCdir)
INTFACE = """+self.mangling+"""
SENDIS = 
BUF = 
TRANSCOMM = """+self.transcomm+"""
WHATMPI = 
SYSERRORS = 
DEBUGLVL = -DBlacsDebugLvl=$(BLACSDBGLVL)
DEFS1 = -DSYSINC $(SYSINC) $(INTFACE) $(DEFBSTOP) $(DEFCOMBTOP) $(DEBUGLVL)
BLACSDEFS = $(DEFS1) $(SENDIS) $(BUFF) $(TRANSCOMM) $(WHATMPI) $(SYSERRORS)

F77            = """+self.mpif77+"""
F77NO_OPTFLAGS = """+self.noopt+"""
F77FLAGS       = $(F77NO_OPTFLAGS) """+self.fcflags+"""
F77LOADER      = $(F77)
F77LOADFLAGS   = """+self.ldflags_f77+"""
CC             = """+self.mpicc+"""
CCFLAGS        = """+self.ccflags+"""
CCLOADER       = $(F77)
CCLOADFLAGS    = """+self.ldflags_c+"""
ARCH      = ar
ARCHFLAGS = r
RANLIB    = """+self.ranlib+"""
      """)
      print('done.')
      return 1;

