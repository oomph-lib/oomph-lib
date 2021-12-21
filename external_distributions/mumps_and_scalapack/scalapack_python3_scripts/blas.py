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


class Blas(Frame):
    """ This class takes care of the BLACS installation. """
    def __init__(self):
        print('\n','='*40)
        print('  BLAS installation/verification')
        print('='*40)

        if(self.blaslib != ''):
            self.check_blas()
        else:
            if (self.downblas):
                self.down_install_blas()
            else:
                print("Need BLAS")
                sys.exit()
            




    def check_blas(self):

        print('Checking if provided BLAS works...', end=' ')
        # This function simply generates a FORTRAN program
        # that contains few calls to BLAS routine and then
        # checks if compilation, linking and execution are succesful

        sys.stdout.flush()
        writefile('tmpf.f',"""
      program ftest
      double precision da, dx(1)
      dx(1)=1
      da = 2
      call dscal(1,da,dx,1)
      stop
      end\n""")

        fcomm = self.mpif77+' -o tmpf '+'tmpf.f '+self.blaslib+' '+self.ldflags_f77
        (output, error, retz) = runShellCommand(fcomm)
        
        if(retz != 0):
            print('\n\nBLAS: provided BLAS cannot be used! aborting...')
            print('error is:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()
        
        comm = './tmpf'
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nBLAS: provided BLAS cannot be used! aborting...')
            print('error is:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()
        
        killfiles(['tmpf.f','tmpf'])
        print('yes')

        return 0;


    def down_install_blas(self):

        print("""
The reference BLAS library is being installed.
Don't expect high performance from this reference library!
If you want performance, you need to use an optimized BLAS library and,
to avoid unnecessary complications, if you need to compile this optimized BLAS
library, use the same compiler you're using here.""")
        sys.stdout.flush()

        savecwd = os.getcwd()

        # creating the build,lib and log dirs if don't exist
        if(not os.path.isdir(os.path.join(os.getcwd(),'build'))):
            os.mkdir(os.path.join(os.getcwd(),'build'))

        if(not os.path.isdir(os.path.join(os.getcwd(),'lib'))):
            os.mkdir(os.path.join(os.getcwd(),'lib'))

        if(not os.path.isdir(os.path.join(os.getcwd(),'log'))):
            os.mkdir(os.path.join(os.getcwd(),'log'))

        # chdir into the build directory            
        os.chdir(os.path.join(os.getcwd(),'build'))

        # Check if blas.tgz is already present in the working dir
        # otherwise download it
        if(not os.path.isfile(os.path.join(os.getcwd(),getURLName(self.blasurl)))):
            print('Downloading reference BLAS...', end=' ')
            downloader(self.blasurl,self.downcmd)
            print('done')

        # unzip and untar      
        print('Unzip and untar reference BLAS...', end=' ')
        comm = 'gunzip -f blas.tgz'
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nBLAS: cannot unzip blas.tgz')
            print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()
        
        comm = 'tar xf blas.tar'
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nBLAS: cannot untar blas.tgz')
            print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()
        os.remove('blas.tar')
        print('done')
     
        # change to BLAS dir
        os.chdir(os.path.join(os.getcwd(),'BLAS'))

        # compile and generate library
        print('Compile and generate reference BLAS...', end=' ')
        sys.stdout.flush()
        comm = self.mpif77+' '+self.fcflags+' -c *.f'
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nBLAS: cannot compile blas')
            print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()

        log = output+error

        comm = 'ar cr librefblas.a *.o'
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nBLAS: cannot create blas library')
            print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()
        print('done')

        # write the log on a file
        log = log+output+error
        fulllog = os.path.join(savecwd,'log/blaslog')
        writefile(fulllog, log)
        print('Installation of reference BLAS successful.')
        print('(log is in ',fulllog,')')

        # move librefblas.a to the lib directory        
        os.rename('librefblas.a',os.path.join(savecwd,'lib/librefblas.a'))

        # set framework variables to point to the freshly installed BLAS library
        self.blaslib  = os.path.join(savecwd,'lib/librefblas.a ')
        Frame.blaslib = os.path.join(savecwd,'lib/librefblas.a ')
        os.chdir(savecwd)
