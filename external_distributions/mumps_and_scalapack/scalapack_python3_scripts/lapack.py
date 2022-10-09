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


class Lapack(Frame):
    """ This class takes care of the LAPACK installation. """
    def __init__(self):
        print('\n','='*40)
        print('  Lapack installation/verification')
        print('='*40)

        if(self.lapacklib != ''):
            self.check_lapack()
        else:
            if (self.downlapack):
                self.down_install_lapack()
            else:
                print("Need LAPACK")
                sys.exit()





    def check_lapack(self):
        """ Checks if the provided LAPACK library is good """
        print('Checking if provided LAPACK works...', end=' ')

        # This function simply generates a FORTRAN program
        # that contains few calls to LAPACK routines and then
        # checks if compilation, linking and execution are succesful

        sys.stdout.flush()
        writefile('tmpf.f',"""
      program ftest
      integer info, ipiv(1)
      double precision A(1)
      A(1)=1
      call dgetrf(1,1,A,1,ipiv,info)
      stop
      end\n""")

        fcomm = self.mpif77+' -o tmpf '+'tmpf.f '+self.lapacklib+' '+self.blaslib+' '+self.ldflags_f77
        (output, error, retz) = runShellCommand(fcomm)

        if(retz != 0):
            print('\n\nLAPACK: provided LAPACK cannot be used! aborting...')
            print('error is:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()
        
        comm = './tmpf'
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nLAPACK: provided LAPACK cannot be used! aborting...')
            print('error is:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()
        
        killfiles(['tmpf.f','tmpf'])
        print('yes')

        return 0;




    def down_install_lapack(self):
        print(""" Download and install LAPACK from netlib.org""")

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

        if(not os.path.isfile(os.path.join(os.getcwd(),getURLName(self.lapackurl)))):
            print("""Download LAPACK...""", end=' ')
            sys.stdout.flush()
            downloader(self.lapackurl, self.downcmd)
            print("done")    
            print("""Unzip and untar LAPACK...""", end=' ')
        comm = 'gunzip -f lapack.tgz'
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nLAPACK: cannot unzip lapack.tgz')
            print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()

        comm = 'tar xf lapack.tar'
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nLAPACK: cannot untar lapack.tgz')
            print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()
        os.remove('lapack.tar')
        print("done")

        # change to LAPACK dir
        os.chdir(os.path.join(os.getcwd(),'lapack-3.2'))

        # compile and generate library
        self.set_etime()
        self.write_makeinc()

        print('Compiling LAPACK (this will take several minutes)...', end=' ')
        sys.stdout.flush()
        comm = self.make+' lib'
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nLAPACK: error building LAPACK')
            print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()

        # write the log on a file
        fulllog = os.path.join(savecwd,'log/laplog')
        writefile(fulllog, output+error)
        print('Installation of LAPACK successful.')
        print('(log is in ',fulllog,')')

        # move librefblas.a to the lib directory        
        os.rename('libreflapack.a',os.path.join(savecwd,'lib/libreflapack.a'))

        # set framework variables to point to the freshly installed BLAS library
        self.lapacklib  = os.path.join(savecwd,'lib/libreflapack.a ')
        Frame.lapacklib = os.path.join(savecwd,'lib/libreflapack.a ')
        os.chdir(savecwd)

        return







    def set_etime(self):
        """ Detects the FORTRAN timer """

        print('Setting ETIME...', end=' ')
        writefile('tmpf.f',"""
      PROGRAM DSECND
      REAL               T1
      REAL               TARRAY( 2 )
      REAL               ETIME_
      EXTERNAL           ETIME_
      T1 = ETIME_( TARRAY )
      END\n""")
        fcomm = self.mpif77+' '+self.fcflags+' -o tmpf tmpf.f'
        (output, error, retz) = runShellCommand(fcomm)

        if(retz == 0):
            self.etime = 'EXT_ETIME_'
            print(self.etime)
            killfiles(['tmpf.f', 'tmpf'])
            return


        writefile('tmpf.f',"""
      PROGRAM DSECND
      REAL               T1
      REAL               TARRAY( 2 )
      REAL               ETIME
      EXTERNAL           ETIME
      T1 = ETIME( TARRAY )
      END\n""")
        
        fcomm = self.mpif77+' '+self.fcflags+' -o tmpf tmpf.f'
        (output, error, retz) = runShellCommand(fcomm)

        if(retz == 0):
            self.etime = 'EXT_ETIME'
            print(self.etime)
            killfiles(['tmpf.f', 'tmpf'])
            return


        writefile('tmpf.f',"""
      PROGRAM DSECND
      REAL               T
      INTRINSIC          ETIME
      CALL CPU_TIME( T )
      END\n""")
        
        fcomm = self.mpif77+' '+self.fcflags+' -o tmpf tmpf.f'
        (output, error, retz) = runShellCommand(fcomm)

        if(retz == 0):
            self.etime = 'INT_CPU_TIME'
            print(self.etime)
            killfiles(['tmpf.f', 'tmpf'])
            return
        

        writefile('tmpf.f',"""
      PROGRAM DSECND
      REAL               T1
      REAL               TARRAY( 2 )
      REAL               ETIME
      INTRINSIC          ETIME
      T1 = ETIME( TARRAY )
      END\n""")
        
        fcomm = self.mpif77+' '+self.fcflags+' -o tmpf tmpf.f'
        (output, error, retz) = runShellCommand(fcomm)

        if(retz == 0):
            self.etime = 'INT_ETIME'
            print(self.etime)
            killfiles(['tmpf.f', 'tmpf'])
            return


        killfiles(['tmpf.f', 'tmpf'])
        self.etime = 'NONE'



    def write_makeinc(self):
        """ Writes the make.inc file for LAPACK installation """

        print('Writing make.inc...', end=' ')
        writefile('make.inc',"""
SHELL = /bin/sh
PLAT =
FORTRAN  = """+self.mpif77+"""
OPTS     = """+self.fcflags+"""
DRVOPTS  = $(OPTS)
NOOPT    = """+self.noopt+"""
LOADER   = """+self.mpif77+"""
LOADOPTS =
TIMER    = """+self.etime+"""
ARCH     = ar
ARCHFLAGS= cr
RANLIB   = """+self.ranlib+"""
BLASLIB      = """+self.blaslib+"""
LAPACKLIB    = libreflapack.a
TMGLIB       = tmglib.a
EIGSRCLIB    = eigsrc.a
LINSRCLIB    = linsrc.a
        """)
        print('done.')
