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


class Scalapack(Frame):
    """ This class takes care of the BLACS installation. """
    def __init__(self):
        print('\n','='*40)
        print('ScaLAPACK installer is starting now. Buckle up!')
        print('='*40)

        self.down_install()






    def write_slmakeinc(self):
        """ Writes the SLmake.inc file for ScaLAPACK installation """

        sdir = os.getcwd()
        print('Writing SLmake.inc...', end=' ')
        sys.stdout.flush()
        writefile('SLmake.inc',"""
SHELL         = /bin/sh

home          = """+sdir+"""

PLAT          = """+self.plat+"""

USEMPI        = -DUsingMpiBlacs

SMPLIB        = 
BLACSFINIT    = """+self.blacsF77lib+"""
BLACSCINIT    = """+self.blacsClib+"""
BLACSLIB      = """+self.blacslib+"""
TESTINGdir    = $(home)/TESTING


CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

PBLASdir      = $(home)/PBLAS
SRCdir        = $(home)/SRC
TESTdir       = $(home)/TESTING
PBLASTSTdir   = $(TESTINGdir)
TOOLSdir      = $(home)/TOOLS
REDISTdir     = $(home)/REDIST
REDISTTSTdir  = $(TESTINGdir)

F77           = """+self.mpif77+"""
CC            = """+self.mpicc+"""
NOOPT         = """+self.noopt+"""
F77FLAGS      =  $(NOOPT) """+self.fcflags+"""
CCFLAGS       = """+self.ccflags+"""
SRCFLAG       =
F77LOADER     = $(F77) 
CCLOADER      = $(F77)
F77LOADFLAGS  = """+self.ldflags_f77+"""
CCLOADFLAGS   = """+self.ldflags_c+"""

CDEFS         = -DNO_IEEE $(USEMPI) """+self.mangling+"""

ARCH          = ar
ARCHFLAGS     = cr
RANLIB        = """+self.ranlib+"""

SCALAPACKLIB  = $(home)/libscalapack.a
BLASLIB       = """+self.blaslib+"""
LAPACKLIB     = """+self.lapacklib+"""

PBLIBS        = $(SCALAPACKLIB) $(FBLACSLIB) $(LAPACKLIB) $(BLASLIB) $(SMPLIB)
PRLIBS        = $(SCALAPACKLIB) $(CBLACSLIB) $(SMPLIB)
RLIBS         = $(SCALAPACKLIB) $(FBLACSLIB) $(CBLACSLIB) $(LAPACKLIB) $(BLASLIB) $(SMPLIB)
LIBS          = $(PBLIBS)
        """)
        print('done.')






    def down_install(self):
        """ Downloads ind installs ScaLAPACK """

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

        if(not os.path.isfile(os.path.join(os.getcwd(),getURLName(self.scalapackurl)))):
            downloader(self.scalapackurl, self.downcmd)

        comm = 'gunzip -f scalapack-1.8.0.tgz'
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nScaLAPACK: cannot unzip scalapack-1.8.0.tgz')
            print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()
        

        comm = 'tar xf scalapack-1.8.0.tar'
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nScaLAPACK: cannot untar scalapack-1.8.0.tar')
            print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
            sys.exit()
        os.remove('scalapack-1.8.0.tar')
        
        os.chdir(os.path.join(os.getcwd(),'scalapack-1.8.0'))

        self.write_slmakeinc()

        print('Compiling ScaLAPACK...', end=' ')
        sys.stdout.flush()
        comm = self.make
        (output, error, retz) = runShellCommand(comm)
        if(retz != 0):
            print('\n\nScaLAPACK: error building ScaLAPACK')
            print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
            writefile(os.path.join(savecwd,'log/scalog'), output+error)
            sys.exit()

        fulllog = os.path.join(savecwd,'log/scalog')
        writefile(fulllog, output+error)
        print('Installation of ScaLAPACK successful..')
        print('(log is in ',fulllog,')')
		
        if(self.testing == 1):
            print('Compiling test routines...')
            sys.stdout.flush()
            comm = self.make+' exe'
            (output, error, retz) = runShellCommand(comm)
            if(retz != 0):
                print('\n\nScaLAPACK: error building ScaLAPACK test routines')
                print('stderr:\n','*'*40,'\n',error,'\n','*'*40)
                writefile(os.path.join(savecwd,'log/scalog'), output+error)
                sys.exit()
            print('done')

        
        os.rename('libscalapack.a',os.path.join(savecwd,'lib/libscalapack.a'))

        self.scalapacklib  = os.path.join(savecwd,'lib/libscalapack.a ')
        Frame.scalapacklib = os.path.join(savecwd,'lib/libscalapack.a ')

        os.chdir(savecwd)
        print("done. ScaLAPACK is installed. Use it in moderation :-)")
