#-----------------------------------------------------------------------
# file: src/compiler.mk
#
# This makefile fragment is included by all other makefiles. It defines
# absolute paths for the root and src/ directories, specifies the choice
# of compiler and various compiler options and (optionally) defines the 
# macro MAKEDEP that enables automatic dependency generation.
# 
# Users may need to modify the following variables:
#
# ROOT_DIR = absolute path to the root (e.g., trunk) directory.
# COMPILER = label for the compiler (gcc, intel, etc.)
#
# Compiler options may also be modified by modifying values of CXXFLAGS
# etc. within the ifeq() ... endif block for the relevant compiler.
# Users may also wish to uncomment the definition of MAKEDEP, as 
# discussed below.
#
# Users may modify the file src/compiler.mk, but should avoid modifying 
# the repository copy src/compiler.mk_r, which is under version control. 
# The operational file src/compiler.mk is created by the configure script,
# by making a copy of src/compiler.mk_r and modifying the value of the
# variable ROOT_DIR (which has a placeholder value in src/compiler.mk_r).
# After the configure script is run from the root simpatico/ directory,
# the value of ROOT_DIR should have been set to the absolute path to 
# that directory on the users machine.
#-----------------------------------------------------------------------
# Absolute directory paths

ROOT_DIR=/home/david/code/pscf++
SRC_DIR=$(ROOT_DIR)/src
BIN_DIR=$(ROOT_DIR)/bin
LIB_DIR=$(ROOT_DIR)/lib
TESTS_DIR=$(ROOT_DIR)/tests

#-----------------------------------------------------------------------
# Automatic dependency generation:
#
# To enable automatic dependency generation: (1) uncomment the line below
# that defines MAKEDEP, and (2) Modify the PYTHONPATH unix environment
# variable to include the $(ROOT_DIR)/tools/python directory. The script 
# makeDep is located in the $(ROOT_DIR)/bin/ directory, and is invoked 
# as a side effect of compilation in order to generate a dependency file.  
# The PYTHONPATH must include the $(ROOT_DIR)/tools/python directory so
# that the python interpreter can find the python modules that are used 
# by the makeDep script.
 
MAKEDEP=$(ROOT_DIR)/bin/makeDep

#-----------------------------------------------------------------------
# COMPILER identifier:
#
# The variable COMPILER is a string that identifies a choice of compiler.
# To select a compiler for which default settings are available, set
# the variable COMPILER by uncommenting one of the following lines
# and commenting out all others. The default choice is gcc, the gnu
# compiler collection.  
#
# String values that contain the string "mpi" set compiler options so 
# as to link to an mpi library. Some these blocks use an mpicxx script
# to invoke the compiler, but users may modify this if needed.

# Choose (uncomment) one value of COMPILER, comment out all others

COMPILER:=gcc
#COMPILER:=gcc_mpicxx
#COMPILER:=gcc_openmpi
#COMPILER:=intel
#COMPILER:=intel_mpicxx
#COMPILER:=pathscale
#COMPILER:=pathscale_mpi

# ------------------------------------------------------------------
# Compiler settings:
#
# Each of the following blocks sets values for the following set
# of makefile variables, using values appropriate to a particular
# compiler:
#
# CXX        - path to C++ compiler
# CXXFLAGS   - flags used to compile source files, without linking
# LDFLAGS    - flags used to compile and link a main program
# TESTFLAGS  - flags usd to compile unit test programs
# AR         - path to archiver, to create library (*.a) files
# ARFLAGS    - flags used by archiver
#
# Blocks for compilers that link to the MPI library also define 
# UTIL_MPI=1. The default value of CXX for gcc_mpicxx and intel_mpicxx
# is CXX=mpicxx. The string mpicxx is the name of a script that 
# invokes an mpi compiler several systems that we have used for
# development.  On these systems, the choices of actual compiler 
# and mpi library are determined by loading appropriate software
# modules (e.g., module load intel-mpi).
# ------------------------------------------------------------------

#-- Gnu GCC compiler (Serial) ------
ifeq ($(COMPILER),gcc)
  CXX=g++
  CXXFLAGS= -O3 -ffast-math -Wall -Winline -std=c++98 -pedantic
  LDFLAGS=
  TESTFLAGS= -Wall -std=c++98 -pedantic
  AR=ar
  ARFLAGS=rcs
endif

# -- Gnu GCC compiler, MPI via mpicxx script -----
ifeq ($(COMPILER),gcc_mpicxx)
  CXX=mpicxx
  CXXFLAGS= -O3 -ffast-math -Wall -Winline -std=c++98
  LDFLAGS=
  TESTFLAGS= -O3 -ffast-math -Wall -std=c++98 -pedantic
  AR=ar
  ARFLAGS=rcs
  UTIL_MPI=1
endif

# -- Gnu GCC compiler, OpenMPI - explicit headers and libs) -----
ifeq ($(COMPILER),gcc_openmpi)
  CXX=g++
  CXXFLAGS= -O3 -ffast-math -Wall -Winline -std=c++98 -I/opt/local/include/openmpi
  LDFLAGS= -L/opt/local/lib -lmpi_cxx -lmpi
  TESTFLAGS= -O3 -ffast-math -Wall -std=c++98 -pedantic
  AR=ar
  ARFLAGS=rcs
  UTIL_MPI=1
endif

# -- Intel ICC compiler (serial) ---
ifeq ($(COMPILER),intel)
  CXX=icpc
  CXXFLAGS= -fast -ansi
  LDFLAGS= -fast
  TESTFLAGS= -ansi
  AR=xiar
  ARFLAGS=rcs
endif

# -- Intel ICC compiler (MPI) ------
ifeq ($(COMPILER),intel_mpicxx)
  CXX=mpicxx
  CXXFLAGS= -ansi
  LDFLAGS=
  TESTFLAGS= -ansi
  AR=xiar
  ARFLAGS=rcs
  UTIL_MPI=1
endif

# -- PathScale compiler ------------
ifeq ($(COMPILER),pathscale)
  CXX=pathCXX
  CXXFLAGS= -Ofast
  LDFLAGS= -Ofast
  TESTFLAGS=
  AR=ar
  ARFLAGS=rcs
endif

# -- PathScale compiler (MPI) ------
ifeq ($(COMPILER),pathscale_mpi)
  CXX=mpicxx
  CXXFLAGS= -Ofast
  LDFLAGS= -Ofast
  TESTFLAGS=
  AR=ar
  ARFLAGS=rcs
  UTIL_MPI=1
endif

