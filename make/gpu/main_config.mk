#=========================================================================
# file: $(BLD_DIR)/config.mk
#
# This makefile fragment is the main configuration file for the pscfpp
# build system.  A copy of this file is included by all other makefiles. 
# This file is created and installed by the "setup" script. One copy of 
# this file is installed in the src/ directory, which is used for in-source 
# compilation. Two others copies are created in the root directories of 
# the bld/parallel and bld/serial build directory trees, which are used for
# out-of-source compilation of parallel and serial programs respectively.
#
#=========================================================================
# This file contains user-modifiable definitions of several types of 
# makefile variables:
# 
#  - Variables ROOT_DIR, SRC_DIR, BLD_DIR and BIN_DIR that contains 
#    absolute paths for the simpatico root directory and some of its
#    subdirectories.
#
#  - Variables UTIL_MPI, UTIL_DEBUG, and UTIL_CXX11 that, if defined
#    enabled conditional compilation of optional features of the code.
#
#  - Environment-dependent variabes that control the command name by
#    which the compiler is invoked and some compiler command line options.
# 
#  - A variable MAKEDEP that enables automatic dependency generation.
#
#=========================================================================
# Variables that define absolute directory paths
#
# In the config.mk file installed in each build directory, correct
# values of these variables should have been set by the setup script,
# and should not need to modified by the user.

# Absolute path to the root simpatico directory
ROOT_DIR=/home/dorfmank/chawl029/GPU_pscfpp/merged/pscfpp

# Path to the build directory (location for intermediate generated files)
# This should also be the directory that contains this script.
BLD_DIR=$(ROOT_DIR)/src

# Path to the source directory (contains C++ source code files)
SRC_DIR=$(ROOT_DIR)/src

# Installation directory for binary executable program files
BIN_DIR=$(ROOT_DIR)/bin

# Directory for shared permanent (read-only) data used by programs.
DAT_DIR=$(ROOT_DIR)/data

#======================================================================
# Variables that (if defined) enable compilation of optional features

# Defining UTIL_MPI enables linking to an MPI library. MPI must be
# enabled by uncommenting this line to build a parallel program.
#UTIL_MPI=1

# Defining UTIL_DEBUG enables a variety of extra sanity checks, at some
# cost in speed. Debugging is disabled (commented out) by default.
#UTIL_DEBUG=1

# Defining UTIL_GPU enables compilation using the nvidia compiler nvcc
#
UTIL_GPU=1

# Defining UTIL_CXX11 enables the use of features of the C++11 standard.
# Use of C++11 is enabled by default. If C++11 is disabled, the code uses 
# only syntax and features defined in the older 1998 C++ standard.
UTIL_CXX11=1

# Comment After setup, the above definitions of UTIL_MPI and UTIL_DEBUG 
# may be uncommented or commented out invoking the "configure" script 
# from the build directory that contains this file, thus enabling or 
# disabling the MPI or and/or debugging features.  Specifically, the 
# configure script may be be invoked with option -m to enable 
# (option -m1) or disable (-m0) MPI, and with option -g to enable 
# (-g1) or disable (-g0) debugging. For example, the command:
#
# ./configure -m1 -g0 
#
# would enable MPI (uncomment the definition of UTIL_MPI) and disable 
# debugging (comment out the definition of UTIL_DEGUB) prior to 
# compilation. 
#
#======================================================================
# Environment-dependent makefile variables.
#
# Variables defined in this block define the names of the commands 
# used to invoke the compiler when compiling and linking files and 
# some of the command line options passed to the compiler to control,
# e.g., optimization level, warnings, and search paths for header files 
# and libraries. Different choices may be used when compiling serial 
# programs (when UTIL_MPI is not defined) or parallel programs (when
# UTIL_MPI is defined). The choice of values of variables defined here 
# may thus be different for different combinations of compiler, MPI
# library and operating system.  See the section of this file 
# entitled "Makefile Patterns and Recipes" for a discussion of how
# these variables are used.
#
# The following block of variable definitions is copied by the setup
# script from a compiler setup file in the make/compiler directory
# and incorporated into this config.mk file. The name of the desired 
# compiler setup file may be specified as an argument to the setup 
# script command (e.g., "> ./setup intel").
#

# ---------------------------------------------------------------
#### Definitions for generic unix environment ###################
#    (gcc compiler, mpicxx mpi wrapper)
#
# The definitions given below work for systems in which:
#
#   - All header files and libraries are in standard locations.
#   - The command g++ is used to invoke the C++ compiler
#   - A compiler wrapper script named "mpicxx" is used for MPI code
#   - Compiler options are compatible with those for the gcc compiler
#
# These definitions work on Mac OSX for serial code, for which 
# g++ invokes the clang compiler, and for MPI code with some MPI
# libraries.
#
# ---------------------------------------------------------------
# General definitions

# Path to search for header files (must include SRC_DIR)
INCLUDES= -I$(SRC_DIR)

# Compiler option to specify ANSI C++ standard
ifdef UTIL_CXX11
   CXX_STD = --std=c++11
else
   CXX_STD = --std=c++98
endif

# ---------------------------------------------------------------
# Compiler and options used for serial programs 

# Command to invoke C++ compiler for serial (non-MPI) code
CXX_SER=g++

# Flags passed to compiler when debugging is enabled
CXXFLAGS_SER_DEBUG= -Wall $(CXX_STD)

# Flags passed to compiler when debugging is disabled (production code)
CXXFLAGS_SER_FAST= -Wall $(CXX_STD) -O3 -ffast-math -Winline

# Compiler flags used in unit tests
TESTFLAGS= -Wall $(CXX_STD)

# ---------------------------------------------------------------
# Compiler, options and execution command for parallel programs 

# Command to invoke the C++ compiler for compiling MPI code.
# Note: This is often name of a wrapper script provided by the 
# MPI library implementation
CXX_PAR=mpicxx

# Flags passed to compiler when debugging is enabled
CXXFLAGS_PAR_DEBUG= -Wall $(CXX_STD)

# Flags passed to compiler when debugging is disabled
CXXFLAGS_PAR_FAST= -Wall $(CXX_STD) -O3 -ffast-math -Winline

# MPI execution command (followed by integer number of processors)
MPIRUN=mpirun -np

# ---------------------------------------------------------------
# Compiler for gpu
CXX_GPU=nvcc

CXXFLAGS_GPU= -O3 -arch=sm_35
# ---------------------------------------------------------------
# Linker / Loader 

# Flags passed to compiler for linking and loading
LDFLAGS=

# ---------------------------------------------------------------
# Archiver

# Library archiver command (for creating static libraries)
AR=ar

# Flags (command line options) passed to archiver
ARFLAGS=rcs

#-----------------------------------------------------------------------
# Choose values of CXX and CXX_FLAGS

# Choose values for CXX (compiler command) and CXX_FLAGS (general 
# compiler options) depending on whether UTIL_MPI and/or UTIL_DEBUG
# are defined (i.e., on whether MPI and/or debugging features are 
# enabled).

ifneq ($(UTIL_MPI),1)
   ifeq ($(UTIL_GPU),1)
      CXX=$(CXX_GPU)
      CXXFLAGS=$(CXXFLAGS_GPU)
   else
      # Serial programs:
      # Serial compiler command (MPI disabled)
      CXX=$(CXX_SER)
      ifdef UTIL_DEBUG
         # Flags for serial programs with debugging
         CXXFLAGS=$(CXXFLAGS_SER_DEBUG)
      else
      # Flags for serial programs with no debugging
         CXXFLAGS=$(CXXFLAGS_SER_FAST)
      endif
   endif
else
   # Parallel programs:
   # Parallel compiler command or wrapper script (MPI enabled)
   CXX=$(CXX_PAR)
   ifdef UTIL_DEBUG
      # Flags for parallel programs with debugging
      CXXFLAGS=$(CXXFLAGS_PAR_DEBUG)
   else
      # Flags for parallel programs with no debugging
      CXXFLAGS=$(CXXFLAGS_PAR_FAST)
   endif
endif

# ======================================================================
# Makefile Patterns and Recipes
#
# The makefile variables defined above are used in the makefile pattern 
# rules and recipes that control compilation of C++ files, creation of 
# libraries, and linking to create executables. The following sections 
# briefly explain these rules, to provide a context for the meaning of 
# the variables defined above.
#
#-----------------------------------------------------------------------
# Compiler Pattern Rules:
#
# The pattern rule for compiling and linking C++ files in a particular 
# namespace is defined in the appropriate namespace level config.mk 
# file. For example, the rule for compiling C++ files in the simp/ 
# directory tree, which contains all classes defined in the Simp C++
# namespace, is given in the file simp/config.mk. The pattern rules for 
# different namespaces are similar except for differences in which 
# preprocessor variable definitions are passed to the compiler. For
# each namespace, the basic compiler pattern rule is of the form:
# 
# $(BLD_DIR)%.o:$(SRC_DIR)/%.cpp
#      $(CXX) $(INCLUDES) $(DEFINES) $(CXXFLAGS) -c -o $@ $<
#
# This pattern compiles a *.cpp file in a subdirectory of the source 
# directory $(SRC_DIR) and creates a *.o object file in a corresponding
# subdirectory of the build directory, $(BLD_DIR). The variables used
# in this pattern are:
#
# CXX         - C++ compiler executable name 
# INCLUDES    - Directories to search for included header files
# DEFINES     - compiler options that define C preprocessor macros
# CXXFLAGS    - compiler options used during compilation
#
# Comments:
# 
# 1) The variable CXX is the name of an executable command that may
# be either the name of the compiler command (e.g., g++) or the name
# of a wrapper script that invokes the compiler (e.g., the mpicxx
# script provided with the OpenMPI MPI library to invoke the compiler
# with appropriate search paths)
#
# 2) The variable INCLUDES is a string that must include the path 
# $(SRC_DIR) to the simpatico/src directory, in order to allow the 
# compiler to find header files that are part of the package.
#
# 3) The variable DEFINES in the above pattern is a stand-in for a 
# variable that specifies a list of C preprocessor macro definitions. 
# This variable is not defined in this main configuration file, and
# is assigned different values for code in different namespace level 
# directories, which are defined in the namespace level patterns.mk 
# files. The value of $(DEFINES) for each namespace contains a string 
# of compiler options that use the compiler "-D" option to define the 
# set of preprocessor macro definitions used to control conditional
# compilation of optional features that are relevant in a particular 
# namespace. Each of these preprocessor variable macros has the
# same name as a corresponding makefile variable that must be defined
# to enable the feature. Thus for, example, when the build system has 
# been configured to enable debugging, the DEFINES string will include
# a substring "-D UTIL_DEBUG" to define the UTIL_DEBUG preprocessor 
# macro and thereby enable conditional compilation of blocks of code
# that contain optional sanity checks.  
#
# 4) The variable $(CXXFLAGS) should specify all flags that are used by 
# the compiler, rather than only the preprocessor, and that are used in
# all namespaces. This string normally contains the $(CXX_STD) string 
# as a substring, as well as options that specify the optimization 
# level (e.g., -O3) and any desired compiler warnings (e.g., "-Wall").
#
#-----------------------------------------------------------------------
# Archiver Recipes:
#
# The simpatico build system creates a static library in each namespace
# level subdirectory of the build directory in which code is compiled.
# The recipe used to compile this library is defined in the sources.mk
# file in the appropriate namespace-level directory. The rule for the
# McMd namespace, as an example, is of the form
#
# $(AR) rcs $(mcMd_LIB) $(mcMd_OBJS)
#
# where $(AR) is the name of archiver command used to create a library,
# $(mcMD_LIB) is an absolute path for the resulting library file and 
# $(mcMd_OBJS) is a string that contains absolute paths for all of the
# *.o object files created by compiling source files in the directory
# src/mcMd. Recipes for other namespaces are analogous.
#
#-----------------------------------------------------------------------
# Linker recipes:
# 
# Executable files are created by linking the compiled main program to
# the required set of static libraries. For example, recipe for creating 
# the mdSim executable is of the form
#
#	$(CXX) -o $(mdSim_EXE) $(mdSim).o $(LIBS) $(LDFLAGS)
#
# Here $(mdSim_EXE) is the path to the executable, which is installed 
# in the bin/ directory by default, $(mdSim).o is the path to the 
# object file created by compiling the src/mcMd/mdSim.cpp source 
# file, $(LIBS) is a list of all required state libraries files, and
# $(LDFLAGS) is a list of flags passed to the linker. 
#
# The variable $(LDFLAGS) is empty by default, but can, if necessary, 
# be used to specify a non-standard path to a directory containing 
# the MPI library when compiling parallel programs. This should not 
# be necessary if the compiler is invoked using the name of a wrapper
# script that sets this path automatically.
 
#=======================================================================
# Automatic dependency generation.
 
# Executable file invoked to compute dependencies among header files.
MAKEDEP=$(BIN_DIR)/../scripts/python/makeDep

# The file $(BIN_DIR)/makeDep is an executable python script that is
# installed in the binary directory specified by the setup script, 
# and that is used during compilation to analyze dependencies among
# C++ files. The makeDep script imports a python module named makeDepend 
# that is located in the $(ROOT_DIR)scripts/python directory. If the 
# python interpreter fails to find the makeDepend module, you may 
# need to add the absolute path to the simpatico/scripts/python directory 
# to the PYTHON_PATH environment variable, as described in the web 
# documentation.
#
