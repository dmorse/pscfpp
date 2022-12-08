
# PSCF - Polymer Self-Consistent Field Theory (C++/CUDA)

PSCF is a package of software for solving the Edwards-Helfand 
self-consistent field theory (SCFT) of polymer liquids. The version 
described here is written primarily in C++, with GPU accelerated code 
in CUDA.  

## History

This C++/CUDA version of PSCF is intended to supersede the older
Fortran PSCF program. The older Fortran program is maintained in a 
separate github.com repository at https://github.com/dmorse/pscf.
This new version currently provides almost all of the capabilities of 
the older Fortran program, and several important new capabilities, as
discussed below.

## Overview

Differences between the C++/CUDA version of PSCF and the older Fortran 
version and expected advantages of the new code include:

   - PSCF (C++/CUDA) is an extensible package of several different 
     programs designed for use with different geometries, boundary 
     conditions, algorithms or hardware, implemented within a common 
     framework. 

   - PSCF (C++/CUDA) enables simulations of mixtures containing arbitrary 
     acyclic branched copolymers, in addition to the linear block 
     copolymers and linear homopolymers allowed by the Fortran PSCF code.

   - Adoption of C/C++ as a base language has simplified implementation 
     of programs that that uses graphical processing units (GPUs).

## Programs

Currently, the package contains the following SCFT programs:

   - **pscf_fd** : A one-dimensional finite difference program for problems 
     that involve variation in a single spatial coordinate in Cartesian, 
     cylindrical or spherical coordinates.

   - **pscf_pc** : CPU-based pseudo-spectral solvers for periodic 
     microstructures that are periodic in 1, 2, or 3 coordinates.

   - **pscf_pg**: Analogous GPU-enabled pseudo-spectral solvers for period 
     microstructures in 1, 2 or 3 dimensions.

*psf_fd* : The one-dimensional finite difference program is useful for 
treating problems involving flat or curved interfaces, as well as cylindrical 
or spherical copolymer micelles. The executable for this program is named 
pscf_fd. The suffix "fd" stands for "finite difference".

*pscf_pc* : The pseudo-spectral CPU programs for periodic microstructures 
are closely analogous to the older PSCF Fortran program, and provide similar 
capabilties and performance.  Different executable files are used in the 
new code to solve 1, 2 and 3 dimensionally periodic structures, which are 
named pscf_pc1, pscf_pc2 and pscf_pc3, respectively. In these names, "pc" 
stands for "periodic CPU".  These programs allow the user to search for a 
solution with any specified crystal system type and space group symmetry, 
and provide efficient algorithms to relax the unit cell parameters so as 
to minimize the free energy. 
 
*pscf_pg* : The GPU-accelerated pseudo-spectral programs for periodic 
structures are based on algorithms similar to those used in the CPU 
pseudo-spectral programs, and provide almost identical capabilities.  
The GPU accelerated programs for solving 1, 2 and 3 dimensionally periodic 
structures are named pscf_pg1, pscf_pg2 and pscf_pg3, respectively, where 
"pg" stands for "periodic GPU". 

## Getting the source code

The PSCF C++/CUDA source code is maintained in the github repository

   <https://github.com/dmorse/pscfpp>.

The source code may be obtained by using a git version control system 
client to clone the repository. To do so on a machine with a git client, 
enter the command:
``` 
git clone --recursive https://github.com/dmorse/pscfpp.git
```
Note the use of the --recursive option to the git clone command. This
is necessary to clone several git submodules that are maintained in 
separate repositories. This command will create a new directory called 
pscfpp/ that contains all of the source code and associated 
documentation, including all required git submodules.

## Documentation

PSCF is distributed with source files for an html web manual. A copy
of the documentation of a recent version is available online at

  <https://dmorse.github.io/pscfpp-man>

After cloning the source code, you can use the doxygen documentation
generation program to generate a local copy of this documentation. 
To do this, doxygen must be installed on your computer, and the 
directory containing the doxygen executable must be in your command 
search PATH. To generate documentation:

   - Change directory (cd) to the pscfpp/ root directory

   - Enter "make html"

This should create many html files in the pscfpp/docs/html directory.  
To begin reading the documentation, point a browser at the file
pscfpp/docs/html/index.html, which is the main page of the manual.

## Dependencies

The PSCF source code is written using C++ as the primary language, with
CUDA used in some programs for GPU acceleration. PSCF is only provided 
in source file format - all programs must be compiled from source.
This package was developed on linux and and Mac OS X operating systems 
using standard unix utilities, and is designed to run on these or other 
unix-like systems. 

To compile and run this or other unix software on Mac OS X, you must 
first install the unix command line tools. These tools are provided as 
part of the much larger XCode Mac development environment, but can also 
be installed separately.

The CPU-based programs within the PSCF package depend on the following 
external libraries:

  - Gnu scientific library (GSL)

  - FFTW fast Fourier transform library

The one-dimensional finite difference program pscf_fd requires only GSL, 
and not FFTW. The pscf_pc CPU-based programs for spatially periodic 
structures require both GSL and FFTW.

The GPU-accelerated pscf_pg programs can only run on a computer with an
appropriate NVIDIA GPU. To compile or run these programs, the system must 
also have an NVIDIA CUDA development environment that provides the CUFFT 
fast Fourier transform library. 

Procedures for installing these dependencies are different for 
different operating system environments and different package managers. 
Instructions for some common environments are given in the web manual.

## Compiling

Complete directions for compiling and installing PSCF are provided in 
section 2 of the html documentation. Short instructions for compiling, 
after cloning the git repository and installing all of the required 
dependencies, are given below:

   - Add the pscfpp/bin directory to your linux command search PATH
     environment variable.
   
   - Add the pscfpp/lib/python directory to your PYTHONPATH 
     environment variable.
   
   - cd to the pscfpp/ root directory
   
   - Run the pscfpp/setup script by entering "./setup" from the pscfpp
     directory, with or without an optional filename argument.
     This script installs default versions of several files that are 
     required by the build system. The setup script usually only needs
     to be run once, before compiling for the first time. In a linux
     environment, it is usually sufficient to run the setup script
     without a filename argument.  See the html documentation for 
     information about how to invoke this script on a Mac. 

   - To compile and install only CPU-based programs in the package 
     (excluding GPU-accelerated programs), enter "make all-cpu"
   
   - To compile all programs, including the GPU-accelerated programs,
     on a machine that has a CUDA development environment install, 
     instead enter "make all". 

The "make all-cpu" and "make all" commands should install executable
programs in the pscfpp/bin directory.

## Command line usage 

PSCF is a package containing several different SCFT programs designed 
for different geometries, different algorithms or different hardware. 
Executable names (given above) are:

   - pscf_fd : 1D finite-difference program 

   - pscf_pcD : CPU based programs for D = 1,2, or 3 dimensional periodic
     structures. 

   - pscf_pgD : GPU based programs for D = 1,2, or 3 dimensional periodic
     structures

In the names pscf_pcD and pscf_pgD, D denotes a dimension dimension
of space that can be D = 1, 2, 3. For example, the CPU program for 
three-dimensionally periodic microstructures is thus pscf_pc3.

Each of these programs reads a parameter file and a command file. The 
parameter file is fixed-format file that contains parameters required 
to initialize the program. The command file is a more flexible script 
containing a sequence of commands that are read and executed sequentially 
to specify a sequence of computational steps. The command file for a 
standard SCFT calculation specifies the name of a file that contains 
an initial guess for monomer chemical potential fields and names of 
files to which final chemical potential and monomer concentration 
fields should be written.

The usual command line syntax for invoking any pscf program is:
```
program -e -p param -c command
```
where "program" denotes the name of the program, "param" denotes the 
path to a parameter file, and "command" denotes the path to a parameter 
file.  For example, one might enter
```
pscf_pc3 -e -p param -c command
```
to run the pscf_pc3 CPU program for three dimensional periodic structures. 
The "-e" command line option causes the program to echo the parameter file 
to standard output as this file is read, which is generally a good idea.

The above form of the command would write log output to the screen.  
Output produced during a computation may also be redirected to a log file 
by using the unix ">" standard output redirection operator. For example,
a command of the form
```
program -e -p param -c command > log
```
would redirect output that is written direct the computation to a file 
named "log".

## Examples

The directory pscfpp/examples contains a set of examples of simple 
calculations. Each example directory contains a parameter file (usually
named "param"), a command file (usually named "command"), and a input 
chemical potential field (w field) file.  Top level subdirectories of 
pscfpp/examples contain examples for different PSCF programs or 
families of closely related programs.

Subdirectory examples/fd subdirectory contains examples for the 1D 
finite-difference program pscf_fd. Top level subdirectories of 
examples/fd contain examples for planar, cylindrical and spherical 
geometries, as indicated by the subdirectory names. One or more 
example is given for each geometry.

Subdirectory examples/pspc contains examples contains examples the CPU 
based pscf_pc programs. Top level subdirectories contain solutions for 
a particular type of physical system, e.g., examples/pspc/diblock 
contains examples for a diblock copolymer melt.  Subdirectories of 
examples/pspc/diblock contain examples for lamellar, hexagonal, and 
BCC structures, among others.

Subdirectory examples/pspg contains examples for the pscf_pg programs
for periodic structures.

Each example directory contains a script named "run" that can used to
run the example using the supplied input files and appropriate options. 
The simplest way to run an example is to change directory (cd) to the
directory containing the example and enter "./run" from that directory.
Users may also inspect the run file to see the command and options
require to run the program. Most example directories also contain a
script named "clean" that can be run to remove all output files that
are created by running the example.

## Contributors

- David Morse
- Guo Kang Cheong
- Anshul Chawla
- Ryan Collanton
- Ben Magruder

## License

The C++/CUDA version of PSCF is free, open source software. It is 
distributed under the terms of the GNU General Public License as (GPL) 
published by the Free Software Foundation, either version 3 of the 
License or (at your option) any later version.  PSCF is distributed 
without any warranty, without even the implied warranty of merchantability 
or fitness for a particular purpose.  See the LICENSE file or the 
<a href=http://www.gnu.org/licenses/>gnu web page</a> for details.

## Support 

Development of PSCF is supported by the National Science Foundation 
program for Cyberinfrastructure for Sustained Scientific Development 
(CSSI) under Grant No. 2103627.

