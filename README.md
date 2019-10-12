
# Pscf++

Pscf++ - C++ programs for polymer self-consistent field theory 

## Overview

Pscf++ is a library of C++ classes and programs for solving the 
Edwards-Helfand self-consistent field theory for polymer liquids. 
Pscf++ is still under development, but is intended to eventually
supersede the existing PSCF Fortran program. 

Differences from the existing PSCF Fortran code and expected advantages 
of the new code include:

   - Pscf++ is an extensible package of several different programs 
     designed for use with different geometries and boundary conditions, 
     different algorithms or different hardware, designed around a 
     common framework. 

   - Pscf++ allows simulations of mixtures containing arbitrary acyclic 
     branched copolymers, in addition to the linear block copolymers and 
     homopolymers allowed by Fortran PSCF code.

   - An object oriented design allows creation of independent objects
     to represent different phases, which will facilitate analysis of 
     phase coexistence.

   - Adoption of C/C++ has simplified implementation of SCFT solvers
     that uses graphical process units (GPUs).

Pscf++ is free, open source software. It is distributed under the terms 
of the GNU General Public License as (GPL) published by the Free Software 
Foundation, either version 3 of the License or (at your option) any later 
version.  Pscf++ is distributed without any warranty, without even the 
implied warranty of merchantability or fitness for a particular purpose. 
See the LICENSE file or the 
<a href=http://www.gnu.org/licenses/> gnu web page </a> for details.

## Programs

Currently, the package contains the following SCFT solvers:

   - A one-dimensional finite difference solver for problems that involve 
     variation in a single spatial coordinate in Cartesian, cylindrical 
     or spherical coordinates.

   - A CPU-based pseudo-spectral solver for 1, 2 , or 3 dimensionally 
     periodic structures.

   - A GPU accelerated pseudo-spectral solver for period structures. 

The one-dimensional finite different solver is useful for treating problems
involving flat or curved interfaces, and cylindrical or spherical micelles.
The executable for the this program is named pscf_fd.

The CPU-based pseudo-spectral solver for periodic microstructures is
similar in most respects to the existing PSCF Fortran program, and 
provides very similar performance. Like the Fortran program, it allows 
the user to search for a solution with a specified space group symmetry.
The new program can read and write the same file formats for representing 
a field in terms of symmetry-adapted basis functions as those used by 
the PSCF Fortran program, and provides very similar perfomance. The
new code does, however, still lack a few features of the original Fortran 
code.  Currently, the most important missing features are the absence of 
the "sweep" continuation feature and the lack of specialized code to 
simulate point-particle solvents. The CPU programs for solving 1, 2 and 
3 dimensionally periodic structures are named pscf_pc1d, pscf_pc2d and 
pscf_pc3d, respectively, where "pc" stands for "periodic CPU".

The GPU-accelerated pseudo-spectral solver for periodic structures is 
based on the same algorithms as CPU pseudo-spectral solver, but is 
less mature. The most important difference in features is that the 
GPU-accelerated code does yet allow the user to use symmetry -adapted
basis functions to constrain the space group symmetry for the solution.
The GPU accelerated programs for solving 1, 2 and 3 dimensionally 
periodic structures are named pscf_pg1d, pscf_pg2d and pscf_pg3d, 
respectively, where "pg" stands for "periodic GPU".

## Getting the Source Code

The pscf++ source code is maintained in the github repository

   <https://github.com/dmorse/pscfpp>.

It may be obtained by using a git version control system client to
clone the repository. To do so, enter the command:

   git clone --recursive https://github.com/dmorse/pscfpp.git

The use of the --recursive option to the git clone command:
This is necessary to clone some git submodules that are maintained
in separate repositories. This command will create a directory 
called pscfpp/ that contains all of the source code and associated
documentation.

## Documentation

Pscf++ is distributed with source files for an html web manual.
After cloning the source code, you can use the doxygen documentation
generator to generate a local copy of this documentation. To do this,
doxygen must be installed on your computer, and the directory 
containing the doxygen executable must be in your command search
PATH. To generate documentation:

   - Change directory (cd) to the pscfpp/ root directory

   - Enter "make html"

This should create many html files in the pscfpp/doc/html directory.
To begin reading the documentation, point a browser at the file
pscfpp/doc/html/index.html, which is the main page of the manual.

## Dependencies

The pscf++ source code is written in a combination of C++ and (for
the GPU accelerated program) Cuda, and must be compiled from source.
The package was developed on linux and and Mac OS X operating systems 
using standard unix utilities, and is designed to run on these 
systems. To compile linux-like software on a Mac, you must first 
install the XCode Mac development environment and the unix command 
line tools.  

The CPU-based programs within the pscf++ package depend on the 
following external libraries:

  - Gnu scientific library (GSL)

  - FFTW fast Fourier transform library

The one-dimensional finite difference program pscf_fd1d requires 
only GSL, and not FFTW. The CPU-based programs for spatially
periodic structures require both GSL and FFTW libraries.

The GPU-accelerated programs can only run on a computer with an
appropriate nVidia graphics card. To compile these programs, the
system must also have an nVidia cuda development environment 
that provides the CUFFT fast Fourier transform library. 

## Compiling

Complete directions for compiling and installing pscf++ are
provided in section 2 of the html documentation. Short instructions
for compiling, after installing all of the required dependencies,
are given below:

- Add the pscfpp/bin directory to your linux command search PATH
  environment variable.

- Add the pscfpp/scripts/python directory to your PYTHONPATH
  environment variable.

- cd to the pscfpp/ root directory

- Enter "./setup" from this root directory to run a setup script
  (you only need to do this once, before compiling the first time).
- Change directory (cd) to the directory pscfpp/bld/.

- To compile and install all CPU-based programs in the package 
  (excluding GPU-accelerated programs), enter "make all-cpu"

- To compile the GPU-accelerated programs on a machine with an
  nVidia GPU, a Cuda compiler and the CUFFT library, enter
  "make pspg". 

The setup script installs a file pscfpp/bld/config.mk that contains
makefile variables that define compiler executable names, compiler 
options and paths to library files. If the default options are not
adequate, may edit this file as needed.

------------------------------------------------------------------------
