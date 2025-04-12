
# PSCF - Polymer Self-Consistent Field (C++/CUDA)

PSCF is a software package for field-theoretic analysis of inhomogeneous 
equilibrium states of polymer materials with components that tend to phase
separate.  PSCF can perform either self-consistent field theory (SCFT) or
some types of stochastic field-theoretic simulation (FTS).

The current version of PSCF is written primarily in C++, supplemented by 
CUDA to enable the use of a graphics processing unit (GPU). See below
for a discussion of the relationship to an older Fortran program of the 
same name. 

## SCFT and FTS

PSCF was originally designed for SCFT calculations, and provides an 
extensive set of tools for this purpose.  The acronym PSCF stands for 
Polymer Self-Consistent Field, reflecting this origin.

The current version of PSCF (v1.2) is the first to also provide tools for
stochastic field theoretic simulations (FTSs) that rely on a partial 
saddle-point approximation (PS-FTS). Such simulations can be performed 
using ether Brownian dynamics (BD) or Monte-Carlo (MC) sampling algorithms. 

The fully-fluctuating formulation of polymer field theory (which is not
yet implemented in PSCF) requires stochastic sampling of complex-valued 
fields. The partial saddle-point approximation for FTS used in this 
version of PSCF is an approximation for incompressible systems in which 
a pressure-like field that imposes a constraint on the total monomer 
concentrated is approximated at a mean-field (or saddle-point) level, 
while field or fields that couple to composition fluctuations are allowed 
to freely fluctuate. The resulting approximate theory involves only 
real-valued fields.

## History

The current C++/CUDA version of PSCF originated as a completely rewritten
version of an older Fortran SCFT program of the same name. The Fortran 
PSCF version is a single program that was designed for SCFT treatment
of systems that can contain any mixture of linear block polymers and 
small-molecule solvents in a domain with periodic boundary conditions. 
The current C++/CUDA version is intended to supersede the Fortran version, 
which is no longer being developed or actively maintained. The Fortran 
PSCF program is still available in a separate github.com repository at 
https://github.com/dmorse/pscf.  The current C++/CUDA version provides 
almost all of the capabilities of the Fortran PSCF program, and several
important new capabilities, as discussed below.

Functional differences between the current C++/CUDA version of PSCF and 
the legacy Fortran version include:

   - The current version is an extensible package of several programs
     that are designed for use with different types of spatial domain,
     and/or different hardware, implemented within a common framework.

   - The current version can perform stochastic PS-FTS calculations as
     well as SCFT.

   - The current version allows analysis of mixtures containing acyclic 
     branched block polymers, while the Fortran program can only treat 
     linear polymers.

   - The current version enables use of a graphics processing unit (GPU)
     to dramatically accelerate some applications.

The only important capability of the Fortran PSCF program that has not 
been ported to the C++/CUDA version is the ability to perform generalized
random-phase approximation (RPA) calculations for periodic ordered phases. 
Such calculations can be used to characterize composition fluctuations 
and SCFT limits of stability for ordered phases.

## Programs

PSCF currently provides three executable programs:

   - **pscf_1d** : The pscf_1d program is is designed to perform SCFT
     calculations for one-dimensional (1D) problems in Cartesian,
     cylindrical or spherical coordinates. A finite difference method is
     used to solve the underlying partial differential equation, known
     as the modified diffusion equation (MDE). This program can be used
     to treat problems involving flat and curved interfaces, as well as
     cylindrical and spherical micelles.

   - **pscf_pc** : The pscf_pc rogram is designed to perform SCFT and
     PS-FTS calculations for systems that are periodic in 1, 2 or 3 
     spatial dimensions, using standard CPU hardware. A pseudo-spectral 
     algorithm is used to solve the MDE. This program provides capabilities
     for SCFT calculations analogous to those of the older Fortran PSCF 
     program, as well as algorithms for PS-FTS simulations.  The suffix 
     "pc" stands for "periodic CPU".

   - **pscf_pg** : The pscf_pg program is a GPU-accelerated version
     of pscf_pc that can also perform SCFT and PS-FTS calculations for
     periodic systems. It is based on algorithms similar to those used 
     in pscf_pc and provides almost identical features, but provides 
     higher performance for large systems. The suffix "pg" stands for 
     "periodic GPU".

## Features

All three PSCF programs are designed to treat an incompressible mixture
containing any number of block polymer, homopolymer and small molecule
(point-like) solvent molecular species. All polymer species are treated 
using the standard Gaussian model of polymer conformations as continuous
random walks.

Features for both SCFT and FTS (all programs):

  - Ability to treat mixtures of any number of block polymer and solvent
    species.

  - Ability to treat acyclic branched block polymers of arbitrary
    topology, in addition to linear polymers

  - Ability to use canonical, grand-canonical or mixed statistical
    ensembles: Users may specify either a volume fraction or a chemical
    potential for each molecular species in a mixture

  - Examples of input files for many types of calculation and structures

  - Python tools for data analysis

  - Thorough user and developer documentation provided as an integrated
    web manual

  - Well documented, open source code written in object oriented C++

Features for SCFT (all programs):

  - Efficient Anderson-mixing SCFT iteration algorithms

  - Efficient SCFT calculations for sequences of parameter choices along
    a path through parameter space (parameter "sweeps")

Features for SCFT or PS-FTS for periodic systems (pscf_pc and pscf_pg):

  - Accurate pseudo-spectral solution of the modified diffusion equation

  - Periodic boundary conditions with 1, 2 or 3 dimensional periodicity

  - Unit cells with any possible 2D and 3D Bravais lattice system (e.g.,
    cubic, orthorhombic, monoclinic, etc.)

  - A companion Matlab package
    [Polymer Visual](<https://github.com/kdorfmanUMN/polymer_visual/>)
    for visualization of periodic structures

Features for SCFT on periodic systems (pscf_pc and pscf_pg): 

  - Automatic optimization of unit cell parameters so as to minimize the 
    SCFT free energy density

  - Imposition of any user-selected space-group symmetry on SCFT solutions

  - Built-in database of symmetry operations for all 230 3D space groups
    and 17 2D plane groups 

  - Inhomogeneous density constraint (a "mask") for SCFT of systems in
    a confined geometry (e.g., thin films)

  - External fields 

  - Tools for treating thin polymer films (implemented using a mask to
    represent confinement and external fields to represent selective
    interactions with surfaces).

Features for PS-FTS (pscf_pc and pscf_pg):

  - Brownian dynamics (BD) and Monte Carlo (MC) sampling 

  - BD algorithms: Leimkuhler-Matthews and predictor-corrector BD step
    algorithms

  - MC move algorithms: real-space and "smart"/force-bias MC moves 

  - Efficient algorithms for adjusting the pressure field so as to 
    find a partial saddle-point (i.e., to impose incompressibility)

  - Tools for calculation of free energy differences by thermodynamic 
    integration, including the Einstein-crystal integration method

  - Parameter "ramps" in which one or more parameters change continuously
    during a simulation, which can be used for continuous thermodynamic 
    integration

  - "Analyzer" classes to compute and analyze quantities of physical
    interest, including the structure factor, order parameters used to
    identify phase transitions, and derivatives of the Hamiltonian
    needed to for thermodynamic integration calculations

  - Tools for performing data analysis either during a simulation or
    by postprocessing field trajectory files that are created during
    a simulation 

## Getting the source code

The PSCF source code is maintained in the github repository

   <https://github.com/dmorse/pscfpp>.

The source code may be obtained by using a git version control system
client to clone this repository. To do so on a machine with a git client,
enter the command:
```
git clone --recursive https://github.com/dmorse/pscfpp.git
```
Note the use of the --recursive option to the git clone command. This
option is necessary to clone two git submodules that are maintained in 
separate repositories. This command will create a new directory called 
pscfpp/ that contains all of the source code and associated documentation, 
including the required git submodules.

We do *not* recommend that users obtain the source code by simply
downloading and unpackng a zip or tar file of a tagged release from the 
PSCF github repository. Doing so would create a directory that does not 
contain source code for two git repositories that are automatically 
downloaded and installed as submodules by the above git command. It is 
possible to install the required submodules after the fact, but simpler 
to follow the instructions given above to clone the repository using 
the --recursive option. 

## Documentation

PSCF is distributed with source files for an html web manual. A copy of
the manual for the most recent numbered release is available online at

  <https://dmorse.github.io/pscfpp-man>

After cloning the source code, users can also use the
[doxygen](<https://www.doxygen.nl/index.html>) documentation generation
program to generate a local copy of the web manual.  To do this, the
doxygen application must be installed on your computer, and the directory
containing the doxygen executable must be in your unix command search
PATH. To generate local documentation:

   - Change directory (cd) to the pscfpp/ root directory

   - Enter "make html"

This should create many html files in the docs/html subdirectory of the
pscfpp/ root directory.  To begin reading the resulting documentation,
point a web browser at the file pscfpp/docs/html/index.html , which is
the main page of the manual.

## Dependencies

PSCF has been developed and tested on both linux and and Mac OS X 
operating systems, and is designed to run on these or other unix-like 
systems. 

The pscf_1d and pscf_pc CPU programs depend on the following external
libraries:

  - [GSL](<https://www.gnu.org/software/gsl/>) - GNU Scientific Library

  - [FFTW](<https://www.fftw.org/>) Fast Fourier Transform library

The pscf_1d one-dimensional program pscf_1d requires only GSL, while
the pscf_pc program relies on both GSL and FFTW.

The GPU-accelerated pscf_pg program can only be compiled and run on a
computer with an appropriate Nvidia GPU and an Nvidia CUDA development
kit. Some GPU-accelerated libraries used by pscf_pg, such as the FFT
library cuFFT, are provided as part of this kit. The pscf_pg program 
cannot be run on an Apple computer, because no Apple computer uses 
the required type of GPU. 

Procedures for installing the required dependencies are different for 
different operating system environments and for different package 
managers.  Instructions for some common environments are given in the 
web manual.

## Compiling

The PSCF source code is written using C++ as the primary language, with
CUDA used in pscf_pg for GPU acceleration. PSCF is only provided in source 
code format, and so must be compiled from source.

The build system used to compile and install PSCF relies on standard unix
utilities that are available in any unix-like command-line environment, and 
also requires access to a Python 3 interpreter.  To compile and run this or
other unix software on Mac OS X, the user must first install the Mac OS X 
unix command line tools.

Complete directions for compiling and installing PSCF are provided in
section 2 of the [web manual](<https://dmorse.github.io/pscfpp-man>).
A brief summary is given below of instructions for steps that must be
taken after cloning the git repository and installing all of the
required dependencies:

   - Add the directory path pscfpp/bin to your unix command search PATH
     environment variable. This is where executables will be installed 
     by default. 

   - Add the directory path pscfpp/lib/python to your PYTHONPATH
     environment variable. This allows a python interpreter to find 
     python modules that are used by the build system.

   - Run the pscfpp/configure script by entering "./configure" from the
     pscfpp/ root directory.  This script installs default versions of 
     several files that are used by the build system. The configure 
     script usually only needs to be run once, before compiling for the 
     first time.  See Section 2.6 of the 
     [web manual](<https://dmorse.github.io/pscfpp-man>) for more
     information about the configure script.

   - To compile and install only the CPU-based programs on a computer 
     that does not have an appropriate NVIDA GPU, simply enter "make" 
     from the pscfpp/ root directory immediately after running the 
     configure script.

   - Some further configuration is required to compile the GPU-enabled
     pscf_pg program. Compilation of pscf_pg is only possible on a
     machine that has an appropriate NVIDA GPU and a CUDA development kit. 
     To enable compilation of CUDA code, first enter "./setopts -c1" from 
     the pscfpp/ directory.  Then use the -a option of the same setopts
     script to set the correct GPU architecture for your system. For
     example, to compile for a V100 GPU with CUDA compute capability 7.0, 
     one would enter "./setopts -a sm_70". Instructions for choosing the 
     correct string parameter for the -a option for a particular GPU can 
     be obtained by entering "./setopts -h" or by consulting the 
     installation section of the web manual.  After these steps, enter 
     "make all" to compile.

The procedure described above for building PSCF is similar to the standard
"configure; make; make install" procedure for GNU open source software. 
The main difference is that PSCF does not require a separate "make install" 
command. This is because the PSCF "make" command is designed to create 
executables in the bin/ subdirectory of the root pscfpp/ directory in 
which the package is configured. The above instructions assume that 
this is their intended final destination.  If the pscfpp/ directory is 
installed within a user's home directory, the above instructions do not 
require the use of sudo or adminstrator permission, but instead do require 
the user to modify their PATH variable to allow the unix shell to find 
executables installed within user's home directory.

## Command line usage

PSCF is a package the provides three executable programs (pscf_1d,
pscf_pc, and pscf_pg) with somewhat different capabilities but very
similar command line interfaces. To perform a calculation, each PSCF
program must read a parameter file and a command file. The parameter
file, which is processed first, is a fixed format file that contains
parameters required to describe the physical system of interest and
initialize the program.  The command file, which is processed after
the parameter file, is a script that contains a sequence of commands
that are read and executed in the order that they appear. Contents and
syntax for PSCF parameter and commands files are discussed in Secs. 
3-5 of the web manual.

The command file usually contains the names of several input and output
data files as arguments to commands that read or write these files.
Specifically, a command file to perform either a SCFT or PS-FTS calculation
normally contains a command to read an input file that contains an
initial guess for the monomer chemical potential fields.  A command
file for either an SCFT or FTS calculation also often contains
commands to output final chemical potential and monomer concentration
fields to output files, i.e., either converged solutions for SCFT or 
the final fields obtained after the last step of a FTS.

Command line usage for different programs is described below, starting
with the simple case of the pscf_1d program for one-dimensional SCFT
calculations:

**pscf_1d** : The usual command line syntax for invoking the pscf_1d 
program is
```
pscf_1d -p param -c command -e
```
where the string "param" denotes the name of a parameter file, and
"command" denotes the name of a command file.  The -p and -c options
are required.

The "-e" command line option in the above example is optional but 
recommended. This option causes the program to "echo" the parameter file 
to standard output as this file is being processed. Use of this option
makes it easier to debug failures that arise from syntax errors in 
the parameter file. 

**pscf_pc and pscf_pg** : The command line interfaces for the pscf_pc and
pscf_pg programs take the same required and optional elements as pscf_1d, 
but also require a value for the dimension of space as an additional 
parameter.  This is an integer parameter of the -d option that must have 
a value 1, 2 or 3, which specifies the number of coordinates along which 
the structure is periodic.  The usual syntax for invoking pscf_pc for an
SCFT calculation of a three dimensionally periodic structure (e.g., a 
network structure or a 3D arrangement of spheres), while also using the 
-e option, is thus 
```
pscf_pc -d 3 -p param -c command -e
```
where "param" and "command" again denote names of parameter and command 
files. The syntax for an SCFT calculation for a one-dimensionally periodic 
lamellar phase would instead contain an option "-d 1", while an SCFT
calculation for a 2D periodic phase of hexagonally packed cylinders would
use an option "-d 2".  The syntax for invoking pscf_pg is the same as that 
for pccf_pc, except for the use of the different program name.

The pscf_pc and pscf_pg programs can also perform PS-FTS simulations. 
Such simulations usually allow field fluctuations vary in all three 
spatial dimensions, and so pscf_pc and pscf_pg are normally invoked 
using an option "-d 3" when used for this purpose.

**Output redirection** : Programs that are invoked as shown in the above 
examples would write log output that is produced during execution to the 
user's screen (i.e., to standard output).  This log output may instead be
redirected to a file by using the unix ">" operator to direct standard 
output to a file.  For example, the command
```
pscf_pc -d 3 -p param -c command -e > log
```
would redirect log output that is written during simulation of a 3D
periodic structure to a file named "log" in the current working 
directory. Standard output is normally redirected to a file for long 
calculations, or whenever a program is run in a queue on a shared 
resource.

## Examples

The quickest way to become familiar with PSCF parameter and command
files is by examining examples.  The directory pscfpp/examples contains
input files for examples of a variety of different types of SCFT and 
PS-FTS calculations.  Top level subdirectories of the examples/ directory
named 1d/, pc/ and pg/ contain examples for different PSCF programs
(pscf_1d, pscf_pc, or pscf_pg).

Subdirectory examples/1d contains examples of SCFT calculations for the
1D finite-difference program pscf_1d. Top level subdirectories of the
directory examples/1d contain examples for planar, cylindrical and
spherical geometries, as indicated by the subdirectory names. One or
more example is given for each geometry.

Subdirectory examples/pc contains examples for the pscf_pc CPU program.
Top level subdirectories of examples/pc named scf/ and fts/ contain
examples of SCFT and FTS calculations, respectively. Top level
subdirectories of examples/pc/scf contain examples of input files for
SCFT calculations performed on different types of physical system.
For example, the directory examples/pc/scf/diblock contains examples
for a diblock copolymer melt.  Subdirectories of examples/pc/scf/diblock
contain examples for lamellar, hexagonal, and BCC structures, among
others.

Subdirectory examples/pg contains examples for the pscf_pg program for
periodic structures. Subdirectories are organized in a manner similar
to that used for the examples/pc directory tree.

We refer to a directory that contains all of the input files for a single
example or a set of closely related examples as an example directory.
Each such example directory contains at least one parameter file (usually
named "param"), at least one command file (usually named "command"), and
an initial chemical potential field (w field) file. Example directories 
that contain input files for two or more examples usually contain a single 
initial chemical potential field, but may contain two or more different 
sets of parameter or command files for use in different examples. 
By convention, input field files are located in a subdirectory of each
example directory named "\in", while output files are created in a 
subdirectory named "\out". 

Example directories that contain files for a single example usually
contain a shell script named "run" that can be executed to run the example
using the supplied input files.  The simplest way to run such an example
is thus to change directory (cd) to the relevant example directory and 
enter "./run" from within that directory. (Note the use of the prefix
prefix "./", which tells the operating system to look for the run script
in the current working directory). Users may also inspect the text of such
a run script to see the command and command line options that are being
used to run the example.  Example directories that contain input files
for two or more closely related examples may contain several such scripts
with names that are variants of "run", each of which can be executed to 
run a specific example. Such directories also usually contain a file 
named CONTENTS that explains the purposes of different run scripts.

Almost every example directory also contains a script named "clean" that 
can be used to remove all of the output files that are created by running 
the example.

## License

This C++/CUDA version of PSCF is free, open source software. It is
distributed under the terms of the GNU General Public License (GPL) as 
published by the Free Software Foundation, either version 3 of the 
License or (at your option) any later version.  PSCF is distributed 
without any warranty, without even the implied warranty of merchantability 
or fitness for a particular purpose.  See the LICENSE file or the
<a href=http://www.gnu.org/licenses/>gnu web page</a> for details.

## Support

Development of PSCF is currently supported by the National Science
Foundation program for Cyberinfrastructure for Sustained Scientific
Development (CSSI) under Grant No. 2103627.

## Contributors

- David Morse
- Guo Kang Cheong
- Anshul Chawla
- Ryan Collanton
- Ben Magruder
- Kexin Chen
- Ying Zheng

