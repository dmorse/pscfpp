
# PSCF - Polymer Self-Consistent Field (C++/CUDA)

PSCF is a software package for field-theoretic analysis of inhomogeneous 
equilibrium states of polymerc materials with constituents that tend to 
phase separate.  PSCF can be used to perform either self-consistent field 
theory (SCFT) calculations or some types of stochastic field-theoretic 
simulation (FTS) calculations.

The current version of PSCF is written primarily in C++, supplemented by 
CUDA to enable the use of a graphics processing unit (GPU). PSCF is 
distributed only in source code form, and so must be compiled by the 
user. 

## Methods: SCFT and PS-FTS

PSCF was originally designed for SCFT calculations, and provides an
extensive set of tools for this purpose .  The acronym PSCF stands for 
Polymer Self-Consistent Field, reflecting this origin.

PSCF can now also perform field theoretic simulation (FTS) calculations
that rely on a partial saddle-point approximation (PS-FTS). Such 
simulations can be performed using ether Brownian dynamics (BD) or 
Monte-Carlo (MC) sampling techniques. 

The fully-fluctuating formulation of polymer field theory (which is not 
implemented in PSCF) requires the stochastic sampling of complex-valued 
fields. The partial saddle-point approximation for FTS used in the current
version of PSCF is an approximation for incompressible models in which a 
Langrange multiplier pressure-like field that imposes a constraint on the 
total monomer density is approximated at a mean-field (or saddle-point) 
level, while fields that couple to composition fluctuations are allowed 
to freely fluctuate. The resulting approximation yields a theory that 
involves only real-valued fields.

## History

The current C++/CUDA version of PSCF originated as a completely rewritten
version of an older Fortran SCFT program of the same name. The Fortran 
PSCF version builds a single program that was designed for SCFT analysis 
of systems that can contain mixtures of linear block polymers and small
molecule solvents in a domain with periodic boundary conditions.  The 
current C++/CUDA version is intended to supersede the Fortran version, 
which is no longer being developed or actively maintained. The Fortran 
PSCF program is still available in a separate github.com repository at 
https://github.com/dmorse/pscf.  The current C++/CUDA version provides 
almost all of the capabilities of the Fortran PSCF program, and many
important new capabilities, as discussed below.

Functional differences between the current C++/CUDA version of PSCF and 
the legacy Fortran version include the following:

   - The current version is an extensible package of several programs
     that are designed for different types of spatial domain and/or 
     different types of hardware, implemented within a common framework.

   - The current version can perform stochastic PS-FTS calculations as
     well as SCFT.

   - The current version can treat systems that contain acyclic branched 
     polymers, while the Fortran program was designed for linear polymers.

   - The current version enables use of a graphics processing unit (GPU)
     to dramatically accelerate some applications.

## Programs

PSCF provides source code for the following three executable programs:

   - **pscf_1d** : The pscf_1d program is is designed to perform SCFT
     calculations for one-dimensional (1D) problems in Cartesian,
     cylindrical or spherical coordinates. A finite difference method is
     used to solve the underlying partial differential equation, known
     as the modified diffusion equation (MDE). This program can be used
     to treat problems involving flat and curved interfaces, as well as
     cylindrical and spherical micelles.

   - **pscf_pc** : The pscf_pc program can be used to perform SCFT and
     PS-FTS calculations for systems that are periodic in 1, 2 or 3 
     spatial dimensions, using standard CPU hardware. A pseudo-spectral 
     algorithm is used to solve the MDE. This program provides capabilities
     for SCFT calculations analogous to those of the older Fortran PSCF 
     program, as well as tools for PS-FTS calculations.  The suffix "pc" 
     stands for "periodic CPU".

   - **pscf_pg** : The pscf_pg program is a GPU-accelerated version
     of pscf_pc that can also perform SCFT and PS-FTS calculations for
     periodic systems. It is based on algorithms analogous to those used 
     in pscf_pc and provides almost identical features, but provides 
     higher performance for large systems. The suffix "pg" stands for 
     "periodic GPU".

## Features

All three PSCF programs are designed to treat an incompressible mixture
containing any number of block polymer, homopolymer and small molecule
(point-like) solvent molecular species. Polymer species may be modelled
using either the standard Gaussian model of polymer conformations as 
continuous random walks or (in pscf_pc and pscf_pg) using a discrete 
bead-spring model with harmonic springs.

Features relevant to both SCFT and PS-FTS (all programs):

  - Ability to treat acyclic branched block polymers

  - Ability to use canonical, grand-canonical or mixed statistical
    ensembles: Users may specify either a volume fraction or a chemical
    potential for each molecular species in a mixture

  - Examples of input files for many types of calculation and structures

  - Python tools for data analysis and file preparation

  - User and developer documentation provided as an integrated web manual

  - Well documented, open source code written in object oriented C++

Features for SCFT (all programs):

  - Efficient Anderson-mixing SCFT iteration algorithms

  - Parameter "sweeps": Continuation algorithms for sequences of SCFT 
    calculations with different parameters along a 1D path through 
    parameter space.

Features for SCFT or PS-FTS for periodic systems (pscf_pc and pscf_pg):

  - Pseudo-spectral solution of the modified diffusion equation

  - Choice of a continuous chain or a discrete bead-spring model for
    polymer conformations.

  - Periodic boundary conditions with 1, 2 or 3 dimensional periodicity

  - Unit cells with any possible 2D and 3D Bravais lattice system (e.g.,
    cubic, orthorhombic, monoclinic, etc.)

  - A companion Matlab package for visualization of periodic structures:
    [Polymer Visual](<https://github.com/kdorfmanUMN/polymer_visual/>)

Features for SCFT on periodic systems (pscf_pc and pscf_pg): 

  - Automatic optimization of unit cell parameters so as to minimize the 
    SCFT free energy density

  - Imposition of any user-selected space-group symmetry on SCFT solutions

  - Built-in database of symmetry operations for all 230 3D space groups
    and 17 2D plane groups 

  - Inhomogeneous density constraint (a "mask") for SCFT of systems in
    a confined geometry.

  - External fields 

  - Tools for thin polymer films (using a mask to represent confinement 
    and external fields for selective surface interactions)

Features for PS-FTS (pscf_pc and pscf_pg):

  - Brownian dynamics (BD) and Monte Carlo (MC) sampling algorithms

  - BD algorithms: Leimkuhler-Matthews and predictor-corrector BD step
    algorithms

  - MC algorithms: real-space and "smart"/force-bias MC moves 

  - Efficient new algorithms for adjusting the pressure field so as to 
    find partial saddle-point states (i.e., to impose incompressibility)

  - Tools for calculation of free energy differences by thermodynamic 
    integration, including the Einstein-crystal method

  - Parameter "ramps" in which one or more parameters change continuously
    during a simulation

  - "Analyzer" classes to compute and analyze quantities of physical
    interest, including the structure factor, order parameters used to
    identify phase transitions, and derivatives of the Hamiltonian
    needed to for thermodynamic integration calculations

  - Tools for performing data analysis either during a simulation or
    by postprocessing field trajectory files that are created during
    a simulation 

Non-features (provided by the Fortran program but not the current version):

  - The only capability of the older Fortran PSCF program that is *not* 
    provided by the current version is the ability to perform SCFT linear 
    susceptibility calculations for periodic ordered phases
    [A. Ranjan, J. Qin, and D. Morse, *Macromolecules* **41**, 942 (2008)].

## Getting the source code

The PSCF source code is maintained in the github repository

   <https://github.com/dmorse/pscfpp>.

The source code should be obtained by using a git version control system
client to clone this repository. To do so from the unix command line on a 
machine with a git client installed, enter the command:
```
git clone --recursive https://github.com/dmorse/pscfpp.git
```
Note the use of the --recursive option of the git clone command. Use of
this option is necessary with PSCF to clone and install two git submodules 
that are maintained in separate github repositories. This command will 
create a new directory called pscfpp that contains all of the source code 
and associated documentation, including the two required git submodules.

We do *not* recommend that users obtain the source code by simply
downloading and unpackng a zip or tar file of a tagged release from the 
PSCF github repository. Doing so would create a pscfpp directory that does
not contain source code for two git repositories that are automatically 
installed as submodules by the above git command. It is possible to install
these submodules after the fact, but simpler to follow the instructions 
given above.

## Documentation

PSCF is distributed with source files for an html web manual. A copy of
the manual for the most recent numbered release is available online at

  <https://dmorse.github.io/pscfpp-man>

All information given in this README file is also given in the web manual,
often in much greater detail. 

If desired, users can use the [doxygen](<https://www.doxygen.nl/index.html>)
documentation generation program to generate a local copy of the web manual.
The resulting html pages will be installed within the docs/html subdirectory
of the pscfpp root directory.  To do this, the doxygen application must be
installed on your computer, and the directory containing the doxygen 
executable must be in the unix command search PATH variable. To generate a 
local copy of the web manual:

   - Change directory (cd) to the pscfpp root directory

   - Enter "make html"

This should create many html files in the pscfpp/docs/html directory. To 
begin reading, point a web browser at the file pscfpp/docs/html/index.html ,
which is the main page of the web manual.

## Dependencies

PSCF has been developed and tested on both linux and and Mac OS X operating
systems, and is designed to be compiled on these or other unix-like systems.

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
cannot be run on Mac OS X, because no Apple computer uses the required 
type of GPU.

Procedures for installing the required dependencies are different for 
different operating system environments and for different package 
managers.  Instructions for some common environments are given in the 
web manual.

## Compiling

The PSCF source code is written using C++ as the primary language, with
CUDA used in pscf_pg for GPU acceleration. PSCF is only provided in source 
code format, and so must be compiled from source.

The build system used to compile and install PSCF uses the unix "make"
utility, which is available in any unix-like command-line environment.
The build system also uses some Python modules that are provided within 
the package, and so requires access to a Python 3 interpreter.  To compile 
and run this or other unix software on Mac OS X, the user must first 
install the Mac OS X unix command line tools.

Complete instructions for compiling PSCF are provided in section 2 of 
the [web manual](<https://dmorse.github.io/pscfpp-man>).  A brief summary 
is given below of instructions for steps that must be taken after cloning 
the pscfpp git repository and installing all of the required dependencies 
(i.e., GSL and FFTW) in a linux or Mac OS X environment.

   - Add the absolute path to the directory pscfpp/bin to your unix 
     command search PATH environment variable. This is where executable
     files will be installed by default. 

   - Add the absolute path to the pscfpp/lib/python to your PYTHONPATH
     environment variable. This allows a python interpreter to find python 
     modules that are used by the build system.

   - Run the pscfpp/configure script by entering "./configure" from the
     pscfpp root directory.  This script installs default versions of 
     several customizable configuration files that are used by the build 
     system. The configure script usually only needs to be run once, 
     before compiling for PSCF the first time.  See Section 2.6 of the 
     [web manual](<https://dmorse.github.io/pscfpp-man>) for more
     information about this script.

   - To compile and install only the CPU-based programs on a computer that
     does not have an appropriate NVIDA GPU, simply enter "make" from the
     pscfpp root directory immediately after running the configure script.

   - Some further configuration is required to compile the GPU-enabled
     pscf_pg program. Compilation of pscf_pg is only possible on a 
     machine that has an appropriate NVIDA GPU and a CUDA development kit. 
     To enable compilation of CUDA code, first enter "./setopts -c1" from 
     the pscfpp directory.  Then use the -a option of the same setopts
     script to set the correct GPU architecture for your system. For
     example, to compile for a V100 GPU with CUDA compute capability 7.0, 
     one would enter "./setopts -a sm_70". Instructions for choosing the 
     correct string parameter for the -a option for a particular GPU can 
     be obtained by entering "./setopts -h", or by consulting the 
     installation section of the web manual.  After taking these steps to 
     enable CUDA compilation, re-enter "make" to re-compile.

The procedure described above for building PSCF is somewhat similar to the 
standard "configure; make; make install" procedure used for GNU software.  
One difference is that PSCF does not require a separate "make install" 
command. The PSCF "make" command creates executables in the bin 
subdirectory of the pscfpp directory tree that also contains source
files.  The above instructions assume that this is the intended final 
destination of the executable files, rather than a standard unix system 
directory such as /user/local in which the unix shell would look for
executable files by default. If the pscfpp directory is created within 
a user's home directory, the above instructions thus do not require the 
use of special sudo or adminstrator permission in order to comple and
install the package. This procedure, however, instead requires the user 
to modify their PATH environment variable to allow the unix shell to find 
executable files that have been placed in a non-standard location.

## Command line usage

PSCF is a package the provides three executable programs (pscf_1d, pscf_pc,
and pscf_pg) with different capabilities but very similar command line 
interfaces. To perform a calculation, each PSCF program must read two files,
which we refer to as a parameter file and a command file. The parameter 
file, which is processed first, is a fixed format file that contains 
physical and algorithmic parameters required to describe the system of 
interest and initialize the program.  The command file, which is processed 
after the parameter file, is a script that contains a sequence of commands 
that are interpreted and executed in the order that they appear. Contents 
and syntax for PSCF parameter and command files are discussed in Secs. 3-5 
of the web manual.

A PSCF command file normally contains the names of several input and 
output data files as arguments to commands that read or write these files.
Specifically, a command file to perform either a SCFT or PS-FTS calculation
normally contains a command to read an input file that contains an initial 
set of monomer chemical potential fields.  This initial state may be used
in SCFT as an initial guess for an iterative solver, or in PS-FTS as the
first state of a Markov chain.  A command file for either an SCFT or FTS 
calculation also often contains commands to output final chemical potential
and monomer concentration fields to output files (i.e., either a converged 
solution for a SCFT calculation or the fields obtained after the last step 
of a stochastic simulation).

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
Here, "param" and "command" again denote names of parameter and command 
files. The syntax for an SCFT calculation for a one-dimensionally periodic 
lamellar phase would instead contain an option "-d 1", while an SCFT
calculation for a 2D periodic phase of hexagonally packed cylinders would
use an option "-d 2".  The syntax for invoking pscf_pg is the same as that 
for pscf_pc, except for the use of the program name pscf_pg rather than
pscf_pc.

The pscf_pc and pscf_pg programs can also both perform PS-FTS simulations. 
Such simulations usually allow field fluctuations to vary in all three 
spatial dimensions, and so pscf_pc and pscf_pg are normally invoked using 
an option "-d 3" when used for this purpose.

**Output redirection** : Programs that are invoked as shown in the above 
examples would write log output that is produced during execution to the 
user's screen (i.e., to standard output).  This log output may instead be
redirected to a file by using the unix ">" redirection operator.  For 
example, the command
```
pscf_pc -d 3 -p param -c command -e > log
```
would redirect log output that is written during a calculation to a file
named "log" in the current working directory. Standard output is normally 
redirected to a file for long calculations, or whenever a program is run 
in a queue on a shared computing cluster.

## Examples

The quickest way to become familiar with PSCF parameter and command files 
is by studying examples.  The directory pscfpp/examples contains input 
files for examples of a variety of different types of SCFT and PS-FTS 
calculations.  Top level subdirectories of the examples directory named 
1d, pc and pg contain examples for the three different PSCF programs,
pscf_1d, pscf_pc, and pscf_pg.

Subdirectory examples/1d contains examples of SCFT calculations for the 
1D finite-difference program pscf_1d. Top level subdirectories of the 
directory examples/1d contain examples for planar, cylindrical and 
spherical geometries, as indicated by the subdirectory names. One or more 
example is given for each geometry.

Subdirectory examples/pc contains examples for the pscf_pc CPU program. 
Top level subdirectories of examples/pc named scf and fts contain examples
of SCFT and PS-FTS calculations, respectively. Top level subdirectories 
of directory examples/pc/scf contain examples of input files for SCFT 
calculations performed on different types of physical system.  For example,
the directory examples/pc/scf/diblock contains SCFT examples for a diblock
copolymer melt.  Subdirectories of examples/pc/scf/diblock contain examples
for lamellar, hexagonal, and BCC structures, among others.

Subdirectory examples/pg contains examples for the pscf_pg program for
periodic structures. Subdirectories are organized in a manner similar to 
that of the examples/pc directory tree.

We refer to a directory that contains all of the input files for a single
example or a set of closely related examples as an example directory.
Each such example directory contains at least one parameter file (usually
named "param"), at least one command file (usually named "command"), and
an initial chemical potential field (w field) file. Example directories 
that contain input files for two or more examples usually contain a single 
initial chemical potential field, but may contain two or more different 
sets of parameter or command files for use in different examples. 
By convention, input field files are located in a subdirectory of each
example directory named "/in", while output files are created in a 
subdirectory named "/out". 

Example directories that contain files for a single example usually
contain a shell script named "run" that can be executed to run the example
using the supplied input files.  The simplest way to run such an example
is thus to change directory (cd) to the relevant example directory and 
enter "./run" from within that directory. (Note the use of the prefix "./",
which tells the operating system to look for the run script in the current 
working directory.) Users may inspect the text of such run scripts to see 
the command and command line options that are being used to run an example.
Example directories that contain input files for two or more closely 
related examples may contain several such scripts with names that are 
variants of "run", each of which can be executed to run a specific 
example. Many such directories also contain a file named CONTENTS that 
explains the purposes of different run scripts.

Almost every example directory also contains a script named "clean" that 
can be used to remove all of the output files that are created by running 
the example. To remove such files, simply enter "./clean". 

## Contributors

- David Morse
- Guo Kang Cheong
- Anshul Chawla
- Ryan Collanton
- Ben Magruder
- Kexin Chen
- Ying Zheng

## License and copyright

This C++/CUDA version of PSCF is free, open source software. It is
distributed under the terms of the GNU General Public License (GPL) as 
published by the Free Software Foundation, either version 3 of the 
License or (at your option) any later version.  PSCF is distributed 
without any warranty, without even the implied warranty of merchantability 
or fitness for a particular purpose.  See the LICENSE file or the
<a href=http://www.gnu.org/licenses/>gnu web page</a> for details.

Copyright 2015-2025, The Regents of the University of Minnesota.

## Support

Development of PSCF is currently supported by the National Science
Foundation program for Cyberinfrastructure for Sustained Scientific
Development (CSSI) under Grant No. 2103627.

