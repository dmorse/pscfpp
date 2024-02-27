
# PSCF - Polymer Self-Consistent Field (C++/CUDA)

PSCF is a package of software for field-theoretic simulation of
inhomogeneous equilibrium states of polymer liquids, with a focus
on microstructures formed by block polymers. PSCF can perform both
self-consistent field theory (SCFT) calculations and some types of
stochastic field-theoretic simulation (FTS).

This current version of PSCF is written primarily in C++, supplemented
by CUDA code to enable the use of a graphics processing unit (GPU).
See below for a discussion of the relationship to an older Fortran
program of the same name. 

## SCFT and FTS

PSCF was originally designed for SCFT calculations, and provides an
extensive set of tools for this.  The acronym PSCF stands for Polymer
Self Consistent Field, reflecting this origin.

The current version of PSCF (v1.2) is the first to also provide a set
of tools for field theoretic simulations (FTS) based on Monte-Carlo (MC)
and Brownian dynamics (BD) sampling algorithms. The FTS algorithms in
the current version all rely the partial saddle-point approximation
(PSPA) to the exact of "fully fluctuating" formulation of the partititon
function as a functional integral.

The fully-fluctuating formulation of polymer field theory (which is
not yet implemented in PSCF) requires the stochastic sampling of
complex-valued fields. The PSPA used here is an approximation for
incompressible models in which the Lagrange multiplier pressure field
that imposes a constraint on the total monomer density is approximated
at a mean-field (i.e., saddle-point) level, while fields that couple
to composition fluctuations are allowed to fluctuate. The resulting
approximation yields a theory that involves only real-valued fields.

## History

The current C++/CUDA version of PSCF is intended to supersede an older
Fortran SCFT program of the same name. The Fortran PSCF program is
maintained in a separate github.com repository at
https://github.com/dmorse/pscf.
The current C++/CUDA version provides almost all of the capabilities of
the older Fortran program, and several important new capabilities, as
discussed below.

Differences between the current C++/CUDA version of PSCF and the legacy
Fortran version include:

   - The current version is an extensible package of several different
     programs designed for use with different types of spatial domain,
     different algorithms or different hardware, implemented within a
     common framework.

   - The current version enables simulations of mixtures containing
     both linear and acyclic branched block polymers, while the Fortran
     PSCF program allowed only linear polymers.

   - The current version enables use of a graphics processing unit (GPU)
     to dramatically accelerate some applications.

   - PSCF can now perform FTS calculations that rely on the PSPA.  A
     variety of MC and BD sampling algorithms have been implemented.

The most important capability of the Fortran PSCF program that has not 
been ported to the current version is the ability to perform generalized
random-phase approximation (RPA) calculations for ordered phases to 
characterize composition fluctuations and stability limits for such 
phases.

## Programs

PSCF currently contains three programs:

   - **pscf_1d** : The program pscf_1d is is designed to perform SCFT
     calculations on one-dimensional (1d) problems in Cartesian,
     cylindrical or spherical coordinates. A finite difference method is
     used to solve the underlying partial differential equation, known
     as the modified diffusion equation (MDE). This program is useful
     for treating problems involving flat and curved interfaces, as well
     as cylindrical or spherical copolymer micelles.

   - **pscf_pc** : The pscf_pc rogram is designed to treat structures
     that are periodic in 1, 2 or 3 dimensions, using conventional
     CPU hardware. A pseudo-spectral algorithm is used to solve the
     MDE. This program provides capabilities for SCFT calculations
     analogous to those of the older PSCF Fortran program, as well as
     code for FTS calculations based on the PSPA.  The suffix "pc"
     stands for "periodic CPU".

   - **pscf_pg** : The pscf_pg program is a GPU-accelerated version
     of pscf_pc that can also perform SCFT and FTS calculations for
     periodic systems. It is based on algorithms similar to those
     used in pscf_pc and provides almost identical features, but
     provides much higher performance for large systems. The suffix
     "pg" stands for "periodic GPU".

## Features

All PSCF programs are designed to treat an incompressible mixture
containing any number of block polymers, homopolymers and small molecule
(point-like) solvent molecular species. All polymeric species are
treated using the standard Gaussian model of each polymer block as a
continuous random walk.

Features that are common to all three PSCF programs, all of which apply
to both SCFT and FTS calculations, include:

  - Ability to treat mixtures of any number of block polymers and
    solvent species.

  - Ability to treat complex acyclic branched block polymers, in addition
    to linear block polymers and homopolymers

  - Ability to use canonical, grand-canonical or mixed statistical
    ensembles - users may specify either a volume fraction or a chemical
    potential for each molecular species

  - Examples of input files for many types of calculation and structures

  - Python tools for data analysis

  - Thorough user and developer documentation provided as an integrated
    web manual

  - Well documented open source code written in object oriented C++

Features used in SCFT calculations that are common to all three PSCF
programs include:

  - Efficient Anderson-mixing iteration algorithms

  - Efficient SCFT calculations for sequences of parameter choices along
    a path in parameter space ("sweeps"), using extrapolation to
    construct initial guesses

Features specific to the pscf_pc and pscf_pg programs for periodic systems
include:

  - Ability to perform both SCFT and FTS calculations

  - Accurate pseudo-spectral solution of the modified diffusion equation

  - Periodic boundary conditions with 1, 2 or 3 dimensional periodicity

  - Unit cells with all possible 2D and 3D Bravais lattice systems (e.g.,
    cubic, orthorhombic, monoclinic, etc.)

  - Automatic optimization of unit parameters in SCFT so as to minimize
    free energy density

  - Imposition of any user-selected space-group symmetry in SCFT

  - Built-in database of symmetry operations for all 230 3D space groups
    and 17 2D plane groups for use in SCFT

  - Availability of a companion package
    [Polymer Visual](<https://github.com/kdorfmanUMN/polymer_visual/>)
    for visualization of periodic structures

Features specific to SCFT calculations performed with the pscf_pc CPU
program are:

  - Inhomogeneous density constraints (a "mask")

  - External fields

  - Thin polymer films

The mask and external field features are used in pscf_pc to implement
simulations of thin films - A mask is used to constrain a polymer
material to a slit within a periodic supercell, and localized external
fields are used to represent selective interactions with the top and
bottom surfaces.

The pscf_pc program also provides a few utilities for manipulating input
field files that have not been reproduced in pscf_pg. Because the two 
programs use the same field file formats, manipulations of input files 
used by pscf_pg can be performed with pscf_pc. 

## Getting the source code

The current PSCF source code is maintained in the github repository

   <https://github.com/dmorse/pscfpp>.

The source code may be obtained by using a git version control system
client to clone the repository. To do so on a machine with a git client,
enter the command:
```
git clone --recursive https://github.com/dmorse/pscfpp.git
```
Note the use of the --recursive option to the git clone command. This
option is necessary to clone several git submodules that are maintained
in separate repositories. This command will create a new directory
called pscfpp/ that contains all of the source code and associated
documentation, including all required git submodules.

We do *not* recommend that users obtain the source code by simply
downloading a zip or tar file of a tagged release from the PSCF github
repository, because this file will not contain source code for the other
git repositories that are automatically downloaded and installed as
submodules by the above git command. It is possible to install the
required submodules after the fact, but simpler to follow the
instructions given above.

## Documentation

PSCF is distributed with source files for an html web manual. A copy of
the manual for the most recent numbered release is available online at

  <https://dmorse.github.io/pscfpp-man>

After cloning the source code, users can also use the
[doxygen](<https://www.doxygen.nl/index.html>) documentation generation
program to generate a local copy of the web manual.  To do this, the
doxygen application must be installed on your computer, and the directory
containing the doxygen executable must be in your unix command search
PATH. To generate documentation:

   - Change directory (cd) to the pscfpp/ root directory

   - Enter "make html"

This should create many html files in the docs/html subdirectory of the
pscfpp/ root directory.  To begin reading the resulting documentation,
point a web browser at the file pscfpp/docs/html/index.html .  This is
the main page of the web manual.

## Dependencies

PSCF has been developed on both linux and and Mac OS X operating systems,
and is designed to run on these or other unix-like systems.

The pscf_1d and pscf_pc CPU programs depend on the following external
libraries:

  - [GSL](<https://www.gnu.org/software/gsl/>) - GNU Scientific Library

  - [FFTW](<https://www.fftw.org/>) Fast Fourier Transform library

The pscf_1d one-dimensional program pscf_1d requires only GSL, while
pscf_pc program relies on both GSL and FFTW.

The GPU-accelerated pscf_pg program can only be compiled and run on a
computer with an appropriate NVIDIA GPU and an NVIDIA CUDA development
kit. Some GPU-accelerated libraries used by pscf_pg, such as the FFT
library cuFFT, are provided as part of this kit.

Procedures for installing these dependencies are different for different
operating system environments and different package managers. Instructions
for some common environments are given in the web manual.

## Compiling

The PSCF source code is written using C++ as the primary language, with
CUDA used in pscf_pg for GPU acceleration. PSCF is only provided in
source file format, and must be compiled from source.

The build system used to compile and install PSCF relies on standard unix
utilities that are available in any unix-like command-line environment,
and requires access to a Python 3 interpreter.  To compile and run this
or any other unix software on Mac OS X, the user must first install the
Mac OS X unix command line tools.

Complete directions for compiling and installing PSCF are provided in
section 2 of the [web manual](<https://dmorse.github.io/pscfpp-man>).
A brief summary is given below of instructions for steps that must be
taken after cloning the git repository and installing all of the
required dependencies:

   - Add the directory path pscfpp/bin to your unix command search PATH
     environment variable. This is where executables will be installed.

   - Add the directory path pscfpp/lib/python to your PYTHONPATH
     environment variable. This is necessary to allow a python interpreter
     to find python modules that are used by the build system.

   - Run the pscfpp/configure script by entering "./configure" from the
     pscfpp/ root directory, with or without an optional filename argument.
     This script installs default versions of several files that are
     required by the build system. The configure script usually only
     needs to be run once, before compiling for the first time.  In a
     linux environment, it is usually sufficient to run the configure
     script without a filename argument.  See Section 2.6 of the
     [web manual](<https://dmorse.github.io/pscfpp-man>) for more
     information about how to invoke the configure script for Apple OS X.

   - To compile and install only CPU-based programs in the package on a
     computer that does not have an appropriate NVIDA GPU, simply enter
     "make" from the pscfpp/ root directory immediately after
     running the configure script.

   - Some further configuration is required to compile the GPU-enabled
     pscf_pg program. Compilation of pscf_pg is only possible on a
     machine that has an
     appropriate NVIDA GPU and a CUDA development kit. To enable
     compilation of CUDA code, first enter "./setopts -c1" from the
     pscfpp/ directory.  Then use the -a option of the same setopts
     script to set the correct GPU architecture for your system. For
     example, to compile for a V100 GPU with CUDA compute capability
     7.0, one would enter "./setopts -a sm_70". Instructions for
     choosing the correct string parameter for the -a option for a
     particular GPU can be obtained by entering "./setopts -h" or by
     consulting the installation section of the web manual.  After
     these steps, enter "make all" to compile.

The procedure described above for building PSCF is similar to the standard
"configure; make; make install" procedure for GNU software, except that
it does not require a separate "make install" command. This is because
the PSCF "make" command is designed to create executables in the bin/
subdirectory of the root pscfpp/ directory in which the package is
configured, and the above instructions assume that this is their intended
final destination.  If the pscfpp/ directory is installed within a user's
home directory, the above instructions do not require the use of sudo or
adminstrator permission, but instead require the user to modify their
PATH variable to allow the unix shell to find executables installed in
psfpp/bin.

## Command line usage

PSCF is a package containing several different SCFT programs (pscf_1d,
pscf_pc, and pscf_pg) with somewhat different capabilities but very
similar command line interfaces. To perform a calculation, each of these
programs must read a parameter file and a command file. The parameter
file, which is processed first, is a fixed format file that contains
parameters required to describe the physical system of interest and
initialize the program.  The command file, which is processed after
the parameter file, is a script that contains a sequence of commands
that are read and executed in the order that they appear. Contents and
syntax for PSCF parameter and commands files are discussed in Sec. 3
of the web manual.

The command file usually contains the names of several input and output
data files as arguments to commands that read or write these files.
Specifically, a command file to perform either a SCFT or FTS calculation
normally contains a command to read an input file that contains an
initial guess for the monomer chemical potential fields.  A command
file for either an SCFT or FTS calculation also normally contains
commands output final chemical potential and monomer concentration
fields to output files, i.e., converted solutions for SCFT or the
field values after the final step of a FTS.

### pscf_1d

The usual command line syntax for invoking the pscf_1d program for
one-dimensional SCFT calculations is
```
pscf_1d -p param -c command -e
```
where the string "param" denotes actual name of a parameter file, and
"command" denotes the name of a command file.  The -p and -c options
are required.

The "-e" command line option in the above example is optional, and
causes the program to "echo" the parameter file to standard output as
this file is being processed. Use of this option makes it easier to
debug any failure arising from a syntax error in the parameter file.
The echoed output also documents allowed optional elements of the
parameter file that were not found in the input file, thus helping
to document the full file format.

### pscf_pc and pscf_pg

The command line interfaces for the pscf_pc and pscf_pg programs take
the same required and optional elements as pscf_1d, but also require a
value for the dimension of space as an additional parameter.  This is
an integer parameter of the -d option that must have a value 1, 2 or 3,
which specifies the number of coordinates along which the structure
is periodic.  The usual syntax for invoking pscf_pc for a three
dimensionally periodic structure (e.g., a network structure or a 3D
arrangement of spheres), while also using the -e option, is thus
```
pscf_pc -d 3 -p param -c command -e
```
where "param" and "command" denote names of parameter and command files.
The syntax for a simulation of one-dimensionally periodic lamellar
phase would instead use 1 as the argument of the -d option. The syntax
for invoking pscf_pg is the same as that for pccf_pc, except for the
use of the different program name.

Programs that are invoked as shown in the above examples would write log
output that is produced during execution to the user's screen (i.e., to
standard output).  This log output may instead be redirected to a file
by using the unix ">" standard output redirection operator. For example,
a command such as
```
pscf_pc -d 3 -p param -c command -e > log
```
would redirect log output that is written during simulation of a
3-dimensionally periodic structure to a file named "log" in the current
working directory. Standard output should normally be redirected to a
file when a program is run in a queue on a shared resource.

## Examples

The quickest way to learn the syntax of the PSCF parameter and command
files is by examining examples.  The directory pscfpp/examples contains
input files for examples of many different types of SCFT and FTS
calculations.  Top level subdirectories of the examples/ directory
named 1d/, pc/ and pg/ contain examples for different PSCF programs
(pscf_1d, pscf_pc, or pscf_pg).

Subdirectory examples/1d contains examples of SCFT calculations for the
1D finite-difference program pscf_1d. Top level subdirectories of the
directory examples/fd contain examples for planar, cylindrical and
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
example or a set of very closely related examples as an example directory.
Each such example directory contains at least one parameter file (usually
named "param"), at least one command file (usually named "command"), and
an input chemical potential field (w field) file. Example directories
that contain input files for two or more examples usually contain a single
file with an initial chemical potential field, but may contain two or
more different sets of parameter or command files for use in different
examples.

Example directories that contain files for a single example usually
contain a shell script named "run" that can be executed to run the
example using the supplied input files.  The simplest way to run such an
example is thus to change directory (cd) to the relevant example directory
and enter "./run" from within that directory. (Note the use of the
prefix "./", which tells the operating system to look for the run script
in the current working directory).  Users may also inspect the text of
such a run script to see the command and command line options used to
run the example.  Example directories that contain input files for two
or more closely related examples may contain several such scripts with
names that are variants of "run", each of which can be execute to run a
specific example.

Almost every example directory also contains a script named "clean"
that can be executed to remove all of the output files that are created
by running the example.

## License

This C++/CUDA version of PSCF is free, open source software. It is
distributed under the terms of the GNU General Public License (GPL)
as published by the Free Software Foundation, either version 3 of the
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

