"""! Module for parsing thermo files. """

##
# Parser and container for PSCF thermo file blocks.
#
#  This module provides tools to parse the PSCF "thermo" file format and
#  store the contents of a such a file block in accessible form. This file
#  reports properties that are outputs of a SCFT calculation, such as 
#  Helmholtz free energy and macroscopic pressure (or grand canonical free 
#  energy density), as well as some quantities that can be either inputs or 
#  outputs depending on the mode of operation (e.g., chemical potentials 
#  and volume fractions). A thermo file block can appear either as the 
#  entire contents of a file (which can be created by the PSCF command
#  WRITE_THERMO) or as a section of a larger data (i.e., as part of the
#  log file that is written to standard output at the end of an SCFT or
#  a section of the type of "state" file created by a Sweep operation). 
#  Users can access and modify values of properties stored in a Thermo 
#  object after it has been used to parse a thermo file block using the 
#  Python dot syntax for object attributes.  
#
#  **Construction and Parsing:**
#
#  Users may parse a thermo file block file by creating a Thermo object, 
#  passing either the name of the thermo file or an open file object as 
#  an argument. This constructor parses the file and returns a Thermo 
#  object that contains its contents.
#
#  For example: To read and parse a thermo file with name 'thermo', one
#  could enter:
#  \code
#      from pscfpp.thermo import *
#      t = Thermo('thermo')
#  \endcode
#
#  **Writing out:**
#
#  A multi-line string representation of a Thermo object in thermo file 
#  format (i.e., the format of the original input file), can be obtained 
#  using the __str__ special method, or written to a file using the write 
#  method. 
#
#  Example: If t is a Thermo object, then
#  \code
#      t.write('thermoOut')
#  \endcode
#  writes the contents of to file named 'thermoOut'.
#
#  **Accessing properties:**
# 
#  After creating a Thermo object, users can retrieve the values of 
#  any property by name, using a dot notation for properties stored
#  within the Thermo object. There are three different types of 
#  properties that are stored within a Thermo object. These are listed 
#  below, along with a summary of what is returned when they are 
#  accessed, and an example Python expression that would access this 
#  type of property in a typical Thermo object:
#
#  *Single variable:* if a property is represented by a single value,
#  such as Helmholtz free energy, it can be accessed using dot notation
#  for a Python attribute. For example, if t is a Thermo object, then
#  the expressions
#  \code
#      t.fHelmholtz       
#      t.pressure
#  \endcode
#  return the Helmholtz free energy per monomer in thermal units 
#  (property fHelmholtz) and the non-dimensionalized pressure 
#  (pressure times monomer reference volume in thermal units).
#
#
#  *Arrays of Species objects:* the Thermo class has two data attributes
#  named Polymers and Solvents, each of which is a list of instances of 
#  class Species. Each Species object represents one polymer or solvent 
#  species, and has two data attributes: phi (volume fraction) and mu 
#  (chemical potential). A specific Species from the list Polymer or 
#  Solvent is accessed by square indexing, which returns the Species 
#  object. Dot notation for attributes of a Species can then be used to 
#  access the value of phi or mu for a specific polymer or solvent.
#
#  Example: If t is a Thermo object, then
#  \code
#     t.Polymers[0].phi  
#     t.Solvents[1].mu
#  \endcode
#  return the volume fraction of polymer 0 (the first polymer) and
#  the chemical potential of solvent 1 (the second solvent). 
#
#  *LatticeParameters:* The attribute LatticeParameter is created only
#  by calculations for a periodic system, and is a list of values for
#  different lattice parameters. The interpretation of the order and
#  meaning of the lattice parameters depends on the lattice system,
#  which is given in the input parameter file. Specific parameters 
#  can be accessed by square bracket indexing, and are float values.
#
#  Example:  The expressions
#  \code
#      t.LatticeParameters   
#      t.LatticeParameters[0]
#  \endcode
#  return the list of lattice parameters and the first lattice parameter,
#  respectively.
#
#  **Modifying properties:**
#  
#  The numerical values of properties stored in a Thermo object are stored 
#  as Python attributes of Thermo or Species, which can be modified as well
#  as accessed for reading.
#
#  Example:  The expressions
#  \code
#     t.fHelmholtz *= 2       
#     t.fHelmholtz = 0.8  
#   \endcode
#   reset the value of fHelmholtz.
#
#  Example: The expressions
#  \code
#     t.Polymers[0].phi += 0.3 
#     t.Solvents[1].mu = 0.2
#  \endcode
#  modify values of phi and mu for specific species.
#
#  Example: The expressions
#  \code
#     t.LatticeParameters = [0.2, 0.2]
#     t.LatticeParameters[0] = 0.2
#  \endcode
#  modify lattice parameters, either by replacing the list in the
#  first case, or by modifying a single element in the second.
#
class Thermo:

   ##
   # Constructor.
   #
   # If the filename parameter is the name of a thermo file, the file
   # is parsed and its contents are stored in attributes of the 
   # resulting Thermo object. If filename is absent or None, then an
   # empty object is created. 
   #
   # \param filename  a filename string, defult to be None.
   #
   def __init__(self, filename=None):
      self.fHelmholtz = None
      self.pressure = None
      self.fIdeal = None
      self.fInter = None
      self.fExt = None
      self.polymers = None
      self.solvents = None
      self.cellParams = None
      self.tableLabel = None

      if filename != None:
         with open(filename) as f:
            self.read(f)

   ##
   # Read the passed-in open-file.
   #
   # This function reads a thermo block from a file that is open for 
   # reading and stores the property values in this Thermo object. 
   #
   # \param file  an open-file object.
   #
   def read(self, file):
      line = self.skipEmptyLine(file)
      l = line.split()
      if l[0] != 'fHelmholtz':
         raise Exception('Not valid Thermo file')
      else:
         self.fHelmholtz = float(l[1])
      line = file.readline()
      l = line.split()
      self.pressure = float(l[1])
      line = self.skipEmptyLine(file)

      while line != '\n':
         l = line.split()
         if l[0] == 'fIdeal':
            self.fIdeal = float(l[1])
         if l[0] == 'fInter':
            self.fInter = float(l[1])
         if l[0] == 'fExt':
            self.fExt = float(l[1])
         line = file.readline()

      line = self.skipEmptyLine(file)

      while line != '\n':
         if line == 'polymers:\n':
            self.polymers = []
            self.tableLabel = file.readline()
            line = file.readline()
            while line != '\n':
               l = line.split()
               self.polymers.append(Species(l))
               line = file.readline()

         if line == 'solvents:\n':
            self.solvents = []
            line = file.readline()
            line = file.readline()
            while line != '\n':
               l = line.split()
               self.solvents.append(Species(l))
               line = file.readline()

         if line == 'cellParams:\n':
            self.cellParams = []
            line = file.readline()
            while line != '\n':
               l = line.split()
               self.cellParams.append(float(l[1]))
               line = file.readline()
         
         line = self.skipEmptyLine(file)

         if line == '':
            break

   ##
   # Skip empty lines in the file.
   # 
   # This function skips empty lines read from the file. It is a helper
   # function used by the read() function.
   #
   # \param file  a file object, open for reading
   #
   def skipEmptyLine(self, file):
      line = file.readline()
      while line == '\n':
         line = file.readline()
      return line

   ##
   # Write the contents of this object in thermo file format to a file.
   #
   # This function writes out the thermo file string to the 
   # specified file with the name of the passed-in parameter, 
   # filename.
   #
   # \param filename  a filename string.
   #
   def write(self, filename):
      with open(filename, 'w') as f:
         f.write(self.__str__())

   ##
   # Return a string representation of this object in thermo file format.
   #
   def __str__(self):
      s = ''
      s += 'fHelmholtz'
      v = f'{self.fHelmholtz:.11e}'
      s += f'{v:>22}\n'
      s += 'pressure  '
      v = f'{self.pressure:.11e}'
      s += f'{v:>22}\n'
      s += '\n'

      if (self.fIdeal != None) or (self.fInter != None) or (self.fExt != None):
         if self.fIdeal != None:
            s += 'fIdeal    '
            v = f'{self.fIdeal:.11e}'
            s += f'{v:>22}\n'
         if self.fInter != None:
            s += 'fInter    '
            v = f'{self.fInter:.11e}'
            s += f'{v:>22}\n'
         if self.fExt != None:
            s += 'fExt      '
            v = f'{self.fExt:.11e}'
            s += f'{v:>22}\n'
         s += '\n'

      s += 'polymers:\n'
      s += self.tableLabel
      for i in range(len(self.polymers)):
         p = f'{self.polymers[i].phi:.11e}'
         m = f'{self.polymers[i].mu:.11e}'
         s += f'{i:>5}{p:>20}{m:>20}\n'
      s += '\n'

      if self.solvents != None:
         s += 'solvents:\n'
         s += self.tableLabel
         for i in range(len(self.solvents)):
            p = f'{self.solvents[i].phi:.11e}'
            m = f'{self.solvents[i].mu:.11e}'
            s += f'{i:>5}{p:>20}{m:>20}\n'
         s += '\n'

      if self.cellParams != None:
         s += 'cellParams:\n'
         for i in range(len(self.cellParams)):
            v = f'{self.cellParams[i]:.11e}'
            s += f'{i:>5}{v:>20}\n'
         s += '\n'

      s += '\n'

      return s

##
# Container for phi and mu for a single species in a Thermo object.
#
# A Species object stores the phi and mu values of a single polymer or
# solvent species.
#
class Species:
   
   ##
   # Constructor.
   #
   # \param l  a list of read values in string type.
   #
   def __init__(self, l):
      if len(l) == 3:
         self.phi = float(l[1])
         self.mu = float(l[2])
      if len(l) == 2:
         self.phi = float(l[0])
         self.mu = float(l[1])

