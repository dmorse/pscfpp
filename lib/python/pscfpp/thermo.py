# -----------------------------------------------------------------------
#   This module provides tools to parse the PSCF "thermo file" and store 
#   all values within it in a single object. A thermo file is output by
#   the PSCF command WRITE_THERMO, and can either be a standalone file or
#   can be a section of a larger data file, as is the case in the output
#   of a Sweep operation, in which the files named *.dat contain both the
#   parameter file data and the thermo data. Users can access and modify 
#   the stored values of the properties after parsing by using specific 
#   statements (commands), and can write the entire thermo object to a 
#   file in proper format. 
#
#   Parsing a thermo file:
#
#       A Thermo object represents a PSCF thermo file that is identified 
#       by a first line that contains the value of fHelmholtz. It stores
#       all the contents within the file in a form that allows each
#       property within it to be accessed and modified.
#
#       A PSCF thermo file always starts with the fHelmholtz value.  
#       Users may parse such a file by creating a Thermo object, 
#       passing the name of thermo file as an argument. This constructor 
#       parses the file and returns a Thermo object that contains its 
#       contents. 
#
#       Example: To read and parse a thermo file with name 'thermo',
#       execute the following code within a python3 interpreter:
#
#           from pscfpp.thermo import *
#           t = Thermo('thermo')
#
#       Alternatively, a Thermo object can be read from an open Python
#       file object. In this case, the object must first be created, 
#       and the open file object must then be passed to the read()
#       method. Note that the file object is not closed by read().
#
#       Example:
#
#           f = open('thermo')
#           t = Thermo()
#           t.read(f)
#           f.close()
#
#       A Thermo object can be written to a file by calling the 
#       writeOut() method, passing in the string of the file name to 
#       which the Thermo should be written.
#
#           Example: t.writeOut('thermoOut')
#
#   Accessing properties:
#
#       After creating a Thermo object, users can retrieve the values of 
#       any property by name, using a dot notation for properties stored
#       within the Thermo object. There are three different types of 
#       properties that are stored within a Thermo object. These are listed 
#       below, along with a summary of what is returned when they are 
#       accessed, and an example Python expression that would access this 
#       type of property in a typical Thermo object:
#
#           1. Single parameter: if a property is a single parameter, it 
#           contains a label followed by one value on a single line.  
#           Accessing a single parameter returns the value of the 
#           parameter, which is a float. 
#           Example: t.fHelmholtz    or     
#                    t.fInter
#
#           2. List of Species: the Thermo class has two members named
#           Polymers and Solvents, which are each a list of Species 
#           objects. Each Species object represents one polymer or 
#           solvent, and has two attributes: phi and mu. A specific
#           Species from the list is accessed by square bracket indexing
#           and returns only the Species object, which is meant to be used
#           only through accessing the phi or mu attributes. To access 
#           the phi or mu value of a Species object, a specific name of 
#           the property is needed and it returns the corresponding value, 
#           which is a float. 
#           Example: t.Polymers[0].phi    or
#                    t.Solvents[1].mu
#
#           3. LatticeParameters: accessing the attribute LatticeParameter 
#           returns a list that stores values of different lattice 
#           parameters. A specific lattice parameter can be accessed by 
#           square bracket indexing, and returns a float.
#           Example: t.LatticeParameters    or
#                    t.LatticeParameters[0]
#
#       The parser also allows users to modify some entries in different 
#       preset formats for particular types of properties with equal sign 
#       operator ('='), which are listed below: 
#
#           1. Single parameter: a single parameter can be modified 
#           by Python arithmetic operators. 
#           Example: t.fHelmholtz *= 2    or    
#                    t.fHelmholtz = 0.8   or 
#
#           2. phi or mu value of Species object: can be modified by 
#           Python arithmetic operators.
#           Example: t.Polymers[0].phi += 0.3    or
#                    t.Solvents[1].mu = 0.2
#
#           3. LatticeParameters: two ways to modify:
#           Example: 1. change the whole list by using a Python list:
#                    t.LatticeParameters = [0.2, 0.2]
#                    2. change a specific value within the list:
#                    t.LatticeParameters[0] = 0.2
#
# Module Contents:
#
#   class Thermo:
#           A Thermo object contains the data output by WRITE_THERMO. It 
#           parses a thermo file and stores the contents within in it in a
#           form that allows all properties to be accessed and modified.
#
#   class Species:
#           A Species object stores the phi and mu values of a single 
#           polymer or solvent. It also stores whether it has index from 
#           the reading.
#
# -----------------------------------------------------------------------


class Thermo:
   '''
   Purpose:
      The class represents the Thermo file
   Instance variables:
      fHelmholtz: Helmholtz free energy per monomer, unit of kT
      pressure: 
         nondimensionalized pressure, Pv/kT (v = monomer reference volume)
      fIdeal: ideal gas component of free energy
      fInter: monomer interaction component of free energy 
      fExt: external field component of free energy
      Polymers: 
         a list of Species objects that represent polymers 
      Solvents:
         a list of Species objects that represent solvents
      LatticeParameters:
         a list of floats represents all Lattice parameters
      tableLabel: 
         a string that represents the label for both Polymers and Solvents
   Methods:
      __init__(self, filename=None):
         constructor, with one argument:
            filename:
               the filename that needed to be read, default to be None
      read(self, openFile):
         method to read the open Thermo file, openFile as the argument, line
         by line and update the read items
      skipEmptyLine(self, openFile):
         method to skip the empty line read from the file, with openFile as
         argument
      writeOut(self, filename):
         method to write out the stored Thermo file to a specific txt file
         with the name of the argument filename
      writeOutStirng(self):
         return the string for writing out
   '''
   def __init__(self, filename=None):
      self.fHelmholtz = None
      self.pressure = None
      self.fIdeal = None
      self.fInter = None
      self.fExt = None
      self.Polymers = None
      self.Solvents = None
      self.LatticeParameters = None
      self.tableLabel = None

      if filename != None:
         with open(filename) as f:
            self.read(f)

   def read(self, openFile):
      line = self.skipEmptyLine(openFile)
      l = line.split()
      if l[0] != 'fHelmholtz':
         raise Exception('Not valid Thermo file')
      else:
         self.fHelmholtz = float(l[1])
      line = openFile.readline()
      l = line.split()
      self.pressure = float(l[1])
      line = self.skipEmptyLine(openFile)

      while line != '\n':
         if line == 'Free energy components:\n':
            line = openFile.readline()
            while line != '\n':
               l = line.split()
               if l[0] == 'fIdeal':
                  self.fIdeal = float(l[1])
               if l[0] == 'fInter':
                  self.fInter = float(l[1])
               if l[0] == 'fExt':
                  self.fExt = float(l[1])
               line = openFile.readline()

         if line == 'Polymers:\n':
            self.Polymers = []
            self.tableLabel = openFile.readline()
            line = openFile.readline()
            while line != '\n':
               l = line.split()
               self.Polymers.append(Species(l))
               line = openFile.readline()

         if line == 'Solvents:\n':
            self.Solvents = []
            line = openFile.readline()
            line = openFile.readline()
            while line != '\n':
               l = line.split()
               self.Solvents.append(Species(l))
               line = openFile.readline()

         if line == 'Lattice parameters:\n':
            self.LatticeParameters = []
            line = openFile.readline()
            l = line.split()
            for i in range(len(l)):
               self.LatticeParameters.append(float(l[i]))

         line = self.skipEmptyLine(openFile)

         if line == '':
            break

   def skipEmptyLine(self, openFile):
      line = openFile.readline()
      while line == '\n':
         line = openFile.readline()
      return line

   def writeOut(self, filename):
      with open(filename, 'w') as f:
         f.write(self.writeOutString())

   def writeOutString(self):
      s = ''
      s += 'fHelmholtz'
      v = f'{self.fHelmholtz:.11e}'
      s += f'{v:>22}\n'
      s += 'pressure  '
      v = f'{self.pressure:.11e}'
      s += f'{v:>22}\n'
      s += '\n'

      if (self.fIdeal != None) or (self.fInter != None) or (self.fExt != None):
         s += 'Free energy components:\n'
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

      s += 'Polymers:\n'
      s += self.tableLabel
      if self.Polymers[0].hasIndex == True:
         for i in range(len(self.Polymers)):
            p = f'{self.Polymers[i].phi:.11e}'
            m = f'{self.Polymers[i].mu:.11e}'
            s += f'{i:>5}{p:>20}{m:>20}\n'
      else:
         for i in range(len(self.Polymers)):
            p = f'{self.Polymers[i].phi:.11e}'
            m = f'{self.Polymers[i].mu:.11e}'
            s += f'     {p:>20}{m:>20}\n'
      s += '\n'

      if self.Solvents != None:
         s += 'Solvents:\n'
         s += self.tableLabel
         if self.Solvents[0].hasIndex == True:
            for i in range(len(self.Polymers)):
               p = f'{self.Solvents[i].phi:.11e}'
               m = f'{self.Solvents[i].mu:.11e}'
               s += f'{i:>5}{p:>20}{m:>20}\n'
         else:
            for i in range(len(self.Polymers)):
               p = f'{self.Solvents[i].phi:.11e}'
               m = f'{self.Solvents[i].mu:.11e}'
               s += f'     {p:>20}{m:>20}\n'
         s += '\n'

      if self.LatticeParameters != None:
         s += 'Lattice parameters:\n'
         s += '       '
         v = f'{self.LatticeParameters[0]:.11e}'
         s += f'{v:>18}'
         if len(self.LatticeParameters) > 1:
            for i in range(1, len(self.LatticeParameters)):
               v = f'{self.LatticeParameters[i]:.11e}'
               s += f'{v:>20}'
         s += '\n'

      s += '\n'

      return s

# End class Thermo ------------------------------------------------------


class Species:
   '''
   Purpose:
      The class represents the object type that stores the phi and mu 
      values for the single species.
   Instance variables:
      hasIndex: boolean to store if the single species has index
      phi: float that stores the phi value of the single species
      mu: float that stores the mu value of the single species
   Methods:
      __init__(self, l):
         constructor, with one argument:
            l: the list of the read line after splitting, required
   '''
   def __init__(self, l):
      if len(l) == 3:
         self.hasIndex = True
         self.phi = float(l[1])
         self.mu = float(l[2])
      else:
         self.hasIndex = False
         self.phi = float(l[0])
         self.mu = float(l[1])

# End class Species -----------------------------------------------------