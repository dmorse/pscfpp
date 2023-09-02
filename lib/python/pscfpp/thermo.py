"""! Module for parsing thermo files. """

##
# Container for data in thermo files.
#
#  This module provides tools to parse the PSCF "thermo file" and store 
#  all values within it in a single object. A thermo file is output by
#  the PSCF command WRITE_THERMO, and can either be a standalone file or
#  can be a section of a larger data file, as is the case in the output
#  of a Sweep operation, in which the files named *.dat contain both the
#  parameter file data and the thermo data. Users can access and modify 
#  the stored values of the properties after parsing by using specific 
#  statements (commands), and can write the entire thermo object to a 
#  file in proper format.
#
#  A Thermo object represents a PSCF thermo file that is identified 
#  by a first line that contains the value of fHelmholtz. It stores
#  all the contents within the file in a form that allows each
#  property within it to be accessed and modified.
#
#  **Constrctions:**
#
#      A PSCF thermo file always starts with the fHelmholtz value.  
#      Users may parse such a file by creating a Thermo object, 
#      passing the name of thermo file as an argument. This constructor 
#      parses the file and returns a Thermo object that contains its 
#      contents.
#
#      Example:
#
#      To read and parse a thermo file with name 'thermo':
#        \code
#         from pscfpp.thermo import *
#         t = Thermo('thermo')
#        \endcode
#
#  **Writing out:**
#
#      A Thermo object can be written to a file by calling the 
#      write() method, passing in the string of the file name to 
#      which the Thermo should be written.
#
#      Example:
#        \code
#           t.write('thermoOut')
#        \endcode
#
#  **Accessing properties:**
# 
#      After creating a Thermo object, users can retrieve the values of 
#      any property by name, using a dot notation for properties stored
#      within the Thermo object. There are three different types of 
#      properties that are stored within a Thermo object. These are listed 
#      below, along with a summary of what is returned when they are 
#      accessed, and an example Python expression that would access this 
#      type of property in a typical Thermo object:
#
#           1. Single parameter: if a property is a single parameter, it 
#           contains a label followed by one value on a single line.  
#           Accessing a single parameter returns the value of the 
#           parameter, which is a float. 
#
#           Example: 
#             \code
#                t.fHelmholtz       
#                t.fInter
#             \endcode
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
#
#           Example: 
#             \code
#                t.Polymers[0].phi  
#                t.Solvents[1].mu
#             \endcode
#
#           3. LatticeParameters: accessing the attribute LatticeParameter 
#           returns a list that stores values of different lattice 
#           parameters. A specific lattice parameter can be accessed by 
#           square bracket indexing, and returns a float.
#
#           Example: 
#             \code
#                t.LatticeParameters   
#                t.LatticeParameters[0]
#             \endcode
#
#  **Modifying properties:**
#
#      The parser also allows users to modify some entries in different 
#      preset formats for particular types of properties with equal sign 
#      operator ('='), which are listed below: 
#
#           1. Single parameter: a single parameter can be modified 
#           by Python arithmetic operators. 
#
#           Example: 
#             \code
#                t.fHelmholtz *= 2       
#                t.fHelmholtz = 0.8  
#             \endcode
#
#           2. phi or mu value of Species object: can be modified by 
#           Python arithmetic operators.
#
#           Example: 
#             \code
#                t.Polymers[0].phi += 0.3 
#                t.Solvents[1].mu = 0.2
#             \endcode
#
#           3. LatticeParameters: two ways to modify:
#
#           Example: 
#
#             1. change the whole list by using a Python list:
#                \code
#                   t.LatticeParameters = [0.2, 0.2]
#                \endcode
#             2. change a specific value within the list:
#                \code
#                   t.LatticeParameters[0] = 0.2
#                \endcode
#
class Thermo:

   ##
   # Constructor.
   #
   # If necessary, a defult Thermo object can be created by 
   # without passing in the filename parameter.
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
   # This function reads the pass in open-file object line
   # by line and update the read items into the instance 
   # variables of the Thermo object. The reading stops when
   # all lines in the file are read. 
   #
   # \param openFile  an open-file object.
   #
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
         l = line.split()
         if l[0] == 'fIdeal':
            self.fIdeal = float(l[1])
         if l[0] == 'fInter':
            self.fInter = float(l[1])
         if l[0] == 'fExt':
            self.fExt = float(l[1])
         line = openFile.readline()

      line = self.skipEmptyLine(openFile)

      while line != '\n':
         if line == 'polymers:\n':
            self.polymers = []
            self.tableLabel = openFile.readline()
            line = openFile.readline()
            while line != '\n':
               l = line.split()
               self.polymers.append(Species(l))
               line = openFile.readline()

         if line == 'solvents:\n':
            self.solvents = []
            line = openFile.readline()
            line = openFile.readline()
            while line != '\n':
               l = line.split()
               self.solvents.append(Species(l))
               line = openFile.readline()

         if line == 'cellParams:\n':
            self.cellParams = []
            line = openFile.readline()
            while line != '\n':
               l = line.split()
               self.cellParams.append(float(l[1]))
               line = openFile.readline()
         
         line = self.skipEmptyLine(openFile)

         if line == '':
            break

   ##
   # Skip empty lines in the file.
   # 
   # This function skips empty lines read from the file,
   # which is a helper function of the read() function.
   #
   # \param openFile  an open-file object.
   #
   def skipEmptyLine(self, openFile):
      line = openFile.readline()
      while line == '\n':
         line = openFile.readline()
      return line

   ##
   # Write out a un-intended thermo file string to a file.
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
   # Return the intended thermo file string.
   #
   # This function returns the intended thermo file string.
   #
   # Return value:
   #
   # The intended thermo file string.
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
# Container for data of a single species.
#
#  A Species object stores the phi and mu values of a single 
#  polymer or solvent.
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

