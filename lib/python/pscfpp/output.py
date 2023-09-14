"""! Module of parsers for PSCF output file formats. """

import pscfpp.param as param
import os

# Contains classes Thermo, Species, State, and Sweep

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
#      from pscfpp.output import *
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



##
# Container for data in state files produced by a sweep.
# 
#  This class is a tool to parse a PSCF "state file" and store all 
#  input parameters and property valiues within it in a single object. 
#  A state file is comprised of a block containing param file block 
#  for a given system, which contains most input parameters, and a
#  subsequent block of text containing the thermodynamic data in 
#  PSCF thermo file format.  These state files are output by PSCF
#  Sweep operation with the extension ".stt". Users can access and 
#  modify the stored values of the parameters after parsing by using 
#  specific statements (commands), and can write the entire object 
#  to a file in proper format.
#
#  **Construction and Parsing:**
#
#  The constructor for class State accepts the name of a state file,
#  and parses that file and return a State object that contains the
#  contents of the file. 
#
#  Example: To read and parse a state file with name 'state':
#  \code
#     from pscfpp.output import *
#     s = State('state')
#  \endcode
#
#  **Acccessing stored variables:**
#
#  A State object contains two data attributes named param and thermo,
#  which corresponding to param and thermo sections of the state file.
#  The param attribute is an instance of class pscfpp.param.Composite 
#  that stores the contents of the parameter file block of the state 
#  file, while the thermo attribute is an instance of class 
#  pscfpp.thermo.Thermo. Users can access either attribute or the
#  contents of either object using dot notation. See the documenation
#  of the param.Composite and thermo.Thermo classes for instructions
#  on accessing values for specific variables.
#
#  In the following examples, suppose that variable s is a State object 
#  that contains the contents of a state file. 
# 
#  Example: The expressions
#  \code
#     s.param   
#     s.param.Mixture
#  \endcode
#  return pscfpp.param.Composite objects that contain the contents of
#  the entire parameter file block and of the Mixture subblock of the
#  parameter block. 
#
#  Example: The expressions
#  \code
#      s.thermo    
#      s.thermo.fHelmholtz
#  \endcode
#  return the pscfpp.thermo.Thermo object that contains the contents
#  of the thermo block, and the value of the Helmholtz free energy
#  per monomer, respectively. 
#
#  **Modifying values:**
# 
#  Users may also use dot notation to modify values of variables stored 
#  in the param and thermo attributes. See documentation of classes
#  pscfpp.param.Composite and pscfpp.thermo.Thermo for details. 
#
class State:

   ##
   # Constructor.
   #
   # \param filename  name of a state file to be parsed (string)
   #
   def __init__(self, filename):
      with open(filename) as f:
         firstline = f.readline()
         fl = firstline.split()
         if fl[0] != 'System{':
            raise Exception('Not valid State file')
         else:
            self.param = param.Composite(f, 'System')
            self.thermo = Thermo()
            self.thermo.read(f)

   ##
   # Return string representation of this object in state file format.
   #
   # This function returns multi-line string containing the contents
   # of this object in file format of a state file (i.e, as a parameter
   # block followed by a thermo block). 
   #
   def __str__(self):
      out = self.param.__str__() + '\n' + self.thermo.__str__()
      return out

   ##
   # Write the contents of this object to file in state file format.
   #
   # This function opens a file with the specified filename and writes
   # the string produced by the __str__ function to that file. 
   #
   # \param filename name of the output file (string)
   #
   def write(self, filename):
      with open(filename) as f:
         f.write(self.__str__())


##
#  Container for data in state files produced by a PSCF sweep.
#
#  A Sweep is a container for all the data contained in the PSCF state
#  files created by a sweep. The data contained in each state file is 
#  stored in an instance of class State (full name pscfpp.state.State). 
#  Individual state objects within a Sweep can be accessed using a 
#  square-bracket syntax, like elements of a python list. Methods
#  summary and summaryString return summary reports containing values 
#  of selected variables.
#
#  **Construction:**
# 
#  The constructor for class Sweep parses and stores the contents 
#  of all of the the state files (which end with '.dat') produced 
#  by a PSCF sweep. The constructor takes a single argument which
#  is the baseFileName prefix string. 
#
#  Example: To read and parse all states files produced by a Sweep with 
#  the baseFileNAme prefix string 'out/':
#  \code
#      from pscfpp.output import *
#      s = Sweep('out/')
#  \endcode
#
#  **Accessing elements:**
#      
#  A Sweep object can be treated as a list of State objects, 
#  each of which contains the contents of a corresponding state 
#  file. Each state object can be accessed by a corresponding 
#  index with square brackets, like an element of a list, such
#  s[0] or s[1]. 
#
#  All elements and properties in each such State object can be
#  accessed using the dot notation for attributes of a State object,
#  such as s[0].param.Mixture.nMonomer or s[1].thermo.fHelmholtz.
#  See the documentation for class pscfpp.state.State for details.
#
#
#  **Summary reports:** 
#
#  See documentation for methods summary and summaryString.
#
class Sweep:

   ##
   # Constructor.
   #
   # The input parameter prefix is the baseFileName prefix string.
   #
   # \param prefix  string that is prefixed to all state file names
   #
   def __init__(self, prefix):
      self.sweep = []
      num = 0
      filename = prefix + str(num) + '.stt'
      while os.path.isfile(filename) == True:
         s = State(filename)
         self.sweep.append(s)
         num += 1
         filename = prefix + str(num) + '.stt'

   ##
   # Make a summary report containing values for selected variables.
   #
   # This function constructs a data structure containing the values 
   # of a user-specified set of variables at each state in the sweep. 
   #
   # **Specifying a list of variable names:**
   #
   # The input parameter vars is a list of strings that specify the 
   # variables of interest. Each element in this list is a string that 
   # gives the the name of an attribute of a single State object, starting 
   # with either 'param' or 'thermo', using a dot notation. For example
   # "param.Mixture.nMonomer" specifies the input parameter nMonomer of 
   # the Mixture block of the parameter section of the state file, while
   # "thermo.fHelmoltz" is the free energy per monomer reference volume.
   #
   # **Return value:**
   # 
   # The quantity returned by this function is a list in which each 
   # element is itself a list of values of values of the specified 
   # variables at a single state point.  If the return value is assigned 
   # to a variable named "summary" then summary[2] is a list of values 
   # of the chosen variables at state 2. 
   #
   # The state index (e.g., 2) may optionally be added as the first 
   # element of this list of variable values, by setting the optional
   # parameter index to True. By default index = False, so no index 
   # numbers are not included when this parameter is absent. 
   #
   # **Example:**
   #
   # If s is a Sweep object representing a mixture that has been 
   # subjected to a sweep that changes composition, then the command
   # \code
   #   vals = s.summary(['param.Mixture.Polymer[0].phi','thermo.pressure'])
   # \endcode
   # constructs list named vals of values for the volume fraction for
   # the first polymer species and the non-dimensionalized pressure of 
   # the system, for each state in the sweep.  In this example, the 
   # resulting quantity vals is a list of lists something like
   # \code
   #            [[0.5, 32.4415250701], [0.55, 30.6782301376], 
   #             [0.6, 28.8344576024], [0.65, 26.8019750534],
   #             [0.7, 24.5205812724]]
   # \endcode
   # In this example, the value of vals[2] is a python list of two floats,
   # given by [0.6, 28.8344576024], in which the first element (0.6) is 
   # the value of phi (volume fraction) for polymer 0 and the second 
   # element (28.8344) is the value of the non-dimensionalized pressure, 
   # both evaluated in state 2 (the third state in the sweep). 
   #
   # **Example:**
   #
   # Adding a final parameter index = True will add the state index
   # as an additional first element of the list of values for each 
   # state. For example, the command
   # \code
   #    s.summary(['param.Interaction.chi[0][1]'], index = True)
   # \endcode
   # applied to a sweep that changes the indicated chi parameter will
   # return an array that may look something like the following:
   # \code
   #   [[0, 12.0], [1, 13.0], [2, 14.0], [3, 15.0], [4, 16.0]]
   # \endcode
   # Here the first value for each state is the integer state index, 
   # and the second is the value of the chi parameter chi[0][1] for 
   # interactions between monomers of types 0 and 1. 
   #
   # \param vars  list of variable name strings in dot notation
   # \param index  if True, add the state index as the first variable
   #
   def summary(self, vars, index = False):
      n = len(vars)

      summary = []

      for i in range(0, len(self.sweep)):
         if index == True:
            s = [i]
         else:
            s = []
         for j in range(0, n):
            a = self.sweep[i]
            string = 'a.' + vars[j]
            try:
               val = eval(string)
            except RecursionError:
               raise Exception('Wrong command or values do not exist')
            s.append(val)
         summary.append(s)

      return summary


   ##
   # Return a summary report as a formatted string suitable for printing.
   #
   # This method produced a formatted string containing the same type of
   # summary report that is generated by the summary method, giving values
   # of selected variables for every state in a sweep. The parameter vars 
   # is a list of variable name strings formatted in dot notation, exactly 
   # as for the corresponding parameter of the Sweep.summary method. The 
   # format includes labels for the variable names and includes a state
   # index for each state.
   #
   # Example: If applied to a sweep with 5 state points, the command
   # \code
   #    report = s.summaryString(['param.Interaction.chi[0][1]',
   #                              'thermo.fHelmholtz']))
   #    print(report)
   # \endcode
   # yields an output that looks somethink like this. 
   # 
   # \code       
   #        step      chi[0][1]     fHelmholtz
   #           0  1.2000000e+01  1.9256750e+00
   #           1  1.3000000e+01  2.1102042e+00
   #           2  1.4000000e+01  2.2716872e+00
   #           3  1.5000000e+01  2.4158122e+00
   #           4  1.6000000e+01  2.5464487e+00
   # \endcode
   # The string returned by summaryString can also be written to a file.
   #
   # \param vars list of variable name strings in dot notation
   #
   def summaryString(self, vars):
      n = len(vars)

      summary = []

      for i in range(0, len(self.sweep)):
         s = [i]
         for j in range(0, n):
            a = self.sweep[i]
            string = 'a.' + vars[j]
            try:
               val = eval(string)
            except RecursionError:
               raise Exception('Wrong command or values do not exist')
            s.append(val)
         summary.append(s)

      nl = [4]
      nameList = ['step']
      for i in range(0, n):
         index = vars[i].rfind('.')
         name = vars[i][index+1:]
         nameList.append(name)
         nl.append(len(name))

      valType = []
      for i in range(0, len(summary)):
         valType.append([])
         for j in range(0, len(summary[0])):
            valType[i].append(type(summary[i][j]))

      for i in range(0, len(summary)):
         for j in range(0, len(summary[0])):
            length = len(str(summary[i][j]))
            if (valType[i][j] == str) and (length > nl[j]):
               nl[j] = length
            if (valType[i][j] == float) and (13 > nl[j]):
               nl[j] = 13
            if (valType[i][j] == int) and (length > nl[j]):
               nl[j] = length


      summaryString = ' '
      for i in range(0, len(nameList)):
         stringFormat = '{:>' + str(nl[i]) + 's}'
         if i != len(nameList)-1:
            summaryString += stringFormat.format(nameList[i]) + '  '
         else:
            summaryString += stringFormat.format(nameList[i]) + '\n '

      for i in range(0, len(summary)):
         for j in range(0, len(summary[0])):
            if valType[i][j] == int:
               stringFormat = '{:>' + str(nl[j]) + '}'
               summaryString += stringFormat.format(summary[i][j])
            if valType[i][j] == float:
               stringFormat = '{:.7e}'
               val = stringFormat.format(summary[i][j])
               summaryString += val
            if valType[i][j] == str:
               stringFormat = '{:>' + str(nl[j]) + 's}'
               summaryString += stringFormat.format(summary[i][j])
            if j == len(summary[0])-1:
               summaryString += '\n '
            else:
               summaryString += '  '

      return summaryString

   ##
   # Get a specific State object, specified by integer index.
   #
   # \param key  integer state index.
   #
   def __getitem__(self, key):
      return self.sweep[key]

   ##
   # Get the number of states in a sweep.
   #
   def __len__(self):
      return len(self.sweep)

