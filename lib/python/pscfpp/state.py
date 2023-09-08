"""! Module for parsing state files. """

import pscfpp.param as param
from pscfpp.thermo import Thermo

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
#     from pscfpp.state import *
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


