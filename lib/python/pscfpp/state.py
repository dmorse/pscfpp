"""! Module for parsing state files. """

import pscfpp.param as param
from pscfpp.thermo import Thermo

##
# Container for data in state files.
# 
#  This class is a tool to parse a PSCF "state file" and store 
#  all values within it in a single object. A state file contains both
#  a block of text representing the param file for a given system
#  and a subsequent block of text containing the thermodynamic data
#  for that state (from the SCFT solution). These state files are output
#  by a PSCF Sweep with the extension ".dat", and can be output manually 
#  by PSCF via the commands WRITE_PARAM and WRITE_THERMO. Users can 
#  access and modify the stored values of the parameters after parsing 
#  by using specific statements (commands), and can write the entire 
#  object to a file in proper format.
#
#  A State object represents a PSCF state file that always contains 
#  two parts, a param part and a thermo part. It stores all the
#  contents within the file in a form that allows each property
#  within it from each part to be accessed and modified.
#
#  **Constrction:**
#
#      A PSCF state file always contains two parts, param and thermo.
#      The param part always contains a main parameter Composite named
#      'System' and the thermo part always starts with the fHelmholtz
#      value. User may parse such a file by creating a State object,
#      passing in the name of state file as an argument. The constructor
#      parses the file and returns a State object that contains its
#      contents.
#
#      Example:
#
#      To read and parse a state file with name 'state':
#        \code
#          from pscfpp.state import *
#          s = State('state')
#        \endcode
#
#  **Acccessing elements:**
#
#      A State object contains two members, param and thermo,
#      corresponding to the different sections of the state file. Users 
#      can retrieve either member by dot notation:
#
#           1. param: the param member can be accessed by calling the
#              'param' attribute, which returns a Composite object. All 
#              stored contents within it can be accessed by the same
#              formats listed in the param module.
#
#              Example: 
#                \code
#                   s.param   
#                   s.param.Mixture
#                \endcode
#
#           2. thermo: the thermo member can be accessed by calling the
#              'thermo' attribute, which returns a Thermo object. All
#              stored contents within it can be accessed by the same
#              formats listed in the thermo module.
#
#              Example: 
#                \code
#                   s.thermo    
#                   s.thermo.fHelmholtz
#                \endcode
#
#  **Modifying elements:**
# 
#      The parser also allows users to modify the entries in different 
#      preset formats for particular types of objects with equal sign 
#      operator ('='). Refer to both the param and thermo modules
#      for details.
#
class State:

   ##
   # Constrctor.
   #
   # \param filename a filename string.
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
   # Return the un-intended string of the State object.
   #
   # This function returns the un-intended string 
   # representation in the state file format of the
   # state file.
   #
   # Return value:
   # 
   # The un-intended string representation in the 
   # state file format.
   #
   def __str__(self):
      out = self.param.__str__() + '\n' + self.thermo.__str__()
      return out

   ##
   # Write out an un-intened state file string to a file.
   #
   # This function writes out un-intended state file string
   # to a file with the name of the passed-in parameter
   # filename.
   # 
   # \param filename  a filename string.
   def write(self, filename):
      with open(filename) as f:
         f.write(self.__str__())


