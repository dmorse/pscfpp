# -----------------------------------------------------------------------
#   This class is a tool to parse a PSCF "state file" and store 
#   all values within it in a single object. A state file contains both
#   a block of text representing the param file for a given system
#   and a subsequent block of text containing the thermodynamic data
#   for that state (from the SCFT solution). These state files are output
#   by a PSCF Sweep with the extension ".dat", and can be output manually 
#   by PSCF via the commands WRITE_PARAM and WRITE_THERMO. Users can 
#   access and modify the stored values of the parameters after parsing 
#   by using specific statements (commands), and can write the entire 
#   object to a file in proper format. 
#
#   Parsing a state file:
#
#       A State object represents a PSCF state file that always contains 
#       two parts, a param part and a thermo part. It stores all the
#       contents within the file in a form that allows each property
#       within it from each part to be accessed and modified.
#
#       A PSCF state file always contains two parts, param and thermo.
#       The param part always contains a main parameter Composite named
#       'System' and the thermo part always starts with the fHelmholtz
#       value. User may parse such a file by creating a State object,
#       passing in the name of state file as an argument. The constructor
#       parses the file and returns a State object that contains its
#       contents.
#
#       Example: To read and parse a state file with name 'state',
#       execute the following code within a python3 interpreter: 
#
#           from pscfpp.state import *
#           s = Composite('state')
#
#       A State object can be written to a file by calling the 
#       writeOut() method, passing in the string of the file name to 
#       which the State should be written.
#
#           Example: s.writeOut('stateOut')
#
#   Accessing elements:
#
#       A State object contains two members, param and thermo,
#       corresponding to the different sections of the state file. Users 
#       can retrieve either member by dot notation:
#       
#           1. param: the param member can be accessed by calling the
#              'param' attribute, which returns a Composite object. All 
#              stored contents within it can be accessed by the same
#              formats listed in the param module.
#              Example: s.param    or
#                       s.param.Mixture
#
#           2. thermo: the thermo member can be accessed by calling the
#              'thermo' attribute, which returns a Thermo object. All
#              stored contents within it can be accessed by the same
#              formats listed in the thermo module.
#              Example: s.thermo    or
#                       s.thermo.fHelmholtz
#
#       The parser also allows users to modify the entries in different 
#       preset formats for particular types of objects with equal sign 
#       operator ('='). Refer to both the param and thermo modules
#       for details.
#
# Module Contents:
#
#   class State:
#           A State object contains all the data from a state file. It
#           has two members, param and thermo, where the param member
#           contains all the data in the 'System' block of the state file,
#           enclosed within curly brackets {}, and the thermo member 
#           contains all the remaining data, starting from the fHelmholtz 
#           value until the end of the state file. This class parses a 
#           state file and stores all contents within it in a form that 
#           allows all properties to be accessed and modified.
#
# -----------------------------------------------------------------------

import pscfpp.param
import pscfpp.thermo

class State:
   '''
   Purpose:
      The class represents the state file.
   Instance variables:
      param: the param part of the state file
      thermo: the thermo part of the state file
   Methods:
      __init__(self, filename):
         the constructor of the State object for initiation, with one 
         argument: 
            filename, the filename that needed to be read, required
      writeOut(self, filename):
         method to write out the State object to a specific txt file 
         with name of the argument filename
      writeOutString(self):
         return the string for writing out
   '''
   def __init__(self, filename):
      with open(filename) as f:
         firstline = f.readline()
         fl = firstline.split()
         if fl[0] != 'System{':
            raise Exception('Not valid State file')
         else:
            self.param = pscfpp.param.Composite(f, 'System')
            self.thermo = pscfpp.thermo.Thermo()
            self.thermo.read(f)

   def writeOut(self, filename):
      with open(filename) as f:
         f.write(self.writeOutString())

   def writeOutString(self):
      out = self.param.writeOutString() + '\n' + self.thermo.writeOutString()
      return out

# End class State -------------------------------------------------------
