"""! Module for processing command scripts. """

from Record import *

##
# Class to parse a PSCF command script.
#
# The Script class can read a command script file and store its
# contents in a form that allows particular commands to be accessed and
# modified. 
# 
# Usage: The command:
#
# commands = readCommandFile("commands")
# 
# parses a command file "commands" and returns a Script object,
# which we have named commands.
# 
# Script:
# 
# A Script object is a container for a list of child Command 
# objects.  Each command corresponds to one line in the command script, 
# in which the first field is the label.  Individual commands can be 
# accessed in either of two ways:
#
#   - Commands may be accessed by integer index using the square bracket
#     notation for elements of a python list. If s is a Script,
#     then s[2] is the 3rd command in the script.
#
#   - Each child may be accessed as an attribute with a name given by
#     the label of the child ParamComposite or Parameter. 
# 
# If two or more lines in a command script have the same label (i.e,.
# the same command name), the attribute associated with that label is
# a list in which the first Command with this label is indexed by 0, 
# the second by 1, etc.
# 
class Script:
 
   ## 
   # Constructor.
   #
   def __init__(self):
      self.label_    = None
      self.commands_ = []
 
   ## 
   # Read and parse a command script. 
   #
   # @param lines  list of strings containing lines of a command script
   #
   def read(self, lines):
      i = 0
      finish = False
      while i < len(lines) and not finish:
         line = lines[i]
         command = Command(line)
         self.commands_.append(command)
         if not hasattr(self, command.label_):
            statement = 'self.' + command.label_ + ' = command'
            exec statement
         else:
            attribute = "self." + command.label_
            statement = "isList = type(" + attribute + ") == type([])"
            exec statement
            if not isList:
               statement = "old = " + attribute
               exec statement
               statement = attribute + " = []"
               exec statement
               statement = attribute + ".append(old)"
               exec statement
               statement = "self." + command.label_ + ".append(command)"
               exec statement
            else:
               statement = "self." + command.label_ + ".append(command)"
         if (command.label_ == 'FINISH'):
            finish = True
         i += 1

   ## 
   # Return a string containing the command script. 
   #
   def __str__(self):
       list = []
       for command in self.commands_:
          list.append(str(command))
       return '\n'.join(list)

   ## 
   # Get a Command using the dot syntax for an attribute.
   #
   # @param key  label used as identifier for a Command.
   #
   def __getitem__(self, key):
      return self.commands_[key]

   ##
   # Get the number of commands in the script.
   #
   def __len__(self):
      return len(self.commands_)

##
# A single command, with a label and zero or more parameters.
#
class Command(Record):

   ##
   # Constructor.
   #
   # @param line  string containing the command line
   # 
   def __init__(self, line):
      Record.__init__(self, line)
      self.label_ = self.fields[0]

   ## 
   # Number of parameter / arguments of the command. 
   #
   def nParam(self):
      return len(self.fields - 1)

   ## 
   # Return parameter number i. 
   #
   # @param i index of command parameter / argument
   #
   def param(self, i = 0):
      return self.fields[i+1]

   ##
   # Set the value of a specific command parameter.
   #
   # @param value  new value for the specified command parameter
   # @param i  index of command parameter 
   #
   def setParam(self, value, i = 0):
      self.fields[i+1] = value


##
# Open and read a command script file, return a Script.
#
# \param filename name of file containing a command script
#
def readCommandFile(filename):
   file  = open(filename, 'r')
   lines = file.readlines() 
   file.close()

   p = Script()
   p.read(lines)

   return p

