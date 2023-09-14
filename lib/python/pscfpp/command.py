"""! Module for processing PSCF command scripts. """

from pscfpp.text import *

##
#  Class to parse a PSCF command script.
#
#  The Script class can read a command script file and store its
#  contents in a form that allows particular commands to be accessed and
#  modified. 
# 
#  **Construction and Parsing:**a
# 
#  The constructor can take the name of a command file or an open file 
#  object as a parameter, and will read, parse and store the contents of
#  the specified file. 
#
#  Example: To parse a command file named 'command', one could enter
#  \code
#     s = Script('command')
#  \endcode
#
#  **Accessing individual commands:**
# 
#  A Script object is a container for a list of child Command 
#  objects.  Each command corresponds to one line in the command script, 
#  in which the first field is the label.  Individual Command objects 
#  within a Script object can be accessed in either of two ways:
#
#   - Command objects may be accessed by integer index using the square 
#     bracket notation for elements of a python list. If s is a Script,
#     then s[2] is the 3rd command in the script. 
#
#   - Commands may also be accessed using square bracket notation
#     using the command name (or label) string as a key. If the script
#     contains several commands with the same name, a list of Command
#     objects is returned, listed in the order in which they appeared.
#
#  String representation of individual Commands objects can be printed 
#  using the str() or print functions.  Values of arguments of 
#  individual Command objects may then be accessed or modified using 
#  the param and setParam methods of class Command (see below).
# 
class Script:
 
   ## 
   # Constructor.
   #
   # If supplied with a parameter "in" that is an open file or a file
   # name string, this function reads and parses that command file and
   # stores the contents in the resulting object. 
   #
   # If the parameter "in" is absent or equal to None, an empty object 
   # is created and returned.
   #
   # \param filename  file name or open file object for input command file
   #
   def __init__(self, filename = None):
      self._commands = []
      self._indices = {}
      if (filename != None):
         self.read(filename)
 
   ## 
   # Read and parse a command script. 
   #
   # \param filename  file name string or open file for command file
   #
   def read(self, filename):

      # Determine type of parameter in
      if isinstance(filename, str):
         file = open(filename, 'r')
      elif isinstance(filename, file):
         file = filename
      else:
         raise Exception("Invalid parameter in")

      # Split file into list of line strings
      lines = file.readlines();

      # Process lines
      i = 0
      finish = False
      while i < len(lines) and not finish:
         line = lines[i]
         if not isinstance(line, str):
            raise Exception('Command line object is not a string')
         command = Command(line)

         # Add command to list self._commands
         self._commands.append(command)

         # Add index i to dictionary self._indices
         label = command.label
         if label in self._indices:
            index = self._indices[label]
            if isinstance(index, int):
               self._indices[label] = [index]
            elif not instance(index, list):
               raise Exception('Values of _indices must be int or list')
            self._indices[label].append(i)
         else:
            self._indices[label] = i

         if (command.label == 'FINISH'):
            finish = True
         i += 1

   ## 
   # Return a string containing the command script. 
   #
   def __str__(self):
      list = []
      for command in self._commands:
         list.append(str(command))
      return '\n'.join(list)

   ## 
   # Get a Command by integer index or command name
   #
   # \param key integer index or command name string
   #
   def __getitem__(self, key):
      if isinstance(key, int):
         return self._commands[key]
      elif isinstance(key, str):
         if key in self._indices:
            index = self._indices[key]
            if isinstance(index, int):
               return self._commands[index]
            elif isinstance(index, list):
               clist = []
               for j in index:
                  clist.append(self._commands[j])
               return clist
      else:
         raise Exception("Key must an int or a command string")

   ##
   # Get the number of commands in the script.
   #
   def __len__(self):
      return len(self._commands)

##
# A single command, with a label and zero or more parameters.
#
# The first string (or field) in a command is the label, and subsequent 
# space separated strings are interpreted as command parameters. All
# fields are stored verbatim as strings.
#
class Command(Record):

   ##
   # Constructor.
   #
   # \param line  string containing the command line
   # 
   def __init__(self, line):
      Record.__init__(self, line)
      self.label = self.fields[0]

   ## 
   # Number of parameter / arguments of the command. 
   #
   def nParam(self):
      return (len(self.fields) - 1)

   ## 
   # Return parameter number i. 
   #
   # \param i index of command parameter / argument
   #
   def param(self, i = 0):
      return self.fields[i+1]

   ##
   # Set the value of a specific command parameter.
   #
   # \param value  new value for the specified command parameter
   # \param i  index of command parameter (0 by default)
   #
   def setParam(self, value, i = 0):
      self.fields[i+1] = str(value)


