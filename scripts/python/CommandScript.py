from Record import *

# --------------------------------------------------------------------------
# The CommandScript class can read a command file block and stores its
# contents in a form that allows particular commands to be accessed and
# modified. 
#
# Usage: The command:
#
# commands = readCommandFile("commands")
#
# parses a command file "commands" and returns a CommandScript object,
# which we have named commands.
#
# CommandScript:
#
# A CommandScript object is a container for a list of child Command objects.
# Each command corresponds to one line in the command script, in which the
# first field is the label.  Individual commands can be accessed in either 
# of two ways:
#
# i) The commands_ attribute is a list containing all of the commands,
# indexed in the order in which they appear. 
#
# ii) Each child is also stored as an attribute with a name given by
# the label of the child ParamComposite or Parameter. 
#
# If two or more lines in a command script have the same label, then
# the attribute associated with that label becomes a list, in which
# the first line with this label is indexed by 0, the second by 1, etc.
# --------------------------------------------------------------------------

class CommandScript:
 
   def __init__(self):
      self.label_    = None
      self.commands_ = []
 
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

   def __str__(self):
       list = []
       for command in self.commands_:
          list.append(str(command))
       return '\n'.join(list)

class Command(Record):

   def __init__(self, line):
      Record.__init__(self, line)
      self.label_ = self.fields[0]

   def nParam(self):
      return len(self.fields - 1)

   def param(self, i = 0):
      return self.fields[i+1]

   def setParam(self, value, i = 0):
      self.fields[i+1] = value

def readCommandFile(filename):

   file  = open(filename, 'r')
   lines = file.readlines() 
   file.close()

   p = CommandScript()
   p.read(lines)

   return p

