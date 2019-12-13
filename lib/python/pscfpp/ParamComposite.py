from Record import *

# -----------------------------------------------------------------------------
# File Contents:
#
#   class ParamComposite
#   class Parameter
#   class Blank
#   class ParamRecord
#   def readParamFile(filename)
#
# A ParamComposite object can parse a parameter file block and stores its
# contents in a form that allows particular parameters to be accessed and
# modified. 
#
# A Parameter object holds the label and value(s) of a particular parameter.
#
#------------------------------------------------------------------------------
# Usage: 
#
# To read a parameter file, use the global function readParamFile().
# This reads the file and returns an instance of the ParamComposite class 
# that contains the contents of the file. For example, to read a parameter
# file named "param":
#
# McSimulation = readParamFile("param")
#
# This returns a ParamComposite object, which we have named McSimulation,
# that contains the contents of the file stored as an object tree.
#
# ------------------------------------------------------------------------------
# ParamComposite:
#
# An initialized ParamComposite object contains a list of children that 
# each can be instances of Parameter or other ParamComposite subobjects. 
# Neither the beginning line, containing the label and opening bracket, 
# nor the closing line that contains a closing bracket are included
# as children. The string representation of a ParamComposite, as defined
# by the __str__() method, reproduces the block of lines from which it 
# was formed.
#
# The children of a ParamComposite can be accessed in either of two ways:
#
# i) The children_ attribute is a list containing all of the children,
# indexed in the order in which they appear. The index starts from 0.
#
# ii) Each child is also stored as an attribute with a name given by
# the label of the child ParamComposite or Parameter. In a parameter
# file for an McSimulation,
#
# parameter = McSimulation.FileMaster.inputPrefix
# 
# will thus the variable parameter a Parameter object containing the 
# label and value for the inputPrefix parameter of the FileMaster 
# within a McSimulation.
#
# ------------------------------------------------------------------------------
# Parameter:
# 
#   A Parameter object stores a parameter label in an attribute label_ and
# stores the text string representations of one or more parameter values. 
# The internal format is flexible enough to store a single value, a 
# single line containing several values, or multiple lines of data.
# All values are stored as strings: The Parameter class has no way 
# to determine if a value should be interpreted as a double, int, 
# string, enumeration, or something else. 
#
# 1) If the parameter contains only one value, it can be returned by 
# the value() method, with no arguments. For example:
#
# print McSimulation.FileMaster.inputPrefix.value()
#
# will print the value of the inputPrefix string.
# 
# 2) If the parameter contains more than one value on a single line,
# the string value is returned by the method value(i), where i is
# the index of a particular value, indexed from 0.
#
# 3) If the parameter contains several lines of data, a particular
# value is obtained by the command parameter.line(i).value(j), where
# i is an index of the line (or row), and j is index for a value 
# within that line (or column). Both i and j are indexed from 0.
#
# -----------------------------------------------------------------------------

class ParamComposite:
 
   def __init__(self, parent = None):
      self.parent_ = parent
      self.label_    = None
      self.indent_   = None
      self.children_ = []
 
   def read(self, records, i):
      record = records[i]
      label  = record.fields[0]
      if label[-1] != '{':
         print label[-1]
         return
      self.label_  = label[:-1]
      self.indent_ = record.spaces[0]
      next = True
      while next and i < len(records):
         i += 1
         record = records[i]
         done   = False
         if len(record.fields) == 0:
            self.children_.append(Blank())
            done = True
         if (not done) and (record.fields[0] == '}'):
            next = False
            done = True
         if not done:
            label = record.fields[0]
            if label[-1] == '{':
               child = ParamComposite()
            else:
               child = Parameter()
            i = child.read(records, i)
            self.children_.append(child)
            if not hasattr(self, child.label_):
               statement = 'self.' + child.label_ + ' = child'
               exec statement
            else:
               attribute = "self." + child.label_
               statement = "isList = type(" + attribute + ") == type([])"
               exec statement
               if not isList:
                  statement = "old = " + attribute
                  exec statement
                  statement = attribute + " = []"
                  exec statement
                  statement = attribute + ".append(old)"
                  exec statement
               statement = "self." + child.label_ + ".append(child)"
               exec statement
            done = True
      return i 

   def write(self, filename):
       file  = open(filename, 'w')
       lines = self.__str__()
       file.write(lines);

   def __str__(self):
       list = []
       list.append(self.indent_ + self.label_ + '{\n')
       for child in self.children_:
          list.append(str(child))
       list.append(self.indent_ + '}\n')
       return ''.join(list)

class Blank:

   def __init__(self, parent = None):
      self.parent_ = parent

   def __str__(self):
      return '\n'



class ParamRecord(Record):

   def __init__(self, line):
      Record.__init__(self, line)

   def n(self):
      return len(self.fields)

   def value(self, i = 0):
      return self.fields[i]

   def setValue(self, value, i = 0):
      self.fields[i] = value


class Parameter:
   
   def __init__(self, parent = None):
      self.parent_ = parent
      self.label_    = None
      self.indent_   = None
      self.records_  = []

   def read(self, records, i):
      ''' Read the line or lines associated with a parameter. '''
      record = records[i]
      self.indent_ = record.spaces[0]
      self.label_  = record.fields[0]
      if self.label_[-1] == '{':
         print "Error: In Parameter::read, label = " + label
         return
      offset = len(record.spaces[0]) + len(record.fields[0])
      line   = str(record)[offset:]
      record = ParamRecord(line)
      next   = True
      while next:
         self.records_.append(record)
         if i == len(records) - 1:
            next = False
         else:
            record = records[i+1]
            if len(record.fields) == 0:
               next = False
            elif len(record.spaces[0]) <= len(self.indent_):
               next = False
            if next:
               line = str(record)
               record = ParamRecord(line)
               i += 1
      return i

   def line(self, i):
      ''' Get line i from a multi-line data format. '''
      return self.records_[i]

   def value(self, i = 0):
      ''' Get value i in line 0. '''
      return self.records_[0].value(i)

   def setValue(self, value, i = 0):
      ''' Set value i in line 0. '''
      return self.records_[0].setValue(value, i)

   def __str__(self):
       list = []
       list.append(self.indent_ + self.label_)
       for record in self.records_:
          list.append(str(record))
          list.append('\n')
       return ''.join(list)


def readParamFile(filename):
   file  = open(filename, 'r')
   lines = file.readlines() 
   records = []
   for line in lines:
      records.append(Record(line))
   p = ParamComposite()
   p.read(records, 0)
   return p

