'''! Utilities for manipulating and searching text strings. '''

# import os
import re
from os.path import *
from string import *
from pscfpp.file import *

##
# Class to wrap line breaks.
#
class Wrapper:

   ## 
   # Constructor.
   #
   # Data attributes:
   #
   #   - eol :     end-of-line character
   #   - nIndent : number of spaces preceding each line
   #   - column :  current column index (cursor)
   #   - limit :   maximum number of columns (78 by default)
   #   - text :    string containing resuling wrapped text 
   #
   # \param eol end-of-line character
   # \param nIndent number of space preceding each line
   #
   def __init__(self, eol = '\n', nIndent = 0):
      self.eol = eol
      self.nIndent = nIndent
      self.column = 0
      self.text = ''
      self.limit = 78 

   ## 
   # Clear mutable attributes (column and text).
   #
   def clear(self):
      self.column = 0
      self.text = ''

   ## 
   # Add contents of string to wrapped text.
   # 
   # \param string  input string
   #
   def append(self, string):
      size = len(string)
      if (self.column + 1 + size > self.limit):
         self.text += self.eol
         self.column = 0
         for i in range(self.nIndent):
            self.text += ' '
         self.column = self.nIndent
      elif (self.column > self.nIndent):
         self.text += ' '
         self.column += 1
      self.text += string
      self.column += size

   ## 
   # Return string containing wrapped text (self.text).
   #
   def __str__(self):
      return self.text

   ## 
   # Return string containing wrapped text (self.text).
   #
   def __repr__(self):
      return self.text

##
# A Record represents a string of fields separated by whitespace.
#
# The constructor of this class takes a line of text containing fields
# separated by white spaces and divides it into fields and white space
# strings.
# 
# Data attributes:
#   line   - original string (usually a complete line)
#   fields - list in which fields[i] is field string number i
#   spaces - list in which spaces[i] is white space preceding field i
#
# The element space[i] is the string of blank characters that precedes 
# field[i], for i >= 0. The space[0] may contain zero characters, if the 
# line has no leading blank spaces, but all other elements must contain 
# one or more blank characters. Trailing white space is disregarded.
#
class Record:

   ##
   # Constructor.
   #
   # \param line string containing fields separated by white space
   #  
   def __init__(self, line):
      if line[-1] == '\n':
         line = line[:-1]
      self.line   = line
      self.spaces  = []
      self.fields  = []
      n = 0
      blank  = True
      gap = ''
      for i in range(len(line)):
         if blank:
            if line[i] != ' ':
               begin = i
               blank = False
            else:
               gap = gap + ' '
         else:
            if line[i] == ' ':
               end = i
               blank = True
               self.fields.append(line[begin:end])
               self.spaces.append(gap)
               gap = ' '
               n += 1
      if not blank:
         end = len(line)
         self.fields.append(line[begin:])
         self.spaces.append(gap)
         n +=1
      self.size = n

   ##
   # String representation - line from which Record was constructed.
   #
   def __str__(self):
      line = ''
      for i in range(self.size):
         line += self.spaces[i]
         line += self.fields[i]
      return line

##
# Class to modify selected Records in a file of records.
#
# The editFile function searches for records in which one field
# matches a regular expression, and resets the value of another.
#
class RecordEditor:

   ##
   # Constructor.
   #
   def __init__(self):
      self.hasFilter = False
      self.hasNew    = False
      self.isReady   = False
      self.isTest    = True

   ##
   # Set filter string.
   #
   # \param filter regular expression string
   # \param field  index of field to be searched for (filter)
   #
   def setFilter(self, filter, field = 0):
      self.filter      = filter
      self.filterField = field
      self.hasFilter   = True
      if self.hasNew:
         self.isReady = True

   ##
   # Set new string.
   #
   # \param new    new value for field
   # \param field  index of field to be modified
   #
   def setNew(self, new, field):
      self.new  = new
      self.dataField = field
      self.hasNew = True
      if self.hasFilter:
         self.isReady = True

   ##
   # Set isTest flag, to decide whether to perform a dry run.
   #
   # \param isTest  if true, perform a dry run (no modificatin)
   #
   def setIsTest(self, isTest):
      self.isTest = isTest

   ##
   # Open and edit the specified file of records.
   #
   # \param filename  name of file
   #
   def editFile(self, filename):
      if (not self.isReady):
         print('RecordEditor is not ready')
         return
 
      oldfile = open(filename, 'r')
      lines   = oldfile.readlines()
      n       = len(lines)
      oldfile.close()
   
      if (not self.isTest): 
         newfile = open (filename, 'w')
      found = False
      for line in lines:
         record = Record(line)
         if (self.filter == record.fields[self.filterField]):
            record.fields[self.dataField] = self.new
         if (not self.isTest):
            newfile.write(str(record) + '\n')
      if (not self.isTest):
         newfile.close()


##
#  Read a line of the form "label = string", return string.
#
#  This function reads a line of the form "label = <string>" containing 
#  a known label, an equal sign, and a string following the label. The
#  must match a known label that is passed as a parameter. If the label
#  does not match (failure), this function prints an error and returns. 
#  If the label does match (success), this function returns the 
#  remaining string verbatim, stripped of any preceding and trailing 
#  white space.
# 
#  \param file open for reading, with cursor at beginning of a line
#  \param label expected label (must match)
#  \return string after label
#   
def readLabelledLine(file, label):
   line = file.readline()
   groups = line.strip().split('=')
   if groups[0].strip() != label:
      print('Error: Expected label = ' + label)
      print('     : Found          = ' + groups[0].strip())
      raise
   if len(groups) > 2:
      print('Error: More than one = sign in line')
      raise
   if len(groups) == 1:
      return ''
   else:
     return groups[1].strip()

##
# Read line of form "label = string", in which string may contains spaces.
#
# This function reads a line of the form "label = <string>" containing a
# known label, an equal sign, and a string following the label, in which
# the string may contain sub-strings spaces. The label must match a known 
# label that is passed as a function parameter. If the label function 
# does not match, this function prints an error and exits. If the label 
# matches (success), this function returns a list of sub-strings 
# separated by spaces. For example, if the string after the equal sign
# is 'a b c', the function returns a list ['a', 'b', 'c']. 
# 
#  \param file open for reading, with cursor at beginning of a line
#  \param label expected label (must match)
#  \return list of space-delimited sub-strings that follow label
#   
def readLabelledList(file, label):
   line = readLabelledLine(file, label)
   list = []
   if line != '':
      groups = line.split(' ')
      for i in range(len(groups)):
         if groups[i] != '':
            list.append(groups[i])
   return list


##
# Class to search for text in a file.
#
# Attributes:
#
#   - filters : list of regular expression description strings
#   - results : lines that match any of the filter strings
#   - nFilter : number of filter strings
#
class Grep:

   ##
   # Constructor.
   #
   def __init__(self):
      self.nFilter = 0
      self.filters = []
      self.results = []

   ##
   # Add a regular expression string to the list.
   #
   # \param filter regular expression string
   #
   def addFilter(self, filter):
      self.filters.append(re.compile(filter))
      self.results.append(0)
      self.nFilter += 1

   ##
   # Clear results list (but not filters)
   #
   def clearResults(self):
      for i in range(self.nFilter):
         self.results[i] = 0

   ##
   # Search for lines in the file that match any filter.
   #
   # \param filename name of file to search 
   #
   def grep(self, filename):
      file = open(filename, 'r')
      lines = file.readlines()
      n     = len(lines)
      file.close()
   
      for line in lines:
         for i in range(self.nFilter):
            if (self.filters[i].search(line)):
               self.results[i] += 1


##
# Class to substitute text in one or more files.
#
# String attributes:
#
#   - filter : a regex used to filter lines for possible editing
#   - old : a regex pattern that should be replaced in lines 
#   - new : the string that should replace the old string
#
# The editFile method checks each line in a specified file to see if it 
# contains a string that matches the regular expression defined by the 
# filter attribute. If a line matches the filter expression, it search 
# for substrings that match the regular expression defined by the "old"
# attribute, and replaces each occurence by the string given by the
# "new" attribute.
# 
# The editFiles method applies the editFile method to every file in a
# specified directory tree with a name that matches a specified filename
# pattern.
#
class FileEditor:

   ##
   # Constructor.
   #
   def __init__(self):
      self.hasFilter = False
      self.hasNew    = False
      self.hasOld    = False
      self.isReady   = False
      self.isTest    = True
      self.blockSize = 1

   ##
   # Set the filter string.
   #
   # Set the regular expression used to identify lines for
   # possible modification - lines that match the filter are
   # checked for a sub-string that matches the "old" string.
   #
   # \param filter  the new "filter" string
   #
   def setFilter(self, filter):
      self.filter    = re.compile(filter)
      self.hasFilter = True
      if (self.hasOld and self.hasNew):
         self.isReady = True

   ##
   # Set the old string (the string to be replaced)
   #
   # The "old" string is replaced by the "new" string in lines that 
   # match the "filter" string.
   #
   # \param old  the "old" string, to be replaced
   #
   def setOld(self, old):
      self.old    = old
      self.hasOld = True
      if (self.hasFilter and self.hasNew):
         self.isReady = True

   ##
   # Set the new string (the replacement string)
   #
   # \param new  the "new" string that replaces the old string
   #
   def setNew(self, new):
      self.new    = new
      self.hasNew = True
      if (self.hasFilter and self.hasOld):
         self.isReady = True

   ##
   # Set the isTest boolean flag.
   #
   # If isTest is true, the edit functions only perform a dry run 
   # in which they report what changes would be made if isTest were
   # false. 
   #
   # \param isTest  perform a dry run (no changes) if true (boolean)
   # 
   def setIsTest(self, isTest):
      self.isTest = isTest

   ##
   # Set the blockSize attribute.
   # 
   # \param blockSize  number of lines in block for editFileBlock
   #
   def setBlockSize(self, blockSize):
      self.blockSize = blockSize

   ##
   # Edit a single file - edits confined to individual lines.
   #
   # \param filename  name of file to be edited
   #
   def editFile(self, filename):
     
      if (not self.isReady):
         print('FileEditor is not ready')
         return
 
      oldfile = open(filename, 'r')
      lines   = oldfile.readlines()
      n       = len(lines)
      oldfile.close()
   
      if (not self.isTest): 
         newfile = open (filename, 'w')
      found = False
      for line in lines:
         if (self.filter.search(line)):
            if (not found):
               print(filename)
               found = True
            print(line)
            line = re.sub(self.old, self.new, line);
            print(line)
         if (not self.isTest):
            newfile.write(line)
      if (not self.isTest): 
         newfile.close()
  
   ##
   # Edit a single file - edits may span several lines.
   # 
   # \param filename  name of file to be edited
   #
   def editFileBlocks(self, filename):
     
      if (not self.isReady):
         print('FileEditor is not ready')
         return

      oldfile = open(filename, 'r')
      lines   = oldfile.readlines()
      n       = len(lines)
      oldfile.close()
   
      if (not self.isTest): 
         newfile = open (filename, 'w')

      found = False
      for i in range(n + 1 - self.blockSize):
         block = join(lines[i:i + self.blockSize], '')
         if (self.filter.search(block)):
            if (not found):
               print(filename)
               found = True
            print(block)
            block = re.sub(self.old, self.new, block)
            print(block)
            if (block[-1] == '\n'):
               block = block[0:-1]
            list  = split(block, '\n')
            for j in range(len(list)):
               list[j] = list[j] + '\n'
            lines[i:i+self.blockSize] = list
      if (not self.isTest):
         for line in lines:
            newfile.write(line)
         newfile.close()
  
   ##
   # Edit all files in specified directory that matches a pattern.
   #
   # \param dirName  name of directory (string)
   # \param pattern  filename pattern (string)
   # 
   def editFiles(self, dirName, pattern):
      if (not self.isReady):
         print('FileEditor is not ready')
         return
      dir = Directory(dirName)
      filenames = dir.filenames(pattern)
      for filename in filenames:
         if isfile(filename): 
            self.editFileBlocks(filename)

