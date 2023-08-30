import os
import re
from os.path import *
from string import *
from file import *

##
# Class substitute text in a files or or multiple files.
#
# String attributes:
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
   # \param filter the new filter string
   #
   def setFilter(self, filter):
      self.filter    = re.compile(filter)
      self.hasFilter = True
      if (self.hasOld and self.hasNew):
         self.isReady = True

   ##
   # Set the old string (the string to be replaced)
   #
   # \param filter the new filter string
   #
   def setOld(self, old):
      self.old    = old
      self.hasOld = True
      if (self.hasFilter and self.hasNew):
         self.isReady = True

   ##
   # Set the new string (the replacement string)
   #
   # \param filter the new filter string
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
   def setIsTest(self, isTest):
      self.isTest = isTest

   ##
   # Set the blockSize attribute.
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
         print "FileEditor is not ready"
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
               print filename
               found = True
            print line
            line = re.sub(self.old, self.new, line);
            print line
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
         print "FileEditor is not ready"
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
               print filename
               found = True
            print block
            block = re.sub(self.old, self.new, block)
            print block
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
   # Edit all files in specified directory that match a filename pattern.
   #
   # \param dirName name of directory (string)
   # \param pattern filename pattern (string)
   # 
   def editFiles(self, dirName, pattern):
      if (not self.isReady):
         print "FileEditor is not ready"
         return
      dir = Directory(dirName)
      filenames = dir.filenames(pattern)
      for filename in filenames:
         if isfile(filename): 
            self.editFileBlocks(filename)

