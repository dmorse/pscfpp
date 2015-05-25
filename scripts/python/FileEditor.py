import os
import re
from os.path   import *
from string    import *
from Directory import *
from File      import *

class FileEditor:

   def __init__(self):
      self.hasFilter = False
      self.hasNew    = False
      self.hasOld    = False
      self.isReady   = False
      self.isTest    = True
      self.blockSize = 1

   def setFilter(self, filter):
      self.filter    = re.compile(filter)
      self.hasFilter = True
      if (self.hasOld and self.hasNew):
         self.isReady = True

   def setOld(self, old):
      self.old    = old
      self.hasOld = True
      if (self.hasFilter and self.hasNew):
         self.isReady = True

   def setNew(self, new):
      self.new    = new
      self.hasNew = True
      if (self.hasFilter and self.hasOld):
         self.isReady = True

   def setIsTest(self, isTest):
      self.isTest = isTest

   def setBlockSize(self, blockSize):
      self.blockSize = blockSize

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
   
   def editFiles(self, dirName, pattern):
      if (not self.isReady):
         print "FileEditor is not ready"
         return
      dir = Directory(dirName)
      filenames = dir.filenames(pattern)
      for filename in filenames:
         if isfile(filename): 
            self.editFileBlocks(filename)
