from Record import *
import os
from os.path import *
from file import *

class RecordEditor:

   def __init__(self):
      self.hasFilter = False
      self.hasNew    = False
      self.isReady   = False
      self.isTest    = True

   def setFilter(self, filter, field = 0):
      self.filter      = filter
      self.filterField = field
      self.hasFilter   = True
      if self.hasNew:
         self.isReady = True

   def setNew(self, new, field):
      self.new  = new
      self.dataField = field
      self.hasNew = True
      if self.hasFilter:
         self.isReady = True

   def setIsTest(self, isTest):
      self.isTest = isTest

   def editFile(self, filename):
     
      if (not self.isReady):
         print "RecordEditor is not ready"
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
   
#  def editFiles(self, dirName, pattern):
#     if (not self.isReady):
#        print "RecordEditor is not ready"
#        return
#     dir = Directory(dirName)
#     filenames = dir.filenames(pattern)
#     for filename in filenames:
#        self.editFileBlocks(filename)
