import os
import re
from os.path   import *
from string    import *

class Grep:

   def __init__(self):
      self.nFilter = 0
      self.filters = []
      self.results = []

   def addFilter(self, filter):
      self.filters.append(re.compile(filter))
      self.results.append(0)
      self.nFilter += 1

   def clearResults(self):
      for i in range(self.nFilter):
         self.results[i] = 0

   def grep(self, filename):
     
      file = open(filename, 'r')
      lines = file.readlines()
      n     = len(lines)
      file.close()
   
      for line in lines:
         for i in range(self.nFilter):
            if (self.filters[i].search(line)):
               self.results[i] += 1

