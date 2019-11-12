#!/usr/bin/env python

from Grep      import *
from Directory import *

x = Grep()
x.addFilter(r'UTIL_THROW')

d = Directory(r'.')
filenames = d.filenames('*.cpp')
for filename in filenames:
   x.grep(filename)
   if x.results[0]:
      print filename
   x.clearResults()


