import os
from os.path import *
from string  import *

class File:

    def __init__(self, path=None, scan=1):
        self.path = path
        self.scan = scan
        if self.path and scan:
            self.mtime = getmtime(path)
            self.size  = getsize(path)

    def __str__(self):
        return self.path

    def __repr__(self):
        if scan:
           return '%-40s  size= %10d  mtime= %10d' \
              % (self.path, self.size, self.mtime)
        else:
           return self.path

    def xml(self,indent=''):
        if scan:
           return indent + '<File path ="' + self.path   + '"' \
                         + ' size ="' + str(self.size)  + '"' \
                         + ' mtime="' + str(self.mtime) + '" />\n'
        else:
           return indent + '<File path ="' + self.path   + '" />\n'

    def open(self, mode):
        return open(self.path, mode)

    def write(self, filename):
        file = open(filename, 'w')
        file.write(self.xml())
        file.close()

    def __eq__(self,other):
        if self.path  != other.path  :  return 0
        if self.size  != other.size  :  return 0
        # if self.mtime != other.mtime :  return 0
        return 1

    def __ne__(self,other):
        if self == other : 
           return 0
        return 1

