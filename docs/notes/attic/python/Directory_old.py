import os
from os.path import *
from string  import *
from glob    import *
from File    import *

class Directory(File):

    def __init__(self, path = None, scan = 1):
        self.path  = path
        self.files = {}
        self.dirs  = {}
        if self.path and scan: 
            for name in os.listdir(self.path):
                if self.path == '.':
                    path = name
                else:
                    path = self.path + os.sep + name
                if isfile(path):
                    self.files[path] = File(path)
                if isdir(path):
                    self.dirs[path]  = Directory(path)

    def scan(self):
        for name in os.listdir(self.path):
            if self.path == '.' :
                path = name
            else:
                path = self.path + os.sep + name
            if isfile(path):
                self.files[path] = File(path)
            if isdir(path):
                self.dirs[path]  = Directory(path)
  
    def clear(self):
        self.files = {}
        self.dirs  = {}

    def filenames(self, pattern = '*', recursive = 1):
        #r = []
        r = glob( self.path + '/' + pattern)
        #for x in self.files.keys() :
        #    r.append(str(self.files[x]))
        if recursive:
           for x in self.dirs.keys():
               r.extend(self.dirs[x].filenames(pattern))
        return r

    def __repr__(self):
        r = []
        r.append( self.path + '\n' ) 
        for x in self.files.keys() :
            r.append( str(self.files[x]) + '\n' )
        for x in self.dirs.keys():
            r.append( repr(self.dirs[x])  )
        return join(r,'')

    def __str__(self):
        r = []
        r.append(self.path + '\n' ) 
        for x in self.files.keys() :
            r.append( str(self.files[x]) + '\n' )
        for x in self.dirs.keys() :
            r.append( str(self.dirs[x]) )
        return join(r,'')

    def xml(self,indent=''):
        r = []
        r.append( indent + '<Directory path="'+ self.path + '" >\n' )
        n_indent = indent + '   '
        for x in self.files.keys():
            r.append( self.files[x].xml(n_indent) )
        for x in self.dirs.keys() :
            r.append( self.dirs[x].xml(n_indent)  )
        r.append( indent + '</Directory>'+'\n')
        return join(r,'')

    def write(self,filename):
        file = open(filename,'w')
        file.write( self.xml() )
        file.close()

    def ls(self):
        for x in self.files.keys() :
            print x
        for x in self.dirs.keys() :
            print x + os.sep

    def __eq__(self,other) :
        # Check files
        for x in self.files.keys():
            if x in other.files.keys():
                if self.files[x] !=  other.files[x]: 
                    return 0
            else:
                return 0
        for x in other.files.keys():
            if not x in self.files.keys(): return 0
        # Check directories
        for x in self.dirs.keys():
            if x in other.dirs.keys():
                if self.dirs[x] !=  other.dirs[x]: 
                    return 0
            else:
                return 0
        for x in other.dirs.keys():
            if not x in self.dirs.keys(): return 0
        # Return true if all tests passed
        return 1

    def __ne__(self,other) :
        if self == other : 
           return 0
        return 1

    def diff(self,other) :
        r = []
        # Check files
        for x in self.files.keys():
            if x in other.files.keys():
                if self.files[x] !=  other.files[x]: 
                    r.append( '> ' + self.files[x].xml() )
                    r.append( '< ' + other.files[x].xml() + '\n')
            else:
                r.append('> ' + self.files[x].xml() + '\n' )
        for x in other.files.keys():
            if not x in self.files.keys(): 
                r.append('< ' + other.files[x].xml() + '\n' )
        # Check directories
        for x in self.dirs.keys():
            if x in other.dirs.keys():
                r.append( self.dirs[x].diff(other.dirs[x]) )
            else:
                r.append('> ' + self.dirs[x].path + '\n' )
        for x in other.dirs.keys():
            if not x in self.dirs.keys(): 
                r.append('< ' + other.dirs[x].path + '\n' )
        # Return true if all tests passed
        return join(r,'')

