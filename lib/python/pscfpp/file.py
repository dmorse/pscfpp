'''!   Utilities for manipulating files and paths. '''

import os
from os.path import *
from string import *
from glob import glob

##
#   Generates relative path of path2 relative to path1. 
#
#   The paths path1 and path2 must be either absolute paths or relative
#   paths defined relative to the same directory.
#
#  \param path2 path of interest
#  \param path1 reference path
#  \param Path for file2 relative to directory containing path1
#
def relative_path(path1,path2):
    root = commonprefix([path1, path2])
    root = dirname(root)
    if (root == '/'):
        raise 'Error in relative_path - common directory cannot be / '
    if root :
        path2 = path2[len(root)+1:]   
    path1 = dirname(path1)
    while not ( root == path1 ):
        path2 = pardir + sep + path2 
        path1 = dirname(path1)
    return path2


##
# Change current working directory and create dir if necessary.
#
# This function attempts to change the current working directory to
# the directory specified by argument dir (like os.chdir), but first
# checks if the directory exists, and creates the directory and any
# missing parent directories if these do not yet exist. 
# 
# \param dir
#
def chdirs(dir):
    if not exists(dir):
       print 'Creating directory ' + dir
       os.makedirs(dir)
    os.chdir(dir)


##            
# Open file with specified path for writing, return file object.
#
# Similar to built-in function open(path,'w'), except that open_w
# will create any non-existent directories in the path as needed.
#
# \param path path to file
# \param file object with specified path, opened for writing
#
def open_w(path):
    dir  = dirname(path)
    if dir:
        if not exists(dir):
            print 'Creating directory ' + dir
            os.makedirs(dir)
    return open(path,'w')

##
# Remove a file if possible (if the path exists and it is a file).
#
# This function exits silently, without throwing an exception or returning
# an error code, if the file does not exist or is not a file.  
#
# \param path path to file
#
def rm(path):
   if os.path.exists(path):
      if os.path.isfile(path):
         os.remove(path)

## 
# Rename a file or directory, creating any necessary parent directories.
#
# Rename a file or directory. Similar to os.rename, except that this
# function will create any non-existent directories in the new path as
# needed.
#
# \param old_path - path of file or directory to be renamed
# \param new_path - new path for file or directory
#
def mv(old_path, new_path):
    if (not isfile(old_path)) and (not isdir(old_path) ):
        print 'Path ' + old_path + ' is not a file or directory'
        return
    new_dir = dirname(new_path)
    if new_dir:
        if not exists(new_dir):
            print 'Creating directory ' + new_dir
            os.makedirs(new_dir)
    return os.rename(old_path, new_path)


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

