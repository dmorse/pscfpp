'''! Python scripts used by the PSCF makefile build system. '''

import os
import os.path
import glob
from string import *

from pscfpp.file import *
from pscfpp.text import *

##
# Create a *.d dependency file for a C/C++ source file.
#
# This command uses a compiler executable to analyze the dependencies of
# a C/C++ source file in the src directory tree and creates a dependency 
# file with an appropriate format for make in the analogous location in 
# the bld directory tree. The function first executes a compiler command
# to compute dependencies and then calls python function editDepend to 
# edit the resulting file. The most important change made during editing 
# is the replacement of relative paths output by the g++ compiler to the
# absolute paths used in the PSCF build system.
#
# \param processor  compiler command to analyze dependencies (i.e., g++)
# \param options  string of options passed to processor
# \param cfile  path to source file
# \param srcdir  path to source directory tree (SRC_DIR in makefiles)
# \param blddir  path to bulkd directory tree (BLD_DIR in makefiles)
# \param extraDependencies  string of additional known dependencies 
#
def createDependencyFileCpp(processor, options, cfile, srcdir, blddir, extraDependencies=''):

   # Isolate base source file name 
   base = os.path.basename(cfile)
   groups = base.split('.')
   base = groups[0]

   # Calculate source file directory path relative to srcdir
   abspath = os.path.abspath(cfile)
   relpath = os.path.relpath(abspath, srcdir)
   reldir  = os.path.dirname(relpath)

   # Construct paths to subdirectories of build and src directories
   if (reldir != '.' and reldir != ''):
      blddir = os.path.normpath(os.path.join(blddir, reldir))
      srcdir = os.path.normpath(os.path.join(srcdir, reldir))

   # Construct paths to dependency files
   # pfile - temporary file created by the compiler (erased)
   # dfile - final dependnecy file
   pfile = os.path.normpath(os.path.join(srcdir, base)) + '.p'
   dfile = os.path.normpath(os.path.join(blddir, base)) + '.d'

   # Create and execute compiler command to calculate dependencies
   command = processor + ' '
   command += pfile + ' '
   command += options + ' '
   command += cfile
   #print('\n' + command + '\n')
   os.system(command)

   #Edit dependency file
   editDepend(pfile, dfile, blddir, extraDependencies)
   os.remove(pfile)


##
# Create a *.d dependency file for a CUDA source file.
#
# This command uses a compiler executable to analyze the dependencies of
# a CUDA source file in the src directory tree and creates a dependency 
# file with an appropriate format for make in the analogous location in 
# the bld directory tree. The function first executes a compiler command
# to compute dependencies and then calls python function editDependLocal 
# to edit the resulting file. The most important change made during editing
# is the replacement of relative paths output by the nvcc compiler to the
# absolute paths used in the PSCF build system.
#
# \param processor  compiler command to analyze dependencies (i.e., g++)
# \param options  string of options passed to processor
# \param cfile  path to source file
# \param srcdir  path to source directory tree (SRC_DIR in makefiles)
# \param blddir  path to build directory tree (BLD_DIR in makefiles)
# \param extraDependencies  string of additional known dependencies 
#
def createDependencyFileCuda(processor, options, cfile, srcdir, blddir, extraDependencies=''):

   # Store unmodified srcdir
   srcroot = srcdir

   # Isolate base source file name 
   base = os.path.basename(cfile)
   groups = base.split('.')
   base = groups[0]

   # Calculate source file directory path relative to srcdir
   abspath = os.path.abspath(cfile)
   relpath = os.path.relpath(abspath, srcdir)
   reldir  = os.path.dirname(relpath)

   # Construct paths to subdirectories of build and src directories
   if (reldir != '.' and reldir != ''):
      blddir = os.path.normpath(os.path.join(blddir, reldir))
      srcdir = os.path.normpath(os.path.join(srcdir, reldir))

   # Construct paths to dependency files
   # pfile - temporary file created by the compiler (erased)
   # dfile - final dependnecy file
   pfile = os.path.normpath(os.path.join(srcdir, base)) + '.p'
   dfile = os.path.normpath(os.path.join(blddir, base)) + '.d'

   # Create and execute compiler command to calculate dependencies
   command = processor + ' '
   command += options + ' '
   command += cfile
   command += ' > ' + pfile
   os.system(command)

   # Edit dependency file and remove temporary pfile
   editDependLocal(pfile, dfile, blddir, srcroot, extraDependencies)
   os.remove(pfile)


##
# Edit the dependency file created for a C/C++ file by the compiler.
#
# This function edits the dependency file created by the g++ or compatible
# compiler so as to use absolute rather than relative paths, and so as to
# include the extraDependencies string.
#
# \param pfile  input dependency file name
# \param dfile  output dependency file name
# \param blddir  path to build directory tree (BLD_DIR in makefiles)
# \param extraDependencies  string of additional known dependencies 
#
def editDepend(pfile, dfile, blddir, extraDependencies):
   file  = open(pfile, 'r')
   lines = file.readlines()
   file.close()

   # Extract target from first line
   groups   = lines[0].split(":")
   target   = groups[0]
   lines[0] = groups[1]

   # Replace target by file in build directory
   base = os.path.basename(target)
   target = os.path.normpath(os.path.join(blddir, base))

   # Replace dependencies by absolute paths
   text = Wrapper('\\\n', 8)
   for i in range(len(lines)):
       line = lines[i]
       if line[-1] == '\n':
           line = line[:-1]
       if line[-1] == '\\':
           line = line[:-1]
       lines[i] = line
       if i == 0:
           text.append(target + ': ')
       deps = line.split()
       for dep in deps:
           path = os.path.abspath(dep)
           text.append(path)

   # Process extraDependencies (if any)
   if extraDependencies:
       deps = extraDependencies.split()
       for dep in deps:
          path = os.path.abspath(dep)
          text.append(path)

   file  = open(dfile, 'w')
   file.write(str(text))
   file.close()

##
# Edit the dependency file created for a CUDA file by the compiler.
#
# This function edits the dependency file created by the nvcc compiler for
# a CUDA file so as to use absolute rather than relative paths, and so as 
# to include the extraDependencies string.
#
# \param pfile  input dependency file name
# \param dfile  output dependency file name
# \param blddir  path to build directory tree (BLD_DIR in makefiles)
# \param srcroot  path to source directory tree (SRC_DIR in makefiles)
# \param extraDependencies  string of additional known dependencies 
#
def editDependLocal(pfile, dfile, blddir, srcroot, extraDependencies):
   file  = open(pfile, 'r')
   lines = file.readlines()
   file.close()

   # Extract target from first line
   groups   = lines[0].split(":")
   target   = groups[0]
   lines[0] = groups[1]

   # Replace target by file in build directory
   base = os.path.basename(target)
   target = os.path.normpath(os.path.join(blddir, base))

   # Replace local dependencies by absolute paths
   text = Wrapper('\\\n', 8)
   for i in range(len(lines)):
       line = lines[i]
       if line[-1] == '\n':
           line = line[:-1]
       if line[-1] == '\\':
           line = line[:-1]
       lines[i] = line
       if i == 0:
           text.append(target + ': ')
       deps = line.split()
       for dep in deps:
           pair = [srcroot, dep]
           if (os.path.commonprefix(pair) == srcroot):
              path = os.path.abspath(dep)
              text.append(path)

   # Process extraDependencies (if any)
   if extraDependencies:
       deps = extraDependencies.split()
       for dep in deps:
          path = os.path.abspath(dep)
          text.append(path)

   file  = open(dfile, 'w')
   file.write(str(text))
   file.close()

##
# Class to construct makefile system for a set of source files.
# 
#
class MakeMaker:

   ##
   # Constructor.
   #
   # \param  path  path from working directory
   # \param  pathFromSrc  path to working directory from src directory
   # \param  pathToSrc  path to src directory 
   # \param  parent  parent MakeMaker object (for parent directory)
   # \param  isTest  boolean flag, perform dry run if true
   #
   def __init__(self, path = '.', pathFromSrc = '.', pathToSrc ='.', parent = None, isTest=False):

      self.path          = path         # path from working directory
      self.pathFromSrc   = pathFromSrc  # path from $(SRC_DIR)
      self.pathToSrc     = pathToSrc    # path to $(SRC_DIR)
      self.parent        = parent
      self.isTest        = isTest
      self.hasLib        = False
      self.globalInclude = ''
      self.defines       = ''

      if parent:
         self.root = parent.root
         self.generation = parent.generation + 1
      else:
         self.root = self
         self.generation = 0

      # Parse directory path 
      if self.pathFromSrc == '.':
         self.dirparts = []
         self.dirname  = 'SRC'
      else:
         self.dirparts = self.pathFromSrc.split(os.sep)
         self.dirname  = '_'.join(self.dirparts)

      # Read dir.txt file
      dirFilePath = self.makePath('dir.txt')
      if os.path.exists(dirFilePath):
         file = open(dirFilePath,'r')
         self.dirList = readLabelledList(file, 'subDirectories')
         file.close()
      else:
         self.dirList = []

      # Create a child MakeMaker for each directory
      self.dirs = []
      for dir in self.dirList:
         path = self.makePath(dir)
         if os.path.isdir(path):
            if self.pathFromSrc == '.':
               pathFromSrc = dir
               pathToSrc   = '..'
            else:
               pathFromSrc = self.pathFromSrc + os.sep + dir
               pathToSrc   = '..' + os.sep + self.pathToSrc
            if dir == 'tests':
               maker = MakeMaker(path, pathFromSrc, pathToSrc, self, True)
            else:
               maker = MakeMaker(path, pathFromSrc, pathToSrc, self, self.isTest)
            self.dirs.append(maker)
         
   def setGlobalInclude(self, text):
      self.globalInclude = text
      for dir in self.dirs:
         dir.setGlobalInclude(text)

   def setDefines(self, text):
      self.defines = text
      for dir in self.dirs:
         dir.setDefines(text)

   def addLibrary(self, libName, libObjs):
      self.hasLib  = True
      self.libName = libName
      self.libObjs = libObjs

   def setLinkObjs(self, text):
      self.linkObjs = text
      for dir in self.dirs:
         dir.setLinkObjs(text)

   def makeInclude(self, base):
      ''' Make include line suitable for a makefile'''
      if self.pathFromSrc == '.':
         return 'include $(SRC_DIR)' + os.sep +  base + '\n'
      else:
         return 'include $(SRC_DIR)' + os.sep + self.pathFromSrc + os.sep + base + '\n'

   def srcSuffix(self):
      ''' Return suffix for source files: cpp or cc '''
      if (self.isTest):
         return 'cc'
      else:
         return 'cpp'

   def find(self):
      ''' 
      Find all header and source files in this directory.
      Note: Does not recursively descend into subdirectories
      '''
      # Find source files
      self.srcs = []
      self.hasSrcs = False
      for name in os.listdir(self.path):
         if name != '.':
            path = self.makePath(name)
            if os.path.isfile(path):
               base = os.path.basename(path)
               groups = base.split('.')
               if len(groups) == 2:
                  base   = groups[0]
                  suffix = groups[1]
                  if suffix == self.srcSuffix():
                     path = self.makePathFromSrc(base)
                     self.srcs.append(base)
                     self.hasSrcs = True

      # Find *.h header files
      self.hdrs = []
      self.hasHdrs = False
      for name in os.listdir(self.path):
         if name != '.':
            path = self.makePath(name)
            if os.path.isfile(path):
               base   = os.path.basename(path)
               groups = base.split('.')
               if len(groups) == 2:
                  base   = groups[0]
                  suffix = groups[1]
                  if suffix == 'h':
                     if base not in self.srcs:
                        self.hdrs.append(base)
                        self.hasHdrs = True

   def make(self):

      # Find all the source and header files in this directory
      self.find()

      # Recursively descend into subdirectories
      for dir in self.dirs:
         dir.make()

      # ---------------------------------------------------------------------
      # Open and write sources.mk file
      file = open_w(self.makePath('sources.mk'))

      # Include subdirectory sources.mk files
      if len(self.dirs) > 0:
         for dir in self.dirs:
            basename = os.path.basename(dir.path)
            if basename != 'tests':
               varpath  = basename + os.sep + 'sources.mk'
               file.write(self.makeInclude(varpath))
         file.write('\n')

      # Construct prefix for aggregate names
      if self.dirname == 'SRC':
         prefix = ''
      else:
         prefix = self.dirname + '_'
      
      wrapper = Wrapper('\\\n', 4)

      # Write aggregate definition for SRCS
      wrapper.clear()
      wrapper.append(prefix + 'SRCS=')
      for dir in self.dirs:
         basename = os.path.basename(dir.path)
         if basename != 'tests':
            wrapper.append('$(' + dir.dirname + '_SRCS) ')
      for base in self.srcs:
         name  = '$(SRC_DIR)' + os.sep + self.makePathFromSrc(base)
         name += '.' + self.srcSuffix() + ' '
         wrapper.append(name)
      file.write(str(wrapper))
      file.write('\n\n')

      # Write aggregate definition for OBJS
      file.write(prefix + 'OBJS=$(')
      file.write(prefix + 'SRCS:.' + self.srcSuffix() + '=.o)')
      file.write('\n\n')
      #file.write(prefix + 'DEPS=$(')
      #file.write(prefix + 'SRCS:.' + self.srcSuffix() + '=.d)')
      #file.write('\n\n')

      # Add library target, if any
      if self.hasLib:
         file.write(self.libName + ': ' + self.libObjs + '\n')
         file.write('\t$(AR) rcs ' + self.libName + ' ' + self.libObjs + '\n')
         file.write('\n')

      file.close()

      # --------------------------------------------------------------------
      # Write makefile

      file = open(self.makePath('makefile'), 'w')
      file.write('SRC_DIR_REL =' + self.pathToSrc + '\n\n')
      file.write(self.globalInclude) 
      if self.isTest:
         file.write('include $(SRC_DIR_REL)/test/sources.mk\n')
         file.write('include sources.mk\n')
      file.write('\n') 
      if self.dirname == 'SRC':
         objs = 'OBJS'
      else:
         objs = self.dirname + '_OBJS'
      targets = '$(' + objs + ')'
      deps = '$(' + objs + ':.o=.d)'
      if self.isTest:
         for base in self.srcs:
            targets += ' ' + base
      file.write('all: ' + targets)
      if self.hasLib:
          file.write(' ' + self.libName)
      file.write('\n\n')
      file.write('clean:\n\trm -f ' + targets + ' ' + deps)
      if self.hasLib:
          file.write(' ' + self.libName)
      file.write('\n\n')
      file.write('clean-deps:\n\trm -f ' + deps)
      file.write('\n\n')
      if self.isTest:
         for base in self.srcs:
            file.write(base + ': ' + base + '.o ' + self.linkObjs + '\n')
            file.write('\t$(CXX) $(LDFLAGS) $(INCLUDES) $(DEFINES)')
            file.write(' -o ' + base + ' ' + base + '.o \\\n')
            file.write('\t       ' + self.linkObjs + '\n\n')
      #file.write('-include $(' + objs + ':.o=.d)\n\n')
      file.write('-include ' + deps + '\n\n')
      file.close()

   def makePath(self, string):
      if self.path and self.path != '.':
         return self.path + os.sep + string
      else:
         return string

   def makePathFromSrc(self, string):
      if self.pathFromSrc and self.pathFromSrc != '.':
         return self.pathFromSrc + os.sep + string
      else:
         return string

   def srcDirPath(self):
      if self.generation == 0:
         return self.path
      else:
         shifts = []
         for i in range(self.generation):
            shifts.append('..')
         return os.sep.join(shifts)

   def clear(self):
      self.hdrs = []
      self.srcs = []
      self.dirs = []

   def filenames(self, pattern = '*', recursive = 1):
      r = glob.glob( self.path + '/' + pattern)
      if recursive:
         for dir in self.dirs:
            r.extend(dir.filenames(pattern))
      return r

   def __repr__(self):
      r = []
      r.append( self.path + '\n' ) 
      for x in self.hdrs:
         r.append( str(x) + '\n' )
      for x in self.srcs.keys() :
         r.append( str(x) + '\n' )
      for x in self.dirs:
         r.append( str(x)  )
      return join(r,'')

   def __str__(self):
      return __repr__()

   def ls(self):
      for x in self.hdrs:
         print(x)
      for x in self.srcs:
         print(x)
      for x in self.dirs:
         print(x + os.sep)

