import os
import os.path
import glob
from string import *

from file_util   import *
from readLabel   import *
from makeDepend  import *
from TextWrapper import *

class MakeMaker:

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
      
      wrapper = TextWrapper('\\\n', 4)

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
         print x
      for x in self.srcs:
         print x
      for x in self.dirs:
         print x + os.sep
