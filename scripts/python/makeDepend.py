import os
import os.path
from string import *
from TextWrapper import *

def makeDepend(cfile, options, srcdir, blddir, extraDependencies=''):

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

   pfile = os.path.normpath(os.path.join(srcdir, base)) + '.p'
   dfile = os.path.normpath(os.path.join(blddir, base)) + '.d'

   # Create compiler command to calculate dependencies
   command  = 'g++ '
   command += options
   command += ' -MM -MF '
   command += pfile + ' '
   command += cfile
   #print command
   os.system(command)

   #Edit dependency file
   editDepend(pfile, dfile, blddir, extraDependencies)
   os.remove(pfile)

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
   text = TextWrapper('\\\n', 8)
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

