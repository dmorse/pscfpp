#****M root/file_util ---------------------------------------------
# MODULE
#   file_util
# PURPOSE
#   Utility functions for manipulating files and paths
#*** --------------------------------------------------------------
import os
from os.path import *

#****f file_util/relative_path  -----------------------------------
# FUNCTION
#   relative_path(path1, path2)
# PURPOSE
#   Generates relative path of path2 relative to path1. 
# ARGUMENTS
#   path1 and path2 are paths relative to a common directory
# RETURN
#   Path for file2 relative to directory containing path1
#*** ---------------------------------------------------------------
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


#****f file_util/chdirs -------------------------------------------
# FUNCTION
#   chdirs(dir)
# PURPOSE
#   Change current working directory to dir (like os.chdir), and
#   create dir and any necessary parents if it does not yet exist.
#*** --------------------------------------------------------------
def chdirs(dir):
    if not exists(dir):
       print 'Creating directory ' + dir
       os.makedirs(dir)
    os.chdir(dir)

            
#****f file_util/open_w -------------------------------------------
# FUNCTION
#   open_w(path)
# PURPOSE
#   Open file with specified path for writing, return file object.
#   Similar to built-in function open(path,'w'), except that open_w
#   will create any non-existent directories in the path.
# RETURN
#   file object with specified path, opened for writing
#*** --------------------------------------------------------------
def open_w(path):
    dir  = dirname(path)
    if dir:
        if not exists(dir):
            print 'Creating directory ' + dir
            os.makedirs(dir)
    return open(path,'w')

#****f file_util/rm -----------------------------------------------
# FUNCTION
#   rm(path)
# PURPOSE
#   Remove a file if possible (if the path exists and is a file).
#   Does not throw an exception if it doesn't exist or isn't a file.
#*** --------------------------------------------------------------
def rm(path):
   if os.path.exists(path):
      if os.path.isfile(path):
         os.remove(path)

#****f file_util/rename_w ------------------------------------------
# FUNCTION
#   rename_w(old_path, new_path)
# ARGUMENTS
#   old_path - path of file or directory to be 
# PURPOSE
#   Rename a file or directory. Similar to os.rename,but will 
#   create any non-existent directories in the path
#*** --------------------------------------------------------------
def rename_w(old_path,new_path):
    if (not isfile(old_path)) and (not isdir(old_path) ):
        print 'Path ' + old_path + ' is not a file or directory'
        return
    new_dir = dirname(new_path)
    if new_dir:
        if not exists(new_dir):
            print 'Creating directory ' + new_dir
            os.makedirs(new_dir)
    return os.rename(old_path, new_path)

