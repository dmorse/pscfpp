# -----------------------------------------------------------------------
#   This class is a tool to store all the PSCF "state files" created by 
#   a Sweep within a single object, where each "state file" can be
#   accessed by corresponding indices, and a "summary" method is provided
#   to aggregate data from the entire sweep into a single python list.
#
#   Creating a Sweep object:
# 
#      A Sweep object represents a collection of all PSCF "state files"
#      created by a Sweep. It stores each 'state file' by a State object 
#      within a Python list, and each 'state file' can be accessed by 
#      the corresponding index. 
#
#      User may collect and store all the 'state files' created by a
#      Sweep by creating a Sweep object, passing a param file name
#      or a directory of the folder with all the 'state files' as an
#      argument. The constructor reads and collects all 'state files'
#      (which end with '.dat') and return a Sweep object.
#
#      Example: 
#      1. By passing in a param file with name 'param':
#      
#         from pscfpp.sweep import *
#         s = Sweep('param')
#
#      2. By passing in a directory of the folder that contains all the 
#         'state files', with name 'out/':
#      
#         from pscfpp.sweep import *
#         s = Sweep('out/')
#
#   Accessing elements:
#      
#      A Sweep object can be treated as a list of 'state files',
#      so each state file store in a Sweep object as a State object
#      can be accessed by corresponding index with square brackets.
#      Example: s[0]    or    
#               s[1]
#
#      All elements and properties in each state file stored in a 
#      Sweep object can be accessed by the same format as the 
#      State object. Refer to the state class for details.
#      Example: s[0].param.Mixture.nMonomer    or    
#               s[1].thermo.fHelmholtz
#
#      Two special functions are built for a Sweep object to access
#      the same types of properties easily:
#
#      1.summary([list of property names, as strings], index = False)
#        User can return a python list containing a set of properties
#        (chosen by the user) at each state point in the Sweep.
#        The method returns a list of lists, where one list is given
#        for each state point, containing the values of all requested
#        properties at that state point. The properties are requested in
#        a list of strings, provided as an input, where each string 
#        gives the command that one would use to access that property
#        from the State class. Optionally, the index of each state point
#        in the sweep can be given as the first property in each list; 
#        this is turned off by default, and can be toggled using the 
#        'index' input parameter. 
#        Examples: 
#          s.summary(['param.Mixture.Polymer[0].phi','thermo.pressure'])
#          
#          The above command will return an array that may look something 
#          like the following:
#            [[0.5, 32.4415250701], [0.55, 30.6782301376], 
#             [0.6, 28.8344576024], [0.65, 26.8019750534],
#             [0.7, 24.5205812724]]
#
#                 or
#
#          s.summary(['param.Interaction.chi[0][1]'], index = True)
#          
#          The above command will return an array that may look something 
#          like the following:
#            [[0, 12.0], [1, 13.0], [2, 14.0], [3, 15.0], [4, 16.0]]
#
#      2.summaryString([list of property names, as strings]):
#        This method performs the same task as the summary method above,
#        but formats the resulting data into an organized data table
#        that is stored as a single string. Printing this string gives
#        a nice visual display of the data.
#        Example: 
#          print(s.summaryString(['param.Interaction.chi[0][1]',\
#                                 'thermo.fHelmholtz']))
#        
#          The above command gives an output that will look something
#          like the following:
#          
#            step      chi[0][1]     fHelmholtz
#               0  1.2000000e+01  1.9256750e+00
#               1  1.3000000e+01  2.1102042e+00
#               2  1.4000000e+01  2.2716872e+00
#               3  1.5000000e+01  2.4158122e+00
#               4  1.6000000e+01  2.5464487e+00
#
# Module Contents:
#      
#  class Sweep:
#          A Sweep object stores a collection of PSCF 'state files' 
#          created by a Sweep by using State object. Each stored 
#          'state file' can be accessed by the corresponding indices.
#          Two special functions are built in for a Sweep object to 
#          manage same types of properties easily.


from param import *
from thermo import *
from state import *
import os

class Sweep:
   '''
   Purpose:
      The class represents a collection of state files created by a 
      Sweep.
   Instance variables:
      sweep: a list to stored all read and parsed state files
   Methods:
      __init__(self, d):
         the constructor of the Sweep object for initiation, with one 
         argument: 
            d:
               the directory of the folder or the name of the param
               file, required
      summary(self, l, index = False):
         return a list of the same types properties from each state file
         according to the passed-in argument, l, list of expected 
         properties names, indices optionally
      summaryString(self,l):
         return the string of a table that shows the same types properties
         from each state file according to the passed-in argument, l,
         list of expected properties names
   '''
   def __init__(self, d):
      old_cwd = os.getcwd()

      if os.path.isdir(d) == True:
         os.chdir(d)
      elif os.path.isfile(d) == True:
         p = Composite(d)
         for x,y in p.children.items():
            if x.rfind('Sweep') != -1:
               nd = old_cwd + '/' + y.baseFileName
         try:
            os.chdir(nd)
         except:
            raise Exception('No Sweep block in the Param file')
      else:
         raise Exception('Not valid argument')

      filenames = []
      for file in os.listdir():
         if file.endswith('.dat'):
            filenames.append(file)
      filenames.sort()

      self.sweep = []
      for i in range(0, len(filenames)):
         s = State(filenames[i])
         self.sweep.append(s)

      os.chdir(old_cwd)

   def summary(self, l, index = False):
      n = len(l)

      summary = []

      for i in range(0, len(self.sweep)):
         if index == True:
            s = [i]
         else:
            s = []
         for j in range(0, n):
            a = self.sweep[i]
            string = 'a.' + l[j]
            try:
               val = eval(string)
            except RecursionError:
               raise Exception('Wrong command or values do not exist')
            s.append(val)
         summary.append(s)

      return summary

   def summaryString(self, l):
      n = len(l)

      summary = []

      for i in range(0, len(self.sweep)):
         s = [i]
         for j in range(0, n):
            a = self.sweep[i]
            string = 'a.' + l[j]
            try:
               val = eval(string)
            except RecursionError:
               raise Exception('Wrong command or values do not exist')
            s.append(val)
         summary.append(s)

      nl = [4]
      nameList = ['step']
      for i in range(0, n):
         index = l[i].rfind('.')
         name = l[i][index+1:]
         nameList.append(name)
         nl.append(len(name))

      valType = []
      for i in range(0, len(summary)):
         valType.append([])
         for j in range(0, len(summary[0])):
            valType[i].append(type(summary[i][j]))

      for i in range(0, len(summary)):
         for j in range(0, len(summary[0])):
            length = len(str(summary[i][j]))
            if (valType[i][j] == str) and (length > nl[j]):
               nl[j] = length
            if (valType[i][j] == float) and (13 > nl[j]):
               nl[j] = 13
            if (valType[i][j] == int) and (length > nl[j]):
               nl[j] = length


      summaryString = ' '
      for i in range(0, len(nameList)):
         stringFormat = '{:>' + str(nl[i]) + 's}'
         if i != len(nameList)-1:
            summaryString += stringFormat.format(nameList[i]) + '  '
         else:
            summaryString += stringFormat.format(nameList[i]) + '\n '

      for i in range(0, len(summary)):
         for j in range(0, len(summary[0])):
            if valType[i][j] == int:
               stringFormat = '{:>' + str(nl[j]) + '}'
               summaryString += stringFormat.format(summary[i][j])
            if valType[i][j] == float:
               stringFormat = '{:.7e}'
               val = stringFormat.format(summary[i][j])
               summaryString += val
            if valType[i][j] == str:
               stringFormat = '{:>' + str(nl[j]) + 's}'
               summaryString += stringFormat.format(summary[i][j])
            if j == len(summary[0])-1:
               summaryString += '\n '
            else:
               summaryString += '  '

      return summaryString

   def __getitem__(self, key):
      return self.sweep[key]

