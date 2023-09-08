"""! Module for parsing field files. """

##
# Container for data in field files.
#
#  This class is a tool to parse a PSCF "field file" and store 
#  all values within it in a single object. A field file contains both
#  a block of text representing the field file header and a subsequent
#  block of text containing all the field data. Field files for all 
#  types of field formats (fd1d, rgrid, kgrid and basis) can be parsed 
#  and stored by the Field object. For more detail on the types of  
#  field formats, please refer to the PSCFPP main documentation. Users 
#  can access and modify the stored values of the parameters in the 
#  header and data part after parsing by using specific statements 
#  (commands), and can write the entire object to a file in proper 
#  format.
#
#  A Field object represents a PSCF field file that always contains 
#  two parts, a header part and a data part. It stores all the
#  contents within the file in a form that allows each value
#  within it from each part to be accessed and modified.
#
#  **Construction:**
#
#      The header part always starts with the format version number
#      and the data part follows right after the header without any
#      empty lines. Users may parse such a file by creating a Field 
#      object, passing in the name of field file as an argument. The
#      constructor parses the file and return a Field object that
#      contains its contents.
#
#      Example:
#
#      To read and parse a field file with name 'field':
#        \code
#           from pscfpp.field import *
#           f = Field('field')
#        \endcode
#
#  **Writing out:**
#
#      A Field object can be written to a file by calling the 
#      writeOut() method, passing in the string of the file name to 
#      which the Field should be written.
#
#      Example:
#        \code
#           f.write('fieldOut')
#        \endcode
#
#  **Accessing elements:**
#
#      A Field object contains three members, header, data and type,
#      corresponding to the different sections of the field file and 
#      its field format. Users can retrieve any member by dot notation:
#
#           1.header: the header member can be accessed by calling the
#             'header' attribute, which returns a dictionary that stores
#             all the properties in the header. All stored contents within 
#             it can be accessed by the formats listed below:
#
#             Example: 
#
#               1. accessing the whole header:
#                  \code
#                     f.header
#                  \endcode 
#               2. accessing single property in the header:
#                  \code
#                     f.header['format']    
#                     f.header['N_monomer']
#                  \endcode
#
#  **Modifying elements:**
#
#      The parser also allows users to modify the properties of the header 
#      in the preset format. Users have to modify the desired properties 
#      with the same types of objects for the original ones, otherwise,
#      the parser will not work as expected.
#
#      Example: 
#        \code
#           f.header['dim'] = 3  
#           f.header['cell_param'] = [2.45]   
#           f.header['group_name'] = 'P_1'
#        \endcode
#
#  **Modyfying data table:**
#
#      See documentationd for method addColumn, deleteColumn, and
#      reorder.
#
class Field:

   ##
   # Constructor.
   #
   # \param filename  a filename string
   #
   def __init__(self, filename):
      self.data = []
      self.header = {}
      self.type = ''
      with open(filename) as f:
         self.read(f)

   ##
   # Read and parse the passed-in file
   # 
   # This function reads the passed-in open-file object line
   # by line and update the read items into instance
   # variables. The reading stops when all lines in file are
   # read.
   #
   # \param openFile  a file object, open for reading
   #
   def read(self, openFile):
      line = openFile.readline()
      l = line.split()
      if l[0] == 'format':
         self.header['format'] = [getValue(l[1]), getValue(l[2])]
      elif l[0] == 'nx':
         self.type ='fd1d'
         self.header['nx'] = getValue(l[1])
      else:
         raise Exception('Not valid field file.')

      if self.type != 'fd1d':
         lineCount = 1
         while lineCount <= 15:
            line = openFile.readline()
            l = line.split()
            lineCount += 1
            if lineCount%2 == 0:
               name = l[0]
               if name == 'N_basis':
                  self.type = 'basis'
            else:
               if name == 'mesh' or name == 'ngrid' or name == 'cell_param':
                  d = []
                  for i in range(0, len(l)):
                     d.append(getValue(l[i]))
                  self.header[name] = d
               else:
                  self.header[name] = getValue(l[0])
      else:
            line = openFile.readline()
            l = line.split()
            self.header[l[0]] = getValue(l[1])
            line = openFile.readline()
            l = line.split()

      while line != '':
         d = []
         for i in range(0, len(l)):
            d.append(getValue(l[i]))
         self.data.append(d)
         line = openFile.readline()
         l = line.split()
         

      if self.type != 'basis' and self.type != 'fd1d':
         if type(self.data[0][0]) is int:
            self.type = 'kgrid'
         else:
            self.type = 'rgrid'

   ##
   # Write out field to a file.
   #
   # This function writes out the field file string to the 
   # specified file with the name of the passed-in parameter, 
   # filename.
   #
   # \param filename  a filename string.
   #
   def write(self, filename):
      with open(filename, 'w') as f:
         f.write(self.__str__())

   ##
   # Return string representation of this Field.
   #
   # This function return the string representation of the Field object 
   # in the appropriate field file format.
   #
   def __str__(self):
      out = ''
      for x, y in self.header.items():
         if self.type == 'fd1d':
            o = f'{x:<7}' + str(y) + '\n'
         else:
            if x == 'format':
               o = x + '   ' + str(y[0]) + '   ' + str(y[1]) + '\n'
            else:
               o = x + '\n'
               if x == 'cell_param':
                  val = f'{y[0]:.10e}'
                  o += f'{val:>20}'
                  if len(y) == 1:
                     o += '\n'
                  for i in range(1, len(y)):
                     val = f'{y[i]:.10e}'
                     o += f'{val:>18}'
                     if i == len(y)-1:
                        o += '\n'
               elif x == 'mesh' or x == 'ngrid':
                  o += f'{y[0]:>20}'
                  if len(y) == 1:
                     o += '\n'
                  for i in range(1, len(y)):
                     o += f'{y[i]:>8}'
                     if i == len(y)-1:
                        o += '\n'
               else:
                  o += f'{str(y):>20s}' + '\n'

         out += o

      for i in range(0, len(self.data)):
         row = ''
         for j in range(0, len(self.data[0])):
            if self.type == 'fd1d':
               if j == 0:
                  row += f'{self.data[i][j]:>6}'
               else:
                  val = f'{self.data[i][j]:.11e}'
                  row += f'{val:>20}'

            if self.type == 'basis':
               if j < self.header['N_monomer']:
                  val = f'{self.data[i][j]:.10e}'
                  row += f'{val:>20}'
               else:
                  if j == self.header['N_monomer']:
                     row += '   '
                  row += f'{self.data[i][j]:>5}'

            if self.type == 'rgrid':
               val = f'{self.data[i][j]:.15e}'
               row += f'{val:>23}'

            if self.type == 'kgrid':
               if j == 0:
                  row += f'{self.data[i][j]:>6}'
               else:
                  val = f'{self.data[i][j]:.12e}'
                  if j%2 != 0:
                     row += f'{val:>22}'
                  else:
                     row += f'{val:>20}'

         out += row + '\n'

      return out

   ##
   # Add a new column to the data list.
   #
   # This function adds a new column to the data list of the 
   # Field object. By passing in two parameters, index, the 
   # desired position to add the new column, and element, a 
   # single value for the whole column a a list of values
   # represents the whole column, the datalist of the Field
   # object can be updated with the desired position and values.
   #
   # Example:
   #   \code
   #      f.addColumn(2, 1.2)    
   #      f.addColumn(2, [list of values for the whole column])
   #   \endcode
   #
   # \param index  an position integer.
   # \param element  a single value or a list of values.
   #
   def addColumn(self, index, element):
      if type(element) is list:
         for i in range(0, len(self.data)):
            self.data[i].insert(index, element[i])
      else:
         for i in range(0, len(self.data)):
            self.data[i].insert(index, element)

   ##
   # Delete an exist column from the sata list.
   #
   # This function deletes an exist column from the data list of
   # the Field object. By passing in one parameter, index, the 
   # desired position of the column to be deleted, the data list
   # of the Field object can be updated without the indicated column.
   #
   # Example:
   #   \code
   #      f.deleteColumn(2)
   #   \endcode
   #
   # \param index  an position integer.
   def deleteColumn(self, index):
      for i in range(0, len(self.data)):
         self.data[i].pop(index)

   ##
   # Reorder the data list.
   #
   # This function reorders the data list of the Field object. By
   # passing in one parameter, order, an integer list that represents
   # the new order of the columns in the data list, the data list of
   # the Field object can be updated with the desired order.
   #
   # Example:
   #   \code
   #      f.reorder([0, 3, 4, 1, 2])
   #   \endcode
   #
   # \param order  a list of position integers.
   #
   def reorder(self, order):
      newData = []
      for i in range(0, len(self.data)):
         d = []
         for j in range(0,len(order)):
            d.append(self.data[i][order[j]])
         newData.append(d)
      self.data = newData

##
# Distinguish the correct type of the value from a string.
#
# This function is a helper function to read a open-file
# object. It distinguishes the correct type of the value
# from a string to either integer, float or string.
#
# Return value:
#
# The exact value with the correct type.
#
# \param v  a string represents a value.
#
def getValue(v):
   if (v.isdigit() == True) or (v[0] == '-' and v[1:].isdigit() == True):
      val = int(v)  
   else:
      try:
         val = float(v)
      except ValueError:
         val = str(v)
   
   return val

