# -----------------------------------------------------------------------
#   This class is a tool to parse a PSCF "field file" and store 
#   all values within it in a single object. A field file contains both
#   a block of text representing the field file header and a subsequent
#   block of text containing all the field data. Field files for all 
#   three types of field formats (rgrid, kgrid and basis) can be parsed 
#   and stored by the Field object. For more detail on the three types of  
#   field formats, please refer to the PSCFPP main documentation. Users 
#   can access and modify the stored values of the parameters in the 
#   header and data part after parsing by using specific statements 
#   (commands), and can write the entire object to a file in proper 
#   format.
#
#   Parsing a field file:
#
#       A Field object represents a PSCF field file that always contains 
#       two parts, a header part and a data part. It stores all the
#       contents within the file in a form that allows each value
#       within it from each part to be accessed and modified.
#
#       The header part always starts with the format version number
#       and the data part follows right after the header without any
#       empty lines. Users may parse such a file by creating a Field 
#       object, passing in the name of field file as an argument. The
#       constructor parses the file and return a Field object that
#       contains its contents.
#
#       Example: To read and parse a field file with name 'field',
#       execute the following code within a python3 interpreter: 
#
#           from pscfpp.field import *
#           f = Field('field')
#
#       A Field object can be written to a file by calling the 
#       writeOut() method, passing in the string of the file name to 
#       which the Field should be written.
#
#           Example: f.writeOut('fieldOut')
#
#   Accessing elements:
#
#       A Field object contains three members, header, data and type,
#       corresponding to the different sections of the field file and 
#       its field format. Users can retrieve any member by dot notation:
#
#           1.header: the header member can be accessed by calling the
#             'header' attribute, which returns a dictionary that stores
#             all the properties in the header. All stored contents within 
#             it can be accessed by the formats listed below:
#             Example: 1. accessing the whole header:
#                         f.header 
#                      2. accessing single property in the header:
#                         f.header['format']    or
#                         f.header['N_monomer']
#
#           2.data: the data member can be accessed by calling the 'data'
#             attribute, which returns a list of lists that stores all 
#             the data of the field. A single row of data can be accessed 
#             by corresponding row number. Any specific data points can be 
#             accessed by corresponding two separate square bracket indices.
#             Example: f.data    or
#                      f.data[4]    or
#                      f.data[0][6]
#
#       The parser also allows users to modify the properties of the header 
#       in the preset format. Users have to modify the desired properties 
#       with the same types of objects for the original ones, otherwise,
#       the parser will not work as expected.
#
#       Example: f.header['dim'] = 3    or
#                f.header['cell_param'] = [2.45]   or
#                f.header['group_name'] = 'P_1'
#
#      Three special functions are built for a Field object to modify the 
#      data stored within it:
#
#      1.addColumn(index, element)
#        User can add a new column to the data list of the Field object by
#        calling this method. By passing in two arguments, index, the 
#        desired position to add the new column, and element, a single 
#        value for the whole column or a list of values represents the whole 
#        column, the data list of the Field object can be updated with 
#        the desired position and values.
#        Example: f.addColumn(2, 1.2)    or
#                 f.addColumn(2, [list of values for the whole column])
#
#      2.deleteColumn(index)
#        User can delete an exist column from the data list of the Field 
#        object by calling this method. By passing in one argument, index, 
#        the desired position of the column to be deleted, the data list 
#        of the Field object can be updated without the indicated column.
#        Example: f.deleteColumn(2)
#
#      3.reorder(order)
#        User can reorder the data list of the Field object with desired 
#        order by calling this method. By passing in one argument, order,
#        an integer list that represents the new order of the columns in
#        the data list, the data list of the Field object can be updated 
#        with the desired order.
#        Example: f.reorder([0, 3, 4, 1, 2])
#
# Module Contents:
#
#   class Field:
#           A Field object contains all the data from a field file. It
#           has three members, header, data and type, where the header member
#           contains all the data in the header of the field file starting
#           with the format version, the data member contains all the data of
#           the field, and the type member stores the field format. This class 
#           parses a field file and stores all contents within it in a form
#           that allows all properties and data to be accessed and modified.
#
#   def getValue(v):
#           A function to distinguish the correct type of the value from the 
#           passed-in string, v, and return it with the correct type. 

class Field:
   '''
   Purpose:
      The class represents a filed file
   Instance variables:
      data: 
         the list stores all the data from the field file
      header: 
         the dictionary stores all the header information from the field
         class
      type:
         a string stores the type of the field file 
   Methods:
   __init__(self, filename):
      constructor, with one argument:
         filename: 
            the filename that needs to be read
   read(self, openFile):
      method to read the open Field file, openFile as the argument, line
      by line and update the read items into the data and header, and
      determine the type of the Field file
   writeOut(self, filename):
      method to write out the stored Field file to a specific named file
      with the name of the argument filename
   writeOutStirng(self):
      return the string for writing out
   addColumn(self, index, element):
      method to add a new column to the data list, with two arguments:
         index, an integer represents desired position to add the new column
         element, an single value for the whole column or a list of values
                  represents the whole column that need to be added
   deleteColumn(self, index):
      method to delete an exist column from the data list, with one argument:
         index, an integer represents desired column to be deleted
   reorder(self, order):
      method to reorder the columns in the data list, with one argument:
         order, a list of integer represents the new order of the exist 
                columns
   '''

   def __init__(self, filename):
      self.data = []
      self.header = {}
      self.type = ''
      with open(filename) as f:
         self.read(f)

   def read(self, openFile):
      line = openFile.readline()
      l = line.split()
      if l[0] == 'format':
         self.header['format'] = [getValue(l[1]), getValue(l[2])]
      else:
         raise Exception('Not valid field file.')

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

      while line != '':
         d = []
         for i in range(0, len(l)):
            d.append(getValue(l[i]))
         self.data.append(d)
         line = openFile.readline()
         l = line.split()

      if self.type != 'basis':
         if type(self.data[0][0]) is int:
            self.type = 'kgrid'
         else:
            self.type = 'rgrid'

   def writeOut(self, filename):
      with open(filename, 'w') as f:
         f.write(self.writeOutString())

   def writeOutString(self):
      out = ''
      for x, y in self.header.items():
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

   def addColumn(self, index, element):
      if type(element) is list:
         for i in range(0, len(self.data)):
            self.data[i].insert(index, element[i])
      else:
         for i in range(0, len(self.data)):
            self.data[i].insert(index, element)

   def deleteColumn(self, index):
      for i in range(0, len(self.data)):
         self.data[i].pop(index)

   def reorder(self, order):
      newData = []
      for i in range(0, len(self.data)):
         d = []
         for j in range(0,len(order)):
            d.append(self.data[i][order[j]])
         newData.append(d)
      self.data = newData

# End class Field  -------------------------------------------------------

def getValue(v):
   '''
   Purpose:
      Distinguishing the correct type of the value form a string
   Argument(s):
      v: string represents a value
   Return:
      return the value with the correct type
   '''
   if (v.isdigit() == True) or (v[0] == '-' and v[1:].isdigit() == True):
      val = int(v)  
   else:
      try:
         val = float(v)
      except ValueError:
         val = str(v)
   
   return val


# End function getValue  -------------------------------------------------
