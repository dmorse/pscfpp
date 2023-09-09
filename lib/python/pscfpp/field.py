"""! Module for parsing field files. """

##
# Container for data in field files.
#
#  A Field object can parse a PSCF "field file" and store its contents.
#  A field file contains both a header section and a data section. 
#  This class can read and parse all of the field file formats used by
#  PSCF programs including the basis, rgrid, and kgrid formats used by
#  the pscf_pc and pscf_pg programs to represent periodic fields, and 
#  the simpler file format used by pscf_fd for one-dimensional fields.
#  For detail on these file formats see the section of the web manual
#  on \ref user_field_page "field files". Users can access and modify
#  stored values of the parameters in the header and values in the data 
#  section after parsing by using specific statements (commands), and 
#  can write the entire object to a file in proper format.
#
#  Every PSCF field file format contains a header section and a data 
#  section. The header contains globale information about the system
#  and spatial discretization. The data section is a table of field
#  components in which each row contains components associated with
#  a single spatial degree of freedom, and each column contains 
#  components associated with a single monomer type. The format 
#  of header is different for different types of field file, which 
#  makes it possible for the Field class to infer information about
#  the type of a field file while parsing the header,
#  
#  A Field object has three attributes named 'header', 'data', and 
#  'type'. The 'header' attribute is a dictionary that contains the
#  contents of the header section, using the variable names or labels
#  as dictionary keys. The 'data' attribute is a list of lists in 
#  which each element of the outer list corresponds to one row of
#  the data section in the field file. When accessing elements of 
#  the data attribute using two subscripts, the first subscript (the 
#  row index) is thus index for a spatial degree of freedom (i.e., 
#  a grid point or basis function) and the second (the column index) 
#  is an index for a monomer type. 
#
#  The 'type' attribute is a string with possible values 'fd1d', 
#  'basis', 'rgrid' or 'kgrid'. 
#
#  **Construction:**
#
#  The Field constructor can parse any PSCF field file, and infer its 
#  type. To read and parse a field file with name 'field', one could 
#  enter:
#  \code
#     from pscfpp.field import *
#     f = Field('field')
#  \endcode
#
#  **Printing and writing to file:**
#
#  The __str__ function can be used to convert a Field object to a
#  multi-line string with a consistent with that of the file from which
#  it was created. The write function can be used to write that string
#  to a specified file. To write a field stored in a Field object named
#  f to a file named 'out', one would enter:
#  \code
#     f.write('out')
#  \endcode
#
#  **Accessing elements:**
#
#  The header attribute is a Python dictionary containing
#  the contents of the header section of the field file. One can
#  access elements using names of header variables as dictionary
#  keys, as in the expressions
#  \code
#     f.header['format']    
#     f.header['N_monomer']
#  \endcode
#  The available keys are different for different file formats.
#  Values for some quantities that are generally represented by more 
#  than one value are stored as lists, such as the list 'cell_param' 
#  of cell parameter values to describe the unit cell in  a periodic 
#  system.
#
#  The data section is a list of lists, in which elements are accessed
#  using square brackets with integer indices.  The expressions
#  \code
#    f.data[10]
#    f.data[10][1]
#  \endcode
#  access an entire row of the data section or a single component,
#  respectively.
# 
#  Elements of the header and data section can also be modified,
#  simply by putting corresponding expressions on the left side of 
#  an equality (assignment) operator. For example:
#  \code
#     f.header['dim'] = 3  
#     f.header['cell_param'] = [2.45]   
#     f.header['group_name'] = 'P_1'
#  \endcode
#
#  **Modifying data table structure:**
#
#  See the documentation for the addColumn, deleteColumn, and
#  reorder methods, which change the structure of the data 
#  section by adding, deleting or reordering columns associated
#  with different data types. 
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
   # Add a new column to the data list.
   #
   # This function adds a new column to the data list of the Field
   # object. The index parameter is the desired position of the new
   # column. The parameter 'element' can be either a single value
   # that is applied to every element of the entire column, or a
   # list of values that represents the entire column of new values.
   #
   # Example:
   # \code
   #    f.addColumn(2, 1.2)    
   #    f.addColumn(2, [list of values for the whole column])
   # \endcode
   #
   # \param index  integer column index 
   # \param element  value (float) or values (list) for elements of column
   #
   def addColumn(self, index, element):
      if isinstance(element, list):
         if not (len(element) == len(self.data)):
             raise Exception('Incorrect dimension for list element')
         for i in range(0, len(self.data)):
            self.data[i].insert(index, element[i])
      else:
         for i in range(0, len(self.data)):
            self.data[i].insert(index, element)

   ##
   # Delete an existing column from the data list.
   #
   # This function deletes an exist column from the data list of the
   # Field object. By passing in one parameter, index, the desired
   # values with that column index is removed from every row of the
   # data table.
   #
   # \param index  column index, numbered from 0
   #
   def deleteColumn(self, index):
      for i in range(0, len(self.data)):
         self.data[i].pop(index)

   ##
   # Reorder the data list.
   #
   # This function reorders the data list of the Field object. By
   # passing in one parameter, order, an integer list that represents
   # the new order of the columns in the data list, the data list of
   # the Field object can be updated with the desired order. The 
   # length of the list 'order' must be equal to the total number 
   # of columns in the data section, including any that do not 
   # contain field component values. 
   #
   # Note: This function treats all columns in the file format 
   # equivalently, whether they contain field component values or
   # other information. Specifically, when treating a field file 
   # in basis file format for a system with C monomer types in a 
   # system of spatial dimension D, the field file format contains 
   # C + D + 1 columns, in which only the first C contain field 
   # components and the remaining D+1 contain Miller indices for
   # a characteristic wave of each star and the number of waves in
   # the star. To re-order columns that contain field components
   # one must enter a permutation of the integers [0,...,D+1] in
   # which only the first C (0,...,C-1) are re-ordered and the
   # remaining D+1 integers appear in their original order.
   #
   # Example:
   #   \code
   #      f.reorder([1, 0, 2, 3, 4, 5])
   #   \endcode
   #
   # \param order  a list of position integers.
   #
   def reorder(self, order):
      nc = len(order)
      newData = []
      for i in range(0, len(self.data)):
         if not (nc == len(self.data[i])):
            raise Exception('Incorrect length for permutation order')
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

