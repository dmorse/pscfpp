"""! Module for parsing param files. """

##
#  Container for data of a Composite in a param file.
# 
#  This module provides tools to parse the PSCF parameter file and store 
#  all values within it in a single object. Users can access and modify 
#  the stored values of the properties after parsing by using specific 
#  statements (commands), and can write the object (or any of its 
#  sub-objects) to a file in the proper parameter file format.
#
#  A Composite object represents any block of the parameter file 
#  that is identified by opening and closing curly brackets ('{' and
#  '}'). The name of a Composite object is given by the identifier 
#  on the first line, preceding the opening curly bracket. A 
#  Composite object stores all the entries within such a block in a 
#  form that allows each element within it to be accessed and 
#  modified. 
#
#  **Construction:**
#
#      A PSCF parameter file always contains a main parameter  
#      Composite named 'System' with nested subelements. Users may  
#      parse such a file by creating a Composite object, passing  
#      the name of parameter file as an argument. This constructor  
#      parses the file and returns a Composite object that contains 
#      its contents.
#     
#      Example:
#
#      To read and parse a parameter file with name 'param':
#        \code
#           from pscfpp.param import *
#           p = Composite('param')
#        \endcode
#
#  **Writing out:**
#
#      A Composite object can be written to a file by calling the
#      write() method, passing in the string of the file name to 
#      which the Composite should be written.
#
#      Example:
#        \code
#           p.write('paramOut')
#        \endcode
#
#  **Accessing elements:**
# 
#      After creating a Composite object, users can retrieve 
#      the values of any element by name, using a dot notation for 
#      subelements of a Composite.  There are four different types of 
#      objects that are stored within a Composite.  These are listed 
#      below, along with a summary of what is returned when they are 
#      accessed, and an example Python expression that would access 
#      this type of object in a typical parameter file:
#
#           1. Composite: if an element of Composite is another 
#           Composite (identified by curly braces '{' and '}') with
#           that is unique within its parent object, accessing this
#           accessing this entry will return the child Composite 
#           itself; use square bracket indexing starting with 0 to 
#           access blocks that have the same name.
#
#           Example: 
#             \code
#                p.Mixture    
#                p.Mixture.Polymer[1]
#             \endcode
#
#           2. Parameter: if an element is a single parameter, it 
#           contains a label followed by one or more values on a single
#           line. Accessing a Parameter returns the value of the 
#           parameter, which can be a string, an integer, a float, or 
#           a Python list. The value of a parameter is stored as a list 
#           if the corresponding line in the parameter file contains
#           multiple values separated by spaces after the label. 
#
#           Example: 
#             \code
#                p.Mixture.nMonomer    
#                p.Domain.mesh
#             \endcode
#
#           3. Array: if an element is an Array, it is identified by 
#           square brackets '[' and ']'. Accessing an Array returns a 
#           Python list of Parameters where each entry of the list 
#           represents one row or element of the Array; a specific 
#           element can be accessed by square bracket indexing.
#
#           Example: 
#             \code
#                p.Mixture.monomers    
#                p.Mixture.monomers[0]
#             \endcode
#
#           4. Matrix: if an element is a Matrix, it is identified by 
#           parentheses '(' and ')'. Accessing a Matrix returns a list 
#           of lists that represents a square, symmetric matrix; specific 
#           values within the Matrix can be accessed by two separate 
#           square bracket indices.
#
#           Example: 
#             \code
#                p.Interaction.chi    
#                p.Interaction.chi[0][1]
#             \endcode
#
#  **Modifying elements:**
# 
#    The parser also allows users to modify the entries in different 
#    preset formats for particular types of objects with equal sign 
#    operator ('='), which are listed below:
#
#           1. Parameter: a Parameter with a single value can be modified 
#           by Python arithmetic operators. A parameter that contains
#           multiple values on a single line is stored as a Python 
#           list, which can only be modified by reassigning a new 
#           Python list to the attribute.
#
#           Example: 
#             \code
#                p.Mixture.Polymer[1].phi *= 2    
#                p.Mixture.Polymer[0].phi = 0.8    
#                p.Domain.mesh = [72, 72, 72]
#             \endcode
#
#           2. Array: change the whole Array by using a Python list.
#
#           Example: 
#             \code
#                p.Mixture.monomers = [2.0, 2.0]
#             \endcode
#
#           3. Matrix: two ways to modify:
#
#           Example: 
#
#             1. change the whole Matrix by using a list of lists that 
#                represents the squared, symmetric Matrix: 
#                \code
#                   p.Interaction.chi = [[0, 1.0], [1.0, 0]]
#                \endcode
#
#             2. change two values that are symmetric at the same time:
#               \code
#                  p.Interaction.chi = [0, 1, 1.0] 
#               \endcode
#                where the first and second value of the list are the 
#                position of the values needs to be changed and the last 
#                value is the new value assigned to the corresponding 
#                positions.
#
class Composite:

   ##
   # Constructor.
   #
   # The input parameter file can be either a string that
   # is the name of the file or a python open file object,
   # and it is defult to be none. The other input parameter
   # label is a string that is the label of the Composite
   # object and it is also defult to be None. If necessary, 
   # a defult Composite object can be created by without
   # passing in neither file nor label parameters.
   #
   # \param file   a filename string or an open-file object.
   # \param label  label string for the Composite, optional.
   #
   def __init__(self, file=None, label=None):
      self.label = label
      self.children = {} 
      
      if file != None:
         if type(file) == str:
            with open(file) as f:
               firstLine = f.readline()
               fl = firstLine.split()
               if fl[0][-1] != '{':
                  raise Exception('This is not a valid parameter file.')
               else:
                  self.label = fl[0][:-1]
                  self.read(f)
         else:
            self.read(file)

   ##
   # Read the passed-in open-file.
   #
   # This function reads the pass in open-file object line
   # by line and update the read items into the children 
   # variable. The reading stops when "}" is read. 
   #
   # \param file  an open-file object.
   #
   def read(self, openFile):
      line = openFile.readline()
      l = line.split()
      while line != '':
         if l[0][-1] == '{':
            if len(l) == 1:
               p = Composite(openFile, l[0][:-1])
            else:
               raise Exception('Not valid syntax for Composite element.')
         elif l[0][-1] == '[':
            if len(l) == 1:
               p = Array(l[0][:-1], openFile, None)
            else:
               val = []
               if l[-1] == ']':
                  for i in range(1, len(l)-1):
                     val.append(getValue(l[i]))
                  p = Array(l[0][:-1], None, val)
               else:
                  for i in range(1,len(l)):
                     val.append(getValue(l[i]))
                  p = Array(l[0][:-1], openFile, val)   
         elif l[0][-1] == '(':
            if len(l) == 1:
               p = Matrix(l[0][:-1], openFile)
            else:
               raise Exception('Not valid syntax for Matrix element.')
         elif l[0] == '}':
            break
         else:
            if len(l) == 2:
               p = Parameter(l[0], l[1])
            else:
               val = []
               for i in range(1, len(l)):
                  val.append(getValue(l[i]))
               p = Parameter(l[0], val)
         self.addChild(p)   
         line = openFile.readline()
         l = line.split()

   ##
   # Add the single item into the children variable.
   #
   # This function adds the pass in child parameter
   # into the children variable of the Composite object.
   # 
   # \param child single item needs to be added
   #
   def addChild(self, child):
      label = child.label
      if label in self.children:
         self.children[label] = [self.children[label]]
         self.children[label].append(child)
      else:
         self.children[label] = child

   ##
   # Get the children variable of the Composite.
   #
   # This function return the children variable of the
   # Composite object, as a Python dictionary.
   #
   # Return value:
   #
   # The children variable of the Composite object which
   # is a Python dictionary.
   #
   def getChildren(self):
      return self.children

   ##
   # Return the un-indented string of the Composite object.
   #
   # This function return the un-indented string
   # representation of the Composite object in the param
   # file format.
   #
   # Return value:
   #
   # The un-indented string representation in the param
   # file.
   #
   def __str__(self):
      return self.getString()

   ##
   # Return the value stored in the Composite object.
   #
   # The function return the value stored in the
   # children variable of the Composite object with the 
   # specific key, attr.
   #
   # Return value:
   #
   # The value stored in the children variable of the 
   # Composite with the specific key.
   #
   # \param attr the string to specify the value returned.
   #
   def __getattr__(self, attr):
      if attr =='children':
         return {}
      if attr in self.children:
         if type(self.children[attr]) is list:
            return self.children[attr]
         else:
            return self.children[attr].returnData()
      else:
         return self.attr

   ##
   # Indented string of the Composite object.
   #
   # This function return the intended string
   # representation of the Composite object in the param
   # file format, with prefix depth as the passed-in 
   # parameter.
   #
   # Return value:
   #
   # The intended string representation of the 
   # Composite object in the param file format with
   # prefix depth.
   #
   # \param depth the spaces string as the prefix depth.
   #
   def getString(self, depth=''):
      s = depth + self.label + '{' + '\n'
      for item in self.children.values():
         if type(item) is list:
            for i in range(len(item)):
               s += item[i].getString(depth+'  ')
         else:
            s += item.getString(depth+'  ')
      s += depth + '}\n'
      return s

   ##
   # Write out an un-intended param file string to a file.
   #
   # This function writes out the un-intended param
   # string to the specified file with the name of the
   # passed-in parameter, filename.
   #
   # \param filename  a filename string.
   #
   def write(self, filename):
      with open(filename, 'w') as f:
         f.write(self.getString())

   ##
   # Return the Coposite object itself.
   #
   # This function returns the Composite object itself.
   #
   def returnData(self):
      return self
      
   ##
   # Set new value to the specific child.
   #
   # This function sets new value, parameter val, to the
   # specified child stored in the children variable of
   # the Composite object with the name of the parameter
   # label.
   #
   # \param label  a string of the child name.
   # \param val  the expected new value.
   #
   def __setattr__(self, label, val):
      if label in self.children:
         self.children[label].setValue(val)
      else:
         self.__dict__[label] = val

##
# Container for data of a parameter in a param file.
#
#  A Parameter object contains a parameter label and its 
#  value for the labeled parameter contained in a single line.
#
class Parameter:
   
   ##
   # Constructor.
   #
   # \param label  label string for the individual parameter
   # \param val  stored individual parameter value
   #
   def __init__(self, label, val):
      self.label = label
      if type(val) is list:
         self.val = val
      else:
         self.val = getValue(val)

   ##
   # Set the new value to the Parameter object.
   #
   # This function sets the new value, parameter val, to the
   # Parameter object.
   #
   # \param val  the expected new value.
   #
   def setValue(self, val):
      if type(val) is list:
         if val[0] is list:
            raise Exception('Not valid input for Parameter.')
         self.val = []
         for i in range(len(val)):
            self.val.append(getValue(str(val[i])))
      else:
         self.val = getValue(str(val))

   ##
   # Un-indented string of the Parameter object.
   #
   # This function returns the un-intended string 
   # representation for the Parameter object.
   #
   # Return value:
   #
   # The un-intended string representation for the
   # Parameter object.
   # 
   def __str__(self):
      return self.getString()

   ##
   # Indented string for the Parameter object.
   #
   # This function returns the intended string
   # representation for the Parameter object in the param
   # file format, with prefix depth as the passed-in 
   # parameter.
   #
   # Return value:
   #
   # The intended string representation of the 
   # Parameter object in the param file format with
   # prefix depth.
   #
   # \param depth the spaces string as the prefix depth.
   #
   def getString(self, depth=''):
      s = ''
      if type(self.val) is list:
         s += depth + f'{self.label:40}'
         s += f'{self.val[0]:>6}'
         for i in range(1, len(self.val)):
            s += f'{self.val[i]:>7}'
         s += '\n'
      else:
         s += depth + f'{self.label:20}'
         if type(self.val) is float:
            v = f'{self.val:.12e}'
            s += f'{v:>20}\n'
         else:
            s += f'{self.val:>20}\n'
      return s

   ##
   # Return the stored value.
   #
   # This function returns the value stored in the
   # Parameter object.
   #
   # Return value:
   #
   # The value stored in the Parameter object.
   #
   def returnData(self):
      return self.val

##
# Container for data of an array in a param file.
#
#  An Array object represents a one-dimensional array of values.
#  An array appears in a parameter file format in multi-line 
#  with one element value per line, delimited by square brackets.
#  The first lines contains the array label followed immediately
#  by an opening bracket, and the last line contains a closing
#  bracket. Elements of the array appear between these lines in
#  order of increasing array index, starting from 0.
#
class Array:

   ##
   # Constructor.
   #
   # \param label  label string for the Array.
   # \param openFile  an open-file object.
   # \param val  stored value for the Array, defult to be None.
   #
   def __init__(self, label, openFile, val=None):
      self.label = label
      if val == None:
         self.val = []
      else:
         self.val = val
      if openFile != None:
         self.read(openFile)

   ## 
   # Read the passed-in open-file.
   #
   # This function reads the passed-in open-file object line
   # by line and update the read items into the val variable. 
   # The reading stops when "]" is read. 
   #
   # \param file  an open-file object.
   #
   def read(self, openFile):
      line = openFile.readline()
      l = line.split()
      while l[0] != ']':
         if len(l) == 1:
            self.val.append(getValue(l[0]))
         else:
            ls = []
            for i in range(len(l)):
               ls.append(getValue(l[i]))
            self.val.append(ls)

         line = openFile.readline()
         l = line.split()

   ##
   # Un-indented string of the Array object.
   #
   # This function return the un-indented string
   # representation for the Array object in the param
   # file format.
   #
   # Return value:
   #
   # The un-indented string representation in the param
   # file.
   #
   def __str__(self):
      return self.getString()

   ##
   # Indented string for the Array object.
   #
   # This function returns the intended string
   # representation for the Array object in the param
   # file format, with prefix depth as the passed-in 
   # parameter.
   #
   # Return value:
   #
   # The intended string representation of the 
   # Array object in the param file format with
   # prefix depth.
   #
   # \param depth the spaces string as the prefix depth.
   #
   def getString(self, depth=''):
      s = ''
      s += depth + self.label + '[' + '\n'
      if type(self.val[0]) != list:
         for i in range(len(self.val)):
            v = f'{self.val[i]:.12e}'
            s += depth + f'{v:>40}\n'
      else:
         if (type(self.val[0][0]) == int) & (len(self.val[0]) == 2):
            for i in range(len(self.val)):
               v = f'{self.val[i][1]:.12e}'
               s += depth + f'{self.val[i][0]:>41}{v:>22}\n'
         else:
            for i in range(len(self.val)):
               s += depth + f'{self.val[i][0]:^20}'
               for j in range(1, len(self.val[0])):
                  if j == (len(self.val[0])-1):
                     if self.val[i][j] < 0:
                        v = f'{self.val[i][j]:.11e}'
                     else:
                        v = f'{self.val[i][j]:.12e}'
                     s += f'{v:>22}\n'
                  elif j == 1:
                     s += f'{self.val[i][j]}'
                  else:
                     s += f'{self.val[i][j]:>5}'
      s += depth + ']\n'
      return s

   ##
   # Return the stored value.
   #
   # This function returns the value stored in the
   # Array object.
   #
   # Return value:
   #
   # The value stored in the Array object.
   #
   def returnData(self):
      return self.val

   ##
   # Set the new value to the Array object.
   #
   # This function sets the new value, parameter val, to the
   # Array object.
   #
   # \param val  the expected new value.
   #
   def setValue(self, val):
      if (type(val) is list) == False:
         raise Exception('Not valid input for Array.')
      else:
         v = []
         if type(val[0]) is list:
            same = True
            for i in range(len(val)):
               if len(val[i]) != 1:
                  same = False
                  break
            if same == True:
               for i in range(len(val)):
                  v.append(getValue(str(val[i][0])))
            else:
               for i in range(len(val)):
                  v.append([])
                  for j in range(len(val[i])):
                     v[i].append(getValue(str(val[i][j])))
         else:
            for i in range(len(val)):
               v.append(getValue(str(val[i])))
         self.val = v

##
# Container for data of a matrix in a param file.
#
#  A Matrix object represents a two-dimensional array or matrix
#  of values. A matrix appears in a parameter file in a multi-line
#  format in which each element appears on a separate line, 
#  delimited by opening and closing parentheses. The first 
#  line contains a label for a matrix immediately followed by 
#  an opening parenthesis, and the last contains a closing 
#  parenthesis.  In between, each line contains a row index, 
#  a column index and value of a single element of the matrix.
#
class Matrix:
   ##
   # Constructor.
   #
   # \param label  label string for the Matrix.
   # \param openFile  the open-file object.
   #
   def __init__(self, label, openFile):
      self.label = label
      self.val = []
      self.read(openFile)

   ## 
   # Read the passed-in open-file.
   #
   # This function reads the passed-in open-file object line
   # by line and update the read items into the val variable. 
   # The reading stops when ")" is read. 
   #
   # \param file  an open-file object.
   #
   def read(self, openFile):
      att = []
      line = openFile.readline()
      l = line.split()
      while l[0] != ')':
         att.append(l)
         line = openFile.readline()
         l = line.split()
      rowMax = att[0][0]
      colMax = att[0][1]
      for i in range(1, len(att)):
         if att[i][0] > rowMax:
            rowMax = att[i][0]
         if att[i][1] > colMax:
            colMax = att[i][1]
      size = int(max(rowMax, colMax))+1
      for i in range(0, size):
         self.val.append([])
         for j in range(0, size):
            self.val[i].append(getValue('0'))
      for i in range(0, len(att)):
         self.val[int(att[i][0])][int(att[i][1])] = getValue(att[i][2])
         self.val[int(att[i][1])][int(att[i][0])] = getValue(att[i][2])

   ##
   # Un-indented string of the Matrix object.
   #
   # This function return the un-indented string
   # representation for the Matrix object in the param
   # file format.
   #
   # Return value:
   #
   # The un-indented string representation in the param
   # file.
   #
   def __str__(self):
      return self.getString()

   ##
   # Indented string for the Matrix object.
   #
   # This function returns the intended string
   # representation for the Matrix object in the param
   # file format, with prefix depth as the passed-in 
   # parameter.
   #
   # Return value:
   #
   # The intended string representation of the 
   # Array object in the param file format with
   # prefix depth.
   #
   # \param depth the spaces string as the prefix depth.
   #
   def getString(self, depth=''):
      s = ''
      s += depth + self.label + '(\n'
      for i in range(len(self.val)):
         for j in range(i+1):
            s += depth + f'{i:>24}{j:>5}   ' + f'{self.val[i][j]:.12e}\n'
      s += depth + ')\n'
      return s

   ##
   # Return the stored value.
   #
   # This function returns the value stored in the
   # Matrix object.
   #
   # Return value:
   #
   # The value stored in the Matrix object.
   #
   def returnData(self):
      return self.val

   ##
   # Set the new value to the Matrix object.
   #
   # This function sets the new value, parameter val, to the
   # Matrix object.
   #
   # \param val  the expected new value.
   #
   def setValue(self, val):
      if type(val) != list:
         raise TypeError('This is not a valid value for Matrix.')
      else:
         if type(val[0]) is list: # input value is a matrix
            if len(val) != len(val[0]):
               raise Exception('Input Matrix should be squared.')
            for i in range(len(val)):
               for j in range(i):
                  if val[i][j] != val[j][i]:
                     raise Exception('This is not a diagonal symmetric squared Matrix.')
            for i in range(len(val)):
               if val[i][i] != 0:
                  raise Exception('The diagonal of the Matrix should all be zero.')
            self.val = []
            for i in range(len(val)):
               self.val.append([])
               for j in range(len(val[0])):
                  self.val[i].append(getValue(str(val[i][j])))
         else: # input value is a list in format: [x-index, y-index, value]
            if len(val) != 3:
               raise Exception('Not valid input format for Matrix modification.')
            elif (type(val[0]) != int) or (type(val[1]) != int) or ((type(val[-1]) != int) and (type(val[-1]) != float)):
               raise Exception('Not valid input format for Matrix modification.')
            self.val[val[0]][val[1]] = getValue(str(val[-1]))
            self.val[val[1]][val[0]] = getValue(str(val[-1]))


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

