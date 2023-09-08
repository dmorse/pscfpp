"""! Module for parsing param files. """

##
#  Container for data of a Composite in a param file.
#
#  The Composite class provides tools to parse and store the contents of a
#  a parameter file block that appears in a parameter file delimited by 
#  curly brackets ('{' and '}').  The class provides special functions
#  that allow users to access and modify stored values of parameters and
#  subblocks within that block using the dot syntax that is normally used 
#  for python object attributes, and to obtain a multi-line string 
#  representation of the associated block in parameter file format.
#
#  **Construction:**
#
#  A PSCF parameter file always consists of a main block with a label
#  'System' that contains nested sub-elements or children. The constructor
#  function for the Composite class accepts the name of such a parameter 
#  file as a function parameter, opens and parses the specified file and 
#  returns a Composite object that  contains the contents of that file.
#
#  Example: To read and parse a parameter file with the file named 
#  'param', one could enter
#  \code
#     from pscfpp.param import *
#     p = Composite('param')
#  \endcode
#  The variable p is then a Composite object that contains the contents
#  of file 'param'. 
#
#  **String representation:**
# 
#  The getString and __str__ methods of a Composite return the string
#  representation of a Composite object in the multi-line, indented 
#  format of a PSCF parameter file. The write method allows a user to
#  write the resulting string to a specified file.
#
#  Example: To write the contents to a file named 'paramOut', enter
#  \code
#     p.write('paramOut')
#  \endcode
#
#  **Accessing elements:**
#
#  After creating a Composite object, users can retrieve the values of
#  any child sub-element by name, using a dot notation for children of
#  a Composite.  Each child of a Composite may be an instance of the
#  class Composite (to represent a subblock), an instance of class
#  Parameter (to represent a single parameter), an Array (for an 
#  array-valued parameter), or a Matrix (a matrix-valued parameter). 
#  Syntax for accessing different types of child elements, and the
#  nature of the return values, are discussed below:
#
#  *Composite:* A nested sub-block that is delimited by curly brackets
#  is represented by a child Composite object. The label of such a
#  subblock is given by the identifier string that precedes the opening
#  curly bracket on the first line.  If this label is unique within its
#  parent object, then child Composite may be accessed by the python
#  dot notation for object attributes, using the unique sub-block label
#  as an attribute name. If there are two or more subblocks of the same
#  parent block with the same label, then using this shared name as an
#  attribute name after a dot returns a list in which each element is
#  a Composite representing one of the blocks with this label, indexed
#  from 0 in the order in which they appear in the parameter file block.
#
#  Example:  If p is an instance of Composite that contains the contents
#  of a PSCF parameter file for a mixture containing two or more polymer
#  species, then the expressions
#  \code
#     p.Mixture
#     p.Mixture.Polymer[1]
#  \endcode
#  return Composite objects that contain the contents of the entire
#  Mixture block and of the second Polymer subblock within the Mixture
#  block, respectively.
#
#  *Parameter:* A single parameter that appears in the parameter file
#  on a single line containing a label and one or more values is
#  represented internally by a child Parameter object. Using dot
#  notation to access a parameter, using the label as an attribute
#  name, returns the value of the parameter, which can be an int,
#  a float, a string, or a Python list.  The value of a parameter
#  is stored and returned as a list if and only if the corresponding 
#  line in the parameter file contains two or more values separated 
#  by spaces after the label.
#
#  Example:  If p is an initialized Composite that contains the
#  contents of a parameter file, then the expressions
#  \code
#     p.Mixture.nMonomer
#     p.Domain.mesh
#  \endcode
#  return the integer value of the number of monomer types (nMonomer)
#  and a list of integer mesh dimensions (mesh).
#
#  *Array:* An array-valued parameter within a parameter file
#  block is delimited by an opening line containing an array label
#  immediately followed by an opening square bracket ('['), and by
#  a closing line containing a closing square bracket (']'),
#  separated by lines that contain element values, with one element
#  value per line. Accessing an Array using dot notation, using
#  the name of the array as an attribute of the parent Composite,
#  returns a Python list in which each element of the list is
#  equal to the value of one element of the Array. Specific
#  elements can then be accessed by square bracket indexing.
#
#  Example: If p is the root Composite for a parameter file, the
#  expressions
#  \code
#     p.Mixture.monomers
#     p.Mixture.monomers[0]
#  \endcode
#  return the list of monomer statistical segment lengths, and the 
#  value for the first monomer type, respectively.
#
#  *Matrix:* A Matrix appears in a parameter file delimited by 
#  parentheses '(' and ')'. Accessing a Matrix returns a list of
#  lists that represents a square, symmetric matrix; specific values
#  within the Matrix can be accessed by two separate square bracket
#  indices.
#
#  Example: If p is the root Composite for a parameter file, then
#  the expressions
#  \code
#     p.Interaction.chi
#     p.Interaction.chi[0][1]
#  \endcode
#  return the entire chi matrix (a list of lists) and the value of
#  the element in row 0 and column 1, respectively. 
#
#  **Modifying elements:**
#
#  The parser also allows users to modify parameter entries that are 
#  accessed by dot notation, as discussed below:
#
#  *Parameter:* A Parameter with a single value can be modified by
#  Python assignment and arithmetic assignment operators. A parameter 
#  that contains multiple values on a single line, which is stored as 
#  a Python list, can can only be modified by assigning a new Python 
#  list to the attribute.
#
#  Example:
#  \code
#     p.Mixture.Polymer[1].phi *= 2
#     p.Mixture.Polymer[0].phi = 0.8
#     p.Domain.mesh = [72, 72, 72]
#  \endcode
#
#  *Array:* The value of an array-valued parameter is stored as a 
#  python list, and may only be modified by assigning a new list to
#  the variable.
#
#  Example: If p is a Composite representing the parameter file for
#  system with two monomer types, we may assign new values for the
#  statistical segment lengths of both monomers by an assignment such
#  as 
#  \code
#     p.Mixture.monomers = [2.0, 2.0]
#  \endcode
#
#  *Matrix:* There are two ways to modify elements of a Matrix
#
#   - Replace the whole Matrix by using a list of lists that
#     represents the squared, symmetric Matrix, like this:
#     \code
#       p.Interaction.chi = [[0, 2.0], [2.0, 0]]
#     \endcode
#
#   - Replace two values that are symmetric at the same time,
#     by passing a list containing a column index, a row index, and
#     a value. For example, the command
#     \code
#       p.Interaction.chi = [0, 1, 2.0]
#     \endcode
#     sets the (0,1) and (1,0) elements to 2.0. 
#
#  **Design notes:**
#
#  The following two attributes are created by the constructor:
#
#    - "children" is a dictionary of child sub-elements corresponding
#       to elements of the parameter file (i.e., parameters or nested
#       subblocks) that appear within the associated block. Values in
#       this dictionary are all instances of Composite, Parameter, 
#       Array or Matrix. 
#
#    - "label" is a string containing the label that precedes the
#      opening curly bracket
#
#  The __getattr__ function is overloaded to treat the attribute
#  name as a key in "children" dictionary and return either a
#  corresponding value. If the value in the dictionary is a 
#  Composite, this function returns the Composite object. If the
#  value is an instance of Parameter, Array or Matrix, it returns
#  the value of object.
#
class Composite:

   ##
   # Constructor.
   #
   # The constructor normally parses a parameter file and returns an
   # object containing the contents of that file. The parameter "file",
   # if present, can be either the name of a parameter file (a string)
   # or a python file object that is open for reading. The parameter
   # "label", if present, is a string that is the label of the Composite
   # object, which is equal to None by default.  Both parameters have
   # default values of None. 
   #
   # If the "file" parameter is absent, the constructor returns an 
   # object with no children. When an object is constructed in this way, 
   # the read() method may be called later in order to parse a parameter
   # file.
   #
   # \param file   a filename string or an file object
   # \param label  label string for the Composite
   #
   def __init__(self, file=None, label=None):
      self.label = label
      self.children = {}

      if file != None:
         if isinstance(file, str):
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
   # Read and parse a parameter block from a file.
   #
   # This function reads a parameter block and stores its contents in
   # this Composite object. The parameter "file" must be a python
   # file object that is open for reading. Processing of the file stops
   # when the closing bracket "}" that matches the opening bracket is
   # encountered.
   #
   # \param file  a python file object that is open for reading
   #
   def read(self, file):
      line = file.readline()
      l = line.split()
      while line != '':
         if l[0][-1] == '{':
            # Composite (nested)
            if len(l) == 1:
               label = l[0][:-1]
               p = Composite(file, label)
               #p = Composite(file, l[0][:-1])
            else:
               raise Exception('Invalid syntax for a Composite element')
         elif l[0][-1] == '[':
            # Array item
            if len(l) == 1:
               # Multi-line Array
               label = l[0][:-1]
               p = Array(label, file, None)
               #p = Array(l[0][:-1], file, None)
            else:
               # Single line Array
               val = []
               if l[-1] == ']':
                  for i in range(1, len(l)-1):
                     val.append(getValue(l[i]))
                  p = Array(l[0][:-1], None, val)
               else:
                  for i in range(1,len(l)):
                     val.append(getValue(l[i]))
                  p = Array(l[0][:-1], file, val)
         elif l[0][-1] == '(':
            # Matrix item
            if len(l) == 1:
               label = l[0][:-1]
               p = Matrix(label, file)
               #p = Matrix(l[0][:-1], file)
            else:
               raise Exception('Invalid syntax for Matrix parameter')
         elif l[0] == '}':
            # End of this Composite (closing bracket)
            break
         else:
            # Parameter item
            if len(l) == 1:
               raise Exception('Invalid parameter file line')
            elif len(l) == 2:
               # Parameter with a single value
               label = l[0]
               p = Parameter(label, l[1])
               #p = Parameter(l[0], l[1])
            else:
               # Parameter with multiple values on a single line
               label = l[0]
               val = []
               for i in range(1, len(l)):
                  val.append(getValue(l[i]))
               p = Parameter(label, val)
               #p = Parameter(l[0], val)
         self.addChild(p)
         line = file.readline()
         l = line.split()
      # end while loop

   ##
   # Add a single child to this Composite.
   #
   # This function adds the child object that is passed as a parameter
   # to the children dictionary of this Composite object. The child
   # object should be an instance of Composite, Parameter, Array, or 
   # Matrix. The dictionary key is set equal to the label of the child
   # which must thus have been assigned before entry to this function.
   #
   # \param child  object that should needs to be added
   #
   def addChild(self, child):
      label = child.label
      if label in self.children:
         current = self.children[label]
         if not isinstance(current, list):
            self.children[label] = [current]
         self.children[label].append(child)
      else:
         self.children[label] = child

   ##
   # Get the dictionary of child items.
   #
   # This function return the python dictionary that holds the children 
   # of this Composite object. Note: A separate function is required to
   # return this attribute because the __getattr__ function has been 
   # redefined so as to cause the dot operator to do something other 
   # than simply return an attribute by name.
   #
   def getChildren(self):
      return self.children

   ##
   # Get the label of this Composite.
   #
   # See documentation for getChildren for an explanation of why this
   # function is necessary.
   #
   def getLabel(self):
      return self.label

   ##
   # Return an indented string representation of this Composite.
   #
   # This function returns a string representation of this Composite in
   # an indented form of the param file format, with no indentation of
   # the first line (i.e., no prefix string of spaces).
   #
   def __str__(self):
      return self.getString()

   ##
   # Return the value of one child of this Composite.
   #
   # If the attribute parameter is the name (or dictionary key) of a
   # child of this this Composite, then this function returns the value
   # of the corresponding item. If the specified child is a single
   # Parameter (rather than an Array or Matrix) with only one value, 
   # it returns the value as a primitive python variable (an int, float 
   # or string).
   #
   # \param attr the string to specify the value returned.
   #
   def __getattr__(self, attr):
      if attr =='children':
         return {}
      if attr in self.children:
         if isinstance(self.children[attr], list):
            return self.children[attr]
         else:
            return self.children[attr].returnData()
      else:
         raise AttributeError     

   ##
   # Get an indented string representation of this Composite.
   #
   # This function return an indented multi-line string representation
   # of this Composite object in the PSCF parameter file format, with
   # neighboring levels of indentation offset by 3 spaces. The optional
   # parameter "depth" is a string of spaces that, which is empty by
   # default, that prepended to each line, and thus is the total 
   # indentation of the first line. 
   #
   # \param depth  string of spaces used as a prefix to all lines
   #
   def getString(self, depth=''):
      s = depth + self.label + '{' + '\n'
      for item in self.children.values():
         if isinstance(item, list):
            for i in range(len(item)):
               s += item[i].getString(depth+'  ')
         else:
            s += item.getString(depth+'  ')
      s += depth + '}\n'
      return s

   ##
   # Write the indented param file string for this to a file.
   #
   # This function writes out the indented multi-line string 
   # representation of this Composite to a specified file, with no 
   # indentation of the first line (i.e., depth == ''). The paramete 
   # "filename" is a the name of the file to which this should be 
   # written,
   #
   # \param filename  name filename string.
   #
   def write(self, filename):
      with open(filename, 'w') as f:
         f.write(self.getString())

   ##
   # Return the Composite object itself.
   #
   # This function returns the Composite object itself.
   #
   def returnData(self):
      return self

   ##
   # Set a new value for a specified child.
   #
   # This function sets new value, given by parameter val, to the
   # specified a child of this Composite with name given by the
   # parameter "label".
   #
   # \param label  name of a child (key string)
   # \param val  the new value to be assigned
   #
   def __setattr__(self, label, val):
      if label in self.children:
         self.children[label].setValue(val)
      else:
         super(Composite, self).__setattr__(label, val)

##
# A Parameter represents a single parameter in a parameter file. 
#
# A Parameter object contains the label and value of a single parameter.
# A single parameter is represented in a parameter file by a single line
# containing a label followed by one or more values separated by spaces.
# 
class Parameter:

   ##
   # Constructor.
   #
   # \param label  label string for the individual parameter
   # \param val  parameter value
   #
   def __init__(self, label, val):
      self.label = label
      if isinstance(val, list):
         self.val = val
      else:
         self.val = getValue(val)

   ##
   # Set the new value to the Parameter object.
   #
   # This function sets the new value, parameter val, for this Parameter
   # object.
   #
   # \param val  the expected new value.
   #
   def setValue(self, val):
      if isinstance(val, list):
         if val[0] is list:
            raise Exception('Not valid input for Parameter.')
         self.val = []
         for i in range(len(val)):
            self.val.append(getValue(str(val[i])))
      else:
         self.val = getValue(str(val))

   ##
   # Unindented string representation of this Parameter object.
   #
   # This function returns the unindented string representation for this
   # Parameter object, containing a label followed by the value.
   #
   def __str__(self):
      return self.getString()

   ##
   # Indented string for the Parameter object.
   #
   # This function returns the indented string
   # representation for the Parameter object in the param
   # file format, with prefix depth as the passed-in
   # parameter.
   #
   # \param depth the spaces string as the prefix depth.
   #
   def getString(self, depth=''):
      s = ''
      if isinstance(self.val, list):
         s += depth + f'{self.label:40}'
         s += f'{self.val[0]:>6}'
         for i in range(1, len(self.val)):
            s += f'{self.val[i]:>7}'
         s += '\n'
      else:
         s += depth + f'{self.label:20}'
         if isinstance(self.val, float):
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
#  An array appears in a parameter file in a multi-line format
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
   # \param label  label string for the Array
   # \param file  a file object, open for reading
   # \param val  stored value for the Array, default to None
   #
   def __init__(self, label, file, val=None):
      self.label = label
      if val == None:
         self.val = []
      else:
         self.val = val
      if file != None:
         self.read(file)

   ##
   # Read the elements of this Array from a file. 
   #
   # This function reads the elements of a multi-line array from the 
   # file object that is passed as a parameter. On entry, the opening
   # line containing the opening bracket must already have been 
   # processed, and the first unread line of the file should contain
   # the value of the first element.  Reading stops stops when the
   # closing "]" is read.
   #
   # \param file  a file object, open for reading
   #
   def read(self, openFile):
      line = openFile.readline()
      l = line.split()
      while l[0] != ']':
         # Process a line
         if len(l) == 1:
            # Single value
            self.val.append(getValue(l[0]))
         else:
            # Multiple space-separated values on one line
            ls = []
            for i in range(len(l)):
               ls.append(getValue(l[i]))
            self.val.append(ls)
         # Read and split the next line
         line = openFile.readline()
         l = line.split()

   ##
   # String representation of this Array object.
   #
   # This function return the string representation of this Array object
   # in the multi-line indented param file format, with no indentation 
   # of the opening line. 
   #
   def __str__(self):
      return self.getString()

   ##
   # Indented string representation of this Array object.
   #
   # This function returns the string representation for this Array
   # object in a multi-line indented parameter file format with 
   # control over indentation of the first line. The parameter 'depth'
   # is a string of spaces that is the indentation of the first line,
   # is prepended to every succeding line. 
   #
   # \param depth  string of spaces prefixed to all lines
   #
   def getString(self, depth=''):
      s = ''
      s += depth + self.label + '[' + '\n'
      if not isinstance(self.val[0], list):
         # Elements with a single value
         for i in range(len(self.val)):
            v = f'{self.val[i]:.12e}'
            s += depth + f'{v:>40}\n'
      else:
         # Elements with a list of values
         if isinstance(self.val[0][0], int) & (len(self.val[0]) == 2):
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
   # Return the value of this Array, as a list of element values.
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
      if not isinstance(val, list):
         raise Exception('Value of an Array must be a list')
      else:
         v = []
         if isinstance(val[0], list):
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
# A Matrix represents a matrix-valued parameter in parameter file.
#
# A Matrix object represents a two-dimensional array or matrix of
# element values. A matrix appears in a parameter file in a multi-line
# format that begins with a line containing a label immediately followed
# by an opening parenthesis and ends with a line containing only a 
# closing parenthesis. In between, each line contains a row index,
# a column index and the value of a single element of the matrix.
#
class Matrix:

   ##
   # Constructor.
   #
   # \param label  label string for the Matrix.
   # \param openFile  python file object, open for reading
   #
   def __init__(self, label, openFile):
      self.label = label
      self.val = []
      self.read(openFile)

   ##
   # Read the passed-in open-file.
   #
   # This function reads the elements of a Matrix from the file 
   # object that is passed as a parameter. On entry, the opening line
   # containing the opening curly bracket ("{") must already have been 
   # processed, and the first unread line of the file should contain
   # the value of the first element.  Reading stops stops when the
   # closing "}" is read.
   #
   # \param file a file object, open for reading
   #
   def read(self, file):
      att = []
      line = file.readline()
      l = line.split()
      while l[0] != ')':
         att.append(l)
         line = file.readline()
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
   # String representation of this Matrix object.
   #
   # This function return a multi-line string representation of this
   # Matrix object in the PSCF param file format, with no indentation 
   # of the first line. 
   #
   def __str__(self):
      return self.getString()

   ##
   # Indented string representation of this Matrix object.
   #
   # This function returns an indented string representation of this
   # Matrix object in PSCF parameter file format, with a controllable
   # indentation of the first line.
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
   # Return the value of this Matrix as a list of lists.
   #
   def returnData(self):
      return self.val

   ##
   # Set a new value for this Matrix.
   #
   # This function sets the new value, parameter val, to the Matrix
   # object.
   #
   # \param val  the new value
   #
   def setValue(self, val):
      if not isinstance(val, list):
         raise TypeError('This is not a valid value for Matrix.')
      else:
         if isinstance(val[0], list): 
            # Input value is a matrix
            if len(val) != len(val[0]):
               raise Exception('Input Matrix should be square')
            for i in range(len(val)):
               for j in range(i):
                  if val[i][j] != val[j][i]:
                     raise Exception('This is not a symmetric Matrix.')
            for i in range(len(val)):
               if val[i][i] != 0:
                  raise Exception('The diagonal of the Matrix should all be zero.')
            self.val = []
            for i in range(len(val)):
               self.val.append([])
               for j in range(len(val[0])):
                  self.val[i].append(getValue(str(val[i][j])))
         else: 
            # Input value is a list in format: [x-index, y-index, value]
            if len(val) != 3:
               raise Exception('Not valid input format for Matrix modification.')
            elif (type(val[0]) != int) or (type(val[1]) != int) or ((type(val[-1]) != int) and (type(val[-1]) != float)):
               raise Exception('Not valid input format for Matrix modification.')
            self.val[val[0]][val[1]] = getValue(str(val[-1]))
            self.val[val[1]][val[0]] = getValue(str(val[-1]))


##
# Distinguish the type of a value from its string representation.
#
# This function infers the type of a variable from its string 
# representation and returns a value of the appropriate type, which 
# may be an int, float or string.  Any input string that can be
# interpreted as an integer is converted into the resulting int value.
# Any string that can be interpreted as a real number but not as an 
# integer is converted into a float value. All other strings are left 
# un-interpreted and returned verbatim as string values.
#
# \param v  a string representation of a parameter value
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

