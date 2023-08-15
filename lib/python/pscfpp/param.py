# -----------------------------------------------------------------------
#   This module provides tools to parse the PSCF parameter file and store 
#   all values within it in a single object. Users can access and modify 
#   the stored values of the properties after parsing by using specific 
#   statements (commands), and can write the object (or any of its 
#   sub-objects) to a file in the proper parameter file format. 
#
#   Parsing a parameter file:
#
#       A Composite object represents any block of the parameter file 
#       that is identified by opening and closing curly brackets ('{' and
#       '}'). The name of a Composite object is given by the identifier 
#       on the first line, preceding the opening curly bracket. A 
#       Composite object stores all the entries within such a block in a 
#       form that allows each element within it to be accessed and 
#       modified. 
#
#       A PSCF parameter file always contains a main parameter Composite 
#       named 'System' with nested subelements. Users may parse such a 
#       file by creating a Composite object, passing the name of 
#       parameter file as an argument. This constructor parses the file 
#       and returns a Composite object that contains its contents. 
#
#       Example: To read and parse a parameter file with name 'param',
#       execute the following code within a python3 interpreter:
#
#           from pscfpp.param import *
#           p = Composite('param')
#
#       Alternatively, a Composite object can be read from an open Python
#       file object. In this case, it is assumed that the first line of
#       the Composite has already been read (in order to identify the 
#       Composite by it's opening curly bracket '{'), and the name of
#       the composite must be passed in along with the open file object.
#       Note that the file object is not closed by this constructor.
#
#       Example:
#
#           f = open('param')
#           line = f.readline().strip()
#           if line[-1] == "{":
#              p = Composite(f,line[:-1])
#           f.close()
#
#       A Composite object can be written to a file by calling the 
#       write() method, passing in the string of the file name to 
#       which the Composite should be written.
#
#           Example: p.write('paramOut')
#
#   Accessing elements:
#
#       After creating a Composite object, users can retrieve 
#       the values of any element by name, using a dot notation for 
#       subelements of a Composite.  There are four different types of 
#       objects that are stored within a Composite.  These are listed 
#       below, along with a summary of what is returned when they are 
#       accessed, and an example Python expression that would access 
#       this type of object in a typical parameter file:
#
#           1. Composite: if an element of Composite is another 
#           Composite (identified by curly braces '{' and '}') with
#           that is unique within its parent object, accessing this
#           accessing this entry will return the child Composite 
#           itself; use square bracket indexing starting with 0 to 
#           access blocks that have the same name.
#           Example: p.Mixture    or    
#                    p.Mixture.Polymer[1]
#
#           2. Parameter: if an element is a single parameter, it 
#           contains a label followed by one or more values on a single
#           line. Accessing a Parameter returns the value of the 
#           parameter, which can be a string, an integer, a float, or 
#           a Python list. The value of a parameter is stored as a list 
#           if the corresponding line in the parameter file contains
#           multiple values separated by spaces after the label. 
#           Example: p.Mixture.nMonomer    or     
#                    p.Domain.mesh
#
#           3. Array: if an element is an Array, it is identified by 
#           square brackets '[' and ']'. Accessing an Array returns a 
#           Python list of Parameters where each entry of the list 
#           represents one row or element of the Array; a specific 
#           element can be accessed by square bracket indexing.
#           Example: p.Mixture.monomers    or
#                    p.Mixture.monomers[0]
#
#           4. Matrix: if an element is a Matrix, it is identified by 
#           parentheses '(' and ')'. Accessing a Matrix returns a list 
#           of lists that represents a square, symmetric matrix; specific 
#           values within the Matrix can be accessed by two separate 
#           square bracket indices.
#           Example: p.Interaction.chi    or
#                    p.Interaction.chi[0][1]
#
#       The parser also allows users to modify the entries in different 
#       preset formats for particular types of objects with equal sign 
#       operator ('='), which are listed below: 
#
#           1. Parameter: a Parameter with a single value can be modified 
#           by Python arithmetic operators. A parameter that contains
#           multiple values on a single line is stored as a Python 
#           list, which can only be modified by reassigning a new 
#           Python list to the attribute.
#           Example: p.Mixture.Polymer[1].phi *= 2    or    
#                    p.Mixture.Polymer[0].phi = 0.8   or    
#                    p.Domain.mesh = [72, 72, 72]
#
#           2. Array: change the whole Array by using a Python list
#           Example: p.Mixture.monomers = [2.0, 2.0]
#
#           3. Matrix: two ways to modify
#           Example: 1. change the whole Matrix by using a list of lists 
#                       that represents the squared, symmetric Matrix: 
#                       p.Interaction.chi = [[0, 1.0], [1.0, 0]]
#                    2. change two values that are symmetric at the 
#                       same time:
#                       p.Interaction.chi = [0, 1, 1.0], where the 
#                       first and second value of the list are the 
#                       position of the values needs to be changed and 
#                       the last value is the new value assigned to the 
#                       corresponding positions
#
# Module Contents:
#
#   class Composite:
#           A Composite object can be identified by opening and closing 
#           curly brackets ('{' and '}'). A composite contains one or
#           more enclosed elements, each of which may be a Composite, 
#           a Parameter, and Array or a Matrix. A Composite can be read 
#           from a file by calling Composite(filename), and can also be a
#           subelement of a parent Composite.
#
#   class Parameter:
#           A Parameter object contains a parameter label and its 
#           value for the labeled parameter contained in a single line.
#
#   class Array:
#           An Array object represents a one-dimensional array of values.
#           An array appears in a parameter file format in multi-line 
#           with one element value per line, delimited by square brackets.
#           The first lines contains the array label followed immediately
#           by an opening bracket, and the last line contains a closing
#           bracket. Elements of the array appear between these lines in
#           order of increasing array index, starting from 0.
#
#   class Matrix:
#           A Matrix object represents a two-dimensional array or matrix
#           of values. A matrix appears in a parameter file in a multi-line
#           format in which each element appears on a separate line, 
#           delimited by opening and closing parentheses. The first 
#           line contains a label for a matrix immediately followed by 
#           an opening parenthesis, and the last contains a closing 
#           parenthesis.  In between, each line contains a row index, 
#           a column index and value of a single element of the matrix.
#
#   def getValue(v):
#           A function to distinguish the correct type of the value from the 
#           passed-in string, v, and return it with the correct type. 
#
# -----------------------------------------------------------------------


class Composite:
   '''
   Purpose:
      The class represents the Composite element of the parameter file
   Instance variables:
      label: the label of the Composite element 
      children: 
         the children items of the Composite in the parameter file, 
         default to be an empty dictionary
   Methods:
      __init__(self, file=None, label=None):
         constructor, with two arguments: 
            file:
               the variable represents a filename or a open file, default
               to be None
            label: the label of the Composite, default to be None
      read(self, openFile):
         method to read the open parameter file, openFile as the argument, 
         line by line and update the read items into the children variable; 
         reading stop when '}' is read
      addChild(self, child):
         method to add the single item, argument child, into the children 
         variable
      getChildren(self):
         return the children variable
      __getattr__(self, attr):
         return the value stored in children with the specific key, 
         argument attr
      __str__(self):
         return un-indented string representation in param file format
      getString(self, depth):
         return indented string in param file format, with prefix depth
      write(self, filename):
         write an un-indented param file string to specified file
      returnData(self):
         return the Composite object itself
      __setattr__(self, label, val):
         set the new value, argument val, to the specific child of the 
         Composite in the children dictionary, with the name of argument 
         label
   '''
   
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

   def addChild(self, child):
      label = child.label
      if label in self.children:
         self.children[label] = [self.children[label]]
         self.children[label].append(child)
      else:
         self.children[label] = child

   def getChildren(self):
      return self.children

   def __str__(self):
      return self.getString()

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

   def write(self, filename):
      with open(filename, 'w') as f:
         f.write(self.getString())

   def returnData(self):
      return self
      
   def __setattr__(self, label, val):
      if label in self.children:
         self.children[label].setValue(val)
      else:
         self.__dict__[label] = val

# End class Composite ---------------------------------------------------


class Parameter:
   '''
   Purpose: 
      An individual parameter element of a parameter composite
   Instance variables:
      label: the label of the individual Parameter element 
      val: the value of the individual Parameter element 
   Methods:
      __init__(self, label, val):
         constructor, with two arguments: 
            label: label string for this Parameter, required
            val: value of this Parameter, required
      setValue(self, val):
         Reset the parameter value to argument val
      __str__(self):
         return un-indented string representation for a parameter
      getString(self, depth):
         return the string representation of this parameter.
         Argument depth, the 
         string of spaces that represents the level of the Parameter 
         element
      returnData(self):
         return the Parameter value as an int, float, string, or list
   '''

   def __init__(self, label, val):
      self.label = label
      if type(val) is list:
         self.val = val
      else:
         self.val = getValue(val)

   def setValue(self, val):
      if type(val) is list:
         if val[0] is list:
            raise Exception('Not valid input for Parameter.')
         self.val = []
         for i in range(len(val)):
            self.val.append(getValue(str(val[i])))
      else:
         self.val = getValue(str(val))

   def __str__(self):
      return self.getString()

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

   def returnData(self):
      return self.val

# End class Parameter ----------------------------------------------------


class Array:
   '''
   Purpose:
      The class represents the Array element of the parameter file
   Instance variables:
      label: the label of the Array element
      value: the value of the Array element, default to be an empty list
   Methods:
      __init__(self, label, openFile, val=None):
         constructor, with three arguments: 
            label: the label of the Array, required
            openFile: the opened parameter file, required
            val: the value of the Array, default to be None
      read(self, openFile):
         method to read the open parameter file, openFile as the argument, 
         line by line and update the value variable according to the read 
         lines; reading stop when ']' is read
      __str__(self):
         return unindented string representation for an array 
      getString(self, depth):
         return the string for writing out with argument depth, the string 
         of spaces that represents the level of the Array element
      returnData(self):
         return the list of exact value of the Array object stored as the 
         Value object
      setValue(self, val):
         set new value to the val variable with argument val
   '''
   def __init__(self, label, openFile, val=None):
      self.label = label
      if val == None:
         self.val = []
      else:
         self.val = val
      if openFile != None:
         self.read(openFile)

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

   def __str__(self):
      return self.getString()

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

   def returnData(self):
      return self.val

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
# End class Array ------------------------------------------------------


class Matrix:
   '''
   Purpose:
      The class represents the Matrix element of the parameter file
   Instance variables:
      label: the label of the Matrix element
      value: the value of the Matrix element, default to be an empty list
   Methods:
      __init__(self, label, openFile):
         constructor, with two arguments: 
            label: the label of the Matrix, required
            openFile: the opened parameter file, required
      read(self, openFile):
         method to read the open parameter file, openFile as the argument, 
         line by line and update the value variable according to the read 
         lines; reading stop when ')' is read
      __str__(self):
         return an unindented string representation
      getString(self, depth):
         return an indented string with indentation prefix depth
      returnData(self):
         return the list of exact value of the Matrix object stored as 
         the Value object
      setValue(self, val):
         set new value to the val variable with argument val
   '''

   def __init__(self, label, openFile):
      self.label = label
      self.val = []
      self.read(openFile)

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

   def __str__(self):
      return self.getString()

   def getString(self, depth=''):
      s = ''
      s += depth + self.label + '(\n'
      for i in range(len(self.val)):
         for j in range(i+1):
            s += depth + f'{i:>24}{j:>5}   ' + f'{self.val[i][j]:.12e}\n'
      s += depth + ')\n'
      return s

   def returnData(self):
      return self.val

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

# End class Matrix ---------------------------------------------------

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