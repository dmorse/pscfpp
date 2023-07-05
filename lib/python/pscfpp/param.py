# ------------------------------------------------------------------------
#   This module provides tools to parse the PSCF parameter file and store 
#   all values within it in a single object. Users can access and modify 
#   the stored values of the parameters after parsing by using specific 
#   statements (commands). 
#
#   Parsing a parameter file:
#
#       A Composite object represents any block of the parameter file that 
#       is identified by opening and closing curly brackets ('{' and '}').
#       The name of a Composite object is given by the identifier on the
#       first line, preceding the opening curly bracket.  A Composite 
#       python object stores all the entries within such a block in a 
#       form that allows each element within it to be accessed and 
#       modified. 
#
#       A PSCF parameter file always contains a main parameter composite 
#       named 'System' with nested sub-elements. Users may parse such a 
#       file by calling the readParamFile function, passing the name of 
#       parameter file as an argument. This function parses the file and
#       returns a Composite object that contains its contents. 
#
#       Example: To read and parse ae parameter file with name 'param',
#       execute the following code within a python3 interpreter:
#
#           from pscfpp.param import *
#           p = readParamFile('param')
#
#       A Composite object can be written to a file by calling the 
#       writeOut() method, passing in the string of the file name to 
#       which the Composite should be written.
#
#           Example: p.writeOut('paramOut')
#
#   Accessing elements:
#
#       After calling the readParamFile() function, users can retrieve 
#       the values of any element by name, using a dot notation for 
#       subelements of a Composite.  There are four different types of 
#       objects that are stored within a Composite.  These are listed 
#       below, along with a summary of what is returned when they are 
#       accessed, and an example python expression that would access 
#       this type of object in a typical parameter file:
#
#           1. Composite: If an element of Composite is another 
#           Composite (identified by curly braces '{' and '}') with
#           that is unique within its parent object, accessing this
#           accessing this entry will return the child Composite 
#           itself; use square bracket 
#           indexing starting with 0 to access blocks that have the 
#           same name.
#           Example: p.Mixture    or    
#                    p.Mixture.Polymer[1]
#
#           2. Parameter: If an element is a single parameter, it 
#           contains a label followed by one or more values on a single
#           line.  Accessing a Parameter returns the value of the 
#           parameter, which can be a string, an integer, a float, or 
#           a Python list. The value of a parameter is stored as a list 
#           if the corresponding line in the parameter file contains
#           multiple values separated by spaces after the label. 
#           Example: p.Mixture.nMonomer    or     
#                    p.Domain.mesh
#
#           3. Array: if an entry is an Array, it is identified by 
#           square brackets '[' and ']'. Accessing an Array returns a 
#           Python list of Parameters where each entry of the list 
#           represents one row or element of the Array; a specific 
#           element can be accessed by square bracket indexing.
#           Example: p.Mixture.monomers    or
#                    p.Mixture.monomers[0]
#
#           4. Matrix: if an entry is a Matrix, it is identified by 
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
#           1. Parameter: A Parameter with a single value can be modified 
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
#   class Element:
#           The Element class is the base class of all types of elements 
#           in the parameter file. It has four subclasses: Composite, 
#           Parameter, Array and Matrix.
#
#   class Composite(Element):
#           A Composite object can be identified by opening and closing 
#           curly brackets ('{' and '}'). A composite contains one or
#           more enclosed elements, each of which may be a Composite, 
#           a Parameter, and Array or a Matrix.
#
#   class Parameter(Element):
#           A Parameter object contains a parameter label and its 
#           value for the labelled parameter contained in a single line.
#
#   class Array(Element):
#           An Array object represents a one-dimensional array of values.
#           An array appears in a parameter file format in multi-line 
#           with one element value per line, delimited by square brackets.
#           The first lines contains the array label followed immediately
#           by an opening bracket, and the last line contains a closing
#           bracket. Elements of the array appear between these lines in
#           order of increasing array index, starting from 0.
#
#   class Matrix(Element):
#           A Matrix object represents a two-dimensional array or matrix
#           of values. A matrix appears in a parameter file in a multi-line
#           format in which each element appears on a separate line, 
#           delimited by opening and closing parentheses. The first 
#           line contains a label for a matrix immediately followed by 
#           an opening parenthesis, and the last contains a closing 
#           parenthesis.  In between, each line contains a row index, 
#           a column index and value of a single element of the matrix.
#
#   class Value:
#           A Value object stores the value of a single primitive 
#           variable (and integer, floating point number or string) that 
#           represents all or part of the value of a parameter. The type 
#           is inferred from the string format - strings that can be 
#           intepreted as integers are stored as integers, others that 
#           can represent floating point numbers are stored as floats,
#           and all others are stored as literal strings. 
#
#
#   function readParamFile(filename)
#           This function takes the name of a parameter file as an
#           argument and returns a ParamComposite object containing
#           the contents of the file.
#
# --------------------------------------------------------------------------

class Element:
   '''
   Purpose: 
      Base class of all types of elements of a parameter block.
      This base class has four concrete subclasses:
         Composite: 
            Parameter block delimited by curly brackets ('{' and '}')
         Parameter: 
            Single value parameter within a single line
         Array: 
            1D array of values delimited by square brackets ('[' and ']')
         Matrix:
            Matrix of values delimited by parentheses
   Instance variables:
      label: label string for the elemente
      parent: The parent of the element, default to None
   Methods:
      __init__(self, label, openFile=None, parent=None, val=None):
         Constructor, with four arguments: 
            label, label string for this Element
            openFile, file name of a parameter file (default to None)
            parent, the parent of the Element (default to None)
            val, the value of the Element (default to None)
      getLabel(self): 
         return the label of the Element
      returnData(self):
         pass for the base class
      writeOutString(self, depth=''):
         pass for the base class
   '''

   def __init__(self, label, openFile=None, parent=None, val=None):
      self.label = label
      self.parent = parent

   def getLabel(self):
      return self.label

   def returnData(self):
      pass

   def writeOutString(self, depth=''):
      pass

# End class Element -------------------------------------------------------


class Parameter(Element):
   '''
   Purpose: 
      An individual parameter element of a parameter composite
   Instance variables:
      label: the label of the individual Parameter element 
      parent: the parent of the current individual Parameter element 
      val: the value of the individual Parameter element 
   Methods:
      __init__(self, label, openFile=None, parent=None, value=None):
         Constructor, with four arguments: 
            label, label string for this Parameter (required)
            openFile, an open parameter file (default to None)
            parent, parent Element of this Parameter (default to None)
            val, value of this Parameter;
      getValue(self): 
         Return the parameter value and print its type
      setValue(self, val):
         Reset the parameter value to argument val
      __repr___(self):
         return the string that represents the stored value
      writeOutString(self, depth):
         return the string representation of this parameter.
         Argument depth, the 
         string of spaces that represents the level of the Parameter 
         element
      returnData(self):
         return the Parameter value as an int, float, string, or list
   '''

   def __init__(self, label, openFile=None, parent=None, val=None):
      super().__init__(label, openFile, parent, val)
      if type(val) is list:
         self.val = val
      else:
         self.val = Value(val)

   def getValue(self):
      print(self.val.getType())
      return self.val.getValue()

   def setValue(self, val):
      if type(val) is list:
         if val[0] is list:
            raise Exception('Not valid input for Parameter.')
         self.val = []
         for i in range(len(val)):
            self.val.append(Value(str(val[i])))
      else:
         self.val = Value(str(val))

   def __repr__(self):
      if type(self.val) is list:
         v = []
         for i in range(len(self.val)):
            v.append(self.val[i].getValue())
         return str(v)
      else:
         return str(self.val.getValue())

   def writeOutString(self, depth=''):
      s = ''
      if type(self.val) is list:
         s += depth + f'{self.label:40}'
         s += f'{self.val[0].getValue():>6}'
         for i in range(1, len(self.val)):
            s += f'{self.val[i].getValue():>7}'
         s += '\n'
      else:
         s += depth + f'{self.label:20}'
         if self.val.getType() is float:
            v = f'{self.val.getValue():.12e}'
            s += f'{v:>20}\n'
         else:
            s += f'{self.val.getValue():>20}\n'
      return s

   def returnData(self):
      if type(self.val) is list:
         v = []
         for i in range(len(self.val)):
            v.append(self.val[i].getValue())
         return v
      else:
         return self.val.getValue()

# End class Parameter ----------------------------------------------------


class Array(Element):
   '''
   Purpose:
      The class represents the Array element of the parameter file, which
      is a subclass of the Element class
   Instance variables:
      label: the label of the Array element
      parent: the parent of the current Array element
      value: the value of the Array element, defult to be an empty list
   Methods:
      __inti__(self, label, openFile=None, parent=None, val=None):
         the constructor of the Array object for initiation, with four 
         arguments: 
            label, the label of the Array; 
            openFile, the opened parameter file;
            parent, the parent of the current Array, defult to be None;
            val, the value of the Array, defult to be None;
         that initiate the label and the parent of the Array object, and 
         pass in the open file for reading
      read(self, openFile):
         method to read the open parameter file, openFile as the argement, line 
         by line and update the value variable according to the read lines;
         reading stop when ']' is read
      __repr__(self):
         return the string of the value variable, in the list format string 
      writeOutString(self, depth):
         return the string for writting out with argument depth, the string of 
         spaces that represents the level of the Array element
      returnData(self):
         return the list of exact value of the Array object stored as the 
         Value object
      setValue(self, val):
         set new value to the val variable with argement val
   '''
   def __init__(self, label, openFile=None, parent=None, val=None):
      super().__init__(label, openFile, parent, val)
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
            self.val.append(Value(l[0]))
         else:
            ls = []
            for i in range(len(l)):
               ls.append(Value(l[i]))
            self.val.append(ls)

         line = openFile.readline()
         l = line.split()

   def __repr__(self):
      v = []
      if type(self.val[0]) is list:
         for i in range(len(self.val)):
            v.append([])
            for j in range(len(self.val[0])):
               v[i].append(self.val[i][j].getValue())
      else:
         for i in range(len(self.val)):
            v.append(self.val[i].getValue())
      return str(v)

   def writeOutString(self, depth=''):
      s = ''
      s += depth + self.label + '[' + '\n'
      if type(self.val[0]) != list:
         for i in range(len(self.val)):
            v = f'{self.val[i].getValue():.12e}'
            s += depth + f'{v:>40}\n'
      else:
         if (self.val[0][0].getType() == int) & (len(self.val[0]) == 2):
            for i in range(len(self.val)):
               v = f'{self.val[i][1].getValue():.12e}'
               s += depth + f'{self.val[i][0].getValue():>41}{v:>22}\n'
         else:
            for i in range(len(self.val)):
               s += depth + f'{self.val[i][0].getValue():^20}'
               for j in range(1, len(self.val[0])):
                  if j == (len(self.val[0])-1):
                     if self.val[i][j].getValue() < 0:
                        v = f'{self.val[i][j].getValue():.11e}'
                     else:
                        v = f'{self.val[i][j].getValue():.12e}'
                     s += f'{v:>22}\n'
                  elif j == 1:
                     s += f'{self.val[i][j].getValue()}'
                  else:
                     s += f'{self.val[i][j].getValue():>5}'
      s += depth + ']\n'
      return s

   def returnData(self):
      v = []
      if type(self.val[0]) is list:
         for i in range(len(self.val)):
            v.append([])
            for j in range(len(self.val[0])):
               v[i].append(self.val[i][j].getValue())
      else:
         for i in range(len(self.val)):
            v.append(self.val[i].getValue())
      return v

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
                  v.append(Value(str(val[i][0])))
            else:
               for i in range(len(val)):
                  v.append([])
                  for j in range(len(val[i])):
                     v[i].append(Value(str(val[i][j])))
         else:
            for i in range(len(val)):
               v.append(Value(str(val[i])))
         self.val = v
# End class Array ------------------------------------------------------


class Matrix(Element):
   '''
   Purpose:
      The class represents the Matrix element of the parameter file, which
      is a subclass of the Element class
   Instance variables:
      label: the label of the Matrix element
      parent: the parent of the current Matrix element
      value: the value of the Matrix element, defult to be an empty list
   Methods:
      __inti__(self, label, openFile=None, parent=None, value=None):
         the constructor of the Matrix object for initiation, with four 
         arguments: 
            label, the label of the Matrix; 
            openFile, the opened parameter file;
            parent, the parent of the current Matrix, defult to be None;
            value, the value of the Matrix, defult to be None;
         that initiate the label and the parent of the Matrix object, and 
         pass in the open file for reading
      read(self, openFile):
         method to read the open parameter file, openFile as the argement, line 
         by line and update the value variable according to the read lines;
         reading stop when ')' is read
      __repr__(self):
         return the string of the value variable, in the list format string 
      writeOutString(self, depth):
         return the string for writting out with argument depth, the string of 
         spaces that represents the level of the Matrix element
      returnData(self):
         return the list of exact value of the Matrix object stored as the 
         Value object
      setValue(self, val):
         set new value to the val variable with argement val
   '''

   def __init__(self, label, openFile=None, parent=None, val=None):
      super().__init__(label, openFile, parent, val)
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
            self.val[i].append(Value('0'))
      for i in range(0, len(att)):
         self.val[int(att[i][0])][int(att[i][1])] = Value(att[i][2])
         self.val[int(att[i][1])][int(att[i][0])] = Value(att[i][2])

   def __repr__(self):
      v = []
      for i in range(len(self.val)):
         v.append([])
         for j in range(len(self.val[0])):
            v[i].append(self.val[i][j].getValue())
      return str(v)

   def writeOutString(self, depth=''):
      s = ''
      s += depth + self.label + '(\n'
      for i in range(len(self.val)):
         for j in range(i+1):
            s += depth + f'{i:>24}{j:>5}   ' + f'{self.val[i][j].getValue():.12e}\n'
      s += depth + ')\n'
      return s

   def returnData(self):
      v = []
      for i in range(len(self.val)):
         v.append([])
         for j in range(len(self.val[0])):
            v[i].append(self.val[i][j].getValue())
      return v

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
                  self.val[i].append(Value(str(val[i][j])))
         else: # input value is a list in format: [x-index, y-index, value]
            if len(val) != 3:
               raise Exception('Not valid input format for Matrix modification.')
            elif (type(val[0]) != int) or (type(val[1]) != int) or ((type(val[-1]) != int) and (type(val[-1]) != float)):
               raise Exception('Not valid input format for Matrix modification.')
            self.val[val[0]][val[1]] = Value(str(val[-1]))
            self.val[val[1]][val[0]] = Value(str(val[-1]))

# End class Matrix ---------------------------------------------------


class Composite(Element):
   '''
   Purpose:
      The class represents the Copmosite element of the parameter file, which
      is a subclass of the Element class
   Instance variables:
      label: the label of the Composite element 
      parent: the parent of the current Composite element 
      children: 
         the children items of the Composite in the parameter file, defult 
         to be an empty dictionary
   Methods:
      __init__(self, label, openFile=None, parent=None, val=None):
         the constructor of the Composite object for initiation, with four 
         arguments: 
            label, the label of the Composite; 
            openFile, the open parameter file;
            parent, the parent of the current Composite, defult to be None;
            val, the value of the Composite, defult to be None;
         that initiate the label and the parent of the Composite object, and 
         pass in the open file for reading
      read(self, openFile):
         method to read the open parameter file, openFile as the argement, line 
         by line and update the read items into the children variable; reading 
         stop when '}' is read
      addChild(self, child):
         method to add the single item, argument child, into the children variable
      getChildren(self):
         return the children variable
      __repr__(self):
         return the string of children variable, in the dictionary format string
      __getattr__(self, attr):
         return the value stored in children with the specific key, argument attr
      writeOut(self, filename):
         method to write out the Composite element to a specific txt file with name
         of the argument filename
      writeOutString(self, depth):
         return the string for writting out with argument depth, the string of 
         spaces that represents the level of the Composite element
      returnData(self):
         return the Composite object itself
      __setattr__(self, label, val):
         set the new value, argument val, to the specific child of the Composite in
         the children dictionary, with the name of argument label
   '''
   def __init__(self, label, openFile, parent=None, val=None):
      super().__init__(label, openFile, parent, val)
      self.children = {} 
      self.read(openFile)

   def read(self, openFile):
      line = openFile.readline()
      l = line.split()
      while line != '':
         if l[0][-1] == '{':
            if len(l) == 1:
               p = Composite(l[0][:-1], openFile, self, None)
            else:
               raise Exception('Not valid syntax for Composite element.')
         elif l[0][-1] == '[':
            if len(l) == 1:
               p = Array(l[0][:-1], openFile, self, None)
            else:
               val = []
               if l[-1] == ']':
                  for i in range(1, len(l)-1):
                     val.append(Value(l[i]))
                  p = Array(l[0][:-1], None, self, val)
               else:
                  for i in range(1,len(l)):
                     val.append(Value(l[i]))
                  p = Array(l[0][:-1], openFile, self, val)   
         elif l[0][-1] == '(':
            if len(l) == 1:
               p = Matrix(l[0][:-1], openFile, self, None)
            else:
               raise Exception('Not valid syntax for Matrix element.')
         elif l[0] == '}':
            break
         else:
            if len(l) == 2:
               p = Parameter(l[0], None, self, l[1])
            else:
               val = []
               for i in range(1, len(l)):
                  val.append(Value(l[i]))
               p = Parameter(l[0], None, self, val)
         self.addChild(p)   
         line = openFile.readline()
         l = line.split()

   def addChild(self, child):
      label = child.getLabel()
      if label in self.children:
         self.children[label] = [self.children[label]]
         self.children[label].append(child)
      else:
         self.children[label] = child

   def getChildren(self):
      return self.children

   def __repr__(self):
      return self.writeOutString()

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

   def writeOut(self, filename):
      with open(filename, 'w') as f:
         f.write(self.writeOutString())

   def writeOutString(self, depth=''):
      s = depth + self.label + '{' + '\n'
      for item in self.children.values():
         if type(item) is list:
            for i in range(len(item)):
               s += item[i].writeOutString(depth+'  ')
         else:
            s += item.writeOutString(depth+'  ')
      s += depth + '}\n'
      return s

   def returnData(self):
      return self
      
   def __setattr__(self, label, val):
      if label in self.children:
         self.children[label].setValue(val)
      else:
         self.__dict__[label] = val

# End class Composite ---------------------------------------------------

class Value:
   '''
   Purpose: 
      A Value represents a primitive parameter value (int, float, or string)
   Instance variables:
      val: variable to store a primitive parameter value
      type: the type of the value, either integer, floating point or string
   Methods:
      __init__(self, val):
         Constructor, with argument val that is the string representation 
         of the value. The constructor automatically determines the type.
      getType(self): return the type of stored value
      getValue(self): return the value of stored value
   Note:
      Debugging by the commented line to check if the constructor has the 
      expected function of distinguishing the exact type of the value
   '''

   def __init__(self, val):
      if val.isdigit() == True:
         # print('int')
         self.val = int(val)
      else:
         try:
            self.val = float(val)
            # print('float')
         except ValueError:
            self.val = val
            # print('string')
      # print(self.value)
      self.type = type(self.val)

   def getType(self):
      return self.type

   def getValue(self):
      return self.val

# End class Value  -------------------------------------------------------


def readParamFile(filename):
   '''
   Argument:
      filename: string, name of the parameter file to be read
   Return:
      the 'System' ParamComposite 
   Note:
      If the first line does not begin with 'System{', this function
      will print the error message 'This is not a valid parameter file'
      and return None.
   '''
   with open(filename) as f:
      firstLine = f.readline()
      # print(firstLine, end='')
      if firstLine != "System{"+'\n':
         p = None
         raise Exception('This is not a valid parameter file.')
      else:
         # line = f.readline()
         # l = line.split()
         # p = Composite(l[0][:-1], f, None, None)
         p = Composite('System', f, None, None)
   return p

# End function readFileName --------------------------------------------
