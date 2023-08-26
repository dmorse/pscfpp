"""! Module for parsing param files. """

##
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
   # This function is 
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

##
# Container for data of a parameter in a parameter file.
#
# 
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

   ##
   # Constructor.
   #
   # \param label  label string for the individual parameter
   # \param val  stored individual parameter value
   def __init__(self, label, val):
      self.label = label
      if type(val) is list:
         self.val = val
      else:
         self.val = getValue(val)

   ##
   # Reset the value of the individual parameter.
   #
   # \param val 
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
   # Unindented string representation.
   #
   # 
   def __str__(self):
      return self.getString()

   ##
   # 
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