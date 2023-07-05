# Composite class ---------------------------------------------------------

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

# readParamFile(filename) Method ----------------------------------------

def readParamFile(filename):
   '''
   Argument:
      filename: string, name of the specific parameter file needed to be read
   Return:
      the 'System' composite of the read parameter file
   Note:
      'This is not a valid parameter file.' will be print if the read file were
      not in the correct format and return None
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

# End Method -----------------------------------------------------------
