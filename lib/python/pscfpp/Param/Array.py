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
