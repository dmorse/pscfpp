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
