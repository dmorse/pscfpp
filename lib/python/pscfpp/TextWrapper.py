class TextWrapper:

   def __init__(self, eol = '\n', nIndent = 0):
      self.eol = eol
      self.nIndent = nIndent
      self.column = 0
      self.text = ''
      self.limit = 78

   def clear(self):
      self.column = 0
      self.text = ''

   def append(self, string):
      size = len(string)
      if (self.column + size > self.limit):
         self.text += self.eol
         for i in range(self.nIndent):
            self.text += ' '
         self.column = self.nIndent
      self.text += string
      self.column += size

   def __str__(self):
      return self.text

   def __repr__(self):
      return self.text
