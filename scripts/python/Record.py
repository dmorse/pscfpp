#---------------------------------------------------------------------------
# A Record represents a string containing fields separated by whitespace.
#
# The constructor takes a line of text and divides into fields and spaces.
# The attribute fields is a list of field strings. The spaces attribute is
# a list containing strings of blank spaces that precede each space. The
# element space[i] is the string of blank characters that precedes field[i],
# for i >= 0. The space[0] may contain zero characters, if the line has no
# leading blank spaces, but all other elements must contain one or more 
# blank characters. Trailing white space is disregarded.
#---------------------------------------------------------------------------

class Record:
  
   def __init__(self, line):
      if line[-1] == '\n':
         line = line[:-1]
      self.line   = line
      self.spaces  = []
      self.fields  = []
      n = 0
      blank  = True
      gap = ''
      for i in range(len(line)):
         if blank:
            if line[i] != ' ':
               begin = i
               blank = False
            else:
               gap = gap + ' '
         else:
            if line[i] == ' ':
               end = i
               blank = True
               self.fields.append(line[begin:end])
               self.spaces.append(gap)
               gap = ' '
               n += 1
      if not blank:
         end = len(line)
         self.fields.append(line[begin:])
         self.spaces.append(gap)
         n +=1
      self.size = n

   def __str__(self):
      line = ''
      for i in range(self.size):
         line += self.spaces[i]
         line += self.fields[i]
      return line
