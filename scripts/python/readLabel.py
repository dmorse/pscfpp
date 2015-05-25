def readLabelledLine(file, label):
   ''' 
   Read a line of the form "label = <string>", return the string that 
   follows the equals sign.
   '''
   line = file.readline()
   groups = line.strip().split('=')
   if groups[0].strip() != label:
      print "Error: Expected label = " + label
      print "     : Found          = " + groups[0].strip()
      raise
   if len(groups) > 2:
      print "Error: More than one = sign in line"
      raise
   if len(groups) == 1:
      return ''
   else:
     return groups[1].strip()

def readLabelledList(file, label):
   ''' 
   Read a line of the form "label = <string>", in which the string after
   the = sign contains several items separated by spaces, such as 
   'a b c', return a list of the items, such as ['a', 'b', 'c'].
   '''
   line = readLabelledLine(file, label)
   list = []
   if line != '':
      groups = line.split(' ')
      for i in range(len(groups)):
         if groups[i] != '':
            list.append(groups[i])
   return list
