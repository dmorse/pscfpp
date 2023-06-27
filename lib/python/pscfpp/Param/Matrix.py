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