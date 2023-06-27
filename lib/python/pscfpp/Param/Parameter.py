class Parameter(Element):
	'''
	Purpose: 
		The class represents the Parameter element of the parameter file, which
		is a subclass of the Element class
	Instance variables:
		label: the label of the individual Parameter element 
		parent: the parent of the current individual Parameter element 
		val: the value of the individual Parameter element 
	Methods:
		__inti__(self, label, openFile=None, parent=None, value=None):
			the constructor of the Parameter object for initiation, with four 
			arguments: 
				label, the label of the Parameter; 
				openFile, the open parameter file, defult to be None;
				parent, the parent of the current Parameter, defult to be None;
				val, the value of the Parameter;
			that initiate the label, parent and the stored value of the Parameter 
			object
		getValue(self): 
			print out the type of the value of the individual parameter and return 
			the stored value, without argument
		setValue(self, val):
			set the new value to the val variable with argument val
		__repr___(self):
			return the string that represents the stored value
		writeOutString(self, depth):
			return the string for writting out with argument depth, the string of 
			spaces that represents the level of the Parameter element
		returnData(self):
			return the exact value of the Parameter object stored as the Value object
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