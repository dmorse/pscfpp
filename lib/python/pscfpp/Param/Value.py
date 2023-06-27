class Value:
	'''
	Purpose: 
		The class represents the object type that store the single value for parameters 
		of the parameter file by distinguishing the exact type of it
	Instance variables:
		val: variable to store a single value of the parameters
		type: the type of the value stored, either integer, floating point or string
	Methods:
		__init__(self, val):
			the constructor of the Value object for initiation, with one argument:
				val, the string represents the value needed to be stored for the Value object
			that initiate the value stored with its correct type
		getType(self): return the type of stored value
		getValue(self): return the value of stored value
	Note:
		Debugging by the commented line to check if the constructor has the expected 
		function of distinguishing the exact type of the value

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