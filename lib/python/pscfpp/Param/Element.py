class Element:
	'''
	Purpose: 
		The base class of all types of elements of the parameter file, with four
		subclasses:
			Composite: 
				Parameter block identified by opening and closing curly brackets 
				('{' and '}')
			Parameter: 
				Single value parameter within a single line
			Array: 
				Parameter has values stores in an one-dimentional array identified 
				by openning and closing square brackets ('[' and ']')
			Matrix:
				Parameter has values stores in a two-dimentional array (matrix)
				identified by openning and closing parenthesis
	Instance variables:
		label: Element label of the parameter file
		parent: The parent of the element, defult to be None
	Methods:
		__init__(self, label, openFile=None, parent=None, val=None):
			the constructor of the Element object for initiation, with four 
			arguments: 
				label, the label of the Element; 
				openFile, file name of the opened parameter file, defult to be None;
				parent, the parent of the Elment, defult to be None;
				val, the value of the Element, defult to be None;
			that initiate the label and the parent of the Element object
		getLabel(self): 
			return the label of the element
		returnData(slef):
			pass for the base class
		writeOutString(self, depth=''):
			pass for the base class
	'''
	def __init__(self, label, openFile=None, parent=None, val=None):
		self.label = label
		self.parent = parent

	def getLabel(self):
		return self.label

	def returnData(self):
		pass

	def writeOutString(self, depth=''):
		pass