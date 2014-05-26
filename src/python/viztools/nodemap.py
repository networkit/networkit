# TODO: check functions and update docstrings if necessary
'''
Created on May 23, 2013

@author: beddig
'''


def NodeMapReader(path):
	"""read nodemap, returns a dictionary

	:param path: Path to the file that contains a nodemap
	:rtype: Dictionary
	"""
	u = 0
	d = {}

	# check if node attributes are int, float or string
	with open(path, "r") as NodeFile:

		string = False
		integer = False
		floating = False

		firstLine = NodeFile.readline()
		firstLine = firstLine.strip()

		if firstLine=="string":
			string = True
		elif firstLine == "int":
			integer = True
		elif firstLine == "float":
			floating = True


	def enterString(line, u):
		"""
		Enters a string to the dictionary

		:param line:
		:param u:
		"""
		d[u] = str(line)

	def enterInteger(line, u):
		"""
		Enters an integer to the dictionary

		:param line:
		:param u:
		"""
		d[u] = int(line)

	def enterFloat(line, u):
		"""
		Enters a float to the dictionary

		:param line:
		:param u:
		"""
		d[u] = float(line)

	if string:
		enterAttribute = enterString
	elif integer:
		enterAttribute = enterInteger
	elif floating:
		enterAttribute = enterFloat

	with open(path, "r") as NodeFile:

		inFirstLine = True

		for line in NodeFile:
			line = line.strip()

			if inFirstLine:
				inFirstLine = False

			else:
				enterAttribute(line, u)
				u += 1

	return d

def NodeMapWriter(d, path):
	"""
	write a nodemap

	:param d: Dictionary which contains the nodemap
	:param path: Path to the file where the nodemap is written to
	"""

	with open(path, "w") as NodeFile:

		# check if node attributes are int, string or float
		if type(d[0]) == str:
			NodeFile.write("string" + " \n")
		elif type(d[0]) == int:
			NodeFile.write("int" + " \n")
		elif type(d[0]) == float:
			NodeFile.write("float" + " \n")

		for u in d:
			NodeFile.write(str(d[u]))
			NodeFile.write("\n")
