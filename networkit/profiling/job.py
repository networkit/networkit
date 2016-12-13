#
# file: multiprocessing.py
# author: Mark Erb
#

class Job:
	""" abstract job class """

	def __init__(self, type, name):
		""" constructor """
		self.__name = name
		self.__type = type
		
		
	def getName(self):
		""" return the name connected with the job """
		return self.__name
		
	def getType(self):
		""" returns the type of the job """
		return self.__type
		
	def run(self):
		""" computation """
		raise Error("abstract class: overwrite this method")