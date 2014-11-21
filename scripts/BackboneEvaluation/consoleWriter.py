import os

class ConsoleWriter:

	def receiveResult(self, taskResult):
		for row in taskResult.data:
			print(row)
