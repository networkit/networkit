import os

class TxtResultWriter:

	def __init__(self, outputDir):
		self.outputDir = outputDir
		self.SEPARATOR = "-------------------------------------"


	def getText(self, result):
		text = ""
		text += "GraphFile=" + result.task.graphPath + "\n"
		text += "GraphLoadingTime=" + str(result.loadingTime) + "\n"
		text += self.SEPARATOR + "\n"
		for bprops in result.backboneProperties:
			text += self.getGraphPropertiesText(bprops)
			text += self.SEPARATOR + "\n"
		return text

	def getGraphPropertiesText(self, props):
		text = ""
		for name in dir(props):
			attr = getattr(props, name)
			if not callable(attr) and not name.startswith('_'):
				text += str(name) + "=" + str(attr) + "\n"
		return text


	def receiveResult(self, taskResult):
		report = self.getText(taskResult)
		filename = os.path.join(self.outputDir, taskResult.task.graphName + "_PROPERTIES.txt")
		with open(filename, 'w') as f:
			f.write(report)
