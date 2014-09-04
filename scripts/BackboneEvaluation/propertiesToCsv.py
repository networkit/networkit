import csv
import glob

def handlePropertiesFile(path):
	print("Opening file ", path)
	graphFile = None
	graphLoadingTime = None
	currentProps = {}
	props = [] # list of dictionaries

	with open(path) as f:
		for line in f:
			if line.startswith("--"):
				if len(currentProps) > 0:
					props.append(currentProps)
				currentProps = {}
			else:
				name,value = line.replace('\n', '').replace('\r','').split("=")
				if name == "GraphFile":
					graphFile = value
				elif name == "GraphLoadingTime":
					graphLoadingTime = value
				else:
					currentProps[name] = value

	writeCsv(path, props)

def writeCsv(propertiesFile, propertiesList):
	if len(propertiesList) == 0 or len(propertiesList[0]) == 0:
		return

	with open(propertiesFile.replace("_PROPERTIES.txt", "_PROPERTIES.csv"), 'w') as csvfile:
		csvwriter = csv.writer(csvfile, delimiter=',')
		#Assume that all property dictionaries contain exactly the same set of keys.
		columnNames = sorted(list(propertiesList[0].keys()))
		csvwriter.writerow(columnNames)
		for properties in propertiesList:
			row = [properties[x] for x in columnNames]
			csvwriter.writerow(row)

inputFiles = glob.glob("./output/*_PROPERTIES.txt")
for inputFile in inputFiles:
	handlePropertiesFile(inputFile)
