# TODO: check functions and update docstrings if necessary
'''
Created on May 23, 2013

@author: beddig
'''

def EdgeMapWriter(d, path):
	"""
	Write an edgemap

	:param path: Path to the file where the edgemap is written to 
	"""
	
	with open(path, "w") as edgemap:
		
		firstLine = True
		
		# write first line
		if firstLine:
			for u in d:
				if type(d[u]) == str:
					
					edgemap.write("string" + " \n")
					break
					
				elif type(d[u]) == int:
					edgemap.write("int" + " \n")
					break
					
				elif type(d[u]) == float:
					edgemap.write("float" + " \n")
					break
				
			firstLine = False
		
		for u in d:
			edgemap.write(str(u[0]) + " " + str(u[1]) + " " + str(d[u]) + " \n")
			
def EdgeMapReader(path):
	"""
	read an edgemap

	:param path: Path to the file that contains the edgemap
	"""
	
	d = {}
	
	# check if edge attributes are int, float or string
	with open(path, "r") as edgeFile:
		
		string = False
		integer = False
		floating = False
		
		firstline = edgeFile.readline()
		firstline = firstline.strip()
		
		if firstline == "string":
			string = True
		elif firstline == "int":
			integer = True
		elif firstline == "float":
			floating = True
			
	def String(line):
		d[int(line[0]), int(line[1])] = str(line[2])
		
	def Integer(line):
		d[int(line[0]), int(line[1])] = int(line[2])	
		
	def Floating(line):
		d[int(line[0]), int(line[1])] = float(line[2])
		
	if string:
		typ = String
	if integer:
		typ = Integer
	if floating:
		typ = Floating	
			
	with open(path, "r") as edgeFile:
		
		inFirstLine = True
		
		for line in edgeFile:
			line = line.strip()
			line = line.split(" ")
			
			# skip first line
			if inFirstLine:
				inFirstLine = False
				
			else:
				typ(line)
					
	return d

def addAttributesFromEdgeMap(G, path):
	"""
	add attributes from edgemap to graph G

	:param G: The graph object
	:param path: Path to the file that contains the edgemap
	:rtype G: The graph object
	"""
	
	from graphdebugger import edgemap
	d = edgemap.EdgeMapReader(path)   
	
	# append attributes
	
	edges = []
	for u in d:
		edges += [u]
	
	# create list of edges and attributes	
	u = 0 
	attributes = []
	while u < len(edges):
		attributes += [(edges[u][0], edges[u][1], {edges[u]:d[edges[u]]})]
		u += 1
	
	G.add_edges_from(attributes)
	
	return G
