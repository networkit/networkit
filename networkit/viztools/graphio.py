'''
Created on 19.03.2013

@author: cls

.. module:: graphio
'''

#TODO: check functions and update docstrings if necessary

# external imports
import sys
import networkx as nx
import logging


class METISGraphReader:
	"""
	Reads a graph from a file in the METIS format.
	In METIS file a graph of n nodes is stored in a file of n+1 lines.
	The first line lists the number of nodes, the number of edges and the weight of edges.
	The remaining lines list the neighbors of the current node.
	Nodes are indexed beginning from zero.
	"""


	def read(self, path):
		"""
		Read networkx graph from METIS file

		:param path: Path to the METIS file of a graph
		:rtype: The graph object
		"""
		G = nx.Graph()
		u = 0

		def parseUnweightedLine(line, u):
			"""
			:param line: current line in METISfile
			:param u: Node-id u
			"""

			for token in line.split(" "):
				if (len(token) > 0):
					v = int(token)

					if not G.has_edge(u, v):
						G.add_edge(u, v)


		def parseWeightedLine(line, u):
			"""
			:param line: current line in METISfile
			:param u: Node-id u
			"""
			ln = line.split(" ")
			i = 0 # i is index for nodes
			j = 1 # j is index for weight

			while j < len(ln):
				v = int(ln[i])
				w = int(ln[j]) #TODO: weights can also be non integer
				if not G.has_edge(u, v):
					G.add_edge(u, v, weight = w)
				i += 2
				j += 2
				
		def comment(path):
			"""returns the line where the graphfile begins"""
			with open(path, "r") as graphFile:
				j = 0
				for line in graphFile:
					line = line.strip()
					fls = line.split()
					if list(fls[0])[0]=="%":
						j += 1
					else:
						break
			return j


		# handle first line
		with open(path, "r") as graphFile:
			begin = comment(path)
			linenumber = 0
			for line in graphFile:
				#skip comments
				if linenumber < begin:
					linenumber += 1
				#read first line
				elif linenumber == begin:
					line = line.strip()
					fls = line.split(" ")
					linenumber += 1

					if len(fls) == 3:
							(n, m, w) = (int(fls[0]), int(fls[1]), int(fls[2]))
							if w == 1:
									weighted = True
							elif w == 0:
									weighted = False
					elif len(fls) == 2:
							(n,m) = (int(fls[0]), int(fls[1]))
							weighted = False
			
		if weighted:
			parseLine = parseWeightedLine
		else:
			parseLine = parseUnweightedLine

		# parse all lines
		with open(path, "r") as graphFile:
			linenumber = 0
			inFirstLine = True
			begin = comment(path)
			
			for line in graphFile:
				line = line.strip()
				#skip comments
				if linenumber < begin:
					linenumber += 1
				# read graphfile
				else:
					if inFirstLine:
						# skip first line
						inFirstLine = False
						u += 1
					else:
						parseLine(line, u)
						u += 1

		logging.debug("n = {0}".format(n))
		logging.debug("number of nodes= {0}".format(G.number_of_nodes()))
		logging.debug("m = {0}".format(m))
		logging.debug("number of edges= {0}".format(G.number_of_edges()))
 
		assert n == G.number_of_nodes()
		assert m == G.number_of_edges()
 
 
		if not n == G.number_of_nodes():
			raise InputError("number of nodes of G is not the same as the number in the file")
		if not m == G.number_of_edges():
			raise InputError("number of edges of G is not the same as the number in the file")

		G = nx.convert_node_labels_to_integers(G, first_label=0)

		return G

class METISGraphWriter:

	def write(self, G, path):
		"""
		Write networkx graph to the METIS file format

		:param G: The graph object
		:param path: Path to the METIS file where graph G is written to
		"""

		G = nx.convert_node_labels_to_integers(G, first_label=1)

		#check if graph is weighted
		weight = True
		h = G.neighbors(1)
		if G.edge[1][h[0]]=={}:
			weight = False

		with open(path, "w") as line:

			#write weighted graph
			if weight:
				firstline = str(G.number_of_nodes())+" "+str(G.number_of_edges())+" "+str(1)
				line.write(firstline+" \n")
				for u in G.nodes():
					h = G.neighbors(u)
					for v in h:
						line.write(str(v)+" "+str(G.edge[u][v]["weight"])+" ")
					line.write("\n")

			#write unweighted graph
			else:
				firstline = str(G.number_of_nodes())+" "+str(G.number_of_edges())+" "+str(0)
				line.write(firstline+" \n")
				for u in G.nodes():
					h = G.neighbors(u)
					for v in h:
						line.write(str(v)+" ")
					line.write("\n")

			
					
def readGraph(path):
	"""
	Read the graph from the given path

	:param path: The path of the graph
	"""
	
	def recognizeFormat(path):
		"""
		Recognize the format of the graphfile

		:param path: The path of the graph
		"""
		parts = path.split(".")
		ending = parts[-1]
		if ending == "graph":
			return "METIS"
		elif ending == "edgelist":
			return "edgelist"
		elif ending == "adjlist":
			return "adjlist"
		else:
			return None		
	
	def readMETIS(path):
		G = METISGraphReader().read(path)
		return G

	def readAdjlist(path):
		G = nx.read_adjlist(path)
		return G
	
	def readEdgelist(path):
		G = nx.read_edgelist(path)
		return G
	
	# TODO: add more formats when needed
		
	format = recognizeFormat(path)
	
	if format == "METIS":
		G = readMETIS(path)
	elif format == "edgelist":
		G = readEdgelist(path)
	elif format == "adjlist":
		G = readAdjlist(path)
	elif format is None:
		raise RuntimeError("format not recognized")
	return G
	

# some helper functions - are they still useful?
def appendNodeAttributes(G, map, attrname):
	""" Read attributes from a map node id -> attribute and append them to the nodes of G

	:rtype : object
	:param G: The graph object
	:param map: 
	:param attrname:
	:rtype: The graph Object 
	"""
	for v in G.nodes():
		if v in map:
			G[v][attrname] = map[v]
	return G
	



if __name__ == "__main__":
	pass
	

