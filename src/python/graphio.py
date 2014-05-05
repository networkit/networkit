from _NetworKit import (Graph, METISGraphReader, METISGraphWriter, DotGraphWriter, EdgeListIO, \
						 VNAGraphWriter, GMLGraphWriter, LineFileReader, SNAPGraphWriter, ClusteringReader, ClusteringWriter, DGSWriter, \
						  DGSStreamParser, GraphUpdater, FastMETISGraphReader)

import os
import logging
import numpy
import scipy.io

class formats:
	metis = "metis"
	fastmetis = "fastmetis"
	#CLUSTERING = "clustering"
	edgelist_tab_one = "edgelist-t1"
	edgelist_tab_zero = "edgelist-t0"
	edgelist_space_one = "edgelist-s1"
	edgelist_space_zero = "edgelist-s0"
	#class input:
		#METIS = "metis"
		#CLUSTERING = "clustering"
		#DGSSTREAM = "dgsstream"
	class output:
		#METIS = "metis"
		#CLUSTERING = "clustering"
		graphviz = "graphviz"
		gml = "gml"
		snap = "snap"
		vna = "vna"

# reading

def readGraph(path, format="fastmetis", **kwargs):
	"""    Read graph file in various formats and return a NetworKit::Graph"""
	
	readers =  {"metis": METISGraphReader(),
			"fastmetis" : FastMETISGraphReader(),
				#"edgelist": EdgeListIO(),
				"edgelist-t1" : EdgeListIO('\t', 1),
				"edgelist-t0": EdgeListIO('\t', 0),
				"edgelist-s1": EdgeListIO(' ', 1), 
				"edgelist-s0": EdgeListIO(' ', 1)
				}

	try:
		reader = readers[format]#(**kwargs)
	except KeyError:
		raise Exception("unrecognized format: {0}".format(format))


	if ("~" in path):
		path = os.path.expanduser(path)
		print("path expanded to: {0}".format(path))
	if not os.path.isfile(path):
		raise IOError("{0} is not a file".format(path))
	else:
		with open(path, "r") as file:    # catch a wrong path before it crashes the interpreter
			G = reader.read(path)
			return G

	return None


def readMat(path):
	""" Reads a Graph from a matlab object file containing an adjacency matrix"""
	matlabObject = scipy.io.loadmat(path)	
	# result is a dictionary of variable names and objects, representing the matlab object
	for (key, value) in matlabObject.items():
		if type(matlabObject[key]) is numpy.ndarray:
			A = matlabObject[key]
			break # found the matrix
			
	(n, n2) = A.shape
	if (n != n2):
		raise Exception("this (%sx%s) matrix is not square".format(n, n2))
	if not ((A.transpose() == A).all()):
		logging.warning("the adjacency matrix is not symmetric")
	G = Graph(n)
	nz = A.nonzero()
	edges = [(u,v) for (u,v) in zip(nz[0], nz[1])]
	for (u,v) in edges:
		G.addEdge(u, v)
	return G

# writing

def writeGraph(G, path, format="metis"):
	""" Write graph to various output formats """
	writers = 	{"metis" : METISGraphWriter(),
				"gexf": None,
				"vna": VNAGraphWriter(),
				"dot": DotGraphWriter(),
				"graphviz": DotGraphWriter(),
				"gml": GMLGraphWriter(),
				"edgelist-t1" : EdgeListIO('\t', 1),
				"edgelist-t0": EdgeListIO('\t', 0),
				"edgelist-s1": EdgeListIO(' ', 1), 
				"edgelist-s0": EdgeListIO(' ', 1)
				}
	try:
		writer = writers[format]
		writer.write(G, path)
		logging.info("wrote graph {0} to file {1}".format(G, path))
	except:
		raise Exception("format {0} currently not supported".format(format))		

class GraphConverter:
	
	def __init__(self, reader, writer):
		self.reader = reader
		self.writer = writer
		
	def convert(self, inPath, outPath):
		G = self.reader.read(inPath)
		self.writer.write(G, outPath)
		
	def __str__(self):
		return "GraphConverter: {0} => {0}".format(self.reader, self.writer)

def getConverter(fromFormat, toFormat):
	
	readers =  {"metis": METISGraphReader(),
				"edgelist-t1" : EdgeListIO('\t', 1),
				"edgelist-t0": EdgeListIO('\t', 0),
				"edgelist-s1": EdgeListIO(' ', 1), 
				"edgelist-s0": EdgeListIO(' ', 1)}    
	writers =  {"metis" : METISGraphWriter(),
				"gexf": None,
				"vna": VNAGraphWriter(),
				"dot": DotGraphWriter(),
				"graphviz": DotGraphWriter(),
				"gml": GMLGraphWriter(),
				"edgelist-t1": EdgeListIO('\t', 1),
				"edgelist-t0": EdgeListIO('\t', 0),
				"edgelist-s1": EdgeListIO(' ', 1), 
				"edgelist-s0": EdgeListIO(' ', 1)
				} 
	
	reader = readers[fromFormat]
	writer = writers[toFormat]
	
	return GraphConverter(reader, writer)


def convertGraph(fromFormat, toFormat, fromPath, toPath=None):
	converter = getConverter(fromFormat, toFormat)
	if toPath is None:
		toPath = "{0}.{1}.graph".format(fromPath.split(".")[0], toFormat)
	converter.convert(fromPath, toPath)
	print("converted {0} to {1}".format(fromPath, toPath))



# dynamic

def readStream(path, mapped=True, baseIndex=0):
	"""
		Read a graph event stream from a file.
	"""
	return DGSStreamParser(path, mapped, baseIndex).getStream()

def writeStream(stream, path):
	"""
		Write a graph event stream to a file.
	"""
	DGSWriter().write(stream, path)


def graphFromStreamFile(path, mapped=True, baseIndex=0):
	stream = readStream(path, mapped, baseIndex)
	G = Graph()
	gu = GraphUpdater(G)
	gu.update(stream)
	return G
