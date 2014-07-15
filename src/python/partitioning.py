import math
import numpy as np

from algebraic import laplacianEigenvectors
from _NetworKit import Partition

class SpectralPartitioner(object):
	"""

	Please note that the code in this class assumes the nodes of a graph to be numbered
	from 0 to n.

	"""
	def __init__(self, graph, depth, balanced=True):
		self.graph = graph
		self.depth = depth

		self.balanced = balanced

	def prepareSpectrum(self):
		spectrum = laplacianEigenvectors(self.graph, cutoff = (self.depth + 1))
		self.eigenvectors = spectrum[1]
		self.eigenvalues = spectrum[0]

	def getMedian(self, eigv, vertices):
		values = [eigv[i] for i in vertices]
		values.sort()
		median = values[math.floor(len(values) / 2)]

		return median

	def getMean(self, eigv, vertices):
		values = [eigv[i] for i in vertices]
		mean = np.mean(values)

		return mean

	def bisect(self, depth, partition=None, iteration=1):
		if partition is None:
			vertices = self.graph.nodes()
		else:
			vertices = self.partitions[partition]


		eigv = self.eigenvectors[iteration]

		if (self.balanced):
			split = self.getMedian(eigv, vertices)
		else:
			split = self.getMean(eigv, vertices)

		partA = self.nextPartition
		partB = self.nextPartition + 1
		self.nextPartition += 2		# TODO this is not thread-safe


		self.partitions[partA] = []
		self.partitions[partB] = []

		for vertex in vertices:
			if (eigv[vertex] < split):
				self.partitions[partA].append(vertex)
			else:
				self.partitions[partB].append(vertex)

		if (not (partition is None)):
			del self.partitions[partition]

		if depth > 1:
			self.bisect(depth - 1, partition = partA, iteration = iteration + 1)
			self.bisect(depth - 1, partition = partB, iteration = iteration + 1)

	def generatePartition(self):
		partition = Partition(size=self.graph.numberOfNodes())

		for partIndex in self.partitions:
			vertices = self.partitions[partIndex]

			if (len(vertices) < 1):
				continue

			firstItem = vertices[0] # TODO come on.. there has to be an easier way of doing this...
			partition.toSingleton(firstItem)
			subsetID = partition[firstItem]

			for vertex in vertices[1:]:
				partition.addToSubset(subsetID, vertex)

		self.partition = partition
		return partition

	def run(self):
		self.nextPartition = 0
		self.partitions = {}
		self.prepareSpectrum()

		self.bisect(self.depth)

		self.generatePartition()

	def getPartition(self):
		return self.partition