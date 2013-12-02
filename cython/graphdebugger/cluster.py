import networkx as nx


class ClusteringGenerator:

	def generateOneClustering(self, G):
		"""
		Generate one-clustering of a networkx graph

		:param G: The graph object
		:rtype: clustering 
		"""

		clustering = {}

		for u in G.nodes():
			clustering[u] = 0

		return clustering

	def generateSingletonClustering(self, G):
		"""
		Generate singleton clustering of a networkx graph

		:param G: The graph object
		:rtype: clustering
		"""

		clust = []
		u = 0
		while u < G.number_of_nodes():
			clust += [u]
			u += 1

		clustering = {}
		u = 0
		for v in G.nodes():
			clustering[v] = clust[u]
			u += 1

		return clustering

class ClusteringReader:

	""" The .clust file format:

			Line i holds the cluster id for node i."""

	

	def read(self, path):

		"""Read a clustering from a .clust file.

		"""

		zeta = {}

		with open(path, "r") as file:

			u = 0

			for line in file:

				zeta[u] = int(line.strip())

				u += 1

		return zeta

	

class ClusteringWriter:

	""" The .clust file format:

		Line i holds the cluster id for node i."""

	

	def write(self, zeta, path):

		"""Write a clustering to .clust file format"""

		maxu = max(zeta.keys())

		

		with open(path, "w") as file:

			for u in range(0, maxu + 1):

				file.write("%d\n".format(u))
