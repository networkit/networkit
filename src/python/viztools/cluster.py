

class ClusteringGenerator:
	"""
	Generation of Clusterings:
	
                Clusterings are saved in dictionaries, so that clustering[i] contains the cluster id for node i.
	"""

	def generateOneClustering(self, G):
		"""
		Generate one-clustering of a networkx graph.
		In a one-clustering every node is in cluster 0.

		:param G: The graph object
		:rtype: dictionary
		"""

		clustering = {}

		for u in G.nodes():
			clustering[u] = 0

		return clustering

	def generateSingletonClustering(self, G):
		"""
                Generate singleton clustering of a networkx graph.
                In a singleton clustering node i is in cluster i.

		:param G: The graph object
		:rtype: dictionary
		"""
		
		clustering = {}
		clust = range(len(G.nodes()))
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

		:param path: Path to the .clust file of a clustering
		:rtype: dictionary
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

		"""Write a clustering to .clust file format

		:param zeta: The clustering
		:param path: path to the .clust file where the clustering shall be written to
		"""


		with open(path, "w") as file:

			for u in range(0, len(zeta)):
				
				file.write(str(zeta[u]) + "\n")
				
