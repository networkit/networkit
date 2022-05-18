from warnings import warn

from .community import EdgeCut, GraphClusteringTools, inspectCommunities
from .community import SpectralPartitioner as CommunitySpectralPartitioner

def computeEdgeCut(partition, graph):
	"""
	computeEdgeCut(partition, graph)

	DEPRECATED. This function (and the networkit.partioning module) will be removed in future updates.
	Use networkit.community.EdgeCut() instead.

	Compute edge cut given by graph and a partition.

	Parameters
	----------
	partition : networkit.Partition
		The input partition.
	graph : networkit.Graph
		The input graph.
	"""
	warn("networkit.partitioning.computeEdgeCut is deprecated, will be removed in future updates. Use networkit.community.EdgeCut instead.")
	eCut = EdgeCut()
	return eCut.getQuality(partition, graph)

def computeImbalance(partition, graph):
	"""
	computeImbalance(partition, graph)

	DEPRECATED. This function (and the networkit.partioning module) will be removed in future updates.
	Use networkit.community.GraphClusteringTools.getImbalance instead.
	
	Compute imbalance given by graph and a partition.

	Parameters
	----------
	partition : networkit.Partition
		The input partition.
	graph : networkit.Graph
		The input graph.
	"""	
	return GraphClusteringTools.getImbalance(partition, graph)

def inspectPartitions(partition, graph):
	"""
	inspectPartitions(partition, graph)

	DEPRECATED. This function (and the networkit.partioning module) will be removed in future updates.
	Use networkit.community.inspectCommunities instead.
	
	Compute and visualize properties of a graph partition.

	Parameters
	----------
	partition : networkit.Partition
		The input partition.
	graph : networkit.Graph
		The input graph.
	"""	
	inspectCommunities(partition, graph)

class SpectralPartitioner(CommunitySpectralPartitioner):
	"""
	SpectralPartitioner(graph, count, balanced=True)

	DEPRECATED. This class (and the networkit.partioning module) will be removed in future updates.
	Use networkit.community.SpectralPartitioner instead.

	Class to do spectral partitioning.

	Please note that the code in this class assumes the nodes of a graph to be numbered
	from 0 to n.

	Parameters
	----------
	graph : networkit.Graph
		The input graph.
	count : int
		The number of partitions to create.
	balanced : bool, optional
		Set this to false if you do not want to enforce balance, possibly increasing quality. Default: True
	"""
	def __init__(self, graph, count, balanced=True):
		warn("networkit.partitioning.SpectralPartitioner is deprecated, will be removed in future updates. Use networkit.community.SpectralPartitioner instead.")
		super().__init__(graph, count, balanced)
