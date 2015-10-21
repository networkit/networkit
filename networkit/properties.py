# other submodules
from .components import ConnectedComponents, ParallelConnectedComponents, StronglyConnectedComponents
from .globals import ClusteringCoefficient
from .distance import Diameter, Eccentricity, EffectiveDiameter
from .correlation import Assortativity
from . import community
from . import centrality
from . import termgraph
from . import auxiliary
from . import nxadapter

# other modules
import textwrap
import collections
import math
import logging

try:
	import powerlaw
except ImportError:
	logging.warning("""WARNING: module 'powerlaw' not installed, which is required by some
						functions.""")

try:
	import networkx as nx
except ImportError:
	logging.warning("""WARNING: module 'networkx' not installed, which is required by some
						functions.""")
try:
	import tabulate
except ImportError:
	logging.warning("""WARNING: module 'tabulate' not installed, which is required by some
						functions. In order to install 'tabulate', Python 3.3 is required""")

try:
	from scipy import stats
except ImportError:
	logging.warning("""WARNING: module 'scipy' not installed, which is required by some
						functions.""")


########  PROPERTIES ########


#
# # TODO: move to profiling module
# def components(G):
# 	""" Find and analyze detected components.
# 		Returns the number of components and the sizes
# 		of each component. For more details use the
# 		ConnectedComponents class.
# 	"""
# 	logging.info("[...] finding connected components....")
# 	if G.isDirected():
# 		cc = StronglyConnectedComponents(G)
# 	else:
# 		cc = ConnectedComponents(G)
# 	cc.run()
# 	components = cc.getPartition()
# 	nComponents = components.numberOfSubsets()
# 	componentSizes = components.subsetSizeMap()
# 	return (nComponents, componentSizes)
#
# def numberOfComponents(G):
# 	""" Find and number of components """
# 	logging.info("[...] finding connected components....")
# 	cc = ConnectedComponents(G)
# 	cc.run()
# 	nComponents = cc.numberOfComponents()
# 	return nComponents




# def degeneracy(G):
# 	""" degeneracy of an undirected graph is defined as the largest k for which
# 	the graph has a non-empty k-core. Degeneracy is only implemented for graphs without
# 	self-loops."""
# 	if G.numberOfSelfLoops() > 0:
# 		raise NotImplementedError("Call Graph.removeSelfLoops() first.")
# 	coreDec = centrality.CoreDecomposition(G)
# 	coreDec.run()
# 	return coreDec.maxCoreNumber()
