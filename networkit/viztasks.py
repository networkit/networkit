# local imports
from . import nxadapter
from . import community
from . import centrality
from .coarsening import ParallelPartitionCoarsening
from .support import MissingDependencyError

from warnings import warn

# external imports
try:
	import networkx
except ImportError:
	have_nx = False
else:
	have_nx = True

def save(name, dir="."):
	""" 
	save(name, dir=".")

	Save a figure.

	DEPRECATED. This function (and the networkit.viztasks module) will be removed in future updates.

	Parameters
	----------
	name : str
		Name of the output file.
	dir : str
		Output directory. Default: "."
	"""
	warn("networkit.viztasks.save is deprecated, will be removed in future updates.")
	savefig(os.path.join(dir, "{0}.pdf".format(name)), bbox_inches="tight", transparent=True)


def coloringToColorList(G, coloring):
	"""
	coloringToColorList(G, coloring)

	Calculate node colors based on an input graph and color dict.

	DEPRECATED. This function (and the networkit.viztasks module) will be removed in future updates.

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	coloring : dict()
		Coloring scheme as dict.

	Returns
	-------
	list(tuple(float, float, float))
		List with color values for each node.
	"""
	warn("networkit.viztasks.coloringToColorList is deprecated, will be removed in future updates.")
	clist = []

	nColors = len(coloring.keys())

	for v in G.iterNodes():
		clist.append(float(coloring[v]) / nColors)

	return clist


def drawGraph(G, **kwargs):
	""" 
	drawGraph(G, **kwargs)
	
	Draws a graph via networkX. Passes additional arguments beyond the graph to networkx.draw(...).
	By default, node sizes are scaled between 30 and 300 by node degree.

	DEPRECATED. This function (and the networkit.viztasks module) will be removed in future updates. Use networkit.vizbridges instead to draw 
	graphs (needs additional plugins).

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	`**kwargs` : dict()
		Additional arguments for networkx.draw function.
	"""
	warn("networkit.viztasks.drawGraph is deprecated, will be removed in future updates. Use networkit.vizbridges instead to draw graphs (needs additional plugins).")
	if not have_nx:
		raise MissingDependencyError("networkx")
	if not G.checkConsistency():
		print("WARNING: Multi-graph has been converted to simple graph for display")
		G.removeMultiEdges()
	nxG = nxadapter.nk2nx(G)
	if not "node_size" in kwargs:
		kwargs["node_size"] = [30+270*s for s in centrality.DegreeCentrality(G,True).run().scores()]
	networkx.draw(nxG, **kwargs)

def drawCommunityGraph(G, zeta, **kwargs):
	""" 
	drawCommunityGraph(G, zeta, **kwargs)	

	Draws the community graph for a given graph and partition. Passes any additional arguments to networkx.draw(...).
	By default, node sizes are scaled between 30 and 500 by community size.

	DEPRECATED. This function (and the networkit.viztasks module) will be removed in future updates. Use networkit.vizbridges instead to draw 
	graphs (needs additional plugins).

	Parameters
	----------
	G : networkit.Graph
		The input graph.
	zeta : networkit.Partition
		The input partition.
	`**kwargs` : dict()
		Additional arguments for networkx.draw function.	
	"""
	warn("networkit.viztasks.drawCommunityGraph is deprecated, will be removed in future updates. Use networkit.vizbridges instead to draw graphs (needs additional plugins).")
	if not have_nx:
		raise MissingDependencyError("networkx")
	cg = ParallelPartitionCoarsening(G,zeta)
	cg.run() # convert communities to nodes
	graph = cg.getCoarseGraph()
	comGraph = nxadapter.nk2nx(graph)
	if not "node_size" in kwargs:
		sizes = list(zeta.subsetSizeMap().values())
		max_size = max(sizes)
		sizes = [elem/max_size for elem in sizes]
		kwargs["node_size"] = [30+470*s for s in sizes]
	networkx.draw(comGraph, **kwargs)
