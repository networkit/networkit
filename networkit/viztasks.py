# local imports
from . import nxadapter
from . import community
from . import centrality
from _NetworKit import ParallelPartitionCoarsening

# external imports
import networkx

def save(name, dir="."):
	""" Save a figure """
	savefig(os.path.join(dir, "{0}.pdf".format(name)), bbox_inches="tight", transparent=True)


def coloringToColorList(G, coloring):
	clist = []

	nColors = len(coloring.keys())

	for v in G.nodes():
		clist.append(float(coloring[v]) / nColors)

	return clist


def drawGraph(G, **kwargs):
	""" Draws a graph via networkX. Passes additional arguments beyond the graph to networkx.draw(...).
	    By default, node sizes are scaled between 30 and 300 by node degree.
	"""
	nxG = nxadapter.nk2nx(G)
	if not "node_size" in kwargs:
		kwargs["node_size"] =[30+270*s for s in centrality.DegreeCentrality(G,True).run().scores()],
	networkx.draw(nxG, **kwargs)

def drawCommunityGraph(G, zeta, **kwargs):
	""" Draws the community graph for a given graph and partition. Passes any additional arguments to networkx.draw(...).
	    By default, node sizes are scaled between 30 and 500 by community size.
	"""
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
