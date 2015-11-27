# local imports
from . import nxadapter
from . import community
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
	""" Draw a graph via networkX. It is possible to pass any optional parameters that networkx.draw(...) takes."""
	nxG = nxadapter.nk2nx(G)
	networkx.draw(nxG, **kwargs)

def drawCommunityGraph(G, zeta, **kwargs):
	""" Draws the community graph of a given graph and partition. Takes the same optional parameters as networkx.draw(...) except node_size."""
	cg = ParallelPartitionCoarsening() 
	graph,_ = cg.run(G,zeta) # convert communities to nodes
	comGraph = nxadapter.nk2nx(graph)
	kwargs["node_size"] = [size*2 for size in list(zeta.subsetSizeMap().values())]
	networkx.draw(comGraph, **kwargs)