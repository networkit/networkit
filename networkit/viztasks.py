# local imports
from . import nxadapter
from . import community
from . import viztools

def save(name, dir="."):
	""" Save a figure """
	savefig(os.path.join(dir, "{0}.pdf".format(name)), bbox_inches="tight", transparent=True)


def coloringToColorList(G, coloring):
	clist = []

	nColors = len(coloring.keys())

	for v in G.nodes():
		clist.append(float(coloring[v]) / nColors)

	return clist


def drawGraph(G, figsize=(7,7), labelled=False, nodeSizes=None, layout=None, coloring=None):
	""" Draw a graph"""
	nxG = nxadapter.nk2nx(G)

	if layout is None:
		pos = None
	else:
		pos = layout.layout(G)

	if not coloring is None:
		colors = coloringToColorList(G, coloring)
	else:
		colors = None

	drawer = viztools.drawing.GraphDrawer(size=figsize, labelled=labelled)
	drawer.draw(nxG, nodeSizes=nodeSizes, pos=pos, nodeColors=colors)


def drawCommunityGraph(G, zeta, labelled=False, figsize=(7,7)):
	cg = community.CommunityGraph() #convert communities to nodes
	cg.run(G, zeta)
	graph = cg.getGraph()
	comGraph = nxadapter.nk2nx(graph) #Convert a NetworKit.Graph to a networkx.Graph
	drawer = viztools.drawing.GraphDrawer(size=figsize, labelled=labelled)
	if graph.isWeighted():
		drawer.drawWeighted(comGraph)
	else:
		drawer.draw(comGraph, nodeSizes=[size*2 for size in list(zeta.subsetSizeMap().values())])

