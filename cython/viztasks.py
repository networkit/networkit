import nxadapter
import community
import viztools

def save(name, dir="."):
	""" Save a figure """
	savefig(os.path.join(dir, "{0}.pdf".format(name)), bbox_inches="tight", transparent=True)


def drawGraph(G, figsize=(7,7), labelled=False, nodeSizes=None):
	""" Draw a graph"""
	nxG = nxadapter.nk2nx(G)
	drawer = viztools.drawing.GraphDrawer(size=figsize, labelled=labelled)
	drawer.draw(nxG, nodeSizes)


def drawCommunityGraph(G, zeta, labelled=False, figsize=(7,7)):
	cg = community.CommunityGraph()
	cg.run(G, zeta)
	comGraph = nxadapter.nk2nx(cg.getGraph())
	drawer = viztools.drawing.GraphDrawer(size=figsize, labelled=labelled)
	drawer.draw(comGraph, nodeSizes=[size*2 for size in list(zeta.clusterSizeMap().values())])