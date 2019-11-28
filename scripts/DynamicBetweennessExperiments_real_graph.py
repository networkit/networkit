from networkit import *
from networkit.dynamic import *
from networkit.centrality import *
import operator
from ..support import MissingDependencyError
try:
	import pandas as pd
except ImportError:
	have_pandas = False
else:
	have_pandas = True
import random

def removeAndAddEdges(G, nEdges, tabu=None):
	if nEdges > G.numberOfEdges() - tabu.numberOfEdges():
		raise Error("G does not have enough edges")

	# select random edges for removal
	removed = set()
	while len(removed) < nEdges:
		(u, v) = G.randomEdge()
		if not tabu.hasEdge(u, v) and not ((u,v) in removed or (v,u) in removed):	# exclude all edges in the tabu graph
			removed.add((u, v))
	print (removed)
	# build event streams
	removeStream = []
	for (u, v) in removed:
		removeStream.append(GraphEvent(GraphEvent.EDGE_REMOVAL, u, v, 0))
	addStream = []
	for (u, v) in removed:
		addStream.append(GraphEvent(GraphEvent.EDGE_ADDITION, u, v, 1.0))

	return (removeStream, addStream)

def extractLargestComponent(G):
	"""
	Extract the subgraph of the largest connected component.

	Parameters
	----------
	G : Graph
		Input graph.
	Returns
	-------
	Graph
		Subgraph of largest component, preserving node ids of orignal graph.
	"""

	cc = properties.ConnectedComponents(G)
	cc.run()
	cSizes = cc.getComponentSizes()
	(largestCompo, size) = max(cSizes.items(), key=operator.itemgetter(1))
	compoNodes = [v for v in G.nodes() if cc.componentOfNode(v) is largestCompo]
	C = graph.Subgraph().fromNodes(G, compoNodes)
	return C


def setRandomWeights(G, mu, sigma):
	"""
	Add random weights, normal distribution with mean mu and standard deviation sigma
	"""
	for (u, v) in G.edges():
		w = random.normalvariate(mu, sigma)
		G.setWeight(u, v, w)
	return G



def test(G, T, nEdges, batchSize, epsilon, delta, size):
	if not have_pandas:
		raise MissingDependencyError("pandas")
	# find a set of nEdges to remove from G
	(removeStream, addStream) = removeAndAddEdges(G, nEdges, tabu=T)
	# remove the edges from G
	updater = dynamic.GraphUpdater(G)
	updater.update(removeStream)
	# run the algorithms on the inital graph
	apprBc = ApproxBetweenness(G, epsilon, delta)
	print("Running approx bc")
	apprBc.run()
	dynApprBc = DynApproxBetweenness(G, epsilon, delta, False)
	print("Running dyn approx bc with predecessors")
	dynApprBc.run()
	# apply the batches
	nExperiments = nEdges // batchSize
	timesApprBc = []
	timesDynApprBc = []
	for i in range(nExperiments):
		batch = addStream[i*batchSize : (i+1)*batchSize]
		# add the edges of batch to the graph
		updater.update(batch)
		# update the betweenness with the static approximated algorithm
		t = stopwatch.Timer()
		apprBc.run()
		x = t.stop()
		timesApprBc.append(x)
		print("ApprBC")
		print(x)
		# update the betweenness with the dynamic approximated algorithm
		t = stopwatch.Timer()
		dynApprBc.update(batch)
		y = t.stop()
		timesDynApprBc.append(y)
		print("Speedup DynApprBC (with preds)")
		print(x/y)

	c = pd.Series(timesApprBc)
	d = pd.Series(timesDynApprBc)
	df1 = pd.DataFrame({"Static approx bc" : c, "Dynamic approx bc" : d})
	return df1


if __name__ == "__main__":
	setNumberOfThreads(1)
	graph_name = "road_usa.metis.graph"
	#G = readGraph("/algoDaten/staudt/Graphs/Collections/DynBC/email-Enron.edgelist-t0.graph", Format.EdgeListTabZero)
	G = graphio.METISGraphReader().read("/algoDaten/staudt/Graphs/Collections/DynBC/"+graph_name)
	size = G.numberOfNodes()
	print("Number of nodes: "+str(size))
	G = extractLargestComponent(G)
	print("Took largest connected component of the graph.")
	T = graph.SpanningForest(G).generate()
	print("Computed spanning forest.")

	for i in range(11):
		batchSize = 2**i
		nEdges = batchSize * 10
		epsilon = 0.2
		delta = 0.1
		df1 = test(G, T, nEdges, batchSize, epsilon, delta, size)
		df1.to_csv("results/times_"+graph_name+"_batch_"+str(batchSize)+"NO_PREDS_epsilon="+str(epsilon)+".csv")
