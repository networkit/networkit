from networkit import *
from networkit.dynamic import *
from networkit.centrality import *
from ..support import MissingDependencyError
try:
	import pandas as pd
except ImportError:
	have_pandas = False
else:
	have_pandas = True
import random

def isConnected(G):
	cc = properties.ConnectedComponents(G)
	cc.run()
	return (cc.numberOfComponents() == 1)

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
		addStream.append(GraphEvent(GraphEvent.EDGE_ADDITION, u, v, G.weight(u, v)))

	return (removeStream, addStream)


def setRandomWeights(G, mu, sigma):
	"""
	Add random weights, normal distribution with mean mu and standard deviation sigma
	"""
	for (u, v) in G.iterEdges():
		w = random.normalvariate(mu, sigma)
		G.setWeight(u, v, w)
	return G



def test(G, nEdges, batchSize, epsilon, delta, size):
	if not have_pandas:
		raise MissingDependencyError("pandas")
	# find a set of nEdges to remove from G
	T = graph.SpanningForest(G).generate()
	(removeStream, addStream) = removeAndAddEdges(G, nEdges, tabu=T)
	# remove the edges from G
	updater = dynamic.GraphUpdater(G)
	updater.update(removeStream)
	# run the algorithms on the inital graph
	print("--- IS G CONNECTED? ")
	print(isConnected(G))
	bc = Betweenness(G)
	print("Running bc")
	bc.run()
	dynBc = DynBetweenness(G, True)
	print("Running dyn bc with predecessors")
	dynBc.run()
	apprBc = ApproxBetweenness(G, epsilon, delta)
	print("Running approx bc")
	apprBc.run()
	dynApprBc = DynApproxBetweenness(G, epsilon, delta, True)
	print("Running dyn approx bc with predecessors")
	dynApprBc.run()
	# apply the batches
	nExperiments = nEdges // batchSize
	timesBc = []
	timesDynBc = []
	timesApprBc = []
	timesDynApprBc = []
	scoresBc = []
	scoresApprBc = []
	for i in range(nExperiments):
		batch = addStream[i*batchSize : (i+1)*batchSize]
		# add the edges of batch to the graph
		print("GRAPH SIZE")
		print(size)
		totalTime = 0.0
		for j in range(0, batchSize):
			updater.update([batch[j]])
			# update the betweenness with the dynamic exact algorithm
			if size <= 2**15:
				t = stopwatch.Timer()
				dynBc.update(batch[j])
				totalTime += t.stop()
			else:
				totalTime = -1
		timesDynBc.append(totalTime)
		# update the betweenness with the static exact algorithm
		t = stopwatch.Timer()
		bc.run()
		x = t.stop()
		timesBc.append(x)
		print("Exact BC")
		print(x)
		print("Speedup Dyn BC (with preds)")
		print(x/totalTime)
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
		bcNormalized = [ k/(size*(size-1)) for k in bc.scores()]
		scoresBc.append(bcNormalized)
		scoresApprBc.append(dynApprBc.scores())

	a = pd.Series(timesBc)
	b = pd.Series(timesDynBc)
	c = pd.Series(timesApprBc)
	d = pd.Series(timesDynApprBc)
	df1 = pd.DataFrame({"Static exact bc": a, "Dynamic exact bc" : b, "Static approx bc" : c, "Dynamic approx bc" : d})
	dic2 = {}
	for experiment in range(nExperiments):
		a = pd.Series(scoresBc[experiment])
		b = pd.Series(scoresApprBc[experiment])
		dic2["Exact scores (exp. "+str(experiment)+")"] = a
		dic2["Approx scores (exp. "+str(experiment)+")"] = b
	df2 = pd.DataFrame(dic2)
	return df1, df2


if __name__ == "__main__":
	setNumberOfThreads(1)
#	setLogLevel("INFO")
	batchSize = 128

	for i in range(10,21):
		size = 2**i
		G = generators.DorogovtsevMendesGenerator(size).generate()
		G1 = Graph(G.numberOfNodes(), True, False)
		for u, v in G.iterEdges():
			G1.addEdge(u, v, 1.0)
		G1 = setRandomWeights(G1, 1, 0.1)
		if (isConnected(G1)) :
			nEdges = batchSize * 5
			epsilon = 0.05
			delta = 0.1
			(df1, df2) = test(G1, nEdges, batchSize, epsilon, delta, size)
			df1.to_csv("results_fixed_batch/times_weighted_size_"+str(size)+"_batch_"+str(batchSize)+".csv")
			df2.to_csv("results_fixed_batch/scores_weighted_size_"+str(size)+"_batch_"+str(batchSize)+".csv")
		else:
			print("The generated graph is not connected.")
