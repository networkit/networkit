from NetworKit import *
from dynamic import *
from centrality import *
import pandas as pd

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


def setRandomWeights(G, mu, sigma):
	"""
	Add random weights, normal distribution with mean mu and standard deviation sigma
	"""
	for (u, v) in G.edges():
		w = random.normalvariate(mu, sigma)
		G.setWeight(u, v, w)
	return G



def test(G, nEdges, batchSize, epsilon, delta):
	# find a set of nEdges to remove from G
	T = graph.SpanningForest(G).generate()
	(removeStream, addStream) = removeAndAddEdges(G, nEdges, tabu=T)
	# remove the edges from G
	updater = dynamic.GraphUpdater(G)
	updater.update(removeStream)
	# run the algorithms on the inital graph
	bc = Betweenness(G)
	print ("Running bc")
	bc.run()
	dynBc = DynBetweenness(G, True)
	print ("Running dyn bc with predecessors")
	dynBc.run()
	apprBc = ApproxBetweenness(G, epsilon, delta)
	print ("Running approx bc")
	apprBc.run()
	print (apprBc.numberOfSamples())
	dynApprBc = DynApproxBetweenness(G, epsilon, delta, True)
	print ("Running dyn approx bc with predecessors")
	dynApprBc.run()
	# apply the batches
	nExperiments = nEdges // batchSize
	timesBc = []
	timesDynBc = []
	timesApprBc = []
	timesDynApprBc = []
	for i in range (0, nExperiments):
		batch = addStream[i*batchSize : (i+1)*batchSize]
		# add the edges of batch to the graph
		totalTime = 0.0
		for j in range (0, batchSize):
			updater.update([batch[j]])
			# update the betweenness with the dynamic exact algorithm
			t = stopwatch.Timer()
			dynBc.update(batch[j])
			totalTime += t.stop()
		timesDynBc.append(totalTime)
		# update the betweenness with the static exact algorithm
		t = stopwatch.Timer()
		bc.run()
		x = t.stop()
		timesApprBc.append(x)
		print ("Exact BC")
		print (x)
		print ("Speedup Dyn BC (with preds)")
		print (x/totalTime)
		# update the betweenness with the static approximated algorithm
		t = stopwatch.Timer()
		apprBc.run()
		x = t.stop()
		timesBc.append(x)
		print ("ApprBC")
		print (x)
		# update the betweenness with the dynamic approximated algorithm
		t = stopwatch.Timer()
		dynApprBc.update(batch)
		y = t.stop()
		timesDynApprBc.append(x)
		print ("Speedup DynApprBC (with preds)")
		print (x/y)
<<<<<<< local
	a = pd.Series(timesBc)
	b = pd.Series(timesDynBc)
	c = pd.Series(timesApprBc)
	d = pd.Series(timesDynApprBc)
	df = pd.DataFrame({"Static exact bc": a, "Dynamic exact bc" : b, "Static approx bc" : c, "Dynamic approx bc" : d})
	return df
=======
		# update the betweenness with the dynamic approximated algorithm
		t = stopwatch.Timer()
		dynApprBc2.update(batch)
		y = t.stop()
		timesDynApprBc2.append(x)
		print ("Speedup DynApprBC (without preds)")
		print (x/y)
	return {
			"bc" : timesBc,
			"dynBc (predecessors stored)" : timesDynBc,
			"dynBc (predecessors not stored)" : timesDynBc2,
			"apprBc" : timesApprBc,
			"dynApprBc (predecessors stored)" : timesDynApprBc,
			"dynApprBc (predecessors not stored)" : timesDynApprBc2
			}
>>>>>>> other


<<<<<<< local
if __name__ == "__main__":
	setNumberOfThreads(1)
	size = 1000
	G = generators.DorogovtsevMendesGenerator(size).generate()
	G1 = Graph(G.numberOfNodes(), True, False)
	for e in G.edges():
		G1.addEdge(e[0], e[1], 1.0)
	cc = properties.ConnectedComponents(G1)
	cc.run()
	if (cc.numberOfComponents() == 1) :
		nEdges = 10
		batchSize = 10
		epsilon = 0.05
		delta = 0.1
		df = test(G1, nEdges, batchSize, epsilon, delta)
		df.to_csv("results/weighted_size_"+str(size)+"_batch_"+str(batchSize)+".csv")
	else:
		print ("The generated graph is not connected.")
=======
setNumberOfThreads(1)
G = generators.DorogovtsevMendesGenerator(10000).generate()
G1 = Graph(G.numberOfNodes(), True, False)
for e in G.edges():
	G1.addEdge(e[0], e[1], 1.0)
cc = properties.ConnectedComponents(G1)
cc.run()
if (cc.numberOfComponents() == 1) :
	nEdges = 10
	batchSize = 1
	epsilon = 0.1
	delta = 0.1
	times = test(G1, nEdges, batchSize, epsilon, delta)
	print (times)
else:
	print ("The generated graph is not connected.")
>>>>>>> other
