from NetworKit import *
from dynamic import *
from centrality import *

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
	print ("Running dyn bc")
	dynBc.run()
	apprBc = ApproxBetweenness(G, epsilon, delta)
	print ("Running approx bc")
	apprBc.run()
	dynApprBc = DynApproxBetweenness(G, epsilon, delta, True)
	print ("Running dyn approx bc with predecessors")
	dynApprBc.run()
	dynApprBc2 = DynApproxBetweenness(G, epsilon, delta, False)
	print ("Running dyn approx bc without predecessors")
	dynApprBc2.run()
	# apply the batches
	nExperiments = nEdges // batchSize
	timesBc = []
	timesDynBc = []
	timesApprBc = []
	timesDynApprBc = []
	timesDynApprBc2 = []
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
		timesBc.append(t.stop())
		# update the betweenness with the static approximated algorithm
		t = stopwatch.Timer()
		apprBc.run()
		timesApprBc.append(t.stop())
		# update the betweenness with the dynamic approximated algorithm
		t = stopwatch.Timer()
		dynApprBc.update(batch)
		timesDynApprBc.append(t.stop())
		# update the betweenness with the dynamic approximated algorithm
		t = stopwatch.Timer()
		dynApprBc2.update(batch)
		timesDynApprBc2.append(t.stop())
	return {
			"bc" : timesBc,
			"dynBc" : timesDynBc,
			"apprBc" : timesApprBc,
			"dynApprBc (predecessors stored)" : timesDynApprBc,
			"dynApprBc (predecessors not stored)" : timesDynApprBc2
			}




setNumberOfThreads(1)
G = generators.DorogovtsevMendesGenerator(100000).generate()
if (properties.numberOfComponents(G) == 1) :
	nEdges = 50
	batchSize = 10
	epsilon = 0.1
	delta = 0.1
	times = test(G, nEdges, batchSize, epsilon, delta)
	print (times)
else:
	print ("The generated graph is not connected.")
