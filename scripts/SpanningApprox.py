from networkit import *
from networkit.dynamic import *
from networkit.centrality import *
import operator
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
	print("Largest component: ", largestCompo, ", size: ", size)
	compoNodes = [v for v in G.nodes() if cc.componentOfNode(v) is largestCompo]
	C = graph.Subgraph().fromNodes(G, compoNodes)
	# I've got to re-map the nodes
	nodes = {}
	i = 0
	for u in range(G.numberOfNodes()):
		if C.hasNode(u):
			nodes[u] = i
			i = i + 1
	print("Len: ", i)
	G1 = Graph(C.numberOfNodes())
	def addToG1(u,v,weight, id):
		uscaled = nodes[u]
		vscaled = nodes[v]
		G1.addEdge(uscaled, vscaled)

	C.forEdges(addToG1)
	print("Done")
	return G1


def setRandomWeights(G, mu, sigma):
	"""
	Add random weights, normal distribution with mean mu and standard deviation sigma
	"""
	for (u, v) in G.iterEdges():
		w = random.normalvariate(mu, sigma)
		G.setWeight(u, v, w)
	return G



def test(G, T, nEdges, batchSize, epsilon, delta, size):
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

def readGraphUnweighted(file, directed_graph):
	print("Opening file")
	f = open(file, "r")
	print("File opened")
	n = 0
	minNode = 0
	i = 0
	#first scan to find out the number of nodes
	for line in f:
		fields = line.strip().split()
		if fields[0].startswith("%") or fields[0].startswith("#"):
			continue
		(u, v) = (int(fields[0]), int(fields[1]))
		if i == 0:
			minNode = u
		if u > n:
			n = u
		if v > n:
			n = v
		if v < minNode:
			minNode = v
		if u < minNode:
			minNode = u
		i = i + 1
	f.seek(0,0)
	n = n + 1
	print("number of nodes: "+str(n))
	print("i : ",i)
	print("min node: ", minNode)
	IDs=[];
	for k in range(n):
		IDs.append(False)
	for line in f:
		fields = line.strip().split()
		if fields[0].startswith("%") or fields[0].startswith("#"):
			continue
		(u, v) = (int(fields[0]), int(fields[1]))
		IDs[u] = True
		IDs[v] = True
		#print(u,v,IDs[u],IDs[v])
	isolatedSoFar = []
	counter = 0
	for k in range(n):
		if (k >= minNode and IDs[k] == False):
			counter = counter + 1
			isolatedSoFar.append(counter)
		else:
			isolatedSoFar.append(counter)
	print("Counter: ", counter, ", n = ", n)
	# we create a (unweighted) graph with the connections and a weighted graph with the timestamps
	if not directed_graph:
		G = Graph(n - counter - minNode)
	else:
		G = Graph(n - counter - minNode, False, True)
	f.seek(0,0)
	for line in f:
		fields = line.strip().split()
		if fields[0].startswith("%") or fields[0].startswith("#"):
			continue
		(u, v) = (int(fields[0]), int(fields[1]))
		u = u - minNode - isolatedSoFar[u]
		v = v - minNode - isolatedSoFar[v]
		if (u < 0 or v < 0):
			print("u = ", u, ", v = ", v, ", isolatedSoFar[u] = ", isolatedSoFar[u], ", isolatedSoFar[v] = ", isolatedSoFar[v], ", minode = ", minNode)
		if (u == v or G.hasEdge(u,v)):
			continue
		if G.hasEdge(v,u) and not directed_graph:
			continue
		G.addEdge(u, v)
	print("number of nodes: ", G.numberOfNodes())
	print ("number of edges: ", G.numberOfEdges())
	print("All the edges have been scanned")
	return G



if __name__ == "__main__":
	setNumberOfThreads(1)

	#graph_name = "road_usa.metis.graph"
	path = "../input/"
#	undirected = ["advogato.txt", "Drosophila_melanogaster.txt", "Caenorhabditis_elegans.txt", "ca-GrQc.txt", "ca-HepTh.txt","dip20090126_MAX.txt",
#	undirected = ["as-skitter_comp.graph", "cit-Patents_comp.graph", "com-youtube.ungraph.graph","oregon1_010526.graph","p2p-Gnutella31_comp.graph","soc-Epinions1_comp.graph",
#"CA-GrQc_comp.graph","com-amazon.ungraph.graph","hollywood2009.graph",
#	undirected = ["roadNet-TX_comp.graph","Wiki-Vote_comp.graph","com-dblp.ungraph.graph","LiveJournal.graph","Slashdot0902.graph","orkut.graph"]
	#undirected = ["PGPgiantcompo.graph","PGP2.graph","PGP3.graph","PGP4.graph","PGP5.graph"]
#	undirected = ["mus_musc2.graph","mus_musc3.graph","mus_musc4.graph","mus_musc5.graph","mus_musc6.graph","mus_musc7.graph","mus_musc8.graph","mus_musc9.graph","mus_musc10.graph"]
	undirected = ["PGPgiantcompo.graph"]
	nsamples = [50, 100, 200, 500, 1000]
	epsilons = [0.1, 0.2, 0.5, 0.05]

	for graph_name in undirected:
		pdDict = {}
		G = readGraph(path+graph_name, fileformat=Format.METIS)
		G = extractLargestComponent(G)
		G.indexEdges()
		print("Extracted largest connected component of the graph.")
		n = G.numberOfNodes()
		m = G.numberOfEdges()
		print("Number of nodes: "+str(n))
		print("Number of edges: ", m)
		# exact SEC
		t = stopwatch.Timer()
		sp = Spanning(G)
		sp.run()
		time_sp_exact = t.stop()
		sp_exact = sp.scores()
		time_sp_jlt = []
		sp_jlt = []
		for epsilon in epsilons:
			# approx SEC - JLT
			t = stopwatch.Timer()
			sp = Spanning(G, epsilon)
			sp.runApproximation()
			x = t.stop()
			time_sp_jlt.append(x)
			pdDict["TimeJLT_"+str(epsilon)] = pd.Series(x)
			sp_jlt.append(sp.scores())
			pdDict["JLT_"+str(epsilon)] = pd.Series(sp.scores())
		# approx SEC - Random
		time_sp_random1 = []
		sp_random1 = []
		time_sp_random2 = []
		sp_random2 = []
		time_sp_random3 = []
		sp_random3 = []
		for samp in nsamples:
			t = stopwatch.Timer()
			sp = Spanning(G)
			sp.runTreeApproximation(samp)
			x = t.stop()
			time_sp_random1.append(x)
			pdDict["TimeRandom1_"+str(samp)] = pd.Series(x)
			sp_random1.append(sp.scores())
			pdDict["Random1_"+str(samp)] = pd.Series(sp.scores())
			# approx SEC - Random
			t = stopwatch.Timer()
			sp = Spanning(G)
			sp.runTreeApproximation2(samp)
			x = t.stop()
			time_sp_random2.append(x)
			pdDict["TimeRandom2_"+str(samp)] = pd.Series(x)
			sp_random2.append(sp.scores())
			pdDict["Random2_"+str(samp)] = pd.Series(sp.scores())
			# approx SEC - Random
			t = stopwatch.Timer()
			sp = Spanning(G)
			sp.runPseudoTreeApproximation(samp)
			x = t.stop()
			time_sp_random3.append(x)
			pdDict["TimeRandom3_"+str(samp)] = pd.Series(x)
			sp_random3.append(sp.scores())
			pdDict["Random3_"+str(samp)] = pd.Series(sp.scores())

		pdDict["Exact"] = pd.Series(sp_exact)
		pdDict["TimeExact"] = pd.Series(time_sp_exact)

		df1 = pd.DataFrame(pdDict)
		df1.to_csv(graph_name+"_approx")
