from NetworKit import *
import urllib.parse
import collections


def analyzeWebCommunities(graphPath, urlsPath):
	print("reading input...")
	G = readGraph(graphPath)
	urls = LineFileReader().read(urlsPath)
	urlmap = retrieveAttributes(G.nodes(), urls)
	print("community detection...")
	zeta = LabelPropagation().run(G)

	communitySizes = zeta.clusterSizes()
	communityIds = [i for i in range(len(communitySizes)) if communitySizes[i] > 0]



	netlocMap = {}	# community -> size distribution of netlocs
	
	for cid in communityIds:
		if communitySizes[cid] > 100: # filter 
		community = zeta.getMembers(cid)
		urllib.parse.urlparse(v).netloc


def getNetloc2nodes(G, urls):
	netlocs = [urllib.parse.urlparse(url).netloc for url in urls]
	netloc2nodes = {}
	for u in range(netlocs):
		uloc = netlocs[u]
		if uloc not in netloc2nodes:
			netlocs2nodes[uloc] = []
		netloc2nodes[uloc].append(u)
	return netloc2nodes