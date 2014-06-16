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
	communityIds = [i for i in range(len(communitySizes)) if communitySizes[i] > 0]



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


def toLocations(urls):
	""" Turn the list of urls into a list of locations"""
	errors = 0	# count the parser errors
	locs = []
	for url in urls:
		try:
			parsed = urllib.parse.urlparse(url)
			locs.append(parsed.netloc)
		except:
			print("URL parser error occurred")
			errors += 1
			locs.append(None)
	return locs

def getLocations(nodes, urls):
	""" Given a collection of nodes, return a set of net locations (domains) to which they belong"""
	theurls = dict((u, urls[u]) for u in nodes)
	loclist = [urllib.parse.urlparse(url).netloc for url in theurls]


def matchAndIndex(substring, ls, exact=False):
	matching = {}
	if exact:
		i = 0
		for s in ls:
			if substring == s:
				matching[i] = s
			i += 1
	else:
		i = 0
		for s in ls:
			if substring in s:
				matching[i] = s
			i += 1
	return matching


def writeSeeds(match, filename):
	file = open(filename, "w")
	string = ",".join([str(s) for s in match.keys()])
	file.write(string)
	file.close()


def matchLocs(substring, ls):
	match = matchAndIndex(substring, ls)
	locs = set(match.values())
	return locs