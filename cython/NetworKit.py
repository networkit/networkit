# NetworKit Python Shell


from _NetworKit import *    # import the extension module

# standard library modules
import math
import textwrap
import collections

# non standard library modules
try:
    import networkx as nx
except ImportError:
    print("""WARNING: module 'networkx' not installed, which is required by some
                        functions.""")
try:
    import tabulate
except ImportError:
    print("""WARNING: module 'tabulate' not installed, which is required by some
                        functions. In order to install 'tabulate', Python 3.3 is required""")

# faulthandler prints debug info on errors like segfaults
try:
    import faulthandler
    faulthandler.enable()
except ImportError:
    print("""" WARNING: module 'faulthandler' not found. If installed, it will 
        provide additional debug info on errors like segmentation faults.""")

# local modules
import termgraph
import stopwatch

#--------- NetworKit Python Shell functions ----------------#

def readGraph(path, format=None):
    """    Read graph and return a NetworKit::Graph"""
    # TODO: detect file format by looking at the file content
    if format is None:    # look at file extension
        if path.endswith(".graph") or path.endswith(".metis"):
            reader = METISGraphReader()
        else:
            raise Exception("unknown graph file format")
    else:
        if (format is "metis"):
            reader = METISGraphReader()
        elif (format is "edgelist"):
            reader = EdgeListIO('\t', 1)

    with open(path, "r") as file:    # catch a wrong path before it crashes the interpreter
        G = reader.read(path)
        return G

    return None

########  CONVERSION ########

def nx2nk(nxG, weightAttr=None):
    """ 
    Convert a networkx.Graph to a NetworKit.Graph
        :param weightAttr: the edge attribute which should be treated as the edge weight
     """
    # TODO: consider weights
    n = nxG.number_of_nodes()
    nkG = Graph(n)
    
    if weightAttr is not None:
        nkG.markAsWeighted()
        for (u, v) in nxG.edges():
            w = nxG.edge[u][v][weightAttr]
            nkG.addEdge(u, v, w)
    else:
        for (u, v) in nxG.edges():
            nkG.addEdge(u, v)
    
    return nkG


def nk2nx(nkG):
    """ Convert a NetworKit.Graph to a networkx.Graph """
    nxG = nx.Graph()
    if nkG.isMarkedAsWeighted():
        for (u, v) in nkG.edges():
            nxG.add_edge(u, v, weight=nkG.weight(u, v))
    else:
        for (u, v) in nkG.edges():
            nxG.add_edge(u, v)
    return nxG
    

########  PROPERTIES ########


def nm(nkG):
    n = nkG.numberOfNodes()
    m = nkG.numberOfEdges()
    return (n, m)

def degrees(nkG):
    minMaxDeg = GraphProperties.minMaxDegree(nkG)
    avgDeg = GraphProperties.averageDegree(nkG)
    return (minMaxDeg[0], minMaxDeg[1], avgDeg)

def components(nxG):
    """Analyze connected components"""
    components = nx.connected_components(nxG)
    nComponents = len(components)
    componentSizes = [len(component) for component in components]
    componentSizes.sort(reverse=True) # sort in descending order
    sizeLargestComponent = componentSizes[0]
    return (nComponents, sizeLargestComponent)


# TODO: core decomposition


def properties(nkG, settings):
    nxG = None
    if settings["networkx"]:
        print("[...] converting to NetworX.Graph for some properties....")
        nxG = nk2nx(nkG)

    print("[...] calculating basic properties")
    
    # size
    n, m = nm(nkG)    # n and m

    # degrees 
    degDist = GraphProperties.degreeDistribution(nkG) 
    minDeg, maxDeg, avgDeg = degrees(nkG)

    # number of isolated nodes
    isolates = degDist[0] if len(degDist > 0) else None
    satellites = degDist[1] if len(dedDist > 1) else None

    # number of cliques
    cliques = len(list(nx.find_cliques(nxG))) if nxG else None


    # number of self-loops
    loops = len(nxG.selfloop_edges()) if nxG else None
    
    # density
    dens = nx.density(nxG) if nxG else None

    # diameter
    dia = None
    if settings["diameter"] and (n < 1000):
        print("calculating diameter...")
        dia = nx.diameter(nxG)


    # calculate eccentricity
    ecc = None
    if settings["eccentricity"] and (n < 1000):
        print("calculating eccentricity...")
        eccentricities = nx.eccentricity(nxG)
        ecc = sum(val for val in eccentricities.values()) / n


    # community detection

    ncomPLP, modPLP = None, None
    ncomPLM, modPLM = None, None
    if settings["communities"]:
        print("[...] detecting communities")
        # perform PLP and PLM community detection
        print("performing community detection: PLP")
        # TODO: avoid printout of status bar
        plp = LabelPropagation()
        zetaPLP = plp.run(nkG)
        ncomPLP = zetaPLP.numberOfClusters()
        modPLP = Modularity().getQuality(zetaPLP, nkG)
        print("performing community detection: PLM")
        PLM = Louvain("balanced")
        zetaPLM = PLM.run(nkG)
        ncomPLM = zetaPLM.numberOfClusters()
        modPLM = Modularity().getQuality(zetaPLM, nkG)

    # degree histogram
    
    labels, histo = None, None
    if settings["degreeDistribution"]:
        print("[...] calculating degree histogram")    
        histo = nx.degree_histogram(nxG)
        (labels, histo) = compressHistogram(histo, nbins=25)

    # connected components
    nComponents, sizeLargestComponent = None, None
    if settings["components"]:
        print("[...] finding connected components")    
        nComponents, sizeLargestComponent = components(nxG)

    # clustering
    avglcc = None
    if settings["clustering"]:
        avglcc = GraphProperties.averageLocalClusteringCoefficient(nkG)

    # degree assortativity
    assort = None
    if settings["assortativity"]:
        assort = nx.degree_assortativity_coefficient(nxG)




    # betweenness centrality
    # TODO: average betweenness centrality?

    props = {
         "name": nkG.getName(),
         "n": n,
         "m": m,
         "minDeg": minDeg,
         "maxDeg": maxDeg,
         "avgDeg": avgDeg,
         "avglcc": avglcc,
         "nComponents": nComponents,
         "sizeLargestComponent": sizeLargestComponent,
         "dia": dia,
         "ecc": ecc,
         "isolates": isolates,
         "loops": loops,
         "ncomPLP": ncomPLP,
         "modPLP": modPLP,
         "ncomPLM": ncomPLM,
         "modPLM": modPLM,
         "dens": dens,
         "assort": assort,
         "cliques": cliques,
         "histo": (labels, histo),
         }

    return props


def showProperties(nkG, settings=collections.defaultdict(lambda: True)):
    props = properties(nkG, settings)
    basicProperties = [
        ["nodes (n)", props["n"]],
        ["edges (m)", props["m"]],
        ["min. degree", props["minDeg"]],
        ["max. degree", props["maxDeg"]],
        ["avg. degree", props["avgDeg"]],
        ["isolated nodes", props["isolates"]],
        ["self-loops", props["loops"]],
        ["density", "{0:.6f}".format(props["dens"]) if props["dens"] else None]
    ]
    pathStructure = [
        ["connected components", props["nComponents"]],
        ["size of largest component", props["sizeLargestComponent"]],
        ["diameter", props["dia"]],
        ["avg. eccentricity", props["ecc"]],
    ]
    
    miscProperties = [
        ["degree assortativity", "{0:.6f}".format(props["assort"]) if props["assort"] else None],
        ["cliques", props["cliques"]]
    ]

    communityStructure = [
        ["avg. local clustering coefficient", "", "{0:.6f}".format(props["avglcc"]) if props["avglcc"] else None],
        ["PLP community detection", "", ""],
        ["", "communities", props["ncomPLP"]],
        ["", "modularity", "{0:.6f}".format(props["modPLP"]) if props["modPLP"] else None],
        ["PLM community detection", "", ""],
        ["", "communities", props["ncomPLM"]],
        ["", "modularity", "{0:.6f}".format(props["modPLM"]) if props["modPLM"] else None],
    ]

    print()
    print("Network Properties")
    print("==================")
    print("Basic Properties")
    print(tabulate.tabulate(basicProperties))
    print("Path Structure")
    print(tabulate.tabulate(pathStructure))
    print("Miscellaneous")
    print(tabulate.tabulate(miscProperties))
    print("Community Structure")
    print(tabulate.tabulate(communityStructure))
    print("Degree Distribution")
    print("-------------------")
    (labels, histo) = props["histo"]
    if labels and histo:
        termgraph.graph(labels, histo)


def compressHistogram(hist, nbins=20):
    """ Compress a histogram to a number of bins"""
    compressed = [None for i in range(nbins)]
    labels = [None for i in range(nbins)]
    nbinsprev = len(hist)
    binwidth = math.ceil(nbinsprev / nbins)
    for i in range(nbins):
        l = i * binwidth
        u = (i+1) * binwidth
        compressed[i] = sum(hist[l : u])
        labels[i] = "{0}-".format(l)
    return (labels, compressed)
        

    
def printDegreeHistogram(nxG, nbins=25):
    """ Prints a degree histogram as a bar chart to the terminal"""
    hist = nx.degree_histogram(nxG)
    
    (labels, hist) = compressHistogram(hist, nbins)
    termgraph.graph(labels, hist)
    


def hpProperties(nkG):
    """ For large graphs: get an overview of some properties"""
    print("min/max degree:")
    

# community detection

def inspectCommunities(zeta, G):
    """ Show information about communities"""
    communitySizes = zeta.clusterSizes()
    mod = Modularity().getQuality(zeta, G)
    commProps = [
        ["# communities", zeta.numberOfClusters()],
        ["min community size", min(communitySizes)],
        ["max community size", max(communitySizes)],
        ["avg. community size", sum(communitySizes) / len(communitySizes)],
        ["imbalance", zeta.getImbalance()],
        ["modularity", mod],
    ]
    print(tabulate.tabulate(commProps))
    

def evalCommunityDetection(algo, G):
    """ Evaluate a community detection algorithm """
    
    t = stopwatch.Timer()
    zeta = algo.run(G)
    t.stop()
    results = [
        ["time [s]", t.elapsed],
        ["# communities", zeta.numberOfClusters()],
        ["modularity", Modularity().getQuality(zeta, G)]
    ]
    print(tabulate.tabulate(results))


# working with attributes

def retrieveAttributes(nodes, attributes):
    """
    For a given collection of nodes and a map (or vector) indexed by node ids
    return a map from node to attribute.
    """
    attrs = dict()
    for u in nodes:
        attrs[u] = attributes[u]
    return attrs
    
    
    
