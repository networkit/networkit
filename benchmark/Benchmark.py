"""
The NetworKit benchmark

This module implements a comprehensive benchmark of NetworKit's analytics algorithms



"""


import pandas
import sys
import warnings
import time
import math

# framework imports
from networkit import *

try:
    import networkx
except ImportError:
    error("networkx not available")

try:
    import igraph
except ImportError:
    error("igraph not available")


# helper function

def loadGraph(key, basePath, framework="networkit"):
    if framework is "networkit":
        (fileName, formatName) = networks[key]
        G = readGraph(os.path.join(basePath, fileName), formatName)
        return G
    else if framework is "networkx":
        pass
    else if framework is "igraph":
        pass
    else:
        raise Exception("unknown framework {0}".format(framework))


class Timer:
    """ Use the Python with-statement to time your code
    with this timer. """

    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.elapsed = round(self.end - self.start, 6)



# what is a  test

# get graph
# init
# multiple times
#       run
#       take time
#


# settings

nRuns = 5   # how many runs for representative results

# collection of networks

networks = {
            "PGPgiantcompo" : ("PGPgiantcompo.metis.graph", Format.METIS),
            "power" : ("power.metis.graph", Format.METIS),
            "caidaRouterLevel" : ("caidaRouterLevel.metis.graph", Format.METIS),
            "as-22july06" : ("as-22july06.metis.graph", Format.METIS),
            "coAuthorsDBLP" : ("coAuthorsDBLP.graph", Format.METIS),
            "uk-2007-05" : ("uk2007-05.metis.graph", Format.METIS),
            "uk-2002" : ("uk-2002.metis.graph", Format.METIS),
            "fb-Texas84" : ("Texas84.edgelist", Format.EdgeListTabZero),
            "fb-Caltech36" : ("Caltech36.edgelist", Format.EdgeListTabZero),
            "fb-MIT8" : ("MIT8.edgelist", Format.EdgeListTabZero),
            "fb-Smith60" : ("Smith60.edgelist", Format.EdgeListTabZero),
            "con-fiber_big" : ("con-fiber_big.metis.graph", Format.METIS),
            }




selected = ["PGPgiantcompo", "power"]
collectionDir = os.path.expanduser("~/workspace/Data/NwkBenchmark")




class nk:
    """networkit namespace"""
# - connected components (properties.ConnectedComponents, properties.ParallelConnectedComponents)

    class Algo:
        """ runner for an algorithm"""
        def run(self, G):
            raise Exception("Not implemented")

        def loadGraph(key, basePath):
            (fileName, formatName) = networks[key]
            G = readGraph(os.path.join(basePath, fileName), formatName)
            return G

    class bConnectedComponents(nk.Algo):
        name = "ConnectedComponents"

        def run(self, G):
            cc = properties.ConnectedComponents(G)
            cc.run()

    class bParallelConnectedComponents(nk.Algo):
        name = "ParallelConnectedComponents"

        def run(self, G):
            cc = properties.ParallelConnectedComponents(G)
            cc.run()

    # - k-core decomposition (properties.CoreDecomposition)

    class bCoreDecomposition(nk.Algo):
        name = "CoreDecomposition"

        def run(self, G):
            cd = properties.CoreDecomposition(G)
            cd.run()

    # - degree distribution power-law estimation (properties.powerLawExponent)

    class bPowerLaw(nk.Algo):
        name = "PowerLaw"

        def run(self, G):
            properties.powerLawExponent(G)

    # - degree assortativity (properties.degreeAssortativity)

    class bDegreeAssortativity(nk.Algo):
        name = "DegreeAssortativity"

        def run(self, G):
            properties.degreeAssortativity(G)


    # - BFS & Dijkstra (graph.BFS, graph.Dijkstra)
    class bBFS(nk.Algo):
        name = "BFS"

        def run(self, G):
            bfs = graph.BFS(G)
            bfs.run()


    # - community detection (community.PLM, community.PLP)

    class bCommunityDetection(nk.Algo):
        name = "CommunityDetection"

        def run(self, G):
            plm = community.PLM(G)
            plm.run()

    # - diameter, exact (properties.Diameter.exactDiameter) and estimate (properties.Diameter.estimatedDiameterRange)

    class bDiameter(nk.Algo):
        name = "Diameter"

        def run(self, G):
            d = properties.Diameter.exactDiameter(G)


    class bDiameterEstimate(nk.Algo):
        name = "Diameter"

        def run(self, G):
            d = properties.Diameter.estimatedDiameterRange(G)

    # - clustering coefficients (average local), exact (properties.ClusteringCoefficient.avgLocal) and approximated (properties.ClusteringCoefficient.approxAvgLocal)

    class bClusteringCoefficient(nk.Algo):
        name = "ClusteringCoefficient"

        def run(self, G):
            c = properties.ClusteringCoefficient.avgLocal(G)

    class bApproxClusteringCoefficient(nk.Algo):
        name = "ApproxClusteringCoefficient"

        def run(self, G):
            # TODO: specify number of trials
            c = properties.ClusteringCoefficient.approxAvgLocal(G)



    # - centrality

    # 	- PageRank (centrality.PageRank, centrality.SciPyPageRank)

    class bPageRank(nk.Algo):
        name = "PageRank"

        def run(self, G):
            pr = centrality.PageRank(G)
            pr.run()

    # 	- Eigenvector centrality (centrality.EigenvectorCentrality, centrality.SciPyEVZ)


    # 	- betweenness,  exact (centrality.Betweenness) and approximated (centrality.ApproxBetweenness, centrality.ApproxBetweenness2)

    class bBetweenness(nk.Algo):
        name = "Betweenness"

        def run(self, G):
            bc = centrality.Betweenness(G)
            bc.run()


    class bApproxBetweenness(nk.Algo):
        name = "ApproxBetweenness"

        def run(self, G):
            bc = centrality.ApproxBetweenness(G, epsilon=0.05, delta=0.1)
            bc.run()

#------------------------------------



class nx:
    """networkx namespace"""

    class Algo:
        """ runner for an algorithm"""
        def run(self, G):
            raise Exception("Not implemented")

# - connected components (properties.ConnectedComponents, properties.ParallelConnectedComponents)

    class bConnectedComponents(nx.Algo):
        name = "ConnectedComponents"

        def run(self, G):
            cc = networkx.connected_components(G)

    # - k-core decomposition (properties.CoreDecomposition)

    class bCoreDecomposition(nx.Algo):
        name = "CoreDecomposition"

        def run(self, G):
            cn = networkx.core_number(G)

    # - degree distribution power-law estimation (properties.powerLawExponent)

    # not available

    # - degree assortativity (properties.degreeAssortativity)

    class bDegreeAssortativity(nx.Algo):
        name = "DegreeAssortativity"

        def run(self, G):
            ac = networkx.degree_assortativity_coefficient(G)


    # - BFS & Dijkstra (graph.BFS, graph.Dijkstra)
    class bBFS(nx.Algo):
        name = "BFS"

        def run(self, G):
            networkx.bfs_predecessors(G)


    # - community detection (community.PLM, community.PLP)
    # not available

    # - diameter, exact (properties.Diameter.exactDiameter) and estimate (properties.Diameter.estimatedDiameterRange)

    class bDiameter(nx.Algo):
        name = "Diameter"

        def run(self, G):
            d = networkx.diameter(G)


    # approximate diameter not available

    # - clustering coefficients (average local), exact (properties.ClusteringCoefficient.avgLocal) and approximated (properties.ClusteringCoefficient.approxAvgLocal)

    class bClusteringCoefficient(nx.Algo):
        name = "ClusteringCoefficient"

        def run(self, G):
            c = properties.ClusteringCoefficient.avgLocal(G)

    class bApproxClusteringCoefficient(nx.Algo):
        name = "ApproxClusteringCoefficient"

        def run(self, G):
            # TODO: number of trials
            c = networkx.average_clustering(G, trials)



    # - centrality

    # 	- PageRank (centrality.PageRank, centrality.SciPyPageRank)

    class bPageRank(nx.Algo):
        name = "PageRank"

        def run(self, G):
            # TODO: check parameters
            pr = networkx.pagerank(G)

    # 	- Eigenvector centrality (centrality.EigenvectorCentrality, centrality.SciPyEVZ)


    # 	- betweenness,  exact (centrality.Betweenness) and approximated (centrality.ApproxBetweenness, centrality.ApproxBetweenness2)

    class bBetweenness(nx.Algo):
        name = "Betweenness"

        def run(self, G):
            networkx.betweenness_centrality(G)


    # approximation not available

#-------------------------



class ig:
    """ igraph namespace"""

    class Algo:
        """ runner for an algorithm"""
        def run(self, G):
            raise Exception("Not implemented")

        def loadGraph(key, basePath):
            raise NotImplementedError("TODO")

    class bConnectedComponents(ig.Algo):
        name = "ConnectedComponents"

        def run(self, G):
            c = igraph.Graph.components(G)

# - k-core decomposition (properties.CoreDecomposition)

class bCoreDecomposition(ig.Algo):
    name = "CoreDecomposition"

    def run(self, G):
        igraph.Graph.coreness(G)

# - degree distribution power-law estimation (properties.powerLawExponent)

class bPowerLaw(ig.Algo):
    name = "PowerLaw"

    def run(self, G):
        raise NotImplementedError("TODO")

# - degree assortativity (properties.degreeAssortativity)

class bDegreeAssortativity(ig.Algo):
    name = "DegreeAssortativity"

    def run(self, G):
        raise NotImplementedError("TODO")


# - BFS & Dijkstra (graph.BFS, graph.Dijkstra)
class bBFS(ig.Algo):
    name = "BFS"

    def run(self, G):
        raise NotImplementedError("TODO")


# - community detection (community.PLM, community.PLP)

class bCommunityDetection(ig.Algo):
    name = "CommunityDetection"

    def run(self, G):
        raise NotImplementedError("TODO")

# - diameter, exact (properties.Diameter.exactDiameter) and estimate (properties.Diameter.estimatedDiameterRange)

class bDiameter(ig.Algo):
    name = "Diameter"

    def run(self, G):
        raise NotImplementedError("TODO")


class bDiameterEstimate(ig.Algo):
    name = "Diameter"

    def run(self, G):
        raise NotImplementedError("TODO")

# - clustering coefficients (average local), exact (properties.ClusteringCoefficient.avgLocal) and approximated (properties.ClusteringCoefficient.approxAvgLocal)

class bClusteringCoefficient(ig.Algo):
    name = "ClusteringCoefficient"

    def run(self, G):
        raise NotImplementedError("TODO")

class bApproxClusteringCoefficient(ig.Algo):
    name = "ApproxClusteringCoefficient"

    def run(self, G):
        # TODO: specify number of trials
        raise NotImplementedError("TODO")



# - centrality

# 	- PageRank (centrality.PageRank, centrality.SciPyPageRank)

class bPageRank(ig.Algo):
    name = "PageRank"

    def run(self, G):
        raise NotImplementedError("TODO")

# 	- Eigenvector centrality (centrality.EigenvectorCentrality, centrality.SciPyEVZ)


# 	- betweenness,  exact (centrality.Betweenness) and approximated (centrality.ApproxBetweenness, centrality.ApproxBetweenness2)

class bBetweenness(ig.Algo):
    name = "Betweenness"

    def run(self, G):
        raise NotImplementedError("TODO")

class bApproxBetweenness(ig.Algo):
    name = "ApproxBetweenness"

    def run(self, G):
        raise NotImplementedError("TODO")


#-------------------------------------------------

class bFail(Algo):
    name = "Fail"

    def run(self, G):
        raise Exception("FAIL!")


# Logging

def info(message):
    print(message)

def error(message):
    print(message)

def debug(message):
    pass
    # print(message)

def main():
    logging.info("start benchmark")

    data = benchmark(ConnectedComponents_(), ["PGPgiantcompo", "power"])
    data = benchmark(CoreDecomposition_(), ["PGPgiantcompo", "power"])
    #data = benchmark(_Fail(), ["PGPgiantcompo", "power    "])


def benchmark(algo, graphs):
    info("benchmarking {algo.name}".format(**locals()))
    table = []  # list of dictionaries, to be converted to a DataFrame

    for graphName in graphs:
        try:
            info("loading {0}".format(graphName))
            G = loadGraph(graphName, basePath=collectionDir)
            try:
                for i in range(nRuns):
                    row = {}    # benchmark data row
                    with Timer() as t:
                        debug("running {algo.name}".format(**locals()))
                        algo.run(G)
                    debug("took {0} s".format(t.elapsed))
                    # store data
                    row["algo"] = algo.name
                    row["graph"] = graphName
                    row["time"] = t.elapsed
                    table.append(row)
            except Exception as ex:
                error("algorithm {algo.name} failed with exception: {ex}".format(**locals()))
        except Exception as ex:
            error("loading graph {graphName} failed with exception: {ex}".format(**locals()))

    return pandas.DataFrame(table)


# - generators
# 	- Erdös-Renyi (generators.ErdosRenyiGenerator)
# 	- Barabasi-Albert (generators.BarabasiAlbertGenerator)
# 	- Chung-Lu (generators.ChungLuGenerator)
# 	- RMAT (generators.RmatGenerator)

if __name__ == "__main__":
    main()
