"""
The NetworKit benchmark

This module implements a comprehensive benchmark of NetworKit's analytics algorithms



"""

from networkit import *
import pandas

# timing
#t = stopwatch.Timer()
#t.elapsed
#t.stop()


class Benchmark_ConnectedComponents(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def run():


# helper function

def loadGraph(key):
    (fileName, formatName) = networks[key]
    G = readGraph(os.path.join(basePath, fileName), formatName)
    return G

class Timing:

    def __enter__(self):
        print("enter called")
        self.watch = stopwatch.Timer()

    def __exit__(self, typ, value, traceback):
        print("exit called")
        self.watch.stop()


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
            "PGPgiantcompo" : ("PGPgiantcompo.graph", Format.METIS),
            "uk2007" : ("uk2007.graph", Format.METIS),
            }



# - connected components (properties.ConnectedComponents, properties.ParallelConnectedComponents)

selected = ["PGPgiantcompo", "uk2007"]



def main():
    for selectedGraphName in selected:
        G = loadGraph(selectedGraphName)
        connectedComponents = properties.ConnectedComponents(G)

# - degree distribution power-law estimation (properties.powerLawExponent)
# - k-core decomposition (properties.CoreDecomposition)
# - BFS & Dijkstra (graph.BFS, graph.Dijkstra)
# - clustering coefficients (average local), exact (properties.ClusteringCoefficient.avgLocal) and approximated (properties.ClusteringCoefficient.approxAvgLocal)
# - degree assortativity (properties.degreeAssortativity)
# - diameter, exact (properties.Diameter.exactDiameter) and estimate (properties.Diameter.estimatedDiameterRange)
# - community detection (community.PLM, community.PLP)
# - centrality
# 	- PageRank (centrality.PageRank, centrality.SciPyPageRank)
# 	- Eigenvector centrality (centrality.EigenvectorCentrality, centrality.SciPyEVZ)
# 	- betweenness,  exact (centrality.Betweenness) and approximated (centrality.ApproxBetweenness, centrality.ApproxBetweenness2)
# - generators
# 	- Erdös-Renyi (generators.ErdosRenyiGenerator)
# 	- Barabasi-Albert (generators.BarabasiAlbertGenerator)
# 	- Chung-Lu (generators.ChungLuGenerator)
# 	- RMAT (generators.RmatGenerator)
