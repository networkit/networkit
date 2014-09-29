"""
The NetworKit benchmark

This module implements a comprehensive benchmark of NetworKit's analytics algorithms



"""


import pandas
import sys
import warnings
import time
import math


import nk
import nx
import ig


# helper function

def loadGraph(key, basePath, framework="networkit"):
    if framework is "networkit":
        (fileName, formatName) = networks[key]
        G = readGraph(os.path.join(basePath, fileName), formatName)
        return G
    elif framework is "networkx":
        pass
    elif framework is "igraph":
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
