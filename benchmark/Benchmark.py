"""
The NetworKit benchmark

This module implements a comprehensive benchmark of NetworKit's analytics algorithms



"""


import pandas
import sys
import warnings
import time
import math
import os
import numpy
import matplotlib.pyplot as plt
import seaborn


import networkit

import nk
import nx
import ig

# helper function



def averageRuns(df, groupby=["graph"]):
    """ Average running time, modularity, edges per second and number of clusters over multiple runs"""
    return df.groupby(groupby, as_index=False).mean()


def graphMeta(graphNames, graphDir):
    meta = []
    for name in graphNames:
        G = networkit.readGraph(os.path.join(graphDir, "{0}.gml.graph".format(name)), networkit.Format.GML)
        (n, m) = networkit.properties.size(G)
        meta.append({"name" : name, "n" : n, "m" : m})
    return pandas.DataFrame(meta, columns=["name", "n", "m"])


class Timer:
    """ Use the Python with-statement to time your code
    with this timer. """

    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.elapsed = round(self.end - self.start, 6)



class Bench:

    def __init__(self):
        self.dataPath = "benchData-{0}".format(time.time()) # TODO: timestamp this
        os.mkdir(self.dataPath)


def saveData(df, name):
    df.to_csv(os.path.join(bench.dataPath, "{0}.csv".format(name)), sep="\t")


# settings

nRuns = 5   # how many runs for representative results

# collection of networks

networks = {
            "PGPgiantcompo" : ("PGPgiantcompo.metis.graph", networkit.Format.METIS),
            "power" : ("power.metis.graph", networkit.Format.METIS),
            "caidaRouterLevel" : ("caidaRouterLevel.metis.graph", networkit.Format.METIS),
            "as-22july06" : ("as-22july06.metis.graph", networkit.Format.METIS),
            "coAuthorsDBLP" : ("coAuthorsDBLP.graph", networkit.Format.METIS),
            "uk-2007-05" : ("uk2007-05.metis.graph", networkit.Format.METIS),
            "uk-2002" : ("uk-2002.metis.graph", networkit.Format.METIS),
            "fb-Texas84" : ("Texas84.edgelist", networkit.Format.EdgeListTabZero),
            "fb-Caltech36" : ("Caltech36.edgelist", networkit.Format.EdgeListTabZero),
            "fb-MIT8" : ("MIT8.edgelist", networkit.Format.EdgeListTabZero),
            "fb-Smith60" : ("Smith60.edgelist", networkit.Format.EdgeListTabZero),
            "con-fiber_big" : ("con-fiber_big.metis.graph", networkit.Format.METIS),
            }




selected = ["PGPgiantcompo", "power"]
collectionDir = os.path.expanduser("~/workspace/Data/NwkBenchmark")




class bFail:
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

# Plots

## plot settings

seaborn.set_style("whitegrid")

### Colors
red = seaborn.xkcd_rgb["crimson"]
orange = seaborn.xkcd_rgb["bright orange"]


def timePlot(data, size=(6,3)):
    pos = numpy.arange(len(data))+.5    # the bar centers on the y axis
    labels = list(data["graph"])
    plt.figure(figsize=size)
    plt.xscale("symlog")
    plt.barh(pos, data["time"], align='center', height=0.25, color=red)    # notice the 'height' argument
    plt.yticks(pos, labels)
    plt.gca().xaxis.set_minor_locator(plt.LogLocator(subs=[0,1,2,3,4,5,6,7,8,9,10]))
    #gca().xaxis.set_minor_formatter(FormatStrFormatter("%.2f"))
    plt.xlabel("time [s]")
    plt.grid(True)


def epsPlot(data, size=(6,3)):
    pos = numpy.arange(len(data))+.5    # the bar centers on the y axis
    labels = list(data["graph"])
    plt.figure(figsize=size)
    plt.xscale("log")
    plt.barh(pos, data["time"], align='center', height=0.25, color=red)    # notice the 'height' argument
    plt.yticks(pos, labels)
    plt.gca().xaxis.set_minor_locator(plt.LogLocator(subs=[0,1,2,3,4,5,6,7,8,9,10]))
    #gca().xaxis.set_minor_formatter(FormatStrFormatter("%.2f"))
    plt.xlabel("time [s]")
    plt.grid(True)



def graphPlot(data, size=None):
    pos = arange(len(data))    # the bar centers on the y axis
    plt.figure(figsize=(5,3.3))
    sizes = data["m"]
    labels = list(data["graph"])
    #barh(pos, data["n"], log=True, align='center', height=0.4, color="darkgrey")    # notice the 'height' argument
    plt.barh(pos, sizes, log=True, align='center', height=0.4, color="lightblue", label="m")    # notice the 'height' argument
    plt.yticks(pos, labels)
    plt.xlabel("number of edges")
    plt.grid(True)


def barPlot(data, labels, x_label="", y_label="", size=None, color="b", transparency=1.0, scale="linear"):
    pos = numpy.arange(len(data))+.5    # the bar centers on the y axis
    plt.figure(figsize=size)
    plt.xscale(scale)
    plt.barh(pos, data, align='center', height=0.25, color=color, alpha=transparency)    # notice the 'height' argumen
    plt.yticks(pos, labels)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.grid(True)


def algoBenchmark(algo, graphs):
    info("benchmarking {algo.name}".format(**locals()))
    table = []  # list of dictionaries, to be converted to a DataFrame

    for graphName in graphs:
        try:
            info("loading {0}".format(graphName))
            G = algo.loadGraph(os.path.join(collectionDir, "{0}.gml.graph".format(graphName)))
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


def generatorBenchmark(generator, argtuples):
    pass


# - generators
# 	- Erdös-Renyi (generators.ErdosRenyiGenerator)
# 	- Barabasi-Albert (generators.BarabasiAlbertGenerator)
# 	- Chung-Lu (generators.ChungLuGenerator)
# 	- RMAT (generators.RmatGenerator)
