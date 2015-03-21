from networkit import *
import time

def printInterval(name, _t1, _t2):
	print("[" + name + " took " + str(_t2 - _t1) + "s]")

t_0 = time.time()
G = readGraph("../../scripts/BackboneEvaluation/input/Tennessee95.graphml", Format.GraphML)
G.indexEdges()
print(properties.size(G))
t_1 = time.time()
printInterval('Loading graph file', t_0, t_1)

#c++
ba = backbones.SimmelianBackboneNonParametric(0.5)
B1 = ba.calculate(G)
print(properties.size(B1))
t_2 = time.time()
printInterval('Direct method', t_1, t_2)

#python
#attributizer = backbones.LocalSimilarityAttributizer()
#attribute = bat_ls.getAttribute(G, [0] * G.upperEdgeIdBound())
#gf = backbones.GlobalThresholdFilter(0.5, False)

chiba = backbones.ChibaNishizekiTriangleCounter()
triangles = chiba.getAttribute(G)
sj = backbones.SimmelianJaccardAttributizer()
attribute = sj.getAttribute(G, triangles)
gf = backbones.GlobalThresholdFilter(0.5, True)

B2 = gf.calculate(G, attribute)
print(properties.size(B2))
t_3 = time.time()
printInterval('Python indirection', t_2, t_3)


