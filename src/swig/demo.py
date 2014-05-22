from NetworKit import *

G = Graph(5)
D = DirectedGraph(5)

G.addEdge(0, 1)
G.addEdge(0, 2)

print('minMaxDegree(G) =', GraphProperties_minMaxDegree(G))
print('minMaxDegree(D) =', GraphProperties_minMaxDegree(D))
