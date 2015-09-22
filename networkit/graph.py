# extension imports
from _NetworKit import Graph, BFS, Dijkstra, Subgraph, DynBFS, DynDijkstra, SpanningForest, GraphTools

from enum import Enum

def findCycle(G):
	Color = Enum("Color", "white gray black")
	color = {}
	for u in G.nodes():
		color[u] = Color.white
	S = list()
	s = G.randomNode()
	S.append(s)
	while len(S) > 0:
		u = S.pop()
		print(u)
		if color[u] is Color.white:
			color[u] = Color.gray
			for v in G.neighbors(u):
				S.append(v)
