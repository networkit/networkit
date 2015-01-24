import networkx as nx
import time
G = nx.read_gml('input/power.gml')
#G = nx.read_gml('input/PGPgiantcompo.gml')
#G = nx.read_gml('input/as-22july06.gml')
#G = nx.read_gml('input/fe_pwt.gml')
#G = nx.read_gml('input/bcsstk30.gml')
start = time.time()
nx.betweenness_centrality(G, None, False, None, False, None)
nodebet = time.time() - start
print(nodebet)
start2 = time.time()
nx.edge_betweenness_centrality(G, None, None)
edgebet = time.time() - start2
print(edgebet)
