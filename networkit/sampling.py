""" Sampling from graphs """

__author__ = "Elisabetta Bergamini"

def bfsSample(G, source=None, k = 50):
    """ Start a BFS from source node, return node-induced subgraph of the first k nodes discovered"""
    if not source:
        source = G.randomNode()
    n = G.numberOfNodes()
    visited = [False]*n
    Q = [source]
    closest = set([source])
    global found
    found = 0
    while len(Q) > 0 and found < k:
        u = Q.pop(0)
        def enqueue(u,v,weight, eid):
            global found
            if not visited[v] and found < k:
                found += 1
                visited[v] = True
                Q.append(v)
                closest.add(v)
        G.forEdgesOf(u, enqueue)
    print("found {0} nodes".format(len(closest)))
    G1 = G.subgraphFromNodes(closest)
    return G1
