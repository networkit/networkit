from _NetworKit import SimmelianBackbone, ChibaNishizekiTriangleCounter

import pygephi

def writeCsv(file, G, B):	
	f = open(file, "w+")
	f.write("Source,Target,Type,{0}\n".format("Backbone"))
	for edge in G.edges():
		value = 1 if edge in B.edges() else 0
		f.write("{0},{1},{2},{3}\n".format(edge[0],edge[1],"Undirected",value))
	f.close()

class NetworKitGephiClient:
    """ A python singleton """

    # storage for the instance reference
    __instance = None

    def __init__(self):
        #Disabling Flushing means a good performance boost.
        if NetworKitGephiClient.__instance is None:
            NetworKitGephiClient.__instance = pygephi.GephiClient('http://localhost:8080/workspace0', autoflush=False)

        # Store instance reference as the only member in the handle
        self.__dict__['_GephiInterface__instance'] = NetworKitGephiClient.__instance
        
        """ Delegate accesseso implementation """
    def __getattr__(self, attr):  
        return getattr(self.__instance, attr)

    def __setattr__(self, attr, value):       
        return setattr(self.__instance, attr, value)

# FOR TESTING PURPOSES ONLY.
def gephi_exportBackbone(G, B):
    g = NetworKitGephiClient()
    g.clean()
        
    nAttrs = {"size":10}
        
    for node in G.nodes():
        g.add_node(str(node), **nAttrs)
        
    g.flush()
    
    edgeId = 1
    eAttrs = {"backbone":0}
    for edge in set(G.edges()) - set(B.edges()):
        g.add_edge(edgeId, edge[0], edge[1], False, **eAttrs)
        edgeId = edgeId + 1
        
    eAttrs = {"backbone":1}
    for edge in B.edges():
        g.add_edge(edgeId, edge[0], edge[1], False, **eAttrs)
        edgeId = edgeId + 1
            
    g.flush()
    
# FOR TESTING PURPOSES ONLY.
def gephi_exportGraph(G):
    g = NetworKitGephiClient()
    g.clean()
    
    nAttrs = {"size":10}
        
    for node in G.nodes():
        g.add_node(str(node), **nAttrs)
    
    edgeId = 1
    for edge in G.edges():
        g.add_edge(edgeId, edge[0], edge[1], False)
        edgeId = edgeId + 1
    
    g.flush()

# FOR TESTING PURPOSES ONLY.
def gephi_clean():
    g = NetworKitGephiClient()
    g.clean()