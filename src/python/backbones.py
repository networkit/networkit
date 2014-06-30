from _NetworKit import SimmelianBackbone, MultiscaleBackbone

import pygephi

def writeCsv(file, G, B):	
	f = open(file, "w+")
	f.write("Source,Target,Type,{0}\n".format("Backbone"))
	for edge in G.edges():
		value = 1 if edge in B.edges() else 0
		f.write("{0},{1},{2},{3}\n".format(edge[0],edge[1],"Undirected",value))
	f.close()

 ----------------------------------------------------------------------------------------
# THE FOLLOWING IS AN EXPERIMENTAL INTERFACE TO GEPHI USING THE GRAPH STREAMING PLUGIN.
# THIS FUNCTIONALITY IS TO BE EXTRACTED INTO A NICE AND GENERAL INTERFACE.
# ----------------------------------------------------------------------------------------
 
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
def gephi_getEdgeId(upperNodeIdBound, source, target):
	return (2*upperNodeIdBound*source) + target

# FOR TESTING PURPOSES ONLY.
def gephi_exportBackbone(G, B):
	g = NetworKitGephiClient()
	g.clean()
        
	nAttrs = {"size":10}
        
	for node in G.nodes():
		g.add_node(str(node), **nAttrs)
        
	g.flush()
    
	eAttrs = {"backbone":0, "Type":"Undirected"}
	for edge in set(G.edges()) - set(B.edges()):
		edgeId = gephi_getEdgeId(G.numberOfNodes(), edge[0], edge[1])
		g.add_edge(edgeId, edge[0], edge[1], False, **eAttrs)
        
	eAttrs = {"backbone":1, "Type":"Undirected"}
	for edge in B.edges():
		edgeId = gephi_getEdgeId(G.numberOfNodes(), edge[0], edge[1])
		g.add_edge(edgeId, edge[0], edge[1], False, **eAttrs)
            
	g.flush()
   
# FOR TESTING PURPOSES ONLY
def gephi_exportBackboneAttribute(G, B, aName):
    g = NetworKitGephiClient()

    eAttrs = {aName:0, "Type":"Undirected"}
    for edge in set(G.edges()) - set(B.edges()):
        edgeId = gephi_getEdgeId(G.numberOfNodes(), edge[0], edge[1])
        g.change_edge(edgeId, edge[0], edge[1], False, **eAttrs)
        
    eAttrs = {aName:1, "Type":"Undirected"}
    for edge in B.edges():
        edgeId = gephi_getEdgeId(G.numberOfNodes(), edge[0], edge[1])
        g.change_edge(edgeId, edge[0], edge[1], False, **eAttrs)
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