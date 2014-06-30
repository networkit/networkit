import gephi.pyclient as _gephipyclient   # we want to hide these from the user
import urllib as _urllib

""" 
This module provides methods to export data from NetworKit directly to Gephi using the
Gephi Streaming Plugin.
"""

class _GephiStreamingClient:
    """ 
     A python singleton providing access to a running Master instance of the Gephi
     Streaming plugin 
    """
    
    __instance = None

    def __init__(self):
        #Disabling Flushing means quite a good performance boost.
        if _GephiStreamingClient.__instance is None:
            _GephiStreamingClient.__instance = _gephipyclient.GephiClient('http://localhost:8080/workspace0', autoflush=False)

        # Store instance reference as the only member in the handle
        self.__dict__['_GephiInterface__instance'] = _GephiStreamingClient.__instance
        
    """ Delegate accesses implementation """
    def __getattr__(self, attr):
        return getattr(self.__instance, attr)

    def __setattr__(self, attr, value):       
        return setattr(self.__instance, attr, value)

def _urlError(e):
    print("Could not connect to the gephi streaming plugin. Did you start the streaming master server in gephi?")

def exportGraph(graph):
    """ Exports the given graph to gephi. No attributes or weights are exported. """
    
    try:
        gClient = _GephiStreamingClient()
        gClient.clean()
    
        nAttrs = {}
    
        for node in graph.nodes():
            gClient.add_node(str(node), **nAttrs)
    
        # TODO: how to determine edge ids in a reasonable way?
        edgeId = 1
        for edge in graph.edges():
            gClient.add_edge(edgeId, edge[0], edge[1], False)
            edgeId = edgeId + 1
    
        gClient.flush()
    except _urllib.error.URLError as e:
        _urlError(e)

def exportNodeValues(graph, values, attribute_name):
    """
     This method exports node values (e.g. community information, betwenness centrality values) 
     to gephi using the Gephi Streaming Plugin. Use exportGraph(Graph) first to export 
     the graph itself.
    
    Parameters:
    - values: python list or Partition object that contains the values to be exported.
    - attribute_name: name of the node attribute in gephi
    """
    
    try:
        gClient = _GephiStreamingClient()
    
        if len(values) != graph.numberOfNodes():
            print("Warning: #Nodes (", graph.numberOfNodes(), ") does not match #Values (", len(values), ").")
    
        for i in range(0, min(len(values), graph.numberOfNodes())):
            nAttrs = {attribute_name:values[i]}
            gClient.change_node(str(graph.nodes()[i]), **nAttrs)
        
        gClient.flush()
    except _urllib.error.URLError as e:
            _urlError(e)
            
def clearGraph():
    """ Cleans the first workspace in Gephi. """
    
    try:
        gClient = _GephiStreamingClient()
        gClient.clean()
        gClient.flush()
        
    except _urrlib.error.URLError as e:
        _urlError(e)
