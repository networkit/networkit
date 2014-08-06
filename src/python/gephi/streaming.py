import gephi.pyclient as _gephipyclient   # we want to hide these from the user
import urllib as _urllib

""" 
This module provides methods to export data from NetworKit directly to Gephi using the
Gephi Streaming Plugin.
"""

class GephiStreamingClient:
    """ 
     A python singleton providing access to a running Master instance of the Gephi
     Streaming plugin 
    """

    def __init__(self, url='http://localhost:8080/workspace0'):
        #Disabling Flushing means quite a good performance boost.
        self._pygephi = _gephipyclient.GephiClient(url, autoflush=False)
    
    def _urlError(self, e):
        print("Could not connect to the gephi streaming plugin. Did you start the streaming master server in gephi?")
    
    def exportGraph(self, graph):
        """ Exports the given graph to gephi. No attributes or weights are exported. """
    
        try:
            self._pygephi.clean()
    
            nAttrs = {}
    
            for node in graph.nodes():
                self._pygephi.add_node(str(node), **nAttrs)
        
            # TODO: how to determine edge ids in a reasonable way?
            edgeId = 1
            for edge in graph.edges():
                self._pygephi.add_edge(edgeId, edge[0], edge[1], False)
                edgeId = edgeId + 1
    
            self._pygephi.flush()
        except _urllib.error.URLError as e:
            self._urlError(e)

    def exportNodeValues(self, graph, values, attribute_name):
        """
        This method exports node values (e.g. community information, betwenness centrality values) 
        to gephi using the Gephi Streaming Plugin. Use exportGraph(Graph) first to export 
        the graph itself.
        
        Parameters:
        - values: python list or Partition object that contains the values to be exported.
        - attribute_name: name of the node attribute in gephi
        """
    
        try:    
            if len(values) != graph.numberOfNodes():
                print("Warning: #Nodes (", graph.numberOfNodes(), ") does not match #Values (", len(values), ").")
    
            for i in range(0, min(len(values), graph.numberOfNodes())):
                nAttrs = {attribute_name:values[i]}
                self._pygephi.change_node(str(graph.nodes()[i]), **nAttrs)
        
            self._pygephi.flush()
        except _urllib.error.URLError as e:
            self._urlError(e)
            
    def clearGraph(self):  
        """ Cleans the first workspace in Gephi. """
    
        try:
            self._pygephi.clean()
            self._pygephi.flush()
        
        except _urllib.error.URLError as e:
            self._urlError(e)
