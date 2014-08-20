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
        self._pygephi = _gephipyclient.GephiClient(url, autoflush=10000)

    def _urlError(self, e):
        print("Could not connect to the gephi streaming plugin. Did you start the streaming master server in gephi?")

    def exportGraph(self, graph):
        """ Exports the given graph to gephi. No attributes or weights are exported.
        Please note that only graphs with indexed edges can be exported, so
        edges will be indexed if neccessary. """

        try:
            self._pygephi.clean()

            #Index edges if neccessary
            graph.indexEdges()

            nAttrs = {}

            for node in graph.nodes():
                self._pygephi.add_node(str(node), **nAttrs)

            for edge in graph.edges():
                edgeId = graph.edgeId(edge[0], edge[1])
                self._pygephi.add_edge(edgeId, edge[0], edge[1], False)

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

    def exportEdgeValues(self, graph, values, attribute_name):
        """
        This method exports an edge attribute to gephi using the Gephi Streaming Plugin.
        Use exportGraph(Graph) first to export graph itself.

        Parameters:
        - values: python list contains the values to be exported.
        - attribute_name: name of the edge attribute in gephi
        """
        try:
            if len(values) != graph.numberOfEdges():
                print("Warning: #Nodes (", graph.numberOfEdges(), ") does not match #Values (", len(values), ").")

            idx = 0
            for edge in graph.edges():
                edgeId = graph.edgeId(edge[0], edge[1])
                eAttrs = {attribute_name:values[edgeId], "Type":"Undirected"}
                self._pygephi.change_edge(edgeId, edge[0], edge[1], False, **eAttrs)
                idx += 1

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
