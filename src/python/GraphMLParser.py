from xml.dom import minidom

from _NetworKit import Graph

class GraphMLParser:
    """
"""

    def __init__(self):
        """
"""

    def write(self, graph, fname):
        """
"""

        doc = minidom.Document()

        root = doc.createElement('graphml')
        doc.appendChild(root)

        # if the graph is weighted, add the attribute
        if graph.isWeighted():
             attr_node = doc.createElement('key')
             attr_node.setAttribute('for','edge')
             attr_node.setAttribute('id', 'd1')
             attr_node.setAttribute('attr.name','weight')
             attr_node.setAttribute('attr.type','double')
             root.appendChild(attr_node)

        # Add attributs
        #for a in graph.get_attributs():
        #    attr_node = doc.createElement('key')
        #    attr_node.setAttribute('id', a.name)
        #    attr_node.setAttribute('attr.name', a.name)
        #    attr_node.setAttribute('attr.type', a.type)
        #    root.appendChild(attr_node)
        
        graph_node = doc.createElement('graph')
        graph_node.setAttribute('id', graph.getName())
        #if graph.directed:
        #    graph_node.setAttribute('edgedefault', 'directed')
        #else:
        #    graph_node.setAttribute('edgedefault', 'undirected')
        graph_node.setAttribute('edgedefault', 'undirected')
        root.appendChild(graph_node)

        # Add nodes
        for n in graph.nodes():
            node = doc.createElement('node')
            node.setAttribute('id', str(n))
            #for a in n.attributes():
            #    if a != 'label':
            #        data = doc.createElement('data')
            #        data.setAttribute('key', a)
            #        data.appendChild(doc.createTextNode(str(n[a])))
            #        node.appendChild(data)
            graph_node.appendChild(node)

        # Add edges
        if graph.isWeighted():
            for e in graph.edges():
                edge = doc.createElement('edge')
                edge.setAttribute('source', str(e[0]))
                edge.setAttribute('target', str(e[1]))
                edge.setAttribute('directed', 'false')
                # add edge weight
                data = doc.createElement('data')
                data.setAttribute('key','d1')
                data.appendChild(doc.createTextNode(str(graph.weight(e[0],e[1]))))
                edge.appendChild(data)
                graph_node.appendChild(edge)
        else:
            for e in graph.edges():
                edge = doc.createElement('edge')
                edge.setAttribute('source', str(e[0]))
                edge.setAttribute('target', str(e[1]))
                edge.setAttribute('directed', 'false')
                graph_node.appendChild(edge)

        f = open(fname, 'w')
        f.write(doc.toprettyxml(indent = ' '))
    
    def read(self, fname):
        """
"""

        dom = minidom.parse(open(fname, 'r'))
        root = dom.getElementsByTagName("graphml")[0]
        graph = root.getElementsByTagName("graph")[0]
        name = graph.getAttribute('id')

        # Get information if the graph is weighted and the attribute identifier
        attributes = []
        weighted = False
        for attr in root.getElementsByTagName("key"):
            if (attr.getAttribute('for') == 'edge' and attr.getAttribute('attr.name') == 'weight'):
                weighted = True
                weightedID = attr.getAttribute('id')
                #print("graph identified as weighted, weighted={0}, weightedID={1}".format(weighted,weightedID))

        g = Graph(0,weighted)
        g.setName(name)

        #print("start to gather node information")
        mapping = dict()
        # Get nodes
        for node in graph.getElementsByTagName("node"):
            val = node.getAttribute('id')
            n = g.addNode()
            mapping[val] = n

	    #ignore node attributes
            #for attr in node.getElementsByTagName("data"):
            #    if attr.firstChild:
            #        n[attr.getAttribute("key")] = attr.firstChild.data
            #    else:
            #        n[attr.getAttribute("key")] = ""
        #print("gathered node information, adding edges")
        # Get edges
        for edge in graph.getElementsByTagName("edge"):
            source = mapping[edge.getAttribute('source')]
            dest = mapping[edge.getAttribute('target')]
            #print("added edge: ({0} - {1})".format(source,dest))
            #g.addEdge(source,dest) #FIXME: no edge weights yet

	    #ignore edge attributes
            edgeweight = 0.0
            for attr in edge.getElementsByTagName("data"):
                #print("found attribute with id={0} and data={1}".format(attr.getAttribute('key'),attr.firstChild.data))
                if attr.getAttribute('key') == weightedID:
                    edgeweight = float(attr.firstChild.data)
                    #print("parsed edgeweight for ({0} - {1}) with {2}".format(source,dest,edgeweight))

            g.addEdge(source,dest,edgeweight)

        return g


if __name__ == '__main__':

    parser = GraphMLParser()
    g = parser.parse('test.graphml')

    g.show(True)
