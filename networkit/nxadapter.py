"""
This module handles compatibility between NetworKit and NetworkX
"""
from typing import Union, Mapping

# local imports
from . import graph
from warnings import warn
from .support import MissingDependencyError

# non standard library modules / external
import numpy as np

try:
    import networkx as nx
except ImportError:
    have_nx = False
else:
    have_nx = True

########  CONVERSION ########


def _inferType(attr: object) -> Union[int, float, str, None]:
    """
    Infer the type of an attribute.
    Integer types are mapped to int,
    Floating point types are mapped to float,
    str is mapped to str,
    any other type is mapped to None.
    """
    if np.issubdtype(type(attr), np.integer):
        return int
    if np.issubdtype(type(attr), np.floating):
        return float
    if type(attr) == str:
        return str
    return None


def nx2nk(
    nxG: Union[nx.Graph, nx.DiGraph],
    weightAttr: Union[str, None] = None,
    data: bool = False,
    typeMap: Union[Mapping[str, type], None] = None,
) -> graph.Graph:
    """
    nx2nk(nxG, weightAttr=None)

    Convert a networkx.Graph to a networkit.Graph.
    If data is true, try to convert networkx.Graph data to networkit.Graph attributes.
    Note that there are limitations to this conversion: networkit only supports int, float, and str attribute types.
    Other types will be converted into their string representation.
    Attribute keys are always converted to strings.
    Optionally, a dictionary that maps attribute names to specific types can be supplied.

    Parameters
    ----------
    nxG : networkx.Graph
            The input networkx graph.
    weightAttr : str, optional
            The edge attribute which should be treated as the edge weight. Default: None
    data : bool, optional
            If true, convert networkx.Graph data into networkit.Graph attributes. Default: False
    typeMap : dict, optional
            Specifies the data type an attribute has. Missing attributes are inferred as described above. Default: None
    """

    if not have_nx:
        raise MissingDependencyError("networkx")
    # map networkx node ids to consecutive numerical node ids
    idmap = dict((id, u) for (id, u) in zip(nxG.nodes(), range(nxG.number_of_nodes())))
    z = max(idmap.values()) + 1
    # print("z = {0}".format(z))

    if weightAttr is not None:
        nkG = graph.Graph(z, weighted=True, directed=nxG.is_directed())
        for u_, v_ in nxG.edges():
            u, v = idmap[u_], idmap[v_]
            w = nxG[u_][v_][weightAttr]
            nkG.addEdge(u, v, w)
    else:
        nkG = graph.Graph(z, directed=nxG.is_directed())
        for u_, v_ in nxG.edges():
            u, v = idmap[u_], idmap[v_]
            assert u < z
            assert v < z
            nkG.addEdge(u, v)

    assert nkG.numberOfNodes() == nxG.number_of_nodes()
    assert nkG.numberOfEdges() == nxG.number_of_edges()

    # convert data
    if data:
        # node attributes
        for node, attributes in nxG.nodes(data=True):
            # when we see a new attr, create/attach to graph. otherwise add to existing (get by name). if type is not compatible, raise exception. type is inferred from the first occurence.
            for key, value in attributes.items():
                valueType = (
                    typeMap.get(key, _inferType(value))
                    if typeMap
                    else _inferType(value)
                )
                try:
                    attribute = nkG.getNodeAttribute(str(key), valueType or str)
                except RuntimeError:  # attribute does not exist or is of different type
                    try:
                        attribute = nkG.attachNodeAttribute(str(key), valueType or str)
                        if valueType is None:
                            warn(
                                f"Info: the node attribute {key} has been converted to its string representation."
                            )
                    except (
                        RuntimeError
                    ):  # attribute exists (with different type because of previous logic)
                        raise RuntimeError(
                            f"Node attribute {key} has multiple data types which is not supported in networkit."
                        )

                if valueType is int:
                    attribute[idmap[node]] = int(value)
                elif valueType is float:
                    attribute[idmap[node]] = float(value)
                else:
                    attribute[idmap[node]] = str(value)

        # edge attributes
        nkG.indexEdges()
        for u, v, attributes in nxG.edges(data=True):
            # when we see a new attr, create/attach to graph. otherwise add to existing (get by name). if type is not compatible, raise exception. type is inferred from the first occurence.
            for key, value in attributes.items():
                if key == weightAttr:
                    continue
                valueType = (
                    typeMap.get(key, _inferType(value))
                    if typeMap
                    else _inferType(value)
                )
                try:
                    attribute = nkG.getEdgeAttribute(str(key), valueType or str)
                except RuntimeError:  # attribute does not exist or is of different type
                    try:
                        attribute = nkG.attachEdgeAttribute(str(key), valueType or str)
                        if valueType is None:
                            warn(
                                f"Info: the edge attribute {key} has been converted to its string representation."
                            )
                    except (
                        RuntimeError
                    ):  # attribute exists (with different type because of previous logic)
                        raise RuntimeError(
                            f"Edge attribute {key} has multiple data types which is not supported in networkit."
                        )

                if valueType is int:
                    attribute[idmap[u], idmap[v]] = int(value)
                elif valueType is float:
                    attribute[idmap[u], idmap[v]] = float(value)
                else:
                    attribute[idmap[u], idmap[v]] = str(value)

    return nkG


def nk2nx(nkG):
	""" 
	nk2nx(nkG)

	Convert a networKit.Graph to a networkx.Graph.

	Parameters
	----------
	nkG : networkit.Graph
		The input NetworKit graph.
	"""

	if not have_nx:
		raise MissingDependencyError("networkx")

	if nkG.isDirected():
		nxG = nx.DiGraph()
	else:
		nxG = nx.Graph()
	nxG.add_nodes_from(nkG.iterNodes())
	if nkG.isWeighted():
		for u, v, w in nkG.iterEdgesWeights():
			nxG.add_edge(u, v, weight=w)
	else:
		nxG.add_edges_from(nkG.iterEdges())

	assert (nkG.numberOfNodes() == nxG.number_of_nodes())
	assert (nkG.numberOfEdges() == nxG.number_of_edges())
	return nxG
