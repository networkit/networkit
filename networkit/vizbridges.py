# local imports
from .graph import Graph
from .support import MissingDependencyError
from .structures import Partition
from .viz import MaxentStress

# external imports
from collections import defaultdict
from enum import Enum
from typing import List, Mapping, Optional, Tuple, Union
from warnings import warn

import numpy as np

try:
    import ipycytoscape
except ImportError:
    hasCyto = False
else:
    hasCyto = True

try:
    import plotly.graph_objs as go
    import plotly.express as px
except ImportError:
    hasPlotly = False
else:
    hasPlotly = True
try:
    import seaborn
except ImportError:
    hasSeaborn = False
else:
    hasSeaborn = True


class Dimension(Enum):
    """
    Supported dimensions for visualization.

    Possible values:

    - networkit.vizbridges.Dimension.Two (visualization in 2D with Cytoscape)
    - networkit.vizbridges.Dimension.TwoForcePlotly (visualization in 2D with Plotly)
    - networkit.vizbridges.Dimension.Three (visualization in 3D with Plotly)
    """

    Two = 0
    TwoForcePlotly = 1
    Three = 2


# returns color palette for nodes and edges
def _getColorPalette(nodePalette=None, nodePartition=None, edgePalette=None):
    # Set color palettes
    if not nodePalette:
        if not hasSeaborn:
            raise MissingDependencyError("seaborn")
        # Partitions and scores have different default color palettes
        if nodePartition is not None:
            nodePalette = seaborn.color_palette("hls", nodePartition.numberOfSubsets())
        else:
            nodePalette = seaborn.color_palette("rocket_r", as_cmap=True).colors
    if not edgePalette:
        edgePalette = seaborn.color_palette("rocket_r", as_cmap=True).colors
    return nodePalette, edgePalette


def _calculateNodeColoring(G, palette, nodeScores=None, nodePartition=None):

    # Color calculation: score = continuous distribution, partition = discrete distribution
    hcColors = []

    # Partition
    if nodePartition is not None:
        if len(palette) < nodePartition.numberOfSubsets():
            raise IndexError(
                "Number of partitions higher than number of colors in provided palette. Provide node_palette with enough colors."
            )

        partitions = nodePartition.getVector()

        if len(palette) < nodePartition.numberOfSubsets():
            raise IndexError(
                "Number of partitions to high for default coloring. Provide node_palette with enough colors."
            )

        for i in range(0, len(partitions)):
            hcColors.append(
                (
                    palette[partitions[i]][0] * 255,
                    palette[partitions[i]][1] * 255,
                    palette[partitions[i]][2] * 255,
                )
            )

    # Score
    elif nodeScores is not None:

        minhc = min(nodeScores)
        maxhc = max(nodeScores)

        # Calculate coloring of nodes
        def getRgb(minimum, maximum, value):
            minimum, maximum, value = float(minimum), float(maximum), float(value)
            ratio = int(
                (len(palette) - 1)
                * (value - minimum)
                / (maximum - minimum)
                * (value - minimum)
                / (maximum - minimum)
            )
            r = int(palette[ratio][0] * 255)
            g = int(palette[ratio][1] * 255)
            b = int(palette[ratio][2] * 255)
            return r, g, b

        if abs(maxhc - minhc) > 0:
            for score in nodeScores:
                hcColors.append(getRgb(minhc, maxhc, score))
        else:
            color = palette[int((len(palette) - 1) / 2)]
            for i in range(0, len(nodeScores)):
                hcColors.append((color[0] * 255, color[1] * 255, color[2] * 255))

    # No node values
    else:
        color = palette[0]
        for i in G.iterNodes():
            hcColors.append((color[0] * 255, color[1] * 255, color[2] * 255))

    return hcColors


def _calculateEdgeColoring(G, palette, edgeScores: List[float] = None):

    # Color calculation: score = continuous distribution, partition = discrete distribution
    hcColors = []

    # Score
    if edgeScores:

        minhc = min(edgeScores)
        maxhc = max(edgeScores)

        # Calculate coloring of nodes
        def getRgb(minimum, maximum, value):
            minimum, maximum, value = float(minimum), float(maximum), float(value)
            ratio = int(
                (len(palette) - 1)
                * (value - minimum)
                / (maximum - minimum)
                * (value - minimum)
                / (maximum - minimum)
            )
            r = int(palette[ratio][0] * 255)
            g = int(palette[ratio][1] * 255)
            b = int(palette[ratio][2] * 255)
            return r, g, b

        if abs(maxhc - minhc) > 0:
            for score in edgeScores:
                hcColors.append(getRgb(minhc, maxhc, score))
        else:
            color = palette[int((len(palette) - 1) / 2)]
            for i in range(0, len(edgeScores)):
                hcColors.append((color[0] * 255, color[1] * 255, color[2] * 255))

    # No edge values
    else:
        color = palette[0]
        hcColors = [
            (color[0] * 255, color[1] * 255, color[2] * 255)
        ] * G.numberOfEdges()

    return hcColors


def _getEdgeScoreDirected(
    edgeScoresDict: Mapping[Tuple[int, int], float], u: int, v: int
) -> float:

    try:
        return edgeScoresDict[u, v]
    except KeyError:
        raise KeyError(
            f"edgeScores should include scores for every edge. ({u},{v}) is missing"
        )


def _getEdgeScoreUndirected(
    edgeScoresDict: Mapping[Tuple[int, int], float], u: int, v: int
) -> float:

    # this try block is requried to handle dict and defaultdict. get() returns None for defaultdict
    try:
        getUV = edgeScoresDict[u, v]
    except KeyError:
        getUV = None
    try:
        getVU = edgeScoresDict[v, u]
    except KeyError:
        getVU = None

    # missing edge score
    if getUV == None and getVU == None:
        raise KeyError(
            f"edgeScores should include scores for every edge. ({u},{v}) is missing"
        )

    # (u,v) is present
    if getUV == None and getVU != None:
        return getVU

    # (v,u) is present
    if getVU == None and getUV != None:
        return getUV

    # both are present and the same
    if getVU == getUV:
        return getUV

    # both directions are present and the score is different (and not None)
    if isinstance(edgeScoresDict, defaultdict):
        defaultValue = edgeScoresDict.default_factory()
        if getUV != defaultValue and getVU != defaultValue:
            raise ValueError(f"edgeScore different for both directions of ({u},{v}).")
        elif getUV != defaultValue:
            return getVU
        else:
            return getUV


def widgetFromGraph(
    G: Graph,
    dimension: Dimension = Dimension.Two,
    nodeScores: List[float] = None,
    nodePartition: Partition = None,
    nodePalette: List[Tuple[float, float, float]] = None,
    showIds: bool = True,
    customSize: int = None,
    edgeScores: Union[List[float], Mapping[Tuple[int, int], float]] = None,
    edgePalette: List[Tuple[float, float, float]] = None,
    edgeAttributes: List[Tuple[str, type]] = None,    
):
    """
    widgetFromGraph(G, dimension=Dimension.Two, nodeScores=None, nodePartition=None, nodePalette=None, showIds=True, customSize=None, edgeScores=None, edgePalette=None)

    Creates a widget with a visualization of a given graph. The widget uses one of the supported
    plugins - either Cytoscape (2D) or Plotly (3D). The returned widget already contains
    all nodes and edges from the graph. The graph is colored using an array of norm. rgb-values
    based on seaborn perceptually uniform color map (rocket) or a user given custom color array.
    See matplotlib color maps for the correct formatting:
    https://matplotlib.org/api/_as_gen/matplotlib.colors.Colormap.html#matplotlib.colors.Colormap


    Parameters
    ----------
    G : networkit.graph.Graph
            The input graph.
    dimension : networkit.vizbridges.Dimension, optional
            Select whether to plot in 2D or 3D. This also influences which plugin is choosen
            for visualization. For 2D Cytoscape with auto-layouting is used, for 3D Plotly with
            a Maxent-Stress layout is used. Option :code:`Dimension.TwoForcePlotly` forces a plot
            in 2D with using Plotly instead of Cytoscape. Default: networkit.vizbridges.Dimension.Two
    nodeScores : list(float) or dict, optional
            List of scores for each node, for example scores from a centrality measure. This is
            used for color-calculation of the nodes (continuous distribution). Provide either
            nodeScores or nodePartition - not both. Default: None
    nodePartition : networkit.structures.Partition, optional
            Partition object. This is used for color-calculation of the nodes (discrete distribution).
            Provide either nodeScores or nodePartition - not both. Default: None
    nodePalette : list(tuple(float, float, float)), optional
            List consisting of normalized rgb-values. If none is given, seaborn.color_palette.colors is used. Default: None
    showIds : boolean, optional
            Set whether node ids should be visible in plot-widget. Default: True
    customSize : int, optional
            If not set, plugins will use a default size. Otherwise the widget will have a certain width and height. Default: None
    edgeScores: list(float) or dict, optional
            List of scores for each edge, for example scores from a centrality measure. This is
            used for color-calculation of the edges (continuous distribution).
            If type is list, indexed by edgeids; if type is dict, indexed by node pair tuples. Default: None
    edgePalette : list(tuple(float, float, float)), optional
            List consisting of normalized rgb-values. If none is given, seaborn.color_palette.colors is used. Default: None
    edgeAttributes: list[tuple[str, type]], optional
            List consisting of all (already attached) edge attributes in the format: Tuple(nameOfAttribute, typeOfAttribute) that will be presented as edge labels. Default: None 
    """
    # Sanity checks
    if dimension == Dimension.Two and not hasCyto:
        raise MissingDependencyError("ipycytoscape")

    if (
        dimension == Dimension.Three or dimension == Dimension.TwoForcePlotly
    ) and not hasPlotly:
        raise MissingDependencyError("Plotly")

    if nodeScores is not None:
        if nodePartition is not None:
            raise ValueError("Provide either nodeScores or nodePartition - not both.")
        if len(nodeScores) != G.upperNodeIdBound():
            raise ValueError("nodeScores should include scores for every node.")

    if edgeScores:
        if isinstance(edgeScores, dict):
            if not G.hasEdgeIds():
                raise ValueError("Edges need to be indexed to draw edge scores.")
            edgeScoresList = [None] * G.numberOfEdges()
            for u, v in G.iterEdges():
                if G.isDirected():
                    edgeScoresList[G.edgeId(u, v)] = _getEdgeScoreDirected(
                        edgeScores, u, v
                    )
                else:
                    edgeScoresList[G.edgeId(u, v)] = _getEdgeScoreUndirected(
                        edgeScores, u, v
                    )
            edgeScores = edgeScoresList
        if None in edgeScores:
            raise ValueError("edgeScores should include scores for every edge.")
        if len(edgeScores) != G.numberOfEdges():
            raise ValueError(
                "edgeScores should include scores for every edge.",
                edgeScores,
                G.numberOfEdges(),
            )

    # Color palette is needed for node coloring with Plotly and Cytoscape
    nodePalette, edgePalette = _getColorPalette(
        nodePalette=nodePalette, nodePartition=nodePartition, edgePalette=edgePalette
    )

    if edgeAttributes:
        edgeAttributeList = []
        for u,v in G.iterEdges():
            att_cur_edge = ""
            for att in edgeAttributes:
                try:
                    cur_attribute = G.getEdgeAttribute(att[0], att[1])
                    if cur_attribute: 
                        att_cur_edge = att_cur_edge + str(att[0]) + ": " + str(cur_attribute[u,v]) + ", "
                except:
                    raise ValueError("Input edge attribute does not exist.")            
            edgeAttributeList.append(att_cur_edge)
            
    if dimension == Dimension.Two:
        # Set styling for cytoscape
        if showIds:
            style = [
                {
                    "selector": "node",
                    "css": {
                        "background-color": "data(color)",
                        "content": "data(id)",
                    },
                },
                {
                    "selector": "edge",
                    "css": {
                        "line-color": "data(color)",
                    },
                },
            ]
        else:
            style = [
                {"selector": "node", "css": {"background-color": "data(color)"}},
                {
                    "selector": "edge",
                    "css": {
                        "line-color": "data(color)",
                    },
                },
            ]
        # Create widget
        graphWidget = ipycytoscape.CytoscapeWidget()

        # Add data
        nodes = []
        edges = []
        if G.isDirected():
            edgeClass = "directed "
        else:
            edgeClass = "undirected "

        # Color list (norm. RGB-values) is needed for node coloring with Cytoscape
        nodeHcColors = _calculateNodeColoring(G, nodePalette, nodeScores, nodePartition)
        if edgeScores:
            edgeHcColors = _calculateEdgeColoring(G, edgePalette, edgeScores)
        else:
            edgeHcColors = None

        for i in G.iterNodes():
            n = ipycytoscape.Node(data={"id": i, "color": nodeHcColors[i]})
            nodes.append(n)

        for u, v in G.iterEdges():
            e = ipycytoscape.Edge(
                data={
                    "source": u,
                    "target": v,
                    "classes": edgeClass,
                    "color": edgeHcColors[G.edgeId(u, v)] if edgeHcColors else None,
                }
            )
            edges.append(e)

        # It is much faster to add edges and nodes in bulk.
        graphWidget.graph.add_nodes(nodes)
        graphWidget.graph.add_edges(edges, G.isDirected())

        # Set layout
        graphWidget.set_style(style)
        graphWidget.set_layout(name="cose")

    else:
        # Create widget
        graphWidget = go.FigureWidget()

        # Set layout
        maxLayout = MaxentStress(G, 3, 3, fastComputation=1, graphDistance=0)
        maxLayout.run()
        coordinates = maxLayout.getCoordinates()

        # Set node coloring
        if nodePartition is not None:
            scores = nodePartition.getVector()
        elif nodeScores is not None:
            scores = nodeScores
        else:
            scores = [0.0] * G.upperNodeIdBound()
        labels = [
            "Node: " + str(id) + "<br>Score: " + str(score)
            for id, score in enumerate(scores)
        ]

        # Initiate widget data
        nodes = [[], [], []]
        nodes[0] = [coordinates[k][0] for k in G.iterNodes()]
        nodes[1] = [coordinates[k][1] for k in G.iterNodes()]
        nodes[2] = [coordinates[k][2] for k in G.iterNodes()]

        index = 0
        if dimension == Dimension.TwoForcePlotly:
            edgeScatter = []

            for e in G.iterEdges():
                endpoint1 = [coordinates[e[0]][0], coordinates[e[1]][0]]
                endpoint2 = [coordinates[e[0]][1], coordinates[e[1]][1]]            
                if edgeAttributes:
                    midpoint1 = [np.mean([coordinates[e[0]][0], coordinates[e[1]][0]])]
                    midpoint2 = [np.mean([coordinates[e[0]][1], coordinates[e[1]][1]])]
                if edgeScores:
                    edgeHcColors = _calculateEdgeColoring(G, edgePalette, edgeScores)
                    colorTuple = edgeHcColors[G.edgeId(e[0], e[1])]
                    rgbString = "rgb(" + ",".join(map(str, colorTuple)) + ")"
                    line = dict(color=rgbString, width=5)
                else:
                    line = dict(color="rgb(180,180,180)", width=2)
                index = index + 1
                edgeScatter.append(
                    go.Scatter(
                        x=endpoint1,
                        y=endpoint2,
                        mode="lines",
                        opacity=0.7,
                        line=line,
                        hoverinfo="none",
                        showlegend=None,
                        name="edges",
                    )
                )
                if edgeAttributes:
                    edgeScatter.append(
                        go.Scatter(
                            x=midpoint1,
                            y=midpoint2,
                            mode="text",
                            opacity=0.7,
                            hoverinfo="none",
                            showlegend=None,
                            name="edges",
                            text=edgeAttributeList[index-1],
                        )
                    )

            nodeScatter = go.Scatter(
                x=nodes[0],
                y=nodes[1],
                mode="markers",
                name="nodes",
                marker=dict(
                    symbol="circle",
                    size=9,
                    colorscale=px.colors.convert_colorscale_to_rgb(
                        px.colors.make_colorscale(nodePalette)
                    ),
                    color=scores,
                    line=dict(color="rgb(50,50,50)", width=0.5),
                ),
                hoverinfo="text",
                text=labels,
            )

        else: # Dimension 3
            if edgeScores:
                edgeScoresMapped = np.zeros(3 * G.numberOfEdges())
            if edgeAttributes:
                edgeAttributePositions = [[],[],[]]
            edges = np.zeros((3, 3 * G.numberOfEdges()))
            for e in G.iterEdges():
                edges[0][index] = coordinates[e[0]][0]
                edges[0][index + 1] = coordinates[e[1]][0]
                edges[0][index + 2] = None
                edges[1][index] = coordinates[e[0]][1]
                edges[1][index + 1] = coordinates[e[1]][1]
                edges[1][index + 2] = None
                edges[2][index] = coordinates[e[0]][2]
                edges[2][index + 1] = coordinates[e[1]][2]
                edges[2][index + 2] = None
                if edgeAttributes: # use mean (mid point) to write edge label next to
                    edgeAttributePositions[0].append(np.mean([coordinates[e[0]][0], coordinates[e[1]][0]]))
                    edgeAttributePositions[1].append(np.mean([coordinates[e[0]][1], coordinates[e[1]][1]]))
                    edgeAttributePositions[2].append(np.mean([coordinates[e[0]][2], coordinates[e[1]][2]]))
                if edgeScores:
                    edgeScoresMapped[index] = edgeScores[G.edgeId(e[0], e[1])]
                    edgeScoresMapped[index + 1] = edgeScores[G.edgeId(e[0], e[1])]
                    edgeScoresMapped[index + 2] = edgeScores[G.edgeId(e[0], e[1])]

                index = index + 3

            nodeScatter = go.Scatter3d(
                x=nodes[0],
                y=nodes[1],
                z=nodes[2],
                mode="markers",
                name="nodes",
                marker=dict(
                    symbol="circle",
                    size=9,
                    colorscale=px.colors.convert_colorscale_to_rgb(
                        px.colors.make_colorscale(nodePalette)
                    ),
                    color=scores,
                    line=dict(color="rgb(50,50,50)", width=0.5),
                ),
                hoverinfo="text",
                text=labels,
            )

            if edgeScores:
                line = dict(
                    color=edgeScoresMapped,
                    colorscale=px.colors.convert_colorscale_to_rgb(
                        px.colors.make_colorscale(edgePalette)
                    ),
                    width=5,
                    autocolorscale=False,
                )
            else:
                line = dict(color="rgb(180,180,180)", width=2)

            edgeScatter = go.Scatter3d(
                x=edges[0],
                y=edges[1],
                z=edges[2],
                mode="lines",
                opacity=0.7,
                line=line,
                hoverinfo="none",
                showlegend=None,
                name="edges",
            )
            if edgeAttributes:
                graphWidget.add_traces(go.Scatter3d(
                    x=edgeAttributePositions[0],
                    y=edgeAttributePositions[1],
                    z=edgeAttributePositions[2],
                    mode="text",
                    opacity=0.7,
                    line=line,
                    hoverinfo="none",
                    showlegend=None,
                    name="edges",
                    text = edgeAttributeList
                    ))

        graphWidget.add_traces(nodeScatter)
        graphWidget.add_traces(edgeScatter)

        # Set layout
        minCoordinate, maxCoordinate = 0.0, 0.0
        for pos in coordinates:
            minCoordinate = min(minCoordinate, min(pos))
            maxCoordinate = max(maxCoordinate, max(pos))

        if customSize:
            width, height = customSize, customSize
        else:
            width, height = 1000, 1000

        axis = dict(
            showline=False,  # hide axis line, grid, ticklabels and  title
            zeroline=False,
            showgrid=True,
            showticklabels=True,
            title="",
            autorange=False,
            range=[1.1 * minCoordinate, 1.1 * maxCoordinate],
        )

        graphWidget.layout = go.Layout(
            font=dict(size=12),
            showlegend=False,
            autosize=True,
            scene_aspectmode="cube",
            width=width,
            height=height,
            scene=dict(
                xaxis=dict(axis),
                yaxis=dict(axis),
                zaxis=dict(axis),
                bgcolor="white",
                camera=dict(eye={"x": 1.3, "y": 1.3, "z": 0.5}),
            ),
            margin=go.layout.Margin(
                l=10,
                r=10,
                b=10,
                t=10,
            ),
            hovermode="closest",
        )
    return graphWidget
