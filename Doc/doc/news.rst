.. |br| raw:: html

   <br />

.. role:: hidden
   :class: hidden

====
News
====

.. just ignore the following header. This is a hack to make the other headings created with ~ smaller.

:hidden:`HiddenBiggerHeadingFont`
---------------------------------

February 25, 2017: **Migration to GitHub**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The NetworKit team is happy to announce that the NetworKit project has been successfully migrated to GitHub. Please join
us on

https://github.com/kit-parco/networkit

We believe the migration will make it easier for developers to contribute to the project and we hope to bring the advantages of efficient large-scale network analysis to even more people.

December 13, 2016: **NetworKit 4.2 released**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Today we announce the next version of NetworKit, the open-source toolkit for large-scale network analysis. NetworKit is a Python package, with performance-critical algorithms implemented in C++/OpenMP.

**Release notes**

Major:

- New graph drawing algorithm for the Maxent-stress model; the algorithm can layout even large graphs quickly. It follows the paper by Gansner et al. with some modifications; the biggest deviation is the use of the LAMG solver for the Laplacian linear systems
- Parallel implementation for the approximation of the neighborhood function; class has been refactored from ApproxNeighborhoodFunction to NeighborhoodFunctionApproximation.
- New heuristic algorithm for the neighborhood function. It is based on sampling and the breadth-first search and offers more flexibility with regards to the tradeoff between running time and accuracy as the number of samples can be specified by the user. It is also much faster than the approximation algorithm for networks with a high diameter (e.g. road networks).
    
Minor:

- Iterative implementation of components.StronglyConnectedComponents, which is now the new default. For graphs where edges have been deleted, it is recommended to use the recursive implementation, which is still available.
- Removed heuristic for vertex diameter estimation from centrality.ApproxBetweenness (now the vertex diameter is estimated as suggested in Riondato, Kornaropoulos: Fast approximation of betweenness centrality through sampling)
- Refactoring of the approximation algorithms in the distance group. ApproxNAME -> NAMEApproximation.
- Simplified installation procedure: Install required dependencies automatically

July 06, 2016: **Publication accepted at CSC 2016**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Our paper on approximating current-flow closeness centrality using NetworKit (authors: Bergamini, Wegner, Lukarski, Meyerhenke) has been accepted for publication at the 7th SIAM Workshop on Combinatorial Scientific Computing (CSC16). |br| |br|


July 05, 2016: **NetworKit 4.1.1 released**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a more of a maintenance release, that fixes the pip package and building with clang is possible again (at least with version 3.8).

Note: You can control which C++ compiler the setup.py of the networkit package is supposed to use with e.g. :code:`CXX=clang++ pip install networkit`. This may be helpful when the setup fails to detect the compiler.


June 23, 2016: **NetworKit 4.1 released**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Today we announce the next version of NetworKit, the open-source toolkit for large-scale network analysis.
NetworKit is a Python package, with performance-critical algorithms implemented in C++/OpenMP.

**Release notes**

Major:

new website

C++ implementation of Lean Algebraic Multigrid (LAMG) by Livne et al.
for solving large Laplacian systems serves as backend for various
network analysis kernels

centrality module

-  centrality.TopCloseness: Implementation of a new algorithm for
   finding the top-k nodes with highest closeness centrality faster than
   computing it for all nodes (E. Bergamini, M. Borassi, P. Crescenzi,
   A. Marino, H. Meyerhenke, "Computing Top-k Closeness Centrality
   Faster in Unweighted Graphs", ALENEX'16)

generator module:

-  generator.HyperbolicGenerator: a fast parallel generator for complex
   netwoks based on hyperbolic geometry (Looz, Meyerhenke, Prutkin '15:
   Random Hyperbolic Graphs in Subquadratic Time)

|  

   
Minor:

re-introduced an overview(G)-function that collects and prints some
infromation about a graph

updated documentation

some IO bugfixes

graph module:

-  Subgraph class has been removed, its functionality is now in
   Graph::subgraphFromNodes(...)

generator module: 

-  Many graph generators now provide fit(G) method that returns an
   instance of the generator such that generated graphs are similar to
   the provided one
-  Improved performance of the BarabasiAlbert generator by implementing
   Batagelj's method

distance module:

-  distance.CommuteTimeDistance: a node distance measure, distance is
   low when there are many short paths connecting two nodes
-  Adapted Diameter class to Algorithm convention; diameter algorithm
   can be chosen via enum in the constructor
-  Adapted EffectiveDiameter class to Algorithm convention resulting in
   the classes ApproxEffectiveDiameter, ApproxHopPlot,
   ApproxNeighborhoodFunction; added exact computation of the
   Neighborhood Function

centrality module:

-  centrality.SpanningEdgeCentraliy: edge centrality measure
   representing the fraction of spanning trees containing the edge
-  centrality.ApproxCloseness: new algorithm for approximating closeness
   centrality based on "Computing Classic Closeness Centrality, at
   Scale", Cohen et al.

|




May 9, 2016: **NetworKit journal paper accepted at Network Science**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Our paper describing NetworKit as a toolkit for large-scale complex network analysis has been accepted by the Cambride University Press journal Network Science. |br| |br|



Apr 12, 2016: **Publication accepted at SNAM**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Our paper on sparsification methods for social networks with NetworKit (authors: Linder, Staudt, Hamann, Meyerhenke, Wagner) has been accepted for publication in Social Network Analysis and Mining. |br| |br|



Apr 12, 2016: **Publication accepted at Internet Mathematics**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Our paper on approximating betweenness centrality in dynamic networks with NetworKit (authors: Bergamini, Meyerhenke) has been accepted for publication in Internet Mathematics. |br| |br|



Nov 16, 2016: **Publication accepted at ALENEX16**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Our paper on finding the top-k nodes with highest closeness centrality with NetworKit (authors: Bergamini, Borassi, Crescenzi, Marino, Meyerhenke) has been accepted at the 18th Meeting on Algorithm Engineering and Experiments, ALENEX 2016. |br| |br|



Nov 10, 2015: **NetworKit 4.0 released**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We have just released NetworKit 4.0. Apart from several improvements to algorithms and architecture, the main feature of this release is a new front end for exploratory network analysis, which is described here:

http://nbviewer.ipython.org/urls/networkit.iti.kit.edu/data/uploads/docs/Profiling.ipynb

The new version is now available from the Python Package index. Try upgrading with
:code:`pip3 install —upgrade networkit` |br| |br|


Aug 19, 2015: **NetworKit 3.6 released**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We have released version 3.6 today. Thank you to all contributors. Here are the release notes.

*Release Notes*

Major Updates:

Link Prediction

Link prediction methods try to predict the likelihood of a future or missing connection between two nodes in a given network. The new module networkit.linkprediction contains various methods from the literature.

Edge Sparsification

Sparsification reduces the size of networks while preserving structural and statistical properties of interest. The module networkit.sparsification provides methods for rating edges by importance and then filtering globally by these scores. The methods are described in http://arxiv.org/abs/1505.00564


Further Updates:

- Improved support for directed graph in analysis algorithms
- Improved support for the Intel compiler
- Reader/writer for the GEXF (Gephi) graph file format
- EdgeListReader now reads edge list with arbitrary node ids (e.g.strings) when continuous=False; getNodeMap() returns a mapping from file node ids to graph node ids
- EdgeListReader/Writer now add weights when reading files/writing graphs to file. |br| |br|


Jun 16, 2015: **Publication accepted at ESA15**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Our paper on the approximation of betweenness centrality in fully-dynamic networks with NetworKit (authors: Bergamini, Meyerhenke) has been accepted at the 23rd European Symposium on Algorithms, ESA 2015. |br| |br|


Jun 9, 2015: **NetworKit 3.5 released**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We have released NetworKit 3.5 a couple days ago. Please upgrade to the latest version to receive a number of improvements. We also appreciate feedback on the new release.

*Release Notes*

This release focused on bugfixes, under-the-hood improvements and refactoring.

- Various bugfixes and stability improvements
- Abort signal handling: developed mechanism to interrupt long-running algorithms via the ctrl+C command -- already supported in community.PLM, centrality.Betweennness, centrality.ApproxBetweenness, centrality.ApproxBetweenness2, centrality.PageRank
- Efficient node and edge iteration on the Python layer: G.forEdges, G.forNodes...
- Constant-time check if a graph has self-loops: Graph.hasSelfLoops()
- networkit.setSeed: set a fixed seed for the random number generator
- Refactoring: CoreDecomposition and LocalClusteringCoefficient now in centrality module
- Refactoring: introduced Python/Cython base classes: Centrality, CommunityDetector
- Removed: CNM community detection algorithm
- The GIL (Global Interpreter Lock) is released for many algorithms in order to make it possible to execute multiple computations in parallel in a single Python process.
- Improved support for directed graphs in many algorithms |br| |br|


Dec 4, 2014: **NetworKit 3.4 released**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Today we have released version 3.4 of NetworKit, the open-source toolkit for high-performance network analysis. This release brings numerous critical bugfixes as well as useful incremental features and performance optimizations. We are also moving towards consistent interfaces for algorithms. We have also further simplified the installation dependencies.

Thank you to the numerous people who have contributed code to this release.

More information can be found on https://networkit.iti.kit.edu/. We welcome user feedback and opportunities for collaboration.

Release Notes

Features

* graph
   * Graph can be copied on Python level
   * spanning tree/forest (graph.SpanningForest)
*  algorithms in general
   * Edmonds-Karp max flow algorithm (flow.EdmondsKarp)
   * core decomposition works for directed graphs (properties.CoreDecomposition)
   * algebraic distance, a structural distance measure in graphs (distance.AlgebraicDistance)
* IO
   * there is no longer a default graph file format
   * read and write the GML graph file format (graphio.GMLGraphReader/Writer)
   * conversion of directed to undirected graph (Graph.toUndirected)
   * reader and writer for the GraphTool binary graph format (graphio.GraphToolBinaryReader)
   * METIS graph reader supports arbitrary edge weights (graphio.METISGraphReader)
* algebraic
   * algebraic backend supports rectangular matrices (Matrix.h)
* community detection
   * turbo mode for PLM community detection algorithm gives a factor 2 speedup at the cost of more memory (community.PLM)
   * Cut Clustering community detection algorithm (community.CutClustering)
* generators
   * Erdös-Renyi generator can generate directed graphs (generators.ErdosRenyiGenerator)
   * configuration model graph generator for generating a random simple graph with exactly the given degree sequence (generators.ConfigurationModelGenerator)
   * generator for power law degree sequences (generators.PowerlawDegreeSequence)

Bugfixes

* GraphMLReader improved (graphio.GraphMLReader)
* ConnectedComponents usability improved
* KONECT reader (graphio.KONECTGraphReader)
* fixed build problem on case-insensitive file systems
* closed memory leaks by adding missing destructors on the Cython
* improved memory management by adding missing move constructors
* DynamicForestFireGenerator fixed

Refactoring

* standardization of analysis algorithm interface: parameters given by constructor, computation triggered in run method, results retrieved via getter methods
* run methods return self to allow chaining
* introducing unit tests on Python layer

Build and Installation

* pip installation does no longer require Cython
* pip installation does no longer require SCons, minimal build system as fallback if SCons is missing |br| |br|



Oct 21, 2014: **Publication accepted at ALENEX15**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Our paper on approximating betweenness centrality in dynamic networks with NetworKit (authors: Bergamini, Meyerhenke, Staudt) has been accepted at the 17th Meeting on Algorithm Engineering and Experiments, ALENEX 2015. |br| |br|



Sep 28, 2014: **NetworKit presented at summer school tutorial on network analysis**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In a joint tutorial on Algorithmic methods for network analysis with Dorothea Wagner for the summer school of the DFG priority programme Algorithm Engineering, Henning Meyerhenke introduced NetworKit to the participants. The PhD students from Germany and other European countries successfully solved various network analysis tasks with NetworKit during the tutorial. |br| |br|



Sep 28, 2014: **Publication accepted**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Our paper on selective community detection with NetworKit (authors: Staudt, Marrakchi, Meyerhenke) has been accepted at the First International Workshop on High Performance Big Graph Data Management, Analysis, and Mining (in Conjunction with IEEE BigData'14). |br| |br|



Aug 22, 2014: **NetworKit 3.3 released**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NetworKit 3.3 has been released, including the following improvements to our network analysis framework:

- renamed package to "networkit" according to Python packaging convention
- restructured package to enable "pip install networkit"
- improved community detection algorithms
- improved diameter algorithms
- added support for efficient, arbitrary edge attributes via edge indexing
- Eigenvector Centrality & PageRank on basis of scipy
- spectral methods for graph partitioning  (partitioning.SpectralPartitioner), drawing  (viztools.layout.SpectralLayout) and coloring  (coloring.SpectralColoring)
- new graph generators: stochastic blockmodel (generators.StochasticBlockmodel), Watts-Strogatz model (generators.WattsStrogatzGenerator) and Forest Fire model (generators.DynamicForestFireGenerator)
- union find data structure (structures/UnionFind)
- simple spanning forest algorithm (graph.SpanningForest)
- fast algorithm for partition intersection (community/PartitionIntersection)
- hub dominance in communities (community.HubDominance)
- reader for Matlab adjacency matrices
- support for reading and writing Covers
- performance improvements in Gephi streaming interface |br| |br|



Jul 1, 2014: **NetworKit 3.2 released**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NetworKit 3.2 has been released, including major improvements to our network analysis framework:

*Critical Bugfixes*

- graph data structure supports directed graphs
- optimized connected components algorithm (properties.ParallelConnectedComponents)
- faster heuristic algorithm for approximating betweenness centrality (centrality.ApproxBetweenness2)
- Gephi support: export of node attributes, Gephi streaming plugin support
- graph generators: Dorogovtsev-Mendes model
- improved portability (Windows)
- overhaul of graph file input |br| |br|



May 15, 2014: **New website online**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NetworKit, our tool suite for high-performance network analysis, has its own website now! |br| |br|



Apr 25, 2014: **Introductory talk**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Christian Staudt gave an introductory talk about the current release of NetworKit. The slides and a video of the talk are available on the Documentation page. |br| |br|



Apr 15, 2014: **NetworKit 3.1 released**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Version 3.1 is an incremental update to our tool suite for high-performance network analysis. Improvements and new features include Eigenvector centrality, PageRank, Betweenness centrality approximation, R-MAT graph generator, BFS/DFS iterators, improved BFS and Dijkstra classes, and improved memory footprint when using large objects on the Python level. More detailed information can be found in the accompanying publication. |br| |br|



Mar 13, 2014: **NetworKit 3.0 released**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NetworKit 3.0 is the next major release of our open-source tookit for high-performance network analysis. Since the last release in November, NetworKit has received several improvements under the hood as well as an extension of the feature set. What started as a testbed for parallel community detection algorithms has evolved into a diverse set of tools that make it easy to characterize complex networks. This has been successfully scaled to large data sets with up to several billions of edges.

This being an open-source project, we are very interested in incorporating feedback from data analysts and algorithm engineers. Feel free to contact us with any question on how NetworKit could be applied in your field of research. |br| |br|



Nov 11, 2013: **NetworKit 2.0 released**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Second major release of NetworKit. The toolkit has been improved by adding several graph algorithms and an interactive shell based on Python/Cython. We begin a more frequent release cycle. |br| |br|



Mar 17, 2013: **NetworKit 1.0 released**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Initial release of the community detection component. With this release of NetworKit, we would like to encourage reproduction of our results, reuse of code and contributions by the community. |br| |br|
