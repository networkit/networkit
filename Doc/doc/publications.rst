.. include:: <isonum.txt>

.. |br| raw:: html

   <br />

.. |separator| raw:: html

	<div style="padding-top: 25px; border-bottom: 1px solid #d4d7d9;"></div>

.. _publications:

============
Publications
============

The following is a list of publications on the basis of NetworKit. We ask you to **cite** the appropriate ones if you found NetworKit useful for your own research.
Also, we would appreciate it if you pointed us to your publications in which you used NetworKit and allowed us to reference them on this page.

Journal Paper on NetworKit as a Software Toolkit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. raw:: html

  <ul>
    <li>
        C. Staudt, A. Sazonovs and H. Meyerhenke: NetworKit: A Tool Suite for Large-scale Complex Network Analysis. To appear in Network Science, Cambridge University Press.
        [<a href="http://arxiv.org/abs/1403.3005">arXiv</a>]
        <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
        <div id="collapseDiv" class="collapse">
          <b>Abstract.</b> We introduce NetworKit, an open-source software package for analyzing the structure of large complex networks. Appropriate algorithmic solutions
          are required to handle increasingly common large graph data sets containing up to billions of connections. We describe the methodology applied to develop scalable
          solutions to network analysis problems, including techniques like parallelization, heuristics for computationally expensive problems, efficient data structures,
          and modular software architecture. Our goal for the software is to package results of our algorithm engineering efforts and put them into the hands of domain experts.
          NetworKit is implemented as a hybrid combining the kernels written in C++ with a Python front end, enabling integration into the Python ecosystem of tested tools for
          data analysis and scientific computing. The package provides a wide range of functionality (including common and novel analytics algorithms and graph generators) and
          does so via a convenient interface. In an experimental comparison with related software, NetworKit shows the best performance on a range of typical analysis tasks.
        </div>
    </li>
  </ul>

|separator|


Publications on Algorithms Available in NetworKit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. raw:: html

  <ul>
  <li>
    E. Bergamini, M. Wegner, D. Lukarski, H. Meyerhenke:  Estimating Current-Flow Closeness Centrality with a Multigrid Laplacian Solver. To appear at <i><a href="http://www.eecs.wsu.edu/~assefaw/CSC16/csc16.html">CSC</a> '16</i>.
    <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
    <div id="collapseDiv" class="collapse">
      <b>Abstract.</b> Matrices associated with graphs, such as the Laplacian, lead to numerous interesting graph problems expressed as linear systems. One field where Laplacian linear systems
      play a role is network analysis, e. g. for certain centrality measures that indicate if a node (or an edge) is important in the network. One such centrality measure is current-flow closeness.
      To allow network analysis workflows to profit from a fast Laplacian solver, we provide an implementation of the LAMG multigrid solver in the NetworKit package, facilitating the computation
      of current-flow closeness values or related quantities. Our main contribution consists of two algorithms that accelerate the current-flow computation for one node or a reasonably small node
      subset significantly. One algorithm is an unbiased estimator using sampling, the other one is based on the Johnson- Lindenstrauss transform. Our inexact algorithms lead to very accurate
      results in practice. Thanks to them one is now able to compute an estimation of current-flow closeness of one node on networks with tens of millions of nodes and edges within seconds or a
      few minutes. From a network analytical point of view, our experiments indicate that current-flow closeness can discriminate among different nodes significantly better than traditional
      shortest-path closeness and is also considerably more resistant to noise – we thus show that two known drawbacks of shortest-path closeness are alleviated by the current-flow variant.
    </div>
  </li>

  <br>

    <li>
      M. von Looz, M. Özdayi, S. Laue, H. Meyerhenke:  Generating massive complex networks with hyperbolic geometry faster in practice. To appear at <i><a href="http://ieee-hpec.org/">HPEC</a> '16</i>. [<a href="http://arxiv.org/abs/1606.09481">arXiv</a>]
      <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
      <div id="collapseDiv" class="collapse">
        <b>Abstract.</b> Generative network models play an important role in algorithm development, scaling studies, network analysis, and realistic system benchmarks for graph data sets.   The commonly used graph-based benchmark model R-MAT has some drawbacks concerning realism and the scaling behavior of network properties.
        A complex network model gaining considerable popularity builds random hyperbolic graphs, generated by distributing points within a disk in the hyperbolic plane and then adding edges between points whose hyperbolic distance is below a threshold.
        We present in this paper a fast generation algorithm for such graphs.
        Our experiments show that our new generator achieves speedup factors of 3-60 over the best previous implementation.
        One billion edges can now be generated in under one minute on a shared-memory workstation.
        Furthermore, we present a dynamic extension to model gradual network change, while preserving at each step the point position probabilities.
      </div>
    </li>

    <br>

    <li>
      M. von Looz, H. Meyerhenke: Querying Probabilistic Neighborhoods in Spatial Data Sets Efficiently. To appear at <i><a href="http://iwoca2016.cs.helsinki.fi/">IWOCA 2016</a></i>. [<a href="http://arxiv.org/abs/1509.01990">arXiv</a>]
      <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
      <div id="collapseDiv" class="collapse">
        <b>Abstract.</b> The probability that two spatial objects establish some kind of mutual connection often depends on their proximity.
        To formalize this concept, we define the notion of a <i>probabilistic neighborhood</i>:
        \(\newcommand{\dist}{\operatorname{dist}}\)
        Let \(P\) be a set of \(n\) points in \(\mathbb{R}^d\),  \(q \in \mathbb{R}^d\) a query point, \(\dist\) a distance metric, and \(f : \mathbb{R}^+ \rightarrow [0,1]\) a monotonically decreasing function.
        Then, the probabilistic neighborhood \(N(q, f)\) of \(q\) with respect to \(f\) is
        a random subset of \(P\) and each point \(p \in P\) belongs to \(N(q,f)\) with probability \(f(\dist(p,q))\).
        Possible applications include query sampling and the simulation of probabilistic spreading phenomena, as well as other scenarios where the probability of a connection between two entities decreases with their distance.
        We present a fast, sublinear-time query algorithm to sample probabilistic neighborhoods from planar point sets.
        For certain distributions of planar \(P\), we prove that our algorithm answers a query in \(O((|N(q,f)| + \sqrt{n})\log n)\) time with high probability.
        In experiments this yields a speedup over pairwise distance probing of at least one order of magnitude, even for rather small data sets with \(n=10^5\) and also for other point distributions not covered by the theoretical results.
      </div>
    </li>

    <br>

    <li>
      E. Bergamini, H. Meyerhenke: Approximating Betweenness Centrality in Fully-dynamic Networks. To appear in <i>Internet Mathematics</i>. Code to appear in NetworKit.
      <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
      <div id="collapseDiv" class="collapse">
        <b>Abstract.</b> Betweenness is a well-known centrality measure that ranks the nodes of a network according to their participation in shortest paths. Since an exact
        computation is prohibitive in large networks, several approximation algorithms have been proposed. Besides that, recent years have seen the publication of dynamic
        algorithms for efficient recomputation of betweenness in networks that change over time. In this paper we propose the first betweenness centrality approximation
        algorithms with a provable guarantee on the maximum approximation error for dynamic networks. Several new intermediate algorithmic results contribute to the
        respective approximation algorithms: (i) new upper bounds on the vertex diameter, (ii) the first fully-dynamic algorithm for updating an approximation of the
        vertex diameter in undirected graphs, and (iii) an algorithm with lower time complexity for updating single-source shortest paths in unweighted graphs after a batch
        of edge actions. Using approximation, our algorithms are the first to make in-memory computation of betweenness in dynamic networks with millions of edges feasible.
        Our experiments show that our algorithms can achieve substantial speedups compared to recomputation, up to several orders of magnitude. Moreover, the approximation
        accuracy is usually significantly better than the theoretical guarantee in terms of absolute error. More importantly, for reasonably small approximation error
        thresholds, the rank of nodes is well preserved, in particular for nodes with high betweenness.
      </div>
    </li>

    <br>

    <li>
      G. Lindner, C. L. Staudt, M. Hamann, H. Meyerhenke, D. Wagner: Structure-Preserving Sparsification Methods for Social Networks. To appear in <i>Social Network Analysis
      and Mining</i>.
      <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
      <div id="collapseDiv" class="collapse">
        <b>Abstract.</b> Sparsification reduces the size of networks while preserving structural and statistical properties of interest. Various sparsifying algorithms have been
        proposed in different contexts. We contribute the first systematic conceptual and experimental comparison of edge sparsification methods on a diverse set of network
        properties. It is shown that they can be understood as methods for rating edges by importance and then filtering globally or locally by these scores. We show
        that applying a local filtering technique improves the preservation of all kinds of properties. In addition, we propose a new sparsifi- cation method (Local Degree)
        which preserves edges leading to local hub nodes. All methods are evaluated on a set of social networks from Facebook, Google+, Twitter and LiveJournal with respect
        to network properties including diameter, connected components, community structure, multiple node centrality measures and the behavior of epidemic simulations.
        In order to assess the preservation of the community structure, we also include experiments on synthetically generated networks with ground truth communities.
        Experiments with our implementations of the sparsification methods (included in the open-source network analysis tool suite NetworKit) show that many network
        properties can be preserved down to about 20% of the original set of edges for sparse graphs with a reasonable density. The experimental results allow us to
        differentiate the behavior of different methods and show which method is suitable with respect to which property. While our Local Degree method is best for
        preserving connectivity and short distances, other newly introduced local variants are best for preserving the community structure.
      </div>
    </li>

    <br>

    <li>
      E. Bergamini, M. Borassi, P. Crescenzi, A. Marino, H. Meyerhenke: Computing Top-k Closeness Centrality Faster in Unweighted Graphs. In <i>Proc. SIAM Algorithm
      Engineering & Experiments</i> (ALENEX 2016). Code to appear in NetworKit.
      <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
      <div id="collapseDiv" class="collapse">
        <b>Abstract.</b> Centrality indices are widely used analytic measures for the importance of nodes in a network. Closeness centrality is very popular among these measures.
        For a single node v, it takes the sum of the distances of v to all other nodes into account. The currently best algorithms in practical applications for
        computing the closeness for all nodes exactly in unweighted graphs are based on breadth-first search (BFS) from every node. Thus, even for sparse graphs,
        these algorithms require quadratic running time in the worst case, which is prohibitive for large networks. <br>
        In many relevant applications, however, it is un- necessary to compute closeness values for all nodes. Instead, one requires only the k nodes with the highest
        closeness values in descending order. Thus, we present a new algorithm for computing this top-k ranking in unweighted graphs. Following the rationale of previous
        work, our algorithm significantly re- duces the number of traversed edges. It does so by computing upper bounds on the closeness and stopping the current BFS search
        when k nodes already have higher closeness than the bounds computed for the other nodes.<br>
        In our experiments with real-world and synthetic instances of various types, one of these new bounds is good for small-world graphs with low diameter (such as
        social networks), while the other one excels for graphs with high diameter (such as road networks). Combining them yields an algorithm that is faster than the state
        of the art for top-k computations for all test instances, by a wide margin for high-diameter.
      </div>
    </li>

    <br>

    <li>
      M. von Looz, R. Prutkin and H. Meyerhenke: Fast Generation of Complex Networks with Underlying Hyperbolic Geometry. In <i>Proc. 26th International Symposium on
      Algorithms and Computation</i> (ISAAC 2015). Code in NetworKit. [<a href="http://arxiv.org/abs/1501.03545">arXiv</a>]
      <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
      <div id="collapseDiv" class="collapse">
        <b>Abstract.</b> Complex networks have become increasingly popular for mod- eling various real-world phenomena. Realistic generative network models are important in
        this context as they avoid privacy concerns of real data and simplify complex network research regarding data sharing, reproducibility, and scalability studies.
        Random hyperbolic graphs are a well-analyzed family of geometric graphs. Previous work provided empir- ical and theoretical evidence that this generative graph
        model creates networks with non-vanishing clustering and other realistic features. How- ever, the investigated networks in previous applied work were small,
        possibly due to the quadratic running time of a previous generator.
        In this work we provide the first generation algorithm for these networks with subquadratic running time. We prove a time complexity of
        \(\mathcal{O}((n^{3/2} + m) \log n)\) with high probability for the generation process. This running time is confirmed by experimental data with our
        implementation. The acceleration stems primarily from the reduction of pairwise distance computations through a polar quadtree, which we adapt to hyperbolic
        space for this purpose. In practice we improve the running time of a previous implementation by at least two orders of magnitude this way. Networks with billions
        of edges can now be generated in a few minutes. Finally, we evaluate the largest networks of this model published so far. Our empirical analysis shows that
        important features are retained over different graph densities and degree distributions.
      </div>
    </li>

    <br>

    <li>
      E. Bergamini and H. Meyerhenke: Fully-dynamic Approximation of Betweenness Centrality. In <i>Proc. 23rd European Symp. on Algorithms</i> (ESA 2015). Code to appear
      in NetworKit. [<a href="http://arxiv.org/abs/1504.07091">arXiv</a>]
      <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
      <div id="collapseDiv" class="collapse">
        <b>Abstract.</b> Betweenness is a well-known centrality measure that ranks the nodes of a network according to their participation in shortest paths. Since an
        exact computation is prohibitive in large networks, several approximation algorithms have been proposed. Besides that, recent years have seen the publication of
        dynamic algorithms for efficient recomputation of betweenness in evolving networks. In previous work we proposed the first semi-dynamic algorithms that recompute
        an approximation of betweenness in connected graphs after batches of edge insertions.In this paper we propose the first fully-dynamic approximation algorithms
        (for weighted and unweighted undirected graphs that need not to be connected) with a provable guarantee on the maximum approxima- tion error. The transfer to
        fully-dynamic and disconnected graphs implies additional algorithmic problems that could be of independent interest. In particular, we propose a new upper bound
        on the vertex diameter for weighted undirected graphs. For both weighted and unweighted graphs, we also propose the first fully-dynamic algorithms that keep
        track of this upper bound. In addition, we extend our former algorithm for semi- dynamic BFS to batches of both edge insertions and deletions. <br>
        Using approximation, our algorithms are the first to make in-memory computation of betweenness in fully-dynamic networks with millions of edges feasible.
        Our experiments show that they can achieve substantial speedups compared to recomputation, up to several orders of magnitude.
      </div>
    </li>

    <br>

    <li>
      E. Bergamini, H. Meyerhenke and  C. Staudt: Approximating Betweenness Centrality in Large Evolving Networks. In <i>Proc. SIAM Algorithm Engineering & Experiments</i>
      (ALENEX 2015). [<a href="http://arxiv.org/abs/1409.6241">arXiv</a>] [<a href="http://dx.doi.org/10.1137/1.9781611973754.12">DOI: 10.1137/1.9781611973754.12</a>]
      <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
      <div id="collapseDiv" class="collapse">
        <b>Abstract.</b> Betweenness centrality ranks the importance of nodes by their participation in all shortest paths of the network. Therefore computing exact
        betweenness values is impractical in large networks. For static networks, approximation based on randomly sampled paths has been shown to be significantly faster
        in practice. However, for dynamic networks, no approximation algorithm for betweenness centrality is known that improves on static recomputation. We address this
        deficit by proposing two incremental approximation algorithms (for weighted and unweighted connected graphs) which provide a provable guarantee on the absolute
        approximation error. Processing batches of edge insertions, our algorithms yield significant speedups up to a factor of 104 compared to restarting the approximation.
        This is enabled by investing memory to store and efficiently update shortest paths. As a building block, we also propose an asymptotically faster algorithm for
        updating the SSSP problem in unweighted graphs. Our experimental study shows that our algorithms are the first to make in-memory computation of a betweenness
        ranking practical for million-edge semi-dynamic networks. Moreover, our results show that the accuracy is even better than the theoretical guarantees in terms of
        absolutes errors and the rank of nodes is well preserved, in particular for those with high betweenness.
      </div>
    </li>

    <br>

    <li>
      C. Staudt, Y. Marrakchi, H. Meyerhenke: Detecting Communities Around Seed Nodes in Complex Networks. In <i>Proc. First International Workshop on High Performance
      Big Graph Data Management, Analysis, and Mining</i>, co-located with the <i>IEEE BigData 2014 conference</i>.
      <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
      <div id="collapseDiv" class="collapse">
        <b>Abstract.</b> The detection of communities (internally dense subgraphs) is a network analysis task with manifold applications. The special task of selective
        community detection is concerned with finding high-quality communities locally around seed nodes. Given the lack of conclusive experimental studies, we perform
        a systematic comparison of different previously published as well as novel methods. In particular we evaluate their performance on large complex networks,
        such as social networks. Algorithms are compared with respect to accuracy in detecting ground truth communities, community quality measures, size of communities
        and running time. We implement a generic greedy algorithm which subsumes several previous efforts in the field. Experimental evaluation of multiple objective
        functions and optimizations shows that the frequently proposed greedy approach is not adequate for large datasets. As a more scalable alternative, we propose
        selSCAN, our adaptation of a global, density-based community detection algorithm. In a novel combination with algebraic distances on graphs, query times can
        be strongly reduced through preprocessing. However, selSCAN is very sensitive to the choice of numeric parameters, limiting its practicality. The
        random-walk-based PageRankNibble emerges from the comparison as the most successful candidate.
      </div>
    </li>

    <br>

    <li>
      C. Staudt and H. Meyerhenke: Engineering Parallel Algorithms for Community Detection in Massive Networks. Accepted by <i>IEEE Transactions on Parallel and
      Distributed Systems</i> (TPDS). [<a href="http://arxiv.org/abs/1304.4453">arXiv</a>]
      [<a href="http://dx.doi.org/10.1109/TPDS.2015.2390633">DOI: 10.1109/TPDS.2015.2390633</a>] &#169; 2015 IEEE
      <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
      <div id="collapseDiv" class="collapse">
        Please note that <a href="https://algohub.iti.kit.edu/parco/NetworKit/NetworKit/archive/9cfacca668d8f4e4740d880877fee34beb276792.zip?subrepos=false">NetworKit 3.5</a>
        is the last version to include an implementation of the EPP algorithm. <br><br>

        <b>Abstract</b> The amount of graph-structured data has recently experienced an enormous growth in many applications. To transform such data into useful
        information, fast analytics algorithms and software tools are necessary. One common graph analytics kernel is disjoint community detection (or graph clustering).
        Despite extensive research on heuristic solvers for this task, only few parallel codes exist, although parallelism will be necessary to scale to the data volume
        of real-world applications. We address the deficit in computing capability by a flexible and extensible community detection framework with shared-memory parallelism.
        Within this framework we design and implement efficient parallel community detection heuristics: A parallel label propagation scheme; the first large-scale
        parallelization of the well-known Louvain method, as well as an extension of the method adding refinement; and an ensemble scheme combining the above.
        In extensive experiments driven by the algorithm engineering paradigm, we identify the most successful parameters and combinations of these algorithms. We also
        compare our implementations with state-of-the-art competitors. The processing rate of our fastest algorithm often reaches 50M edges/second. We recommend the
        parallel Louvain method and our variant with refinement as both qualitatively strong and fast. Our methods are suitable for massive data sets with billions of
        edges.
      </div>
    </li>

    <br>

    <li>
      C. Staudt and H. Meyerhenke: Engineering High-Performance Community Detection Heuristics for Massive Graphs. In: <i>Proceedings of the 2013 International Conference on
      Parallel Processing</i>. [<a href="http://arxiv.org/abs/1304.4453">updated and extended version on arXiv</a>,
      <a href="https://networkit.iti.kit.edu/data/uploads/publications/sm2013ehpcdh.bib">bibtex</a>]
      <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
      <div id="collapseDiv" class="collapse">
        <b>Abstract</b> The amount of graph-structured data has recently experienced an enormous growth in many applications. To transform such data into useful information,
        high-performance analytics algorithms and software tools are necessary. One common graph analytics kernel is community detection (or graph clustering). Despite
        extensive research on heuristic solvers for this task, only few parallel codes exist, although parallelism will be necessary to scale to the data volume of
        real-world applications. We address the deficit in computing capability by a flexible and extensible community detection framework with shared-memory parallelism.
        Within this framework we design and implement efficient parallel community detection heuristics: A parallel label propagation scheme; the first large-scale
        parallelization of the well-known Louvain method, as well as an extension of the method adding refinement; and an ensemble scheme combining the strengths of the
        above. In extensive experiments driven by the algorithm engineering paradigm, we identify the most successful parameters and combinations of these algorithms.
        We also compare our implementations with state of the art competitors. The processing rate of our fastest algorithm often exceeds 10M edges/second, making it
        suitable for massive data streams. We recommend the parallel Louvain method and our variant with refinement as both qualitatively strong and relatively fast.
        Moreover, our fast ensemble algorithm yields a good tradeoff between quality and speed for community detection in very large networks.
      </div>
    </li>

  </ul>


|separator|


Publications Using NetworKit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. raw:: html

  <ul>
    <li>
      M. Lozano, C. García-Martínez, F. J. Rodríguez, H. M. Trujillo: Optimizing network attacks by artificial bee colony.
      To appear in <i>Information Sciences, Volume 377</i>, pp. 30-50, January 2017.
       <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
       <div id="collapseDiv" class="collapse">
         <b>Abstract</b> Over the past few years, the task of conceiving effective attacks to complex networks has arisen as an optimization problem. Attacks are modelled as the process of removing a number k of vertices, from the graph that represents the network, and the goal is to maximise or minimise the value of a predefined metric over the graph. In this work, we present an optimization problem that concerns the selection of nodes to be removed to minimise the maximum betweenness
         centrality value of the residual graph. This metric evaluates the participation of the nodes in the communications through the shortest paths of the network.
        <br>
         To address the problem we propose an artificial bee colony algorithm, which is a swarm intelligence approach inspired in the foraging behaviour of honeybees. In this framework, bees produce new candidate solutions for the problem by exploring the vicinity of previous ones, called food sources. The proposed method exploits useful problem knowledge in this neighbourhood exploration by considering the partial destruction and heuristic reconstruction of selected solutions. The
         performance of the method, with respect to other models from the literature that can be adapted to face this problem, such as sequential centrality-based attacks, module-based attacks, a genetic algorithm, a simulated annealing approach, and a variable neighbourhood search, is empirically shown.
       </div>
    </li>

    <br>

    <li>
      M. Riondato, E. Upfal: ABRA: Approximating Betweenness Centrality in Static and Dynamic Graphs with Rademacher Averages.
      To appear in <i>Proc. 22nd ACM SIGKDD Conference on Knowledge Discovery and Data Mining</i> (KDD 2016), August 2016. [<a href="http://arxiv.org/abs/1602.05866">arXiv</a>]
       <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
       <div id="collapseDiv" class="collapse">
         <b>Abstract</b> We present ABRA, a suite of algorithms that compute and maintain probabilistically-guaranteed, high-quality, approximations of the betweenness centrality of all nodes
          (or edges) on both static and fully dynamic graphs. Our algorithms rely on random sampling and their analysis leverages on Rademacher averages and pseudodimension, fundamental
          concepts from statistical learning theory. To our knowledge, this is the first application of these concepts to the field of graph analysis. The results of our experimental evaluation
           show that our approach is much faster than exact methods, and vastly outperforms, in both speed and number of samples, current state-of-the-art algorithms with the same quality guarantees.
       </div>
    </li>

    <br>

    <li>
      M. von Looz, M. Wolter, C. Jacob, H. Meyerhenke: Better partitions of protein graphs for subsystem quantum chemistry. To appear in <i>Proc. 15th Intl. Symp. on Experimental
      Algorithms</i> (SEA 2016), June 2016.
      <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
      <div id="collapseDiv" class="collapse">
        <b>Abstract</b> Determining the interaction strength between proteins and small molecules is key to analyzing their biological function. Quantum-mechanical calculations
        such as Density Functional Theory (DFT) give accurate and theoretically well-founded results. With common implementations the running time of DFT calculations increases
        quadratically with molecule size. Thus, numerous subsystem-based approaches have been developed to accelerate quantum-chemical calculations. These approaches partition
        the protein into different fragments, which are treated separately. Interactions between different fragments are approximated and introduce inaccuracies in the
        calculated interaction energies. <br>
        To minimize these inaccuracies, we represent the amino acids and their interactions as a weighted graph in order to apply graph partitioning. None of the existing graph
        partitioning work can be directly used, though, due to the unique constraints in partitioning such protein graphs. We therefore present and evaluate several algorithms,
        partially building upon established concepts, but adapted to handle the new constraints. For the special case of partitioning a protein along the main chain, we
        also present an efficient dynamic programming algorithm that yields provably optimal results. In the general scenario our algorithms usually improve the previous
        approach significantly and take at most a few seconds.
      </div>
    </li>

    <br>

    <li>
      P. Crescenzi, G. D’Angelo, L. Severini, Y. Velaj: Greedily Improving Our Own Centrality in A Network. In <i>Proc. 14th Intl. Symp. on Experimental Algorithms</i> (SEA 2015).
       LNCS 9125, pp. 43-55. Springer International Publishing, 2015.
       <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
       <div id="collapseDiv" class="collapse">
         <b>Abstract</b> The closeness and the betweenness centralities are two well-known measures of importance of a vertex within a given complex network. Having high
         closeness or betweenness centrality can have positive impact on the vertex itself: hence, in this paper we consider the problem of determining how much a vertex
         can increase its centrality by creating a limited amount of new edges incident to it. We first prove that this problem does not admit a polynomial-time approximation
         scheme (unless \(P=NP\)), and we then propose a simple greedy approximation algorithm (with an almost tight approximation ratio), whose performance is then tested
         on synthetic graphs and real-world networks.
       </div>
    </li>

    <br>

    <li>
      D. Hoske, D. Lukarski, H. Meyerhenke, M. Wegner: Is Nearly-linear the same in Theory and Practice? A Case Study with a Combinatorial Laplacian Solver. In <i>Proc. 14th Intl.
      Symp. on Experimental Algorithms</i> (SEA 2015). LNCS 9125, pp. 205-218. Springer International Publishing, 2015. [<a href="http://arxiv.org/abs/1502.07888">arXiv</a>]
       <button type="button" class="btn-link collapsed" data-toggle="collapse" data-target="#collapseDiv"></button>
       <div id="collapseDiv" class="collapse">
       For the paper follow the arXiv link above. If you are interested in the implementation, see ParCo's <a href="http://parco.iti.kit.edu/software-en.shtml" >software page</a>.
       </div>
    </li>

  </ul>


|separator|

Projects Using NetworKit
~~~~~~~~~~~~~~~~~~~~~~~~

Further projects using NetworKit can be found `here <projects.html>`_.


|separator|

Student Theses Using NetworKit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A list of student theses based on NetworKit can be found `here <student_theses.html>`_.
