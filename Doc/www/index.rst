.. title:: NetworKit

.. toctree::
   :hidden:
   :maxdepth: 2

   Developer Guide <api/DevGuide>
   IPython Notebook <api/notebooks>
   Python Documentation <api/modules>
   C++ Documentation <api/cppdoc>

.. raw:: html

   <section class="Top_Section" style="clear: both; border-bottom: 1px solid #d4d7d9;">
      <div class="Top_Section" style="padding-top: 30px; padding-bottom: 30px">
        <div class="Introduction_Text" style="border-right: 1px solid #d4d7d9; display: table-cell; width: 66.66%; padding-right: 30px">

**NetworKit** is a growing open-source toolkit for high-performance network analysis. Its aim is to provide tools for the analysis of large networks in the size range from thousands to billions of edges. For this purpose, it implements efficient graph algorithms, many of them parallel to utilize multicore architectures. These are meant to compute standard measures of network analysis, such as degree sequences, clustering coefficients and centrality. In this respect, NetworKit is comparable to packages such as NetworkX, albeit with a focus on parallelism and scalability. NetworKit is also a testbed for algorithm engineering and contains novel algorithms from recently published research, especially in the area of community detection (see list of :ref:`publications`).

**NetworKit** is a Python module. High-performance algorithms are written in C++ and exposed to Python via the Cython toolchain. Python in turn gives us the ability to work interactively and with a rich environment of tools for data analysis and scientific computing. Furthermore, NetworKit's core can be built and used as a native library if needed.


.. raw:: html

        </div>
        <div class="Downloads" style="display: table-cell; width: 33.33%; padding-left: 30px">
          <div>Clone from Mercurial</div>
          <span style="display: block;overflow: hidden;"><input onClick="this.setSelectionRange(0, this.value.length)" style="width: 100%" type="text" value="hg clone https://algohub.iti.kit.edu/parco/NetworKit/NetworKit" readonly=""/></span>

          <div style="padding-top: 15px">Install via pip</div>
          <span style="display: block;overflow: hidden;"><input onClick="this.setSelectionRange(0, this.value.length)" style="width: 100%" type="text" value="pip install networkit" readonly=""/></span>

          <div style="padding-top: 15px">Download from <a href="https://algohub.iti.kit.edu/parco/NetworKit/NetworKit">Algohub</a> or as <a href="https://networkit.iti.kit.edu/data/uploads/networkit.zip">Zip file</a></div>

          <div style="padding-top: 15px">Download the <a href="https://networkit.iti.kit.edu/data/uploads/NetworKit-Doc.zip">Documentation</a> or <a href="https://networkit.iti.kit.edu/data/uploads/publications/staudt2014networkit.pdf">Technical Report</a></div>

          <div style="padding-top: 15px;"> <div style="float: left;">Mailing List</div> <div><a style="padding-left: 10px" href="https://lists.ira.uni-karlsruhe.de/mailman/listinfo/networkit"><img style="padding-bottom:2px" src="_static/mailinglist.png"></a> </div> </div>

          <div style="padding-top: 15px">View the <a href="https://lists.ira.uni-karlsruhe.de/pipermail/networkit/">mailing list archive</a></div>

          <div style="padding-top: 15px"><a href="http://nbviewer.ipython.org/urls/networkit.iti.kit.edu/data/uploads/docs/NetworKit_UserGuide.ipynb">NetworKit UserGuide</a></div>

        </div>
      </div>
    </section>

    <section class="MainFeatures" style="clear: both;border-bottom: 1px solid #d4d7d9; padding-top: 20px; padding-bottom: 20px;">
      <div class="FeatureTable" >
        <div style="text-align: center; font-size:16pt; font-weight: bold;">Main Features</div>
        <div style="border-right: 1px solid #d4d7d9; display: table-cell; width: 33.33%; padding: 30px; padding-bottom: 0px;">
          <p style="text-align: center; font-size:14pt">High Performance</p>
          <p style="word-break: normal; text-align:justify;">
            In NetworKit, algorithms and data structures are selected and implemented with high performance and parallelism in mind. Some implementations are among the fastest in published research. For example, community detection in a 3 billion edge web graph can be performed on a 16-core server in a matter of minutes.
          </p>
        </div>

        <div style="border-right: 1px solid #d4d7d9; display: table-cell; width: 33.33%; padding: 30px; padding-bottom: 0px;">
          <p style="text-align: center; font-size:14pt">Interactive Workflow</p>
          <p style="word-break: normal; text-align:justify;">
            NetworKit takes inspiration from other software like R, MATLAB or Mathematica and provides an interactive shell via Python. This allows users to freely combine functions from NetworKit and even use the results ith other popular Python packages. In combination with IPython Notebook, NetworKit provides an easy computing environment for scientific workflows, even on a remote compute server.
          </p>
        </div>

        <div style="display: table-cell; width: 33.33%; padding: 30px; padding-bottom: 0px;">
          <p style="text-align: center; font-size:14pt">Easy Integration</p>
          <p style="word-break: normal; text-align:justify;">
            As a Python module, NetworKit enables seamless integration with Python libraries for scientific computing and data analysis, e.g. pandas for data framework processing and analytics, matplotlib for plotting, networkx for additional network analysis tasks, or numpy and scipy for numerical and scientific computing. Furthermore, NetworKit aims to support a variety of input/output formats.
          </p>
        </div>
      </div>
    </section>

    <section class="NewsVideo" style="clear: both; padding-top: 20px; padding-bottom: 20px">
      <div class="NewsVideo" >
        <div style="float: left; display: table-cell; width: 33.33%; padding-right: 30px">
          <p style="text-align: left; font-size:14pt">Latest News</p>

:ref:`news-1`

:ref:`news-2`

:ref:`news-3`

:ref:`news-4`

:ref:`news-5`

.. raw:: html

          <p style="word-break: normal; padding-top: 15px">
            <a href="news.html">All News</a>
          </p>
        </div>

        <div style="display: table-cell; border-left: 1px solid #d4d7d9; width: 66.66%; padding-left: 30px"> 
          <iframe class="youtube-player" type="text/html" width="640" height="385" src="https://www.youtube.com/embed/RtZyHCGyeIk?html5=1" frameborder="0"></iframe>
        </div>
      </div>
    </section>
