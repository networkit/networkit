.. |separator| raw:: html

	<div style="padding-top: 25px; border-bottom: 1px solid #d4d7d9;"></div>

========
Projects
========

On this page we present projects that use our NetworKit tool suite.


Image Segmentation
~~~~~~~~~~~~~~~~~~

There are serveral approaches to handle segmenting images into its main parts. One approach is to represent the image as a graph and apply graph clustering algorithms to compute a segmentation. This project is based on NetworKit and demonstrates how the framework can be used to segment images. A detailed project description with a the basic idea behind the implemented algorithms can be found in this `pdf-file <https://networkit.iti.kit.edu/uploads/projects/networkit-imagesegmentation.pdf>`_. Furthermore, an interactive `iPython Notebook <http://nbviewer.ipython.org/urls/networkit.iti.kit.edu/uploads/projects/graph-based-segmentation.ipynb>`_ is also available.

The project can be found on `AlgoHub <https://algohub.iti.kit.edu/parco/NetworKit/NetworKit-ImageSegmentation>`_.


|separator|


Protein Interaction Networks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In his Master thesis, Patrick Flick analyzed protein-protein interaction (PPI) networks of human cells. This work specifically looks at PPI networks of various cell and tissue types, here called tissue-specific PPIs (TSPPIs).

This work follows the goal to gain insights into the structure of interactions as well as into the properties of specific groups of proteins inside the TSPPI networks. To that end, an analysis pipeline was implemented and efficient analysis algorithms were developed, which operate on a sub-graph representation for TSPPI networks. The graph properties of TSPPI networks and properties of certain classes of proteins in the network were investigated. This work then re-evaluated prior research on a large set of TSPPIs, and demonstrated that some previous conclusions have to be reconsidered. Finally, NetworKit community-detection algorithms were employed in order to identify tissue-specific functional modules within TSPPIs.

The code, the thesis and more information is available on `github <https://github.com/r4d2/tsppi>`_.
