# Welcome to the CSE 21 Minitutorial!

## Topic: Large-Scale Network Analysis using NetworKit

##### __Session 1__: Friday, March 5, 8:30 a.m. - 10.10 a.m. CST
##### __Session 2__: Friday, March 5, 10:20 a.m. - 12:00 p.m. CST

### Speakers

- Eugenio Angriman, Humboldt University of Berlin, Germany
- Henning Meyerhenke, Humboldt University of Berlin, Germany
- Alexander van der Grinten, Humboldt University of Berlin, Germany

### Jupyter Notebooks

In the following there are links to the notebooks, we will use today for the Minitutorial-session. You can either click on the links or navigate on the left side in the two folders: __session-1__ and __session-2__.

- Session 1 (by Eugenio Angriman):
  - [Graph IO](session-1/graphio.ipynb)
  - [Centrality](session-1/centrality.ipynb)
  - [Graph Generators](session-1/generators.ipynb)


- Session 2 (by Alexander van der Grinten)
  - [Community Detection (Label Propagation)](session-2/community-labelprop.ipynb)
  - [Community Detection (Louvain Method)](session-2/community-louvain.ipynb)
  - [Dynamic Graphs](session-2/dynamic-graphs.ipynb)
  - [Embeddings](session-2/embeddings.ipynb)


#### __Important Note!__

NetworKit is designed to work with billions of nodes and edges. Therefore you can easily produce examples, which take up a lot of memory. To ensure a good Minitutorial-experience for everyone involved in this cloud environment, please only use smaller graphs found in __input__. Also if you generate graphs on your own, try to be sensible with sizes. Since resource usage on your instance is capped, computation may not finish if it tries to allocate too much memory. 

In case this or something else unexpected happens, be sure to stop the running kernel (normally Python 3). This can be done by clicking te circle on the left side (icon directly underneath the folder-symbol). In the following overview either click "Shut Down All" besides the "KERNELS" tab or select the appropriate crashed notebook and click the "x" besides it. After you shut down the kernel, also restart the crashed notebook. If the system seems slow, it can also help to close notebooks (since every notebook is a running Python-kernel) and kernels.

#### Folder structure

Besides the two folders we are covering in our session, another folder called __base-notebooks__ consists of all other notebooks NetworKit provides for users and developers. The __input__-folder consists of example instances used for computation. Feel free to play around with everything :)
