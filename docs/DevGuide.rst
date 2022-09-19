.. _devGuide:

NetworKit Development Guide
===========================

This guide is designed to give potential contributors a quick understanding of the structure of NetworKit as well as what they need to look out for when adding new functionality to NetworKit.

The following text assumes some basic familiarity with the Git
version control software. It is not a Git tutorial, because you
will find a good one `here <https://try.github.io/>`__. Rather, it
explains concepts and workflows for the development of this project.

If you want to contribute, you should consider reading the `technical
report <https://arxiv.org/pdf/1403.3005.pdf>`__ on NetworKit to get
familiar with the architecture.

If you use NetworKit in your research publications, please cite the
mentioned technical report or the specific algorithm. A list of
publications is available on the
`website <https://networkit.github.io/publications.html>`__.

Contribution Workflow
---------------------

The development of the NetworKit project takes place at our `official repository on GitHub <https://github.com/networkit/networkit>`__.

In order to successfully contribute to the project and add your feature to the master branch, there are a couple of steps you need to follow:

1. Fork the NetworKit repository
2. Create a feature branch
3. Add your code
4. Send a pull request (against the master branch)
5. Participate in the pull request discussion

In the following part we will take a closer look at each step of the workflow outlined above.
To make it more tangible, let's assume you want to contribute to NetworKit by adding a new algorithm to the existing suite.

1. Fork the NetworKit repository
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Visit our official repository on GitHub and click the 'Fork' button in the top right to create a personal fork in your account.

Next up, you want to clone your fork to your local machine:

::

   git clone git@github.com:[YOUR-USERNAME]/[FORKED-NETWORKIT].git

Since the algorithm you're working on will take some time to implement, you want to make sure to stay up to date with the original NetworKit repository. This can be accomplished by tracking the original ``upstream`` branch:

::

   git remote add upstream https://github.com/networkit/networkit.git

Now you can update your fork on a regular basis with changes from the original repository:

::

   # Fetch changes from original repository
   git fetch upstream

   # Merge changes into your local branch
   git checkout --track origin/master
   git merge upstream/master

This will for example allow you to keep the ``master`` branch up to date. Doing this will prevent an unpleasent surprise once you're ready to submit your pull request since there will be less or no merge conflicts.

2. Create a Feature Branch
~~~~~~~~~~~~~~~~~~~~~~~~~~

To make sure that the work on your algorithm does not disrupt the development process, be sure to create your own branch where you can add new code.
In order to have a coherent naming scheme that allows for easier communication, please name your branch in the following way: ``feature/descriptive-name-of-feature``.
If we're assuming you're creating an implementation of the shortest path algorithm from Dijkstra, this could be ``feature/dijkstra-shortest-path`` for example. Make sure to pick a meaningful and short description that is easy to understand.

New feature branches should be based off the ``master`` branch where they will be merged back into at a later stage.

::

   # Make sure your branch is based off the master branch
   git checkout master

   # Create your new branch
   git checkout -b feature/[my-awesome-feature-name]

   # You're now on your new feature branch

3. Add Your Code
~~~~~~~~~~~~~~~~

In this step you're going to make and commit the changes needed for your new feature. Please make sure to write clean code that adheres to the style guide outlined further below. It is also important that each feature has appropriate unit tests that cover all of the expected behaviour of the code. Please see the Test-driven development section below for details.

::

   git checkout feature/[my-awesome-feature-name]
   git add [files]
   git commit -m "[descriptive message about the changes you made]"
   git push

Also, from time to time, you should make sure to keep your feature branch up to date with the changes on the ``master`` branch in the main repository. If you followed step 1 from above, this can be easily accomplished:

::

   git checkout feature/[my-awesome-feature-name]
   git fetch upstream
   git merge upstream/master


4. Send a Pull Request
~~~~~~~~~~~~~~~~~~~~~~

Once you finished the development and testing of your new feature, it is time to create a pull request to get your changes merged into the master branch of the NetworKit repository.

This can be done by visiting the **Pull requests page** (https://github.com/[YOUR-USERNAME]/[FORKED-NETWORKIT]/pulls) of your NetworKit fork on GitHub and clicking on the green **New pull request** button at the top right side of the page.

Here the ``base fork`` at the top should point to ``networkit/networkit`` and the base should be ``master``. The ``head fork`` should point to your fork of networkit and the ``compare`` branch to the right should point to the feature branch (``feature/[my-awesome-feature-name]``) you would like to create the pull request for.

Once you've reviewed all changes, click the green **Create pull request** button and your pull request will be created.


5. Participate in the Pull Request Discussion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once your pull request has been submitted, other developers of NetworKit, including the core development team, will take a look at your code and explanation.
In this process there are oftentimes questions that arise or small adjustments that need to be incorporated into the pull request. For this reason, it is important that you actively participate in the discussions around your pull request. This ensures your new feature will eventually make it to the next release.

In case a developer points out a potential issue that needs to be resolved, please make the appropriate changes to your code and push these changes to your feature branch:

::

   git checkout feature/[my-awesome-feature-name]
   # Make appropriate changes to files
   git add [files]
   git commit -m "[message-about-the-resolved-issue]"
   git push


The pull request will automatically show your newest changes and developers will know that you resolved the issue. Once all issues have been resolved and your code is accepted, the pull request will be closed and your feature will be merged into the master branch. In the next release, all users of NetworKit will have access to your awesome feature. Hooray!

Style guide
-----------

We want to ensure that code across NetworKit is easy to understand for existing as well as new developers. This is why new code added to the project should adhere to the existing code style. At this point in time, there is no comprehensive documentation about the code style being used in NetworKit but there are a few things to look out for:

-  Compiler warnings are likely to turn into future errors. Try to fix
   them as soon as they appear. Use the ``-Wall`` flag when compiling C++ code.
-  Read some code to get used to the code style and try to adopt it.
-  Document classes, methods and attributes in Doxygen style.
-  Use the ``count`` and ``index`` integer types for non-negative
   integer quantities and indices.
-  All member variables should be pointers and not references.
-  In most cases, objects are passed by reference. New objects are
   stack-allocated and returned by value. Avoid ``new``
   where possible.
-  Use the ``override`` keyword to indicate that a method overrides a
   virtual method in the superclass.
-  A class should be declared ``final`` unless it is a superclass.
-  ``virtual`` methods should only be declared in superclasses.
-  In Python, indent using tabs, not spaces.

In order to maintain the same standard of code across the entire NetworKit code base, some coding standards are enforced. However, there is some automation to help developers with this. Below is a list of these standards and instructions on how to use the available automation tools that ensure your code adheres to them.

-  ``CppClangFormat`` applies clang-format to all C++ files.
-  ``CppIndentation`` checks that all C++ code is indented with spaces
   and not tabs.
-  ``CppIncludeGuards`` ensures that the header files contain an include guard and
   that it follows the following naming convention: ``NETWORKIT_MODULENAME_HEADERFILE_HPP_``.

The executable file ``check_code.sh`` in NetworKit's root directory carries out all checks in read-only mode and reports if errors are found. Running ``./check_code.sh -w`` will fix these errors. Run this script before commiting your files to make sure your changes are in complaince with the guidelines. The script is executed during CI and will cause your pull request to fail if your code does not conform to the style guide.

On top of the aforementioned mentioned points concerning style, the NetworKit C++ code base also complies to a selection of ``clang-tidy`` static-code analysis checks.
New code must pass these tests before being merged into the development branch. The list of checks can be found in the ``.clang-tidy`` file.
In order to run the ``clang-tidy`` checks while building NetworKit, set the ``CMake`` flag ``-NETWORKIT_CLANG_TIDY`` to ``ON`` in addition to the other compile flags, e.g.
::

    cmake -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Debug -DNETWORKIT_WARNINGS_AS_ERRORS=ON -DNETWORKIT_CLANG_TIDY=ON ..

In a nutshell, new developers should familiarise themselves with the existing code base and adapt the existing style in the C++ as well as Python code base when contributing to NetworKit. Always ensure your code is easy to understand and properly documented.

Report bugs
-----------

Please report any bugs on the `issues page <https://github.com/networkit/networkit/issues>`__ of the official NetworKit repository on GitHub.
In very urgent cases it might also make sense to write on the `mailing list <https://sympa.cms.hu-berlin.de/sympa/subscribe/networkit>`__.
Please provide a minimal example so that others can reproduce that bug.

Tags
----

A tag is nothing more than a “symbolic name” for a revision. In
NetworKit tags are used to mark release versions in the ``master``
branch, with a ``MAJOR.MINOR`` version name scheme.

Using NetworKit
---------------

If you want to build and use NetworKit and do not plan to contribute
changes, simply clone the repository. By default, you will be on the
``master`` branch, which represents the current release. Follow the
setup instructions in the ``README`` file.

Student Projects
~~~~~~~~~~~~~~~~

Students with long-term projects like Bachelor's or Master's theses
should familiarize themselves with the guidelines and select a
forking/branching model with their advisor.

Branching Cheat Sheet
---------------------

-  list all available branches with highlight for the current branch: ``git branch``
-  switch to a specific branch: ``git checkout <branchname>``
-  start a new branch: ``git checkout -b <branchname>``
-  merge ``branch-y`` into ``branch-x``: ``git checkout branch-x``, then
   ``git merge branch-y``
-  see heads (most recent commits) of all branches: ``git show-ref --heads``

Conventions
-----------

The following general conventions apply to all NetworKit developers.

Versioning
~~~~~~~~~~

-  Before you commit, make sure your code compiles and run the unit
   tests. Never push code which breaks the build for others.
-  Commit regularly and often to your local repository.
-  Use meaningful commit messages.
-  Get the newest changes from the repository regularly and merge them
   into your local repository.
-  Make sure that you merged correctly and did not break other people's
   work.
-  Push correct code early if possible. Merging is easier if all
   developers are up to date.
-  Never ``push --force`` to the main repository.

.. _devGuide-unitTests:

Unit Tests and Testing
----------------------

Every new feature must be covered by a unit test. Omitting unit tests
makes it very likely that your feature will break silently as the
project develops, leading to unneccessary work in tracing back the
source of the error. Also your pull request for this feature will most
likely not be accepted.

Unit tests for the C++ part of NetworKit are based on the ``googletest``
library. For more information read the `googletest
primer <https://google.github.io/googletest/primer.html>`__. The Python
test framework currently relies on ``nose`` to collect the tests.

-  Each source folder contains a ``test`` folder with ``googletest``
   classes. Create the unit tests for each feature in the appropriate
   ``test/*GTest`` class by adding a ``TEST_F`` function.
-  Prefix standard unit tests with ``test`` and experimental feature
   tests with ``try``. A ``test*`` must pass when pushed to the main
   repository, a ``try*`` is allowed to fail.
-  Keep the running time of test functions to the minimum needed for
   testing functionality. Testing should be fast, long-running unit
   tests look like infinite loops.
-  If the unit test requires a data set, add the file to the ``input/``
   folder. Only small data sets (a few kilobytes maximum) are acceptable
   in the repository.
-  Any output files produced by unit tests must be written to the
   ``output/`` folder.

To build and run the tests you need the `gtest
library <https://github.com/google/googletest>`__. Assuming, gtest is
successfully installed and you add the paths to your build.conf, the
unit tests should be compiled with:

::

   cd build/
   cmake -DNETWORKIT_BUILD_TESTS=ON ..
   make -jX # To speed up the compilation with make a multi-core machine, you can append `-jX` where X denotes the number of threads to compile with.

To verify that the code was built correctly: Run all unit tests with

::

   ctest -V

To select only a subset of tests, you can run instead

::

   cd .. # Navigate to the project root directory
   build/networkit_tests [options]

Here's a rundown of the available options: Non-performance tests can be selected with

::

   build/networkit_tests --tests/-t

while performance tests are called with

::

   build/networkit_tests --benchmarks/-b

Further options are:

::
   -r, --run         Run unit tests which don't use assertions
   -d, --debug       Run tests to debug some algorithms
       --threads     set the maximum number of threads; 0 (=default) uses OMP 
                     default
       --loglevel    set the log level (TRACE|DEBUG|INFO|WARN|ERROR|FATAL)
       --srcloc      print source location of log messages

To run only specific unit tests, you can also add a filter expression,
e. g.:

::

   build/networkit_tests --gtest_filter=*PartitionGTest*/-f*PartitionGTest*

initiates unit tests only for the Partition data structure.

For the **Python** unit tests, run:

::

   python3 setup.py test [--cpp-tests/-c]

This command will compile the \_NetworKit extension and then run all
test cases on the Python layer. If you append ``--cpp-tests/-c``, the
unit tests of the c++ side will be compiled and run before the Python
test cases.

Test-driven development
~~~~~~~~~~~~~~~~~~~~~~~

If you implement a new feature for NetworKit, we encourage you to adapt
your development process to test driven development. This means that you
start with a one or ideally several test-cases for your feature and then
write the feature for the test case(s). If your feature is mostly
implemented in C++, you should write your test cases there. If you
expose your feature to Python, you should also write a test case for the
extension module on the Python layer. The same applies for features in
Pyton.


Algorithm interface and class hierarchy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We use the possibilities provided through inheritance to generalize the
common behaviour of algorithm implementations:

-  Data and paramters should be passed in the constructor.
-  A void run()-method that takes no parameter triggers the execution.
-  To retrieve the result(s), getter-functions() may be defined.

The ``Algorithm`` base class also defines a few other other functions to
query whether the algorithm can be run in parallel or to retrieve a
string representation.

There may be more levels in the class hierarchy between an algorithm
implementation and the base class, e.g. a single-source shortest-path
class ``SSSP`` that generalizes the behaviour of BFS and Dijkstra
implementations or the ``Centrality`` base class. When implementing new
features or algorithms, make sure to adapt to the existing class
hierarchies. The least thing to do is to inherit from the ``Algorithm``
base class. Changes to existing interfaces or suggestions for new
interfaces should be discussed through the `mailing
list <https://sympa.cms.hu-berlin.de/sympa/subscribe/networkit>`__.

Exposing C++ Code to Python
---------------------------

Assuming the unit tests for the new feature you implemented are correct
and successful, you need to make your feature available to Python in
order to use it. NetworKit uses Cython to bridge C++ and Python. All of
this bridge code is contained in the ``networkit/`` directory. The Cython
files in this directory correspond to the C++ modules. Files with a ``.pxd``
extension declare C++ data types, functions and variables that are imported
by other files. Therefore, if the new code does not introduce new C++ types
or functions that are needed elsewhere, the code should only be added to the
correct ``.pyx`` file. The content is automatically translated into C++ and
then compiled to a Python extension module.

Cython syntax is a superset of Python that knows about static type
declarations and other things from the C/C++ world. The best way to
getting used to it is working on examples. Take the most common case of
exposing a C++ class as a Python class. Care for the following example
that exposes the class ``NetworKit::Dijkstra`` in ``distance.pyx``:

::

        [...]
        from .base cimport _Algorithm, Algorithm
        from .graph cimport _Graph, Graph
        [...]

In order to inherit from ``Algorithm`` and use the ``Graph`` data structure,
we must import the C++ and Python types like is done above.

::

        cdef extern from <networkit/distance/Dijkstra.hpp>:
            cdef cppclass _Dijkstra "NetworKit::Dijkstra"(_SSSP):
                _Dijkstra(_Graph G, node source, bool_t storePaths, bool_t storeNodesSortedByDistance, node target) except +

The code above exposes the C++ class definition to Cython - but not yet
to Python. First of all, Cython needs to know which C++ declarations to
use so the the first line directs Cython to place an ``#include``
statement. The second line defines a class that is only accessible in
the Cython world. Our convention is that the name of the new class is
the name of the referenced C++ class with a prepended underscore to
avoid namespace conflicts. What follows is the "real" C++ name of the
class. After that, the declarations of the methods you want to make
available for Python are needed. The ``except +`` statement is necessary
for exceptions thrown by the C++ code to be rethrown as Python
exceptions rather than causing a crash. Also, take care that the Cython
declarations match the declarations from the referenced header file.

::

        cdef extern from <networkit/distance/_SSSP.hpp>:
            cdef cppclass _SSSP "NetworKit::SSSP"(_Algorithm):
                _SSSP(_Graph G, node source, bool_t storePaths, bool_t storeNodesSortedByDistance, node target) except +
                vector[edgeweight] getDistances(bool moveOut) except +
                [...]

        cdef class SSSP(Algorithm):
            """ Base class for single source shortest path algorithms. """

            cdef Graph _G

            def __init__(self, *args, **namedargs):
                if type(self) == SSSP:
                    raise RuntimeError("Error, you may not use SSSP directly, use a sub-class instead")

            def __dealloc__(self):
                self._G = None # just to be sure the graph is deleted

            def getDistances(self, moveOut=True):
                """
                Returns a vector of weighted distances from the source node, i.e. the
                length of the shortest path from the source node to any other node.

                Returns
                -------
                vector
                    The weighted distances from the source node to any other node in the graph.
                """
                return (<_SSSP*>(self._this)).getDistances(moveOut)
            [...]

We mirror the class hierarchy of the C++ world also in Cython and
Python. This also saves some boiler plate wrapping code as the functions
shared by Dijkstra and BFS only need to be wrapped through SSSP.

::

        cdef class Dijkstra(SSSP):
            """ Dijkstra's SSSP algorithm.

            Returns list of weighted distances from node source, i.e. the length of the shortest path from source to
            any other node.

            Dijkstra(G, source, [storePaths], [storeStack], target)

            Creates Dijkstra for `G` and source node `source`.

            Parameters
            ----------
            G : networkit.Graph
                The graph.
            source : node
                The source node.
            storePaths : bool
                store paths and number of paths?
            storeStack : bool
                maintain a stack of nodes in order of decreasing distance?
            target : node
                target node. Search ends when target node is reached. t is set to None by default.
            """
            def __cinit__(self, Graph G, source, storePaths=True, storeStack=False, node target=none):
                self._G = G
                self._this = new _Dijkstra(G._this, source, storePaths, storeStack, target)

For the class to be accessible from the Python world, you need to define
a Python wrapper class which delegates method calls to the native class.
The Python class variable ``_this`` holds a pointer to an instance of
the native class. Please note that the parameters are now Python
objects. Method wrappers take these Python objects as parameters and
pass the internal native objects to the actuall C++ method call. The
constructor of such a wrapper class is called ``__cinit__``, and it
creates an instance of the native object.

The docstring between the triple quotation marks can be accessed through
Python's ``help(...)`` function and are the main documentation of
NetworKit. Always provide at least a short and precise docstring so the
user can get in idea of the functionality of the class. For C++ types
available to Python and further examples, see through the various
Cython files. The whole process has certainly some
intricacies, e.g. some tricks are needed to avoid memory waste when
passing around large objects such as graphs. When in doubt, look at
examples of similar classes already exposed. Listen to the Cython
compiler - coming from C++, its error messages are in general pleasantly
human-readable.

Adding new C++ Modules to Python
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After having added a new C++ module to NetwoKit, you need to do the same
in Cython before the new code is available in Python.
According to the NetwoKit naming conventions, the corresponding Cython file
should have the same name as your C++ module. For illustration purposes,
assume the new C++ module is called "features". The new the Cython file is,
therefore, also named ``features.pyx``. In order to expose "features" to Python,
the following must be done:

    #. Create the ``features.pyx`` file in the ``networkit/networkit`` directory and
       add the following line:

        ::

                # distutils: language=c++
    #. Add the string "features" to the CMakeLists file in the root directory where 
       the other Python extension modules are listed. At this point the new file
       should be found and compiled.
    #. Export the C++ code as explained in the previous section to the new ``features.pyx`` file.
       - Note that you will need to import all the NetworKit and C++ data types
       needed in your code.
       - ``graph.pxd`` and ``generators.pyx`` are a good starting point to find various
       types as they import several data types.
    #. Build using the usual commands, i.e.

          ::

            python3 setup.py build_ext [-jX]
            pip3 install -e .
    #. Once the features shared object is successfully created, import it in the
       ``networkit/__init__.py`` file as is done with the other modules.
       The ``features`` module has now been exposed to Python.

Make algorithms interruptable with CTRL+C/SIGINT
------------------------------------------------

When an algorithms takes too long to produce a result, it can be
interrupted with a SIGINT signal triggered by CTRL+C. When triggering
from the Python shell while the runtime is in the C++ domain, execution
is aborted and even terminates the Python shell. Therefore, we
implemented a signal handler infrastructure in C++ that raises a special
exception instead of aborting. When implementing an algorithm, it is
strongly encouraged to integrate the signal handler into the
implementation. There are many examples of how to use it, e.g.
``networkit/cpp/centrality/Betweenness.cpp`` or
``networkit/cpp/community/PartitionFragmentation.cpp``

Contact
-------

To discuss important changes to NetworKit, use the `e-mail
list <https://sympa.cms.hu-berlin.de/sympa/subscribe/networkit>`__
(``networkit@lists.hu-berlin.de``).

We also appreciate new issues or pull requests on the GitHub repository.

Building the documentation
--------------------------

The documentation can be automatically generated with sphinx. You will need the following software to generate the documentation:

-  `Sphinx <http://www.sphinx-doc.org>`__ (e.g. via ``pip3 install sphinx``)
-  `Pandoc <http://pandoc.org>`__
-  `Doxygen <http://www.stack.nl/~dimitri/doxygen/>`__

After you installed the above mentioned software, you can build the
class documentation by calling ``./make_doc.sh`` in the folder
``Doc/doc``. This will generate the class documentation for C++ and
Python in ``Doc/Documentation``.

Further Reading
---------------

-  `Interactive Git tutorial <https://try.github.io/>`__
-  `Working with named
   branches <http://humblecoder.co.uk/blog/2010/02/24/working-with-named-branches-in-mercurial/>`__
-  `Managing releases and branchy
   development <http://hgbook.red-bean.com/read/managing-releases-and-branchy-development.html>`__
-  `Cython Documentation <http://docs.cython.org/index.html>`__
