#!/bin/bash
#
# Use this script to create the documentation for NetworKit

rm -rf ../../../NetworKit-Doc
doxygen c++/Doxyfile

# convert markdowns to rsts
pandoc --from=markdown --to=rst --output=python/source/Readme.rst ../../Readme.mdown
pandoc --from=markdown --to=rst --output=python/source/DevGuide.rst ../DevGuide.mdown

make -C python html

# leave the repository in a clean state
rm -f debug.txt
rm python/source/Readme.rst
rm python/source/DevGuide.rst
echo ""
echo "Finished building documentation to ../../../Networkit-Doc"
echo "Open the documentation for c++ with 'open ../../../Networkit-Doc/doxyhtml/index.html'"
echo "Open the documentation for python with 'open ../../../Networkit-Doc/python/html/index.html'"
