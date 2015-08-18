#!/bin/bash
#
# Use this script to create the documentation for NetworKit

rm -rf ../../../Networkit-Doc/*
doxygen c++/Doxyfile
make -C python html
rm -f debug.txt
echo ""
echo "Finished building documentation to ../../../Networkit-Doc"
echo "Open the documentation for c++ with 'open ../../../Networkit-Doc/doxyhtml/index.html'"
echo "Open the documentation for python with 'open ../../../Networkit-Doc/python/html/index.html'"
