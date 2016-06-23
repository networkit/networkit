#!/bin/bash
#
# Use this script to create the website for NetworKit

rm -rf ../Documentation

make latexpdf

# Call doxygen to produce the C++ documentation.
doxygen Doxyfile_Latex
cd ../Documentation/c++
pdflatex -synctex=1 -interaction=nonstopmode refman.tex
pdflatex -synctex=1 -interaction=nonstopmode refman.tex

cd ../
mv python/NetworKit-Python.pdf NetworKit-Python.pdf
mv c++/refman.pdf NetworKit-C++.pdf

rm -rf python
rm -rf c++
rm -rf __pycache__


echo ""
echo "Finished building Documentation to ../Documentation"
