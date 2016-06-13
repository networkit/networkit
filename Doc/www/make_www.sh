#!/bin/bash
#
# Use this script to create the website for NetworKit

rm -rf ../Website

# Do some preparations in rst files prior to the build.
python3 -c'import sphinxPreparation; sphinxPreparation.prepareBuild()'

# Call sphinx to generate the website. This includes a build of NetworKit as python module to get the latest docstrings.
make html-networkit

# sphinx is run twice to make sure references are used correctly
make html

mkdir ../Website/latex

# Call doxygen to produce the C++ documentation.
doxygen Doxyfile

# Clean up the modifications in the rst files modified during the preparation.
python3 -c'import sphinxPreparation; sphinxPreparation.cleanUp()'
rm -rf __pycache__

echo ""
echo "Finished building Website to ../Website"
