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

# Call doxygen to produce the C++ documentation.
doxygen Doxyfile

# Clean up the modifications in the rst files modified during the preparation.
python3 -c'import sphinxPreparation; sphinxPreparation.cleanUp()'
rm -rf __pycache__

# remove doctrees (not needed for html) and move html one up
rm -rf ../Website/doctrees
mv ../Website/html/* ../Website/
rm -rd ../Website/html

# copy userguides from Notebooks/
cp ../Notebooks/NetworKit_UserGuide.ipynb ../uploads/docs/NetworKit_UserGuide.ipynb
cp ../Notebooks/GephiStreaming_UserGuide.ipynb ../uploads/docs/GephiStreaming_UserGuide.ipynb
cp ../Notebooks/SpectralCentrality.ipynb ../uploads/docs/SpectralCentrality.ipynb
cp ../Notebooks/SpectralCentralityWithPandas.ipynb ../uploads/docs/SpectralCentralityWithPandas.ipynb

# copy uploads folder
mkdir ../Website/uploads/
cp -r ../uploads ../Website/

# create documentation
./make_doc.sh

# zip documentation and repository and move it to uploads/
zip -r ../Website/uploads/Documentation.zip ../Documentation/
zip -r ../Website/uploads/NetworKit.zip ../../../
hg archive -t zip ../Website/uploads/NetworKit.zip






echo ""
echo "Finished building Website to ../Website"
