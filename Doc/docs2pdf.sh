#!/bin/bash
# convert documentation in markdown format to PDF
# get pandoc at http://code.google.com/p/pandoc/downloads/list

OPTIONS="-V geometry:margin=2.5cm"
pandoc ${OPTIONS} ../Readme.mdown -o Readme.pdf
pandoc ${OPTIONS} DevGuide.mdown -o DevGuide.pdf
