//============================================================================
// Name        : EnsembleClustering.cpp
// Author      : Christian Staudt
// Version     :
// Copyright   : © 2012, Christian Staudt
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>

#include "METISGraphParser.h"

int main() {


	EnsembleClustering::METISGraphParser *parser = new EnsembleClustering::METISGraphParser;
	parser->parse("/Users/cls/workspace/Data/DIMACS/example.graph");

	return 0;
}
