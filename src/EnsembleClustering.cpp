//============================================================================
// Name        : EnsembleClustering.cpp
// Author      : Christian Staudt
// Version     :
// Copyright   : © 2012, Christian Staudt
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <fstream>

#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"


#include "input/METISGraphParser.h"


int main() {


	// configure logging
	log4cxx::BasicConfigurator::configure();

	// start
	EnsembleClustering::METISGraphParser *parser = new EnsembleClustering::METISGraphParser;
	parser->parse("/Users/cls/workspace/Data/DIMACS/example.graph");

	return 0;
}
