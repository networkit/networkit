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

#include "globals.h"
#include "input/METISGraphParser.h"



int main() {


	// configure logging
	log4cxx::BasicConfigurator::configure();
	log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getInfo());


	// start
	EnsembleClustering::METISGraphParser *parser = new EnsembleClustering::METISGraphParser;

	// available graphs
	//	- kron_g500-simple-logn16.graph

	parser->parse("/Users/cls/workspace/Data/DIMACS/kron_g500-simple-logn16.graph");

	return 0;
}
