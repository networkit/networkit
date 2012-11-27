//============================================================================
// Name        : EnsembleClustering.cpp
// Author      : Christian Staudt
// Version     :
// Copyright   : © 2012, Christian Staudt
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <utility>

#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"

#include "input/METISParser.h"

void testMETISParser() {

	EnsembleClustering::METISParser* parser = new EnsembleClustering::METISParser();

	parser->open("/Users/cls/workspace/Data/DIMACS/kron_g500-simple-logn16.graph");
	std::pair<int, int> header = parser->getHeader();

	int lc = 0;
	while (parser->hasNext()) {
		parser->getNext();
		lc++;
	}

	LOG4CXX_INFO(log4cxx::Logger::getRootLogger(), "parsed " << lc << " lines");

	parser->close();
}


int main() {

	std::cout << "running EnsembleClustering" << std::endl;

	// configure logging
	log4cxx::BasicConfigurator::configure();
	log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getInfo());


	testMETISParser();

	// start
	//EnsembleClustering::METISGraphParser *parser = new EnsembleClustering::METISGraphParser;

	// available graphs
	//	- kron_g500-simple-logn16.graph

	//parser->parse("/Users/cls/workspace/Data/DIMACS/kron_g500-simple-logn16.graph");

	return 0;
}
