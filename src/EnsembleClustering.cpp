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

#include "aux/log.h"
#include "graph/Graph.h"
#include "input/METISParser.h"
#include "input/METIStoSTINGER.h"

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


void testMETIStoSTINGER() {

	std::string graphPath = "/Users/cls/workspace/Data/DIMACS/kron_g500-simple-logn16.graph";
	LOG4CXX_INFO(log4cxx::Logger::getRootLogger(), "trying to read from graph file " << graphPath);


	EnsembleClustering::Graph* G;
	EnsembleClustering::METIStoSTINGER* m2s = new EnsembleClustering::METIStoSTINGER();
	G = m2s->read(graphPath);

	LOG4CXX_INFO(log4cxx::Logger::getRootLogger(), "read graph " << G << " from file " << graphPath);

}


int main() {

	std::cout << "running EnsembleClustering" << std::endl;



	// configure logging
	log4cxx::BasicConfigurator::configure();
	log4cxx::Logger::getRootLogger()->setLevel(log4cxx::Level::getInfo());

	debug("test debug macro");

	testMETIStoSTINGER();

	return 0;
}
