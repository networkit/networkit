/*
 * METIStoSTINGER.cpp
 *
 *  Created on: 27.11.2012
 *      Author: cls
 */

#include "METIStoSTINGER.h"

#include <string>
#include <vector>

#include "STINGERFromAdjacencies.h"
#include "METISParser.h"

namespace EnsembleClustering {

METIStoSTINGER::METIStoSTINGER() {
	// TODO Auto-generated constructor stub

}

METIStoSTINGER::~METIStoSTINGER() {
	// TODO Auto-generated destructor stub
}

graph* METIStoSTINGER::read(std::string graphPath) {

	METISParser* parser = new METISParser();
	STINGERFromAdjacencies* builder = new STINGERFromAdjacencies();

	builder->createGraph();

	parser->open(graphPath);
	std::pair<int, int> header = parser->getHeader();



	while (parser->hasNext()) {
		builder->addAdjacencies(parser->getNext());
	}

	parser->close();

	return builder->getGraph();
}

} /* namespace EnsembleClustering */
