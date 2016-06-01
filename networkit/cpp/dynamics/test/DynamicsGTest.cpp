/*
 * DynamicsGTest.cpp
 *
 *  Created on: 24.12.2013
 *      Author: cls
 */

#include "DynamicsGTest.h"

#include "../DGSStreamParser.h"
#include "../../auxiliary/Log.h"
#include "../GraphEvent.h"
#include "../GraphUpdater.h"

namespace NetworKit {

TEST_F(DynamicsGTest, testDGSStreamParser) {
	DGSStreamParser parser("input/example2.dgs");
	auto stream = parser.getStream();
	for (auto event : stream) {
		INFO(event.toString(), " ");
	}
	INFO("\n");
}


TEST_F(DynamicsGTest, tryDGSStreamParserOnRealGraph) {
	std::string path;
	std::cout << "enter .dgs file path: ";
	std::cin >> path;
	DGSStreamParser parser(path);
	auto stream = parser.getStream();
}

TEST_F(DynamicsGTest, testGraphEventIncrement) {
	Graph G(2, true, false); //undirected
	Graph H(2, true, true); //directed
	G.addEdge(0, 1, 3.14);
	H.addEdge(0, 1, 3.14);
	GraphEvent event(GraphEvent::EDGE_WEIGHT_INCREMENT, 0, 1, 2.1);
	std::vector<GraphEvent> eventstream(1);
	eventstream.push_back(event);
	GraphUpdater Gupdater(G);
	GraphUpdater Hupdater(H);
	Gupdater.update(eventstream);
	Hupdater.update(eventstream);
	EXPECT_EQ(G.weight(0,1), 5.24);
	EXPECT_EQ(H.weight(0,1), 5.24);



}

} /* namespace NetworKit */
