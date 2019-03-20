/*
 * BiconnectedComponentsGTest.cpp
 *
 * Created on: March 2018
 * 		 Author: Eugenio Angriman
 */

#include <gtest/gtest.h>

#include "../../../include/networkit/auxiliary/Log.hpp"
#include "../../../include/networkit/generators/ErdosRenyiGenerator.hpp"
#include "../../../include/networkit/components/BiconnectedComponents.hpp"
#include "../../../include/networkit/components/ConnectedComponents.hpp"

namespace NetworKit {

class BiconnectedComponentsGTest : public testing::Test {};

TEST_F(BiconnectedComponentsGTest, testBiconnectedComponentsTiny) {
	Graph G(9, false, false);
	G.addEdge(0, 1);
	G.addEdge(1, 2);
	G.addEdge(1, 3);
	G.addEdge(1, 4);
	G.addEdge(0, 5);
	G.addEdge(0, 6);
	G.addEdge(4, 5);
	G.addEdge(2, 3);
	G.addEdge(6, 8);
	G.addEdge(6, 7);
	G.addEdge(7, 8);
	BiconnectedComponents bc(G);
	bc.run();

	EXPECT_EQ(bc.numberOfComponents(), 4);
}

TEST_F(BiconnectedComponentsGTest, testBiconnectedComponents) {
	Graph G = ErdosRenyiGenerator(200, 0.2, false).generate();

	BiconnectedComponents bc(G);
	bc.run();

	auto components = bc.getComponents();

	for (auto component : components) {
		ConnectedComponents cc(G);
		cc.run();
		count nComps = cc.numberOfComponents();
		for (node v : component) {

			std::vector<node> neighbors(G.neighborRange(v).begin(),
										G.neighborRange(v).end());
			// Simulating node removal
			for (node w : neighbors) {
				G.removeEdge(v, w);
			}
			ConnectedComponents cc1(G);
			cc1.run();
			count nComps1 = cc1.numberOfComponents();
			// Isolated node counts as a new component
			EXPECT_EQ(nComps, nComps1 - 1);

			// Simulating node re-insertion
			for (node w : neighbors) {
				G.addEdge(v, w);
			}
		}
	}
}

} // namespace NetworKit
