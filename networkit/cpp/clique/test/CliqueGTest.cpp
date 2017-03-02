/*
 * CliqueGTest.cpp
 *
 *  Created on: 08.12.2014
 *      Author: henningm
 */

#include "CliqueGTest.h"
#include "../MaxClique.h"
#include "../../io/METISGraphReader.h"
#include "../../io/SNAPGraphReader.h"
#include "../../auxiliary/Log.h"
#include "../../io/EdgeListReader.h"


namespace NetworKit {


TEST_F(CliqueGTest, testMaxCliqueOnSmallerGraphs) {
	EdgeListReader r(' ',1,"%");
	Graph gJohnson = r.read("input/johnson8-4-4.edgelist");
	Graph gHamming = r.read("input/hamming6-4.edgelist");

	MaxClique mcJohnson(gJohnson);
	MaxClique mcHamming(gHamming);

	mcJohnson.run();
	count maxCliqueSizeJohnson = mcJohnson.getMaxCliqueSize();
	mcHamming.run();
	count maxCliqueSizeHamming = mcHamming.getMaxCliqueSize();

	EXPECT_EQ(14u,maxCliqueSizeJohnson) << "maximum clique size on graph johnson8-4-4 is not correct";
	EXPECT_EQ(4u,maxCliqueSizeHamming) << "maximum clique size on graph hamming6-4 is not correct";

	std::unordered_set<node> cliqueJohnson = mcJohnson.getMaxClique();
	std::unordered_set<node> cliqueHamming = mcHamming.getMaxClique();
	EXPECT_EQ(14u, cliqueJohnson.size());
	EXPECT_EQ(4u, cliqueHamming.size());
}


} /* namespace NetworKit */
