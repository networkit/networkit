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

//TEST_F(CliqueGTest, testMaxClique1) {
//	SNAPGraphReader sreader;
//	Graph G1 = sreader.read("input/email-Enron.txt");
//	MaxClique mc1(G1);
//	count maxCliqueSize1 = mc1.run();
//	EXPECT_EQ(20, maxCliqueSize1);

//	METISGraphReader mreader;
//	Graph G2 = mreader.read("input/kkt_power.graph");
//	MaxClique mc2(G2);
//	count maxCliqueSize2 = mc2.run();
//	EXPECT_EQ(11, maxCliqueSize2);
//}

TEST_F(CliqueGTest, testMaxCliqueOnSmallerGraphs) {
	EdgeListReader r(' ',1,"%");
	Graph gKeller = r.read("input/keller4.edgelist");
	Graph gJohnson = r.read("input/johnson8-4-4.edgelist");
	Graph gHamming = r.read("input/hamming6-4.edgelist");

	MaxClique mcKeller(gKeller);
	MaxClique mcJohnson(gJohnson);
	MaxClique mcHamming(gHamming);

	//count maxCliqueSizeKeller = mcKeller.run();
	mcJohnson.run();
	count maxCliqueSizeJohnson = mcJohnson.getMaxCliqueSize();
	mcHamming.run();
	count maxCliqueSizeHamming = mcHamming.getMaxCliqueSize();

	//EXPECT_EQ(11,maxCliqueSizeKeller) << "maximum clique size on graph keller4 is not correct";
	EXPECT_EQ(14u,maxCliqueSizeJohnson) << "maximum clique size on graph johnson8-4-4 is not correct";
	EXPECT_EQ(4u,maxCliqueSizeHamming) << "maximum clique size on graph hamming6-4 is not correct";

	std::unordered_set<node> cliqueJohnson = mcJohnson.getMaxClique();
	std::unordered_set<node> cliqueHamming = mcHamming.getMaxClique();
	EXPECT_EQ(14u, cliqueJohnson.size());
	EXPECT_EQ(4u, cliqueHamming.size());
}


} /* namespace NetworKit */
