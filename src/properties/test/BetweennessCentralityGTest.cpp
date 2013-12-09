/*
 * BetweennessCentralityGTest.cpp
 *
 *  Created on: 09.12.2013
 *      Author: Lukas Barth, David Weiß
 */
#ifndef NOGTEST

#include <iostream>

#include "BetweennessCentralityGTest.h"
#include "../../graph/GraphGenerator.h"
#include "../../io/METISGraphReader.h"

#define DELTA 0.0001

namespace GrauBart {

using namespace NetworKit;

BetweennessCentralityGTest::BetweennessCentralityGTest() {
	// TODO Auto-generated constructor stub
}

BetweennessCentralityGTest::~BetweennessCentralityGTest() {
	// TODO Auto-generated destructor stub
}

TEST_F(BetweennessCentralityGTest, testStaticGraphs) {
	// construct graph
	Graph clique = GraphGenerator().makeCompleteGraph(5);

	BetweennessCentrality completeBC(clique);
	completeBC.run();

	std::hash_map<node, double> completeBCs = completeBC.getCentrality();

	double centrality = completeBCs.begin()->second;
	clique.forNodes([&](node v) {
		EXPECT_TRUE(completeBCs[v] - DELTA < centrality);
		EXPECT_TRUE(completeBCs[v] + DELTA > centrality);
	});

}

TEST_F(BetweennessCentralityGTest, testOnGraphFiles) {
  METISGraphReader input;
  
  Graph G;
  std::ofstream output;
  
  std::list<std::string> graphList = {"celegans_metabolic", "polblogs", "hep-th"};
  
  std::for_each(graphList.begin(), graphList.end(), [&](std::string filename){
    G = input.read(std::string("input/") + filename + ".graph");

    BetweennessCentrality BC (G);
    BC.run();

    std::hash_map<node, double> centralities = BC.getCentrality();

    output.open(std::string("output/") + filename + ".sol");

    double max = 0.0;
    int maxAt = 0;
    int i = 0;
    for(auto it = centralities.begin(); it != centralities.end(); ++it) {
    	// sorted by vertex
    	std::cout << it->first << ": " << it->second << std::endl;
    	output << it->second << std::endl;

    	if (it->second > max) {
    		max = it->second;
    		maxAt = i;
    	}

    	i++;
    }

    std::cout << "Maximum centrality: " << max << " at vertex " << maxAt << std::endl;

    output.close();
  });
}


} /* namespace GrauBart */

#endif /*NOGTEST */

