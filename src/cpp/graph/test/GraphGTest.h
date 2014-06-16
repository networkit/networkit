/*
 * GraphGTest.h
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#ifndef GRAPHGTEST_H_
#define GRAPHGTEST_H_

#include <tuple>
#include <gtest/gtest.h>

#include "../Graph.h"

namespace NetworKit {

class GraphGTest: public testing::TestWithParam< std::tuple<bool, bool> > {
public:
	virtual void SetUp();

protected:
	Graph Ghouse;
	std::vector< std::pair<node, node> > houseEdgesOut;
	std::vector< std::vector<edgeweight> > Ahouse;
	count n_house;
	count m_house;

	bool isWeightedParameterized() const;
	bool isDirectedParameterized() const;
	Graph createParameterizedGraph(count n = 0) const;
};

} /* namespace NetworKit */

#endif /* GRAPHGTEST_H_ */

#endif /*NOGTEST */
