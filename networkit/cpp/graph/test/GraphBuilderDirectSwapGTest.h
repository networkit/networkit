/*
 * GraphBuilderDirectSwapGTest.h
 *
 *  Created on: 14.08.2014
 *      Author: Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#ifndef GRAPH_BUILDER_DIRECT_SWAP_GTEST_H_
#define GRAPH_BUILDER_DIRECT_SWAP_GTEST_H_

#include <tuple>
#include <gtest/gtest.h>

#include "../Graph.h"
#include "../GraphBuilder.h"

namespace NetworKit {

class GraphBuilderDirectSwapGTest: public testing::TestWithParam< std::tuple<bool, bool> > {
public:
	virtual void SetUp();

protected:
	GraphBuilder bHouse;
	std::vector< std::pair<node, node> > houseEdgesOut;
	std::vector< std::vector<edgeweight> > Ahouse;
	count n_house;
	count m_house;

	bool isGraph() const { return !isWeighted() && !isDirected(); }
	bool isWeightedGraph() const { return isWeighted() && !isDirected(); }
	bool isDirectedGraph() const { return !isWeighted() && isDirected(); }
	bool isWeightedDirectedGraph() const { return isWeighted() && isDirected(); }

	bool isWeighted() const;
	bool isDirected() const;

	GraphBuilder createGraphBuilder(count n = 0) const;
	Graph toGraph(GraphBuilder& b) const;
};

} /* namespace NetworKit */

#endif /* GRAPH_BUILDER_DIRECT_SWAP_GTEST_H_ */

#endif /* NOGTEST */
