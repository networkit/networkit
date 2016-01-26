/*
 * SpanningEdgeCentralityGTest.h
 *
 *  Created on: Jan 17, 2016
 *      Author: Michael
 */

#ifndef NETWORKIT_CPP_CENTRALITY_TEST_SPANNINGEDGECENTRALITYGTEST_H_
#define NETWORKIT_CPP_CENTRALITY_TEST_SPANNINGEDGECENTRALITYGTEST_H_

#include "gtest/gtest.h"
#include "../Spanning.h"

#include <vector>
#include <string>

namespace NetworKit {

using namespace std;

class SpanningEdgeCentralityGTest : public testing::Test {
public:
	SpanningEdgeCentralityGTest() = default;
	virtual ~SpanningEdgeCentralityGTest() = default;

protected:
	vector<string> instances = {/*"instances/SEC/as-skitter_comp.graph", "instances/SEC/CA-GrQc_comp.graph", "instances/SEC/CA-HepTh_comp.graph", "instances/SEC/cit-Patents_comp.graph",
								"instances/SEC/com-amazon.ungraph.graph", "instances/SEC/com-dblp.ungraph.graph", "instances/SEC/com-youtube.ungraph.graph", "instances/SEC/hollywood2009.graph",
								"instances/SEC/LiveJournal.graph", "instances/SEC/oregon1_010526.graph", "instances/SEC/orkut.graph", "instances/SEC/p2p-Gnutella08_comp.graph",
								"instances/SEC/p2p-Gnutella31_comp.graph", "instances/SEC/roadNet-TX_comp.graph", "instances/SEC/Slashdot0902.graph", "instances/SEC/soc-Epinions1_comp.graph",
								"instances/SEC/Wiki-Vote_comp.graph"*/ "/Users/Michael/Downloads/SNAP/CA-GrQc_comp.graph"};

};

} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_CENTRALITY_TEST_SPANNINGEDGECENTRALITYGTEST_H_ */
