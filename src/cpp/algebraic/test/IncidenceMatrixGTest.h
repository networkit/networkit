/*
 * IncidenceMatrixGTest.h
 *
 *  Created on: 01.04.2014
 *      Author: Michael
 */

#ifndef INCIDENCEMATRIXGTEST_H_
#define INCIDENCEMATRIXGTEST_H_

#include "gtest/gtest.h"
#include "../IncidenceMatrix.h"
#include "../Vector.h"
#include "../../graph/Graph.h"

namespace NetworKit {

class IncidenceMatrixGTest : public testing::Test {
protected:
	virtual void SetUp() {
		graph = NetworKit::Graph(5);
		graph.addEdge(0,1);
		graph.addEdge(0,2);
		graph.addEdge(0,3);
		graph.addEdge(2,3);
		graph.addEdge(4,1);
		graph.addEdge(4,4);
	}

	NetworKit::Graph graph;

public:
	IncidenceMatrixGTest();
	virtual ~IncidenceMatrixGTest();
};


} /* namespace NetworKit */

#endif /* INCIDENCEMATRIXGTEST_H_ */
