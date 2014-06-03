/*
 * IDirectedGraphGTest.h
 *
 *  Created on: 10.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#ifndef IDIRECTEDGRAPHGTEST_H_
#define IDIRECTEDGRAPHGTEST_H_

#include <vector>
#include <gtest/gtest.h>

#include "../BasicGraph.h"

namespace NetworKit {

template <typename T>
class IDirectedGraphGTest : public testing::Test {
public:
	virtual void SetUp();

protected:
	T Ghouse;
	std::vector< std::pair<node, node> > houseEdgesOut;
	count m_house;
};

} /* namespace NetworKit */

#endif /* IDIRECTEDGRAPHGTEST_H_ */

#endif /*NOGTEST */
