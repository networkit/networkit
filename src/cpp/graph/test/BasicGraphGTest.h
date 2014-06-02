/*
 * BasicGraph.h
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#ifndef BASICGRAPHGTEST_H_
#define BASICGRAPHGTEST_H_

#include <gtest/gtest.h>

#include "../BasicGraph.h"

namespace NetworKit {

template <typename T>
class BasicGraphGTest: public testing::Test {
	
	virtual void SetUp();

protected:
	T Ghouse;
	std::vector< std::pair<node, node> > houseEdgesOut;
	count n_house;
	count m_house;
};

} /* namespace NetworKit */

#endif /* BASICGRAPHGTEST_H_ */

#endif /*NOGTEST */
