/*
 * IGraphGTest.h
 *
 *  Created on: 10.05.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 */

#ifndef NOGTEST

#ifndef IGRAPHGTEST_H_
#define IGRAPHGTEST_H_

#include <gtest/gtest.h>

#include "../IGraph.h"
#include "../../auxiliary/Random.h"

namespace NetworKit {

template <typename T>
class IGraphGTest : public testing::Test {

};

} /* namespace NetworKit */

#endif /* IGRAPHGTEST_H_ */

#endif /*NOGTEST */
