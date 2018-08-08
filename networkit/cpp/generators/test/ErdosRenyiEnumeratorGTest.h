/*
 * ErdosReyniEnumeratorGTest.h
 *
 *  Created on: 09. Aug. 2018
 *      Author: Manuel Penschuck
 */

#ifndef ERDOSRENYIENUMERATORGTEST_H_
#define ERDOSRENYIENUMERATORGTEST_H_

#include <gtest/gtest.h>
#include <tuple>
#include "../../Globals.h"

namespace NetworKit {

class ErdosRenyiEnumeratorGTest : public testing::TestWithParam<std::tuple<bool, node, double> > {};

} /* namespace NetworKit */
#endif /* ERDOSRENYIENUMERATORGTEST_H_ */

