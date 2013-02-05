/*
 * STINGERGTest.h
 *
 *  Created on: 10.01.2013
 *      Author: cls
 */

#ifndef STINGERGTEST_H_
#define STINGERGTEST_H_

#include <gtest/gtest.h>

extern "C" {
#include "stinger.h"
}
#include "../Graph.h"


namespace EnsembleClustering {

/**
 * Diagnose STINGER problems.
 */
class STINGERGTest: public testing::Test {
};



} /* namespace EnsembleClustering */
#endif /* STINGERGTEST_H_ */
