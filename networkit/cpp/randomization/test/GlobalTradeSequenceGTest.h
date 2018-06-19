/*
 * GlobalTradeSequenceGTest.h
 *
 *  Created on: 20.06.2018
 *      Author: Manuel Penschuck (networkit@manuel.jetzt)
 */

#ifndef RANDOMIZATION_TEST_GLOBAL_TRADE_SEQUENCE_GTEST_H
#define RANDOMIZATION_TEST_GLOBAL_TRADE_SEQUENCE_GTEST_H

#include <gtest/gtest.h>

namespace NetworKit {

template <typename T>
class GlobalTradeSequenceGTest : public ::testing::Test  {
public:
    GlobalTradeSequenceGTest() = default;
    virtual ~GlobalTradeSequenceGTest() = default;
};

}

#endif // RANDOMIZATION_TEST_GLOBAL_TRADE_SEQUENCE_GTEST_H
