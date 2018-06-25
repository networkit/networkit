/*
 * CurveballGlobalTradeGeneratorGTest.cpp
 *
 *  Created on: Jul 11, 2017
 *      Author: Hung Tran <htran@ae.cs.uni-frankfurt.de>
 */

#include "CurveballGlobalTradeGeneratorGTest.h"
#include "../CurveballGlobalTradeGenerator.h"

namespace NetworKit {

TEST_F(CurveballGlobalTradeGeneratorGTest, testTradeNumber) {
    CurveballGlobalTradeGenerator genEven(1, 10);
    auto tradesEven = genEven.generate();
    ASSERT_EQ(tradesEven.size(), 5u);

    CurveballGlobalTradeGenerator genOdd(1, 15);
    auto tradesOdd = genOdd.generate();
    ASSERT_EQ(tradesOdd.size(), 7u);

    CurveballGlobalTradeGenerator genOddTimesFive(5, 15);
    auto tradesOddTimesFive = genOddTimesFive.generate();
    ASSERT_EQ(tradesOddTimesFive.size(), 35u);
}

}
