/*
 * SortedListGTest.cpp
 *
 *  Created on: 17.07.2025
 *      Author: Andreas Scharf (andreas.b.scharf@gmail.com)
 */

#include <gtest/gtest.h>
#include <networkit/auxiliary/SortedList.hpp>
namespace Networkit {

class SortedListGTest : public testing::Test {};

TEST_F(SortedListGTest, testBasicInsertionSortedOrder) {
    Aux::SortedList list(3, 10);

    list.insert(1, 5.0);
    list.insert(2, 7.0);
    list.insert(3, 6.0);

    EXPECT_EQ(list.getSize(), 3);
    EXPECT_EQ(list.getElement(0), 2);
    EXPECT_EQ(list.getValue(0), 7.0);
    EXPECT_EQ(list.getElement(1), 3);
    EXPECT_EQ(list.getValue(1), 6.0);
    EXPECT_EQ(list.getElement(2), 1);
    EXPECT_EQ(list.getValue(2), 5.0);
}

TEST_F(SortedListGTest, testInsertExistingElementWithHigherValue) {
    Aux::SortedList list(3, 10);

    list.insert(4, 3.0);
    list.insert(5, 4.0);
    list.insert(6, 5.0);

    list.insert(4, 6.0);

    EXPECT_EQ(list.getSize(), 3);

    EXPECT_EQ(list.getElement(0), 4);
    EXPECT_EQ(list.getValue(0), 6.0);
    EXPECT_EQ(list.getElement(1), 6);
    EXPECT_EQ(list.getValue(1), 5.0);
    EXPECT_EQ(list.getElement(2), 5);
    EXPECT_EQ(list.getValue(2), 4.0);
}

TEST_F(SortedListGTest, testExceedingCapacityKeepsTopK) {
    Aux::SortedList list(3, 10);

    list.insert(1, 1.0);
    list.insert(2, 2.0);
    list.insert(3, 3.0);
    EXPECT_EQ(list.getElement(2), 1);
    EXPECT_EQ(list.getSize(), 3);

    list.insert(4, 4.0); // Should evict (1, 1.0)
    EXPECT_EQ(list.getSize(), 3);
    EXPECT_NE(list.getElement(2), 1);
    EXPECT_EQ(list.getElement(2), 2);
}

TEST_F(SortedListGTest, testDuplicateInsertSameElement) {
    Aux::SortedList list(3, 10);

    list.insert(1, 5.0);
    list.insert(2, 6.0);
    list.insert(1, 7.0);

    EXPECT_EQ(list.getSize(), 3);
    EXPECT_EQ(list.getElement(0), 1);
    EXPECT_EQ(list.getValue(0), 7.0);
    EXPECT_EQ(list.getElement(1), 2);
    EXPECT_EQ(list.getValue(1), 6.0);
    EXPECT_EQ(list.getElement(2), 1);
    EXPECT_EQ(list.getValue(2), 5.0);
}

TEST_F(SortedListGTest, testClear) {
    Aux::SortedList list(3, 10);
    list.insert(1, 10.0);
    list.insert(2, 9.0);
    list.clear();

    EXPECT_EQ(list.getSize(), 0);
    list.insert(3, 8.0);
    EXPECT_EQ(list.getSize(), 1);
    EXPECT_EQ(list.getElement(0), 3);
    EXPECT_EQ(list.getValue(0), 8.0);
}

TEST_F(SortedListGTest, testInsertingOutOfRangeKeyThrows) {
    EXPECT_THROW(Aux::SortedList(5, 3), std::runtime_error);
}

} // namespace Networkit
