/*
* SortedListGTest.cpp
 *
 *  Created on: 17.07.2025
 *      Author: Andreas Scharf (andreas.b.scharf@gmail.com)
 */

#include <gtest/gtest.h>
#include <networkit/auxiliary/SortedList.hpp>
namespace Networkit {

TEST(SortedListTest, testBasicInsertionSortedOrder) {
    Aux::SortedList sl(3, 10);

    sl.insert(1, 5.0);
    sl.insert(2, 7.0);
    sl.insert(3, 6.0);

    EXPECT_EQ(sl.getSize(), 3u);
    EXPECT_EQ(sl.getElement(0), 2u);
    EXPECT_EQ(sl.getValue(0), 7.0);
    EXPECT_EQ(sl.getElement(1), 3u);
    EXPECT_EQ(sl.getValue(1), 6.0);
    EXPECT_EQ(sl.getElement(2), 1u);
    EXPECT_EQ(sl.getValue(2), 5.0);
}

TEST(SortedListTest, UpdateExistingElementHigherValue) {
    Aux::SortedList sl(3, 10);

    sl.insert(4, 3.0);
    sl.insert(5, 4.0);
    sl.insert(6, 5.0);

    sl.insert(4, 6.0);

    EXPECT_EQ(sl.getElement(0), 4u);
    EXPECT_EQ(sl.getValue(0), 6.0);
    EXPECT_EQ(sl.getSize(), 3u);
}

TEST(SortedListTest, ExceedingCapacityOnlyKeepsTopK) {
    Aux::SortedList sl(3, 10);

    sl.insert(1, 1.0);
    sl.insert(2, 2.0);
    sl.insert(3, 3.0);
    EXPECT_EQ(sl.getElement(2), 1u);
    EXPECT_EQ(sl.getSize(), 3u);

    sl.insert(4, 4.0);  // Should evict (1, 1.0)
    EXPECT_EQ(sl.getSize(), 3u);
    EXPECT_NE(sl.getElement(2), 1u);
    EXPECT_EQ(sl.getElement(2), 2u);
}

TEST(SortedListTest, DuplicateInsertSameValueIgnoredPosition) {
    Aux::SortedList sl(3, 10);

    sl.insert(1, 5.0);
    sl.insert(2, 6.0);
    sl.insert(1, 7.0);

    EXPECT_EQ(sl.getSize(), 3u);
    EXPECT_EQ(sl.getElement(0), 1u);
    EXPECT_EQ(sl.getValue(0), 7.0);
    EXPECT_EQ(sl.getElement(1), 2u);
    EXPECT_EQ(sl.getValue(1), 6.0);
    EXPECT_EQ(sl.getElement(2), 1u);
    EXPECT_EQ(sl.getValue(2), 5.0);
}

TEST(SortedListTest, ClearEmptiesAll) {
    Aux::SortedList sl(3, 10);
    sl.insert(1, 10.0);
    sl.insert(2, 9.0);
    sl.clear();

    EXPECT_EQ(sl.getSize(), 0u);
    sl.insert(3, 8.0);
    EXPECT_EQ(sl.getSize(), 1u);
    EXPECT_EQ(sl.getElement(0), 3u);
}

TEST(SortedListTest, InsertingOutOfRangeKeyThrows) {
    EXPECT_THROW(Aux::SortedList(5, 3), std::runtime_error);
}

TEST(SortedListTest, InsertingNonMonotonicValueIgnored) {
    Aux::SortedList sl(3, 10);

    sl.insert(2, 5.0);
    sl.insert(2, 3.0);
    EXPECT_EQ(sl.getElement(0), 2u);
    EXPECT_EQ(sl.getValue(0), 5.0);
}
}