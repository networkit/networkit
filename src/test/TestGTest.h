/*
 * TestGTest.h
 *
 *  Created on: 11.12.2012
 *      Author: cls
 */

#ifndef TESTGTEST_H_
#define TESTGTEST_H_

#include "gtest/gtest.h"

class GTestTest : public ::testing::Test {
 protected:

  virtual void SetUp() {
  }

  // virtual void TearDown() {}
};


TEST_F(GTestTest, myFirstTest) {
  EXPECT_EQ(0, 0);
}




#endif /* TESTGTEST_H_ */
