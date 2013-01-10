/*
 * AuxGTest.cpp
 *
 *  Created on: 10.01.2013
 *      Author: cls
 */

#include "AuxGTest.h"

namespace EnsembleClustering {

TEST_F(AuxGTest, testRandomInteger) {
	int64_t l = 0; 	// lower bound
	int64_t u = 10;	// upper bound
	RandomInteger randInt(l, u);
	std::vector<int64_t> rVector;
	int n = 1000;
	for (int i = 0; i < n; ++i) {
		int64_t r = randInt.generate();
		assert(l <= r <= u);
		rVector.push_back(r);
	}

	int64_t minR = *(min_element(rVector.begin(), rVector.end()));
	int64_t maxR = *(max_element(rVector.begin(), rVector.end()));

	EXPECT_EQ(minR, l);
	EXPECT_EQ(maxR, u);

	double sum = 0.0;
	for (int64_t r : rVector) {
		sum += r;
	}
	double avg = sum / n;


	DEBUG("avg rand integer: " << avg);
	EXPECT_LE(avg, 6.0);
	EXPECT_GE(avg, 4.0);
}


TEST_F(AuxGTest, testRandomProbability) {
	RandomProbability randPr;
	std::vector<double> rVector;
	int n = 1000;
	for (int i = 0; i < n; ++i) {
		double r = randPr.generate();
		assert(0.0 <= r <= 1.0);
		rVector.push_back(r);
	}

	double sum = 0.0;
	for (double r : rVector) {
		sum += r;
	}
	double avg = sum / n;


	DEBUG("avg rand probability: " << avg);
	EXPECT_LE(avg, 0.6);
	EXPECT_GE(avg, 0.4);
}

} /* namespace EnsembleClustering */
