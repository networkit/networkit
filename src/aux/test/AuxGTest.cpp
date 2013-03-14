/*
 * AuxGTest.cpp
 *
 *  Created on: 10.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "AuxGTest.h"

TEST_F(AuxGTest, produceRandomIntegers) {
	int64_t l = 0; 	// lower bound
	int64_t u = 100;	// upper bound
	Aux::RandomInteger randInt(l, u);

	for (int i = 0; i < 100; ++i) {
		TRACE(randInt.generate());
	}
}


TEST_F(AuxGTest, testRandomInteger) {
	int64_t l = 0; 	// lower bound
	int64_t u = 10;	// upper bound
	Aux::RandomInteger randInt(l, u);
	std::vector<int64_t> rVector;
	int n = 1000;
	for (int i = 0; i < n; ++i) {
		int64_t r = randInt.generate();
		assert(l <= r);
		assert(r <= u);
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
	Aux::RandomProbability randPr;
	std::vector<double> rVector;
	int n = 1000;
	for (int i = 0; i < n; ++i) {
		double r = randPr.generate();
		assert(0.0 <= r);
		assert(r <= 1.0);
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




TEST_F(AuxGTest, testTimer) {
	int64_t sleepTime = 1000; // sleep time in ms
	int64_t tolerance = 10;

	Aux::Timer timer;
	DEBUG("sleeping for " << sleepTime << " ms");
	timer.start();
	std::this_thread::sleep_for(std::chrono::milliseconds(sleepTime));
	timer.stop();
	std::chrono::duration<int64_t, std::milli> elapsed = timer.elapsed();
	int64_t ec = elapsed.count();

	EXPECT_LE(sleepTime - tolerance, ec) << "elapsed time should be roughly equal to sleep time";
	EXPECT_GE(sleepTime + tolerance, ec) << "elapsed time should be roughly to sleep time";
}

