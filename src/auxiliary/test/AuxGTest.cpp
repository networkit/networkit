#ifndef NOGTEST

/*
 * AuxGTest.cpp
 *
 *  Created on: 10.01.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include "AuxGTest.h"

#include <iostream>
#include <algorithm>
#include <chrono>
#include <thread>

#include "../Log.h"
#include "../Random.h"
#include "../RandomInteger.h"
#include "../RandomProbability.h"
#include "../Timer.h"
#include "../MissingMath.h"
#include "../Debug.h"
#include "../PriorityQueue.h"

TEST_F(AuxGTest, produceRandomIntegers) {
	int64_t l = 0; 	// lower bound
	int64_t u = 100;	// upper bound

	for (int i = 0; i < 100; ++i) {
		TRACE(Aux::RandomInteger::generate(l, u));
	}
}

TEST_F(AuxGTest, produceRandomIntegersNew) {
	int64_t l = 0; 	// lower bound
	int64_t u = 100;	// upper bound
	Aux::Random rand;

	for (int i = 0; i < 100; ++i) {
		TRACE(rand.integer(l, u));
	}
}


TEST_F(AuxGTest, testRandomInteger) {
	int64_t l = 0; 	// lower bound
	int64_t u = 10;	// upper bound
	std::vector<int64_t> rVector;
	int n = 1000;
	for (int i = 0; i < n; ++i) {
		int64_t r = Aux::RandomInteger::generate(l, u);
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


TEST_F(AuxGTest, testRandomIntegerNew) {
	int64_t l = 0; 	// lower bound
	int64_t u = 10;	// upper bound
	Aux::Random rand;
	std::vector<int64_t> rVector;
	int n = 1000;
	for (int i = 0; i < n; ++i) {
		int64_t r = rand.integer(l, u);
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


TEST_F(AuxGTest, testRandomIntegerFaster) {
	int64_t l = 0; 	// lower bound
	int64_t u = 10;	// upper bound
	std::vector<int64_t> rVector;
	int n = 1000;
	for (int i = 0; i < n; ++i) {
		int64_t r = Aux::RandomInteger::generateFaster(l, u);
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
	std::vector<double> rVector;
	int n = 1000;
	for (int i = 0; i < n; ++i) {
		double r = Aux::RandomProbability::generate();
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


TEST_F(AuxGTest, testRandomProbabilityNew) {
	Aux::Random rand;
	std::vector<double> rVector;
	int n = 1000;
	for (int i = 0; i < n; ++i) {
		double r = rand.probability();
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
	int64_t tolerance = 20;

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


TEST_F(AuxGTest, testBinomial) {
	EXPECT_EQ(1, Aux::MissingMath::binomial(7,0));
	EXPECT_EQ(7, Aux::MissingMath::binomial(7,1));
	EXPECT_EQ(21, Aux::MissingMath::binomial(7,2));
	EXPECT_EQ(35, Aux::MissingMath::binomial(7,3));
	EXPECT_EQ(35, Aux::MissingMath::binomial(7,4));
	EXPECT_EQ(21, Aux::MissingMath::binomial(7,5));
	EXPECT_EQ(7, Aux::MissingMath::binomial(7,6));
	EXPECT_EQ(1, Aux::MissingMath::binomial(7,7));


}


TEST_F(AuxGTest, benchmarkBinomial) {
	Aux::Timer timer;
	INFO("starting calculation");
	timer.start();

	int64_t n = 500;
	for (int64_t k = 0; k < n; ++k) {
		Aux::MissingMath::binomial(n, k);
	}

	timer.stop();

	INFO("calculation finished after " << timer.elapsedTag());

}


TEST_F(AuxGTest, testVectorDebug) {
	std::vector<int> vec(10, 42);
	std::cout << Aux::vectorToString(vec) << std::endl;
}


TEST_F(AuxGTest, testPriorityQueue) {
	typedef std::pair<double, uint64_t> ElemType;

	// fill vector
	std::vector<ElemType> vec;
	vec.push_back(std::make_pair(0.5, 0));
	vec.push_back(std::make_pair(3.5, 1));
	vec.push_back(std::make_pair(4.5, 2));
	vec.push_back(std::make_pair(2.5, 3));
	vec.push_back(std::make_pair(0.75, 4));
	vec.push_back(std::make_pair(1.5, 5));
	vec.push_back(std::make_pair(8.5, 6));
	vec.push_back(std::make_pair(3.25, 7));
	vec.push_back(std::make_pair(4.75, 8));
	vec.push_back(std::make_pair(5.0, 9));
	vec.push_back(std::make_pair(11.5, 10));
	vec.push_back(std::make_pair(0.25, 11));

	// construct pq from vector
	Aux::PriorityQueue<double, uint64_t> pq(vec);
	EXPECT_EQ(pq.size(), vec.size());

	ElemType elem = pq.extractMin();
	EXPECT_EQ(elem.first, 0.25);
	EXPECT_EQ(elem.second, 11);
	EXPECT_EQ(pq.size(), vec.size() - 1);

	elem = pq.extractMin();
	EXPECT_EQ(elem.first, 0.5);
	EXPECT_EQ(elem.second, 0);
	EXPECT_EQ(pq.size(), vec.size() - 2);

	elem = pq.extractMin();
	EXPECT_EQ(elem.first, 0.75);
	EXPECT_EQ(elem.second, 4);
	EXPECT_EQ(pq.size(), vec.size() - 3);

	elem = pq.extractMin();
	EXPECT_EQ(elem.first, 1.5);
	EXPECT_EQ(elem.second, 5);
	EXPECT_EQ(pq.size(), vec.size() - 4);

	elem = pq.extractMin();
	EXPECT_EQ(elem.first, 2.5);
	EXPECT_EQ(elem.second, 3);
	EXPECT_EQ(pq.size(), vec.size() - 5);
}


#endif /*NOGTEST */
