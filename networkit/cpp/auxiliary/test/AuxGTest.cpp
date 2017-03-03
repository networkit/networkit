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
#include <numeric>
#include <chrono>
#include <thread>
#include <fstream>
#include <set>

#include "../Log.h"
#include "../Random.h"
#include "../Timer.h"
#include "../MissingMath.h"
#include "../PrioQueue.h"
#include "../PrioQueueForInts.h"
#include "../BucketPQ.h"
#include "../StringTools.h"
#include "../SetIntersector.h"
#include "../Enforce.h"
#include "../NumberParsing.h"
#include "../Enforce.h"
#include "../BloomFilter.h"


TEST_F(AuxGTest, produceRandomIntegers) {
	Aux::Random::setSeed(1, false);
#if (LOG_LEVEL == LOG_LEVEL_TRACE)
	int64_t l = 0; 	// lower bound
	int64_t u = 100;	// upper bound

	for (int i = 0; i < 100; ++i) {
		TRACE(Aux::Random::integer(l, u));
	}
#else
	for (int i = 0; i < 100; ++i) {
		TRACE(Aux::Random::integer(1, 100));
	}
#endif
}

TEST_F(AuxGTest, testRandomInteger) {
	Aux::Random::setSeed(1, false);
	int64_t l = 0; 	// lower bound
	int64_t u = 10;	// upper bound
	std::vector<int64_t> rVector;
	int n = 1000;
	for (int i = 0; i < n; ++i) {
		int64_t r = Aux::Random::integer(l, u);
		assert(l <= r);
		assert(r <= u);
		rVector.push_back(r);
	}

	int64_t minR = *(min_element(rVector.begin(), rVector.end()));
	int64_t maxR = *(max_element(rVector.begin(), rVector.end()));

	EXPECT_EQ(minR, l);
	EXPECT_EQ(maxR, u);

	double sum = std::accumulate(rVector.begin(), rVector.end(), uint64_t{0});
	double avg = sum / n;


	DEBUG("avg rand integer: ", avg);
	EXPECT_LE(avg, 6.0);
	EXPECT_GE(avg, 4.0);
}

TEST_F(AuxGTest, testRandomReal) {
	Aux::Random::setSeed(1, false);
	std::vector<double> rVector;
	int n = 1000;
	for (int i = 0; i < n; ++i) {
		double r = Aux::Random::real();
		assert(0.0 <= r);
		assert(r < 1.0);
		rVector.push_back(r);
	}

	double sum = std::accumulate(rVector.begin(), rVector.end(), 0.0);
	double avg = sum / n;


	DEBUG("avg rand probability: ", avg);
	EXPECT_LE(avg, 0.6);
	EXPECT_GE(avg, 0.4);
}

TEST_F(AuxGTest, testRandomProbability) {
	Aux::Random::setSeed(1, false);
	std::vector<double> rVector;
	int n = 1000;
	for (int i = 0; i < n; ++i) {
		double r = Aux::Random::probability();
		assert(0.0 <= r);
		assert(r <= 1.0);
		rVector.push_back(r);
	}

	double sum = std::accumulate(rVector.begin(), rVector.end(), 0.0);
	double avg = sum / n;


	DEBUG("avg rand probability: ", avg);
	EXPECT_LE(avg, 0.6);
	EXPECT_GE(avg, 0.4);
}


TEST_F(AuxGTest, testTimer) {
	int64_t sleepTime = 1000; // sleep time in ms
	int64_t tolerance = 20;

	Aux::Timer timer;
	DEBUG("sleeping for ", sleepTime, " ms");
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

	INFO("calculation finished after ", timer.elapsedTag());

}


TEST_F(AuxGTest, testVectorDebug) {
	std::vector<int> vec(10, 42);
	INFO(vec);
}


TEST_F(AuxGTest, testPriorityQueue) {
	typedef std::pair<double, uint64_t> ElemType;

	// fill vector with keys
	std::vector<double> vec = {0.5, 3.5, 4.5, 2.5, 0.75, 1.5, 8.5, 3.25, 4.75, 5.0, 11.5, 0.25};

	// construct pq from vector
	Aux::PrioQueue<double, uint64_t> pq(vec);
	EXPECT_EQ(pq.size(), vec.size());

	ElemType elem = pq.extractMin();
	EXPECT_EQ(0.25, elem.first);
	EXPECT_EQ(11u, elem.second);
	EXPECT_EQ(pq.size(), vec.size() - 1);

	elem = pq.extractMin();
	EXPECT_EQ(0.5, elem.first);
	EXPECT_EQ(0u, elem.second);
	EXPECT_EQ(pq.size(), vec.size() - 2);

	elem = pq.extractMin();
	EXPECT_EQ(0.75, elem.first);
	EXPECT_EQ(4u, elem.second);
	EXPECT_EQ(pq.size(), vec.size() - 3);

	elem = pq.extractMin();
	EXPECT_EQ(1.5, elem.first);
	EXPECT_EQ(5u, elem.second);
	EXPECT_EQ(pq.size(), vec.size() - 4);

	elem = pq.extractMin();
	EXPECT_EQ(2.5, elem.first);
	EXPECT_EQ(3u, elem.second);
	EXPECT_EQ(pq.size(), vec.size() - 5);
}

TEST_F(AuxGTest, testPrioQueueForIntsWithEmptiness) {
	// fill vector with priorities
	std::vector<int64_t> vec = {17, 4, 1, 5, 3, 11, 9, 19, -9, 1, 4, 20, 8, 8};

	Aux::BucketPQ pq(vec, -20, 20);
	EXPECT_EQ(pq.size(), vec.size());

	// delete everything
	while (pq.size() > 0) {
		pq.extractMin();
	}

	// reinsert entries
	for (uint64_t i = 0; i < vec.size(); ++i) {
		pq.insert(vec[i], i);
	}
	EXPECT_EQ(pq.size(), vec.size());

	// check top
	std::pair<int64_t, uint64_t> mini = pq.extractMin();
	EXPECT_EQ(mini.first, -9);
	EXPECT_EQ(mini.second, 8u);
}

TEST_F(AuxGTest, testPrioQueueForInts) {
	// fill vector with priorities
	std::vector<int64_t> vec;

	// 0-4
	vec.push_back(17);
	vec.push_back(4);
	vec.push_back(1);
	vec.push_back(5);
	vec.push_back(3);

	// 5-9
	vec.push_back(11);
	vec.push_back(9);
	vec.push_back(19);
	vec.push_back(-9);
	vec.push_back(1);

	// 10-14
	vec.push_back(4);
	vec.push_back(17);
	vec.push_back(8);
	vec.push_back(8);
	vec.push_back(12);

	// 15-19
	vec.push_back(16);
	vec.push_back(14);
	vec.push_back(11);
	vec.push_back(7);
	vec.push_back(7);

	// 20-23
	vec.push_back(7);
	vec.push_back(-4);
	vec.push_back(8);
	vec.push_back(0);

	// construct pq from vector
	Aux::BucketPQ pq(vec, -100, 100);

	// check op: extractMin
	std::pair<int64_t, uint64_t> mini = pq.extractMin();
	EXPECT_EQ(mini.first, -9);
	EXPECT_EQ(mini.second, 8u);

	// check op: changeKey
	pq.changeKey(-20, 0);
	mini = pq.extractMin();
	EXPECT_EQ(mini.first, -20);
	EXPECT_EQ(mini.second, 0u);

	// multiply vec by -1 and try again
	for (int64_t& currkey : vec) {
		currkey *= -1;
	}
	Aux::BucketPQ pq2(vec, -100, 100);
	mini = pq2.extractMin();
	EXPECT_EQ(mini.first, -19);
	EXPECT_EQ(mini.second, 7u);

	// check op: changeKey
	pq.changeKey(-20, 0);
	mini = pq.extractMin();
	EXPECT_EQ(mini.first, -20);
	EXPECT_EQ(mini.second, 0u);
}

//FIXME make this working again
/*TEST_F(AuxGTest, testLogging) {
	std::string cl = Aux::currentLogLevel();
	//FATAL("FATAL ERROR");
	//ERROR("This may be here");
	Aux::setLoglevel("ERROR");
	//FATAL("fatal error"<<3<<2.f);
	//ERROR("normal error"<<3<<2.f);
	//WARN("just a warning"<<3<<2.f);
	//INFO("for the sake of information", 3, 2.f);
	//DEBUG("important debug outputs", 3, 2.f);
	//TRACE("trace", 3, 2.f);
	Aux::setLoglevel("INFO");
	EXPECT_EQ("INFO",Aux::currentLogLevel());
	//FATAL("fatal error", 3, 2.f);
	//ERROR("normal error", 3, 2.f);
	//WARN("just a warning", 3, 2.f);
	//INFO("for the sake of information", 3, 2.f);
	//DEBUG("important debug outputs", 3, 2.f);
	//TRACE("trace", 3, 2.f);
	Aux::setLoglevel("TRACE");
	EXPECT_EQ("TRACE",Aux::currentLogLevel());
	//FATAL("fatal error", 3, 2.f);
	//ERROR("normal error", 3, 2.f);
	//WARN("just a warning", 3, 2.f);
	//INFO("for the sake of information", 3, 2.f);
	//DEBUG("important debug outputs", 3, 2.f);
	//TRACE("trace", 3, 2.f);
	Aux::setLoglevel(cl);	
}*/

TEST_F(AuxGTest, testFormatting) {
	using Aux::toStringF;
	EXPECT_EQ("",toStringF(""));
	EXPECT_EQ("",toStringF("%s", ""));
	EXPECT_EQ("%",toStringF("%%", ""));
	EXPECT_EQ("aaabbbaaa",toStringF("aaa%saaa", "bbb"));
	EXPECT_THROW(toStringF("%"), std::invalid_argument);
	EXPECT_THROW(toStringF("%i"), std::invalid_argument);
	EXPECT_THROW(toStringF("%s"), std::invalid_argument);
}

TEST_F(AuxGTest, testRandomChoice) {
	std::vector<uint64_t> data = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	for (uint64_t i = 0; i < 1000; ++i) {
		std::ignore = Aux::Random::choice(data);
	}
}


TEST_F(AuxGTest, testRandomWeightedChoice) {
	std::vector<std::pair<uint64_t, double> > data = {{0, 1.0}, {1, 0.0}};
	for (uint64_t i = 0; i < 100; ++i) {
		auto element = Aux::Random::weightedChoice(data);
		EXPECT_EQ(0u, element);
	}
}

TEST_F(AuxGTest, testRandomIndex) {
	using namespace Aux::Random;
	setSeed(1, false);
	
	for (unsigned i = 0; i < 10; i++) {
		EXPECT_EQ(0u, index(1));
	}
	
	for (unsigned i = 0; i < 100; i++) {
		auto tmp = index(10);
		EXPECT_LE(tmp, 9u);
		EXPECT_GE(tmp, 0u);
	}
}

TEST_F(AuxGTest, testSplit) {
	using Vec = std::vector<std::string>;
	using namespace Aux::StringTools;
	
	EXPECT_EQ(Vec{}, split(""));
	EXPECT_EQ(Vec{""}, split(" "));
	
	{
		auto expected = Vec{"", ""};
		EXPECT_EQ(expected, split("  "));
	}
	{
		auto expected = Vec{"", "a"};
		EXPECT_EQ(expected, split(" a"));
	}
	{
		auto expected = Vec{"a"};
		EXPECT_EQ(expected, split("a "));
	}
	{
		auto expected = Vec{"a"};
		EXPECT_EQ(expected, split("a"));
	}
	{
		auto expected = Vec{"a", "b"};
		EXPECT_EQ(expected, split("a b"));
	}
	{
		auto expected = Vec{"", "a", "b"};
		EXPECT_EQ(expected, split(" a b "));
	}
	{
		auto expected = Vec{"abc", "def", "ghi"};
		EXPECT_EQ(expected, split("abc def ghi"));
	}
	{
		auto expected = Vec{"abc", "def", "ghi"};
		EXPECT_EQ(expected, split("abc def ghi "));
	}
}


TEST_F(AuxGTest, testSetIntersector) {
	uint64_t n = 50;
	std::vector<uint64_t> A = {0, 5, 7, 3, 11, 22, 45, 17, 23, 11, 45};
	std::vector<uint64_t> B = {27, 11, 13, 17, 44, 1, 2, 3, 9, 44, 11, 20};

	Aux::SetIntersector<uint64_t> si(n);
	std::set<uint64_t> intersection = si.intersect(A, B);

	std::set<uint64_t> expectedResult = {3, 11, 17};
	EXPECT_EQ(expectedResult, intersection);
}

TEST_F(AuxGTest, testEnforce) {
	using Aux::enforce;
	EXPECT_THROW(enforce(false), std::runtime_error);
	EXPECT_NO_THROW(enforce(true));
	
	EXPECT_THROW(enforce(false, "foo"), std::runtime_error);
	EXPECT_NO_THROW(enforce(true, "bar"));
	
	EXPECT_THROW(enforce<std::logic_error>(false, "foo"), std::logic_error);
	EXPECT_NO_THROW(enforce<std::logic_error>(true, "foo"));
	
	std::string msg = "some message in a std::string";
	
	EXPECT_THROW(enforce(false, msg), std::runtime_error);
	EXPECT_NO_THROW(enforce(true, msg));
	
	EXPECT_THROW(enforce<std::logic_error>(false, msg), std::logic_error);
	EXPECT_NO_THROW(enforce<std::logic_error>(true, msg));
}

TEST_F(AuxGTest, testEnforceOpened) {
	using Aux::enforceOpened;
	std::ifstream stream;
	EXPECT_THROW(enforceOpened(stream), std::runtime_error);
}

TEST_F(AuxGTest, testNumberParsingInteger) {
	using namespace Aux::Parsing;
	const std::string str = "0 00 1 123 001 1200 12345678    ";
	std::vector<unsigned> expectedValues = {0, 0, 1, 123, 1, 1200, 12345678};
	auto it = str.begin();
	auto end = str.end();
	std::size_t i = 0;
	while(it != end) {
		unsigned result;
		std::tie(result, it) = strTo<unsigned>(it, end);
		EXPECT_EQ(expectedValues[i], result);
		++i;
	}
	
	EXPECT_EQ(it, end);
}

TEST_F(AuxGTest, testNumberParsingSignedInteger) {
	using namespace Aux::Parsing;
	const std::string str = "-0 -00 -1 -123 -001 -1200 -12345678    ";
	std::vector<int> expectedValues = {0, 0, -1, -123, -1, -1200, -12345678};
	auto it = str.begin();
	auto end = str.end();
	std::size_t i = 0;
	while(it != end) {
		int result;
		std::tie(result, it) = strTo<int>(it, end);
		EXPECT_EQ(expectedValues[i], result);
		++i;
	}
	
	EXPECT_EQ(it, end);
}

TEST_F(AuxGTest, testOverflowCatching) {
	using namespace Aux::Parsing;
	const std::string str = "1000";
	EXPECT_THROW(
			(strTo<uint8_t, std::string::const_iterator, Aux::Checkers::Enforcer>(
				str.begin(), str.end())),
			std::runtime_error
			);
	EXPECT_THROW(
			(strTo<int8_t, std::string::const_iterator, Aux::Checkers::Enforcer>(
				str.begin(), str.end())),
			std::runtime_error
			);
}

TEST_F(AuxGTest, testNumberParsingBasicReal) {
	using namespace Aux::Parsing;
	const std::string str = 
		"0 00 1 123 001 1200 12345678    "
		"0.00000 -0000.000 -0000.000e-100"
		;
	std::vector<double> expectedValues = {
		0, 0, 1, 123, 1, 1200, 12345678,
		0, 0, 0
	};
	auto it = str.begin();
	auto end = str.end();
	std::size_t i = 0;
	while(it != end) {
		double result;
		std::tie(result, it) = strTo<double>(it, end);
		EXPECT_EQ(expectedValues[i], result);
		++i;
	}
	
	EXPECT_EQ(it, end);
}


TEST_F(AuxGTest, testNumberParsingAdvancedReal) {
	using Aux::Parsing::strTo;
	auto helper = [](const std::string& str, double expected) {
		auto result = std::get<0>(strTo<double>(str.begin(), str.end()));
		EXPECT_DOUBLE_EQ(result, expected);
	};
#define TEST_CASE_REAL(number) helper(#number, number)
	TEST_CASE_REAL(0);
	TEST_CASE_REAL(1);
	TEST_CASE_REAL(8.76464e+23);
	TEST_CASE_REAL(5.56628e+12);
	TEST_CASE_REAL(9.563574701896000000e+23);
	TEST_CASE_REAL(9.563574701896270940e+23);
	TEST_CASE_REAL(-669585235219883571019776.);
	TEST_CASE_REAL(9.563574701896270946e+23);
	TEST_CASE_REAL(141265720771594800000000.);
	TEST_CASE_REAL(141265720771594850000000.);
	TEST_CASE_REAL(1.4126572077159485e+23);
	TEST_CASE_REAL(1.4126572077159480e+23);
	TEST_CASE_REAL(-784008854951861199831040.);
	TEST_CASE_REAL(-78400885495186119983104.0);
#undef TEST_CASE_REAL
}

TEST_F(AuxGTest, testBloomFilter) {
	Aux::Random::setSeed(1, false);
	Aux::BloomFilter bf(5);
	uint64_t size = 750;
	std::set<uint64_t> randomKeys;

	while (randomKeys.size() < size) {
		uint64_t key = Aux::Random::integer();
		bf.insert(key);
		randomKeys.insert(key);
	}

	for (auto key: randomKeys) {
		EXPECT_TRUE(bf.isMember(key));
	}

	for (uint64_t i = 0; i < size; ++i) {
		uint64_t key = Aux::Random::integer();
		if (randomKeys.count(key) == 0) {
			if (bf.isMember(key)) {
				WARN("Bloom filter: Maybe only false positive, maybe indication of bug!");
			}
		}
	}
}

#endif /*NOGTEST */
