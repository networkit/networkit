/*
 * BenchmarkGTest.cpp
 *
 *  Created on: 30.01.2013
 *      Author: cls
 */

#include "BenchmarkGTest.h"

namespace EnsembleClustering {

BenchmarkGTest::BenchmarkGTest() {
	// TODO Auto-generated constructor stub

}

BenchmarkGTest::~BenchmarkGTest() {
	// TODO Auto-generated destructor stub
}


TEST_F(BenchmarkGTest, sequentialSum) {
	Aux::Timer runtime;

	int64_t n = 1e+9;
	double sum = 0.0;
	runtime.start();
	for (int64_t i = 0; i < n; ++i) {
		sum += i;
	}
	runtime.stop();

	INFO("sum = " << sum << " [" << runtime.elapsed().count() << " ms ]");
}

TEST_F(BenchmarkGTest, parallelSumIncorrect) {
	Aux::Timer runtime;

	int64_t n = 1e+9;
	double sum = 0.0;
	runtime.start();
	#pragma omp parallel for
	for (int64_t i = 0; i < n; ++i) {
		sum += i;
	}
	runtime.stop();

	INFO("sum = " << sum << " [" << runtime.elapsed().count() << " ms ]");
}

TEST_F(BenchmarkGTest, parallelSumAtomicUpdate) {
	Aux::Timer runtime;

	int64_t n = 1e+9;
	double sum = 0.0;
	runtime.start();
	#pragma omp parallel for
	for (int64_t i = 0; i < n; ++i) {
		#pragma omp atomic update
		sum += i;
	}
	runtime.stop();

	INFO("sum = " << sum << " [" << runtime.elapsed().count() << " ms ]");

}

//
//TEST_F(BenchmarkGTest, parallelSumCritical) {
//	Aux::Timer runtime;
//
//	int64_t n = 1e+9;
//	double sum = 0.0;
//	runtime.start();
//	#pragma omp parallel for
//	for (int64_t i = 0; i < n; ++i) {
//		#pragma omp critical
//		sum += i;
//	}
//	runtime.stop();
//
//	INFO("sum = " << sum << " [" << runtime.elapsed().count() << " ms ]");
//
//}


TEST_F(BenchmarkGTest, seqVectorWrite) {
	Aux::Timer runtime;
	int64_t n = 1e+8;

	std::vector<int64_t> vec;
	vec.resize(n);

	runtime.start();
	for (int64_t i = 0; i < n; ++i) {
		vec[i] = i;
	}
	runtime.stop();

	INFO("vector written in [" << runtime.elapsed().count() << " ms ]");
}



TEST_F(BenchmarkGTest, parVectorWrite) {
	Aux::Timer runtime;
	int64_t n = 1e+8;

	std::vector<int64_t> vec;
	vec.resize(n);

	runtime.start();
	#pragma omp parallel for
	for (int64_t i = 0; i < n; ++i) {
		vec[i] = i;
	}
	runtime.stop();

	INFO("vector written in [" << runtime.elapsed().count() << " ms ]");
}




TEST_F(BenchmarkGTest, lambdaSummation_seq) {
	Aux::Timer runtime;
	int64_t n = 1e+9;
	double sum = 0.0;

	auto func = [&](int64_t x) {
		sum += x;
	};

	runtime.start();
	for (int64_t i = 0; i < n; ++i) {
		func(i);
	}
	runtime.stop();

	INFO("sum = " << sum << " [" << runtime.elapsed().count() << " ms ]");
}

TEST_F(BenchmarkGTest, lambdaSummation_parWrong) {
	Aux::Timer runtime;
	int64_t n = 1e+9;
	double sum = 0.0;

	auto func = [&](int64_t x) {
		sum += x;
	};

	runtime.start();
	#pragma omp parallel for
	for (int64_t i = 0; i < n; ++i) {
		func(i);
	}
	runtime.stop();

	INFO("sum = " << sum << " [" << runtime.elapsed().count() << " ms ]");
}


TEST_F(BenchmarkGTest, lambdaSummation_par) {
	Aux::Timer runtime;
	int64_t n = 1e+9;
	double sum = 0.0;

	auto func = [&](int64_t x) {
		#pragma omp atomic update
		sum += x;
	};

	runtime.start();
	#pragma omp parallel for
	for (int64_t i = 0; i < n; ++i) {
		func(i);
	}
	runtime.stop();

	INFO("sum = " << sum << " [" << runtime.elapsed().count() << " ms ]");
}



TEST_F(BenchmarkGTest, lambdaVectorWrite_seq) {
	Aux::Timer runtime;
	int64_t n = 1e+8;

	std::vector<int64_t> vec;
	vec.resize(n);

	auto insert = [&](int64_t i, int64_t x) {
		vec[i] = x;
	};

	runtime.start();
	for (int64_t i = 0; i < n; ++i) {
		insert(i, i);
	}
	runtime.stop();

	INFO("vector written in [" << runtime.elapsed().count() << " ms ]");
}

TEST_F(BenchmarkGTest, lambdaVectorWrite_par) {
	Aux::Timer runtime;
	int64_t n = 1e+8;

	std::vector<int64_t> vec;
	vec.resize(n);

	auto insert = [&](int64_t i, int64_t x) {
		vec[i] = x;
	};

	runtime.start();
	#pragma omp parallel for
	for (int64_t i = 0; i < n; ++i) {
		insert(i, i);
	}
	runtime.stop();

	INFO("vector written in [" << runtime.elapsed().count() << " ms ]");
}

} /* namespace EnsembleClustering */
