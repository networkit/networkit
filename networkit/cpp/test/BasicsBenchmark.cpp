/*
 * BasicsBenchmark.cpp
 *
 *  Created on: 01.02.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NOGTEST

#include "BasicsBenchmark.h"

namespace NetworKit {

TEST_F(BasicsBenchmark, sequentialSum) {
	Aux::Timer runtime;

	int64_t n = 1e+9;
	double sum = 0.0;
	runtime.start();
	for (int64_t i = 0; i < n; ++i) {
		sum += i;
	}
	runtime.stop();

	INFO("sum = " , sum , " [" , runtime.elapsed().count() , " ms ]");
}

TEST_F(BasicsBenchmark, parallelSumIncorrect) {
	Aux::Timer runtime;

	int64_t n = 1e+9;
	double sum = 0.0;
	runtime.start();
	#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
		sum += i;
	}
	runtime.stop();

	INFO("sum = " , sum , " [" , runtime.elapsed().count() , " ms ]");
}

TEST_F(BasicsBenchmark, parallelSumAtomicUpdate) {
	Aux::Timer runtime;

	int64_t n = 1e+9;
	double sum = 0.0;
	runtime.start();
	#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
		#pragma omp atomic update
		sum += i;
	}
	runtime.stop();

	INFO("sum = " , sum , " [" , runtime.elapsed().count() , " ms ]");

}


TEST_F(BasicsBenchmark, parallelSumReduction) {
	Aux::Timer runtime;

	int64_t n = 1e+9;
	double sum = 0.0;
	runtime.start();
	#pragma omp parallel for reduction(+:sum)
	for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
		sum += i;
	}
	runtime.stop();

	INFO("sum = " , sum , " [" , runtime.elapsed().count() , " ms ]");

}

//
//TEST_F(BasicsBenchmark, parallelSumCritical) {
//	Aux::Timer runtime;
//
//	int64_t n = 1e+9;
//	double sum = 0.0;
//	runtime.start();
//	#pragma omp parallel for
//	for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
//		#pragma omp critical
//		sum += i;
//	}
//	runtime.stop();
//
//	INFO("sum = " , sum , " [" , runtime.elapsed().count() , " ms ]");
//
//}


TEST_F(BasicsBenchmark, seqVectorWrite) {
	Aux::Timer runtime;
	int64_t n = 1e+8;

	std::vector<int64_t> vec;
	vec.resize(n);

	runtime.start();
	for (int64_t i = 0; i < n; ++i) {
		vec[i] = i;
	}
	runtime.stop();

	INFO("vector written in [" , runtime.elapsed().count() , " ms ]");
}



TEST_F(BasicsBenchmark, parVectorWrite) {
	Aux::Timer runtime;
	int64_t n = 1e+8;

	std::vector<int64_t> vec;
	vec.resize(n);

	runtime.start();
	#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
		vec[i] = i;
	}
	runtime.stop();

	INFO("vector written in [" , runtime.elapsed().count() , " ms ]");
}




TEST_F(BasicsBenchmark, lambdaSummation_seq) {
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

	INFO("sum = " , sum , " [" , runtime.elapsed().count() , " ms ]");
}

TEST_F(BasicsBenchmark, lambdaSummation_parWrong) {
	Aux::Timer runtime;
	int64_t n = 1e+9;
	double sum = 0.0;

	auto func = [&](int64_t x) {
		sum += x;
	};

	runtime.start();
	#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
		func(i);
	}
	runtime.stop();

	INFO("sum = " , sum , " [" , runtime.elapsed().count() , " ms ]");
}


TEST_F(BasicsBenchmark, lambdaSummation_par_atomic) {
	Aux::Timer runtime;
	int64_t n = 1e+9;
	double sum = 0.0;

	auto func = [&](int64_t x) {
		#pragma omp atomic update
		sum += x;
	};

	runtime.start();
	#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
		func(i);
	}
	runtime.stop();

	INFO("sum = " , sum , " [" , runtime.elapsed().count() , " ms ]");
}

TEST_F(BasicsBenchmark, lambdaSummation_par_reduction) {
	Aux::Timer runtime;
	int64_t n = 1e+9;
	double sum = 0.0;

	auto func = [&](int64_t x) {
		sum += x;
	};

	runtime.start();
	#pragma omp parallel for reduction(+:sum)
	for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
		func(i);
	}
	runtime.stop();

	INFO("sum = " , sum , " [" , runtime.elapsed().count() , " ms ]");
}



TEST_F(BasicsBenchmark, lambdaVectorWrite_seq) {
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

	INFO("vector written in [" , runtime.elapsed().count() , " ms ]");
}

TEST_F(BasicsBenchmark, lambdaVectorWrite_par) {
	Aux::Timer runtime;
	int64_t n = 1e+8;

	std::vector<int64_t> vec;
	vec.resize(n);

	auto insert = [&](int64_t i, int64_t x) {
		vec[i] = x;
	};

	runtime.start();
	#pragma omp parallel for
	for (omp_index i = 0; i < static_cast<omp_index>(n); ++i) {
		insert(i, i);
	}
	runtime.stop();

	INFO("vector written in [" , runtime.elapsed().count() , " ms ]");
}

} /* namespace NetworKit */

#endif /* NOGTEST */
