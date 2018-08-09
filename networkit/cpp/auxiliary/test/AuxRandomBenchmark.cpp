#include "AuxRandomBenchmark.h"

#include "../Random.h"
#include "../Timer.h"

#include <atomic>
#include <random>
#include "omp.h"

namespace NetworKit {
template<typename F>
static double measure(F f, const size_t iterations = 50000000) {
	Aux::Timer timer;
	timer.start();

	for (size_t i = 0; i < iterations; i++) {
		f();
	}

	timer.stop();
	const auto ms = timer.elapsedMilliseconds();
	return (1.0e6 * ms / iterations);
}

TEST_F(AuxRandomBenchmark, benchmarkInteger) {
	uint64_t tmp = 0;
	auto atime = measure([&] {
		tmp += Aux::Random::integer();
	});
	volatile auto dummy = tmp;
	std::cout << "Average time of operation: " << atime << "ns\n";
}

TEST_F(AuxRandomBenchmark, benchmarkLocalDistrInteger) {
	uint64_t tmp = 0;

	auto &prng = Aux::Random::getURNG();
	std::uniform_int_distribution<uint64_t> distr;

	auto atime = measure([&] {
		tmp += distr(prng);
	});

	volatile auto dummy = tmp;
	std::cout << "Average time of operation: " << atime << "ns\n";
}

TEST_F(AuxRandomBenchmark, benchmarkIntegerParallel) {
	double atime{0};

	#pragma omp parallel
	{
		uint64_t tmp = 0;
		auto local_time = measure([&] {
			tmp += Aux::Random::integer();
		}) / omp_get_num_threads();
		volatile auto dummy = tmp;

		#pragma omp atomic
		atime += local_time;
	}

	std::cout << "Average time of operation: " << atime << "ns\n";
}

TEST_F(AuxRandomBenchmark, benchmarkLocalDistrIntegerParallel) {
	double atime{0};

	#pragma omp parallel
	{
		uint64_t tmp = 0;
		auto &prng = Aux::Random::getURNG();
		std::uniform_int_distribution<uint64_t> distr;

		auto local_time = measure([&] {
			tmp += distr(prng);
		}) / omp_get_num_threads();

		volatile auto dummy = tmp;

		#pragma omp atomic
		atime += local_time;
	}

	std::cout << "Average time of operation: " << atime << "ns\n";
}


TEST_F(AuxRandomBenchmark, benchmarkProb) {
	double tmp = 0.0;
	auto atime = measure([&] {
		tmp += Aux::Random::probability();
	});
	volatile auto dummy = tmp;
	std::cout << "Average time of operation: " << atime << "ns\n";
}

TEST_F(AuxRandomBenchmark, benchmarkLocalDistrProb) {
	double tmp = 0.0;

	auto &prng = Aux::Random::getURNG();
	std::uniform_real_distribution<double> distr{0.0, 1.0};

	auto atime = measure([&] {
		tmp += distr(prng);
	});

	volatile auto dummy = tmp;
	std::cout << "Average time of operation: " << atime << "ns\n";
}

TEST_F(AuxRandomBenchmark, benchmarkProbParallel) {
	double atime{0};

	#pragma omp parallel
	{
		double tmp = 0.0;
		auto local_time = measure([&] {
			tmp += Aux::Random::probability();
		}) / omp_get_num_threads();
		volatile auto dummy = tmp;

		#pragma omp atomic
		atime += local_time;
	}

	std::cout << "Average time of operation: " << atime << "ns\n";
}

TEST_F(AuxRandomBenchmark, benchmarkLocalDistrProbParallel) {
	double atime{0};

	#pragma omp parallel
	{
		double tmp = 0.0;
		auto &prng = Aux::Random::getURNG();
		std::uniform_real_distribution<double> distr{0.0, 1.0};

		auto local_time = measure([&] {
			tmp += distr(prng);
		}) / omp_get_num_threads();

		volatile auto dummy = tmp;

		#pragma omp atomic
		atime += local_time;
	}

	std::cout << "Average time of operation: " << atime << "ns\n";
}

} // ! namespace NetworKit
