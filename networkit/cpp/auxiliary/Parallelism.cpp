#include "Parallelism.h"
#include "Log.h"

 // OpenMP
#include <omp.h>



void Aux::setNumberOfThreads(int nThreads) {
#ifdef _OPENMP
		omp_set_num_threads(nThreads);
#else
		ERROR("Thread option ignored since OpenMP is deactivated.");
#endif
}


int Aux::getCurrentNumberOfThreads() {
	#ifdef _OPENMP
		return omp_get_num_threads();
	#else
		return 1;
		ERROR("OpenMP is not available");
	#endif
}


int Aux::getMaxNumberOfThreads() {
	#ifdef _OPENMP
		return omp_get_max_threads();
	#else
		return 1;
		ERROR("OpenMP is not available");
	#endif
}

void Aux::enableNestedParallelism() {
	#ifdef _OPENMP
	omp_set_nested(1); // enable nested parallelism
	#else
	ERROR("OpenMP is not available");
	#endif
}
