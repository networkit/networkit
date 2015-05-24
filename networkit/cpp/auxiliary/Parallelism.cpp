#include "Parallelism.h"
#include "Log.h"

 // OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif


void Aux::setNumberOfThreads(int nThreads) {
	#ifdef _OPENMP
		omp_set_num_threads(nThreads);
	#endif
}


int Aux::getCurrentNumberOfThreads() {
	#ifdef _OPENMP
		return omp_get_num_threads();
	#else
		return 1;
	#endif
}


int Aux::getMaxNumberOfThreads() {
	#ifdef _OPENMP
		return omp_get_max_threads();
	#else
		return 1;
	#endif
}

void Aux::enableNestedParallelism() {
	#ifdef _OPENMP
		omp_set_nested(1); // enable nested parallelism
	#endif
}

