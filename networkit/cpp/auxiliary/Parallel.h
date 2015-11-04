/*
 * Header which provides parallel STL implementations when available.
 */
#ifndef PARALLEL_H
#define PARALLEL_H

#if defined __GNUC__ && defined _OPENMP
#include <parallel/algorithm>
#else
#define NETWORKIT_NO_PARALLEL_STL
#include <algorithm>
#endif

namespace Aux {
	namespace Parallel {
		#ifdef NETWORKIT_NO_PARALLEL_STL
		using std::sort;
		using std::max_element;
		#else
		using __gnu_parallel::sort;
		using __gnu_parallel::max_element;
		#endif
	}
}

#endif