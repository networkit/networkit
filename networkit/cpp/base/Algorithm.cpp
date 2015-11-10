#include "Algorithm.h"
#include "../auxiliary/SignalHandling.h"
#include "../auxiliary/Log.h"
#include <exception>

namespace NetworKit {
	Algorithm::Algorithm() : hasRun(false) {

	}

	std::string Algorithm::toString() const {
		throw std::runtime_error("TODO: implement in subclass and return string representation");
	}

	bool Algorithm::isParallel() const {
		throw std::runtime_error("TODO: Implement in subclass");
		return false;
	}

} /* NetworKit */
