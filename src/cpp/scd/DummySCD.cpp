#include "DummySCD.h"

namespace NetworKit {

DummySCD::DummySCD(const Graph& G) : SelectiveCommunityDetector(G) {

}


void DummySCD::run(std::set<unsigned int> seeds)  {
	for (auto seed : seeds) {
		result[seed] = {seed};
	}

}


}


