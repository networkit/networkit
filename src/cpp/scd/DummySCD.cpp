#include "DummySCD.h"

namespace NetworKit {

DummySCD::DummySCD(const Graph& G) : SelectiveCommunityDetector(G) {

}


void DummySCD::run(std::unordered_set<node> seeds)  {
	
}

std::unordered_map<node, std::unordered_set<node> > DummySCD::getResult() {

}

std::unordered_map<node, double> DummySCD::getTimings() {

}


}


