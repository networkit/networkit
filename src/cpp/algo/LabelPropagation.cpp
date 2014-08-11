#include "LabelPropagation.h"

namespace NetworKit {

void LabelPropagation::run() {


	auto propagateLabels = [&](node u){
		//
		std::set<label> local;
		G.forNeighborsOf(u, [&](node v){
			// 
		});
	};

	G.balancedParallelForNodes(propagateLabels);
}

}
