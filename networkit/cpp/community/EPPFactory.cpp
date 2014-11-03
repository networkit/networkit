#include "EPPFactory.h"
#include "PLM.h"
#include "PLP.h"
#include "../overlap/HashingOverlapper.h"
#include "PLM.h"

namespace NetworKit {

EPP EPPFactory::make(const Graph& G, count ensembleSize, std::string baseAlgorithm, std::string finalAlgorithm) {
	EPP ensemble(G);

	for (count i = 0; i < ensembleSize; ++i) {
		if (baseAlgorithm == "PLP") {
			ensemble.addBaseClusterer(*(new PLP(G)));
		} else if (baseAlgorithm == "PLM") {
			ensemble.addBaseClusterer(*(new PLM(G,false)));
		} else {
			throw std::runtime_error("unknown base algorithm name");
		}
	}

	if (finalAlgorithm == "PLM") {
		ensemble.setFinalClusterer(*(new PLM(G,false)));
	} else if (finalAlgorithm == "PLP") {
		ensemble.setFinalClusterer(*(new PLP(G)));
	} else if (finalAlgorithm == "PLMR") {
		ensemble.setFinalClusterer(*(new PLM(G,true)));
	} else throw std::runtime_error("unknown final algorithm name");

	ensemble.setOverlapper(*(new HashingOverlapper));

	return ensemble;
}

} /* namespace NetworKit */
