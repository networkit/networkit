#include "EPPFactory.h"
#include "PLM.h"
#include "PLP.h"
#include "../overlap/HashingOverlapper.h"

namespace NetworKit {

EPP EPPFactory::make(count ensembleSize, std::string baseAlgorithm, std::string finalAlgorithm) {
	EPP ensemble;

	for (count i = 0; i < ensembleSize; ++i) {
		if (baseAlgorithm == "PLP") {
			ensemble.addBaseClusterer(*(new PLP()));
		} else if (baseAlgorithm == "PLM") {
			ensemble.addBaseClusterer(*(new PLM()));
		} else {
			throw std::runtime_error("unknown base algorithm name");
		}
	}

	if (finalAlgorithm == "PLM") {
		ensemble.setFinalClusterer(*(new PLM()));
	} else if (finalAlgorithm == "PLP") {
		ensemble.setFinalClusterer(*(new PLM()));
	}

	ensemble.setOverlapper(*(new HashingOverlapper));

	return ensemble;
}

} /* namespace NetworKit */
