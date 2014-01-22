#include "EPPFactory.h"
#include "PLM.h"
#include "PLP.h"
#include "../overlap/HashingOverlapper.h"
#include "PLM2.h"

namespace NetworKit {

EPP EPPFactory::make(count ensembleSize, std::string baseAlgorithm, std::string finalAlgorithm) {
	EPP ensemble;

	for (count i = 0; i < ensembleSize; ++i) {
		if (baseAlgorithm == "PLP") {
			ensemble.addBaseClusterer(*(new PLP()));
		} else if (baseAlgorithm == "PLM") {
			ensemble.addBaseClusterer(*(new PLM2(false)));
		} else {
			throw std::runtime_error("unknown base algorithm name");
		}
	}

	if (finalAlgorithm == "PLM") {
		ensemble.setFinalClusterer(*(new PLM2(false)));
	} else if (finalAlgorithm == "PLP") {
		ensemble.setFinalClusterer(*(new PLP()));
	} else if (finalAlgorithm == "PLMR") {
		ensemble.setFinalClusterer(*(new PLM2(true)));
	} else throw std::runtime_error("unknown final algorithm name");

	ensemble.setOverlapper(*(new HashingOverlapper));

	return ensemble;
}

} /* namespace NetworKit */
