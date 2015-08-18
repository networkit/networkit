#include "EPPFactory.h"
#include "PLM.h"
#include "PLP.h"
#include "../overlap/HashingOverlapper.h"
#include <memory>

namespace NetworKit {

namespace EPPFactory {
	EPP make(const Graph& G, count ensembleSize, std::string baseAlgorithm, std::string finalAlgorithm) {
		EPP ensemble(G);

		for (count i = 0; i < ensembleSize; ++i) {
			CommunityDetectionAlgorithm* p;
			if (baseAlgorithm == "PLP") {
				p = new PLP(G);
			} else if (baseAlgorithm == "PLM") {
				p = new PLM(G,false);
			} else {
				throw std::runtime_error("unknown base algorithm name");
			}
			std::unique_ptr<CommunityDetectionAlgorithm> base;
			base.reset(p);
			ensemble.addBaseClusterer(base);
		}

		CommunityDetectionAlgorithm* p;
		if (finalAlgorithm == "PLM") {
			p = new PLM(G,false);
		} else if (finalAlgorithm == "PLP") {
			p = new PLP(G);
		} else if (finalAlgorithm == "PLMR") {
			p = new PLM(G,true);
		} else throw std::runtime_error("unknown final algorithm name");
		std::unique_ptr<CommunityDetectionAlgorithm> final;
		final.reset(p);
		ensemble.setFinalClusterer(final);


		auto overlap = new HashingOverlapper();
		std::unique_ptr<Overlapper> overlap_ptr(overlap);
		ensemble.setOverlapper(overlap_ptr);

		return std::move(ensemble);
	}

	EPP* makePtr(const Graph& G, count ensembleSize, std::string baseAlgorithm, std::string finalAlgorithm) {
		EPP* ensemble = new EPP(G);

		for (count i = 0; i < ensembleSize; ++i) {
			CommunityDetectionAlgorithm* p;
			if (baseAlgorithm == "PLP") {
				p = new PLP(G);
			} else if (baseAlgorithm == "PLM") {
				p = new PLM(G,false);
			} else {
				throw std::runtime_error("unknown base algorithm name");
			}
			std::unique_ptr<CommunityDetectionAlgorithm> base;
			base.reset(p);
			ensemble->addBaseClusterer(base);
		}

		CommunityDetectionAlgorithm* p;
		if (finalAlgorithm == "PLM") {
			p = new PLM(G,false);
		} else if (finalAlgorithm == "PLP") {
			p = new PLP(G);
		} else if (finalAlgorithm == "PLMR") {
			p = new PLM(G,true);
		} else throw std::runtime_error("unknown final algorithm name");
		std::unique_ptr<CommunityDetectionAlgorithm> final;
		final.reset(p);
		ensemble->setFinalClusterer(final);


		auto overlap = new HashingOverlapper();
		std::unique_ptr<Overlapper> overlap_ptr(overlap);
		ensemble->setOverlapper(overlap_ptr);

		return ensemble;
	}

}

}
