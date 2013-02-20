/*
 * EnsembleClusterer.cpp
 *
 *  Created on: 17.12.2012
 *      Author: cls
 */

#include "EnsembleClusterer.h"

namespace EnsembleClustering {

EnsembleClusterer::EnsembleClusterer() :
		Clusterer() {
	this->finalClusterer = NULL;
	this->qm = NULL;
	this->h = -1;	// start with 0 in loop
}

EnsembleClusterer::~EnsembleClusterer() {
	// TODO Auto-generated destructor stub
}

void EnsembleClusterer::setQualityMeasure(QualityMeasure& qm) {
	this->qm = &qm;
}

void EnsembleClusterer::addBaseClusterer(Clusterer& base) {
	this->baseClusterers.push_back(&base);
}

void EnsembleClusterer::setFinalClusterer(Clusterer& final) {
	this->finalClusterer = &final;
}

void EnsembleClusterer::setOverlapper(Overlapper& overlap) {
	this->overlap = &overlap;
}

Clustering EnsembleClusterer::run(Graph& G) {
	// DEBUG
	INFO("STARTING EnsembleClusterer on G=" << G.toString());
	// DEBUG

	// config flags
	bool calcBaseClusteringDissimilarity = false;

	// TODO: add setter methods
	// sub-algorithms
	ClusterContracter contract;
	ClusteringProjector project;

	// hierarchies
	std::vector<Graph> graph; // hierarchy of graphs G^{i}
	std::vector<Clustering> clustering; // hierarchy of core clusterings \zeta^{i}
	std::vector<Clustering> clusteringBack; // hierarchy of core clusterings projected back to the original graph
	std::vector<NodeMap<node> > map; // hierarchy of maps M^{i->i+1}
	std::vector<double> quality; // hierarchy of clustering quality values q^{i} = q(\zeta^{i}, G^{0})

	// other data collections
	std::vector<Clustering> baseClustering(baseClusterers.size(),
			Clustering(0)); // collection of base clusterings - fill with empty clustering

	// DEBUG
	GraphIO graphio;
	// DEBUG

	bool repeat;
	int h = -1; // iteration counter, starts with 0 in loop

	graph.push_back(G); // store G^{0}
	Clustering empty(0);
	clusteringBack.push_back(empty); // push a dummy clustering so that clusteringBack[i] contains the clustering projected from G^{i} to G^{0}  (there is no clusteringBack[0])

	Clustering best(0); // best clustering seen so far
	double bestQuality = -2.0; // quality of best clustering
	int bestLevel = -1; // level of best clustering
	Modularity modularity;

	int iterationsWithoutImprovement = 0;
	const int maxIterationsWithoutImprovement = 3;

	do {
		h += 1; 	// increment iteration/hierarchy counter
		INFO("EnsembleClusterer *** ITERATION " << h << " ***");

		// *** base clusterers calculate base clusterings ***
#pragma omp parallel for
		for (int b = 0; b < baseClusterers.size(); b += 1) {
			try {
				baseClustering.at(b) = baseClusterers.at(b)->run(graph.at(h));
				// DEBUG
				DEBUG("created base clustering: k=" << baseClustering.at(b).numberOfClusters());
				if (baseClustering.at(b).isOneClustering(graph.at(h))) {
					WARN("base clusterer created 1-clustering");
				}
				// DEBUG
			} catch (...) {
				ERROR("base clusterer failed with exception.");
				throw std::runtime_error("base clusterer failed.");
			}
		}

		for (int b = 0; b < baseClusterers.size(); b += 1) {
			double qual = modularity.getQuality(baseClustering[b], graph.at(h));
			if (qual > bestQuality) {
				bestQuality = qual;
				best = baseClustering[b];
				bestLevel = h;
			}
		}

		// ANALYSIS
		if (calcBaseClusteringDissimilarity) {
			JaccardMeasure dm;
			for (int b = 0; b < baseClustering.size(); b += 1) {
				for (int c = b; c < baseClustering.size(); c += 1) {
					double d = dm.getDissimilarity(graph.at(h), baseClustering.at(b), baseClustering.at(c));
					INFO("dm " << b << " <-> " << c << ": " << d);
				}
			}
		}
		//

		// *** overlap clusters to create core clustering ***
		INFO("[BEGIN] finding core clustering");
		clustering.push_back(this->overlap->run(graph.at(h), baseClustering));
		DEBUG("created core clustering: k=" << clustering.at(h).numberOfClusters());
		INFO("[DONE] finding core clustering");

		if (h == 0) {			// first iteration
			// *** calculate quality of first core clustering with respect to first graph ***
			quality.push_back(this->qm->getQuality(clustering.at(h), graph.at(h)));
			DEBUG("pushed quality: " << quality.back());

			if (quality.back() > bestQuality) {
				bestQuality = quality.back();
				best = clustering.back();
				bestLevel = h;
			}

			// *** contract the graph according to core clustering **
			INFO("[BEGIN] contracting graph");
			auto con = contract.run(graph.at(h), clustering.at(h));	// returns pair (G^{i+1}, M^{i->i+1})
			graph.push_back(con.first);		// store G^{i+1}
			map.push_back(con.second);		// store M^{i->i+1}

			//DEBUG
			DEBUG("contracted graph G^" << (h+1) << " created: " << graph.back().toString());
			//DEBUG

			// new graph created => repeat
			repeat = true;
		} else { 	// other iterations
			clusteringBack.push_back(project.projectBackToFinest(clustering.at(h), map, G));
			assert (clustering.at(h).numberOfClusters() == clusteringBack.at(h).numberOfClusters());
			// DEBUG
			DEBUG("created projected clustering: k=" << clusteringBack.at(h).numberOfClusters());
			// DEBUG

			quality.push_back(this->qm->getQuality(clusteringBack.at(h), graph.at(h)));
			DEBUG("pushed quality: " << quality.back());

			if (quality.back() > bestQuality) {
				bestQuality = quality.back();
				best = clustering.at(h);
				bestLevel = h;
			}

			// *** test if new core clustering is better than previous one **
			if (quality.at(h) > quality.at(h-1)) {
				DEBUG("quality[" << h << "] = " << quality.at(h) << " > quality[" << (h-1) << "] = " << quality.at(h-1));

				// better quality => repeat
				repeat = true;
				iterationsWithoutImprovement = 0;
			} else {
				DEBUG("quality[" << h << "] = " << quality.at(h) << " <= quality[" << (h-1) << "] = " << quality.at(h-1));

				++iterationsWithoutImprovement;
				if (iterationsWithoutImprovement
						>= maxIterationsWithoutImprovement) {
					// new core clustering is not better for some time => do not contract according to new core clustering and do not repeat
					repeat = false;
				}
			}

			if (repeat) {
				INFO("[BEGIN] contracting graph");
				auto con = contract.run(graph.at(h), clustering.at(h)); // returns pair (G^{i+1}, M^{i->i+1})
				graph.push_back(con.first); // store G^{i+1}
				map.push_back(con.second); // store M^{i->i+1}
				INFO("[DONE] contracting graph");
			}
		}

	} while (repeat);

	Clustering zetaFine(0);
	if (bestLevel > 0) {
		Clustering zetaCoarse = this->finalClusterer->run(graph.at(h)); // TODO: check: index correct?
//	Clustering zetaFine = project.projectBackToFinest(zetaCoarse, map, G);
		zetaFine = project.projectBackToFinest(best, map, graph.at(bestLevel));
	} else {
		zetaFine = best;
	}

	// FIXME: Rueckprojektion scheint noch nicht zu stimmen, denn wenn man von Level 0 auf Level 0
	// zurueckwirft, kommt ein anderes Clustering raus

	INFO("Best clustering was found on level " << bestLevel << ", quality: " << modularity.getQuality(zetaFine, G) << " vs " << bestQuality << " vs " << modularity.getQuality(best, graph.at(bestLevel)));

	return zetaFine;
}

std::string EnsembleClusterer::toString() const {
	std::stringstream strm;
	strm << "EnsembleClusterer(" << "base=" << this->baseClusterers.front()->toString() << ",ensemble=" << this->baseClusterers.size() << ",final=" << this->finalClusterer->toString() << ")(h=" << this->h << ")";
	return strm.str();
}

} /* namespace EnsembleClustering */
