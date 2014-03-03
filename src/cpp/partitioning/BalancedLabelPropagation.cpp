/*
 * BalancedLabelPropagation.cpp
 *
 *  Created on: Jun 19, 2013
 *      Author: Henning
 */

#include "../auxiliary/Random.h"

#include "BalancedLabelPropagation.h"
#include "../community/GraphClusteringTools.h"

namespace NetworKit {

BalancedLabelPropagation::BalancedLabelPropagation(double myexponent): exponent(myexponent) {

}

BalancedLabelPropagation::~BalancedLabelPropagation() {

}

Partition BalancedLabelPropagation::run(Graph& graph, count numParts) {
	ClusteringGenerator gen;
	// FIXME: change to region growing from random
	Partition partition = gen.makeContinuousBalancedClustering(graph, numParts);
	partition = rerun(graph, numParts, partition);
	return partition;
}


Partition& BalancedLabelPropagation::rerun(Graph& graph, count numParts, Partition& partition) {
	float avg = ceil((float) graph.numberOfNodes() / (float) numParts);
	std::vector<float> adjustByFactor(numParts);
	EdgeCut edgeCut;


	auto isMoveAccepted([&](edgeweight gain, count t) {
		if (gain >= 0) {
			return true;
		}
		else {
			double prob = Aux::Random::real();
			return (prob <= exp((double) gain / (double) t));
		}
	});

	count numIters = 35;
	if (exponent >= 4.0) {
		numIters = 3;
	}
	DEBUG("cut/balance before loop: " , edgeCut.getQuality(partition, graph) , ", " , GraphClusteringTools::getImbalance(partition));

	for (index i = 0; i < numIters; ++i) { // FIXME: different termination criterion
		// read cluster sizes and compute scale values
		TRACE("compute cluster sizes... ");
		std::vector<count> clusterSizes = partition.subsetSizes();
		TRACE("done");
		for (index p = 0; p < numParts; ++p) {
			TRACE("size of cluster " , p , ": " , clusterSizes[p]);
			adjustByFactor[p] = avg / clusterSizes[p];
		}


		TRACE("start loop over nodes... ");
		// perform LP
		graph.forNodes([&](node v) {
#if 0
			node neighbor = graph.randomNeighbor(v);

			if (neighbor != none) {
				index vBlock = partition[v];

				// *** compute gain of relabeling according to neighbor
				index neighBlock = partition[neighbor];

				// what is the weighted degree with neighCluster?
				edgeweight wdegNeigh = GraphClusteringTools::weightedDegreeWithCluster(graph, partition, v, neighBlock);

				// what is the weighted degree to the current cluster?
				edgeweight wdegCurrent = GraphClusteringTools::weightedDegreeWithCluster(graph, partition, v, vBlock);

				edgeweight gain = adjustByFactor[neighBlock] * wdegNeigh - adjustByFactor[vBlock] * wdegCurrent;
				// ***


				// accept label based on simulated annealing acceptance
				if (isMoveAccepted(gain, i, probGen)) {
					partition.moveToSubset(neighBlock, v);
				}
			}
#else

			std::map<index, double> labelWeights; // neighborLabelCounts maps label -> frequency in the neighbors
			index heaviest = none;

			// weigh the labels in the neighborhood of v
			graph.forWeightedNeighborsOf(v, [&](node w, edgeweight weight) {
				index lw = partition[w];
				labelWeights[lw] += weight * adjustByFactor[lw];
				heaviest = lw; // init heaviest by some neighbor's cluster
			});

//			heaviest = std::max_element(labelWeights.begin(),
//							labelWeights.end(),
//							[](const std::pair<index, edgeweight>& p1, const std::pair<index, edgeweight>& p2) {
//								return p1.second < p2.second;})->first;

			if (labelWeights.size() > 1) {

				std::vector<std::pair<index, double> > neighborhood;
				for (std::map<index, double>::iterator iter = labelWeights.begin();
						iter != labelWeights.end(); ++iter) {
					neighborhood.push_back(std::make_pair(iter->first, iter->second));
				}

				double sum = 0.0;
				for (std::vector<std::pair<index, double>>::iterator iter = neighborhood.begin();
						iter != neighborhood.end(); ++iter) {
					sum += pow(iter->second, exponent); // favor large int-/ext-degrees
				}

				// normalize to sum 1 and prefix sum
				for (std::vector<std::pair<index, double>>::iterator iter = neighborhood.begin();
						iter != neighborhood.end(); ++iter) {
					iter->second = pow(iter->second, exponent) / sum;
				}

				double prefixSum = 0.0;
				double temp = 0.0;
				for (std::vector<std::pair<index, double>>::iterator iter = neighborhood.begin();
						iter != neighborhood.end(); ++iter) {
					temp = iter->second;
					iter->second += prefixSum;
					prefixSum += temp;
				}

				double prob = Aux::Random::real();

				for (std::vector<std::pair<index, double>>::iterator iter = neighborhood.begin();
						iter != neighborhood.end(); ++iter) {
					if (iter->second >= prob) {
						heaviest = iter->first;
						break;
					}
				}
			}

			if (partition[v] != heaviest) { // UPDATE
				partition[v] = heaviest;
			}
#endif


//			DEBUG("partition[" , v , "]: " , partition[v]);
		});

		DEBUG("cut/balance in iter " , i , ": " , edgeCut.getQuality(partition, graph) , ", " , GraphClusteringTools::getImbalance(partition));
	}

	std::vector<count> clusterSizes = partition.subsetSizes();
	for (index p = 0; p < clusterSizes.size(); ++p) {
		DEBUG("after: size of cluster " , p , ": " , clusterSizes[p]);
	}

	return partition;
}

Partition& BalancedLabelPropagation::postsmooth(Graph& graph, count numBlocks,
		Partition& partition) {
	double tmpExp = this->exponent;
	this->setExponent(8.0); // quasi-deterministic
	partition = this->rerun(graph, numBlocks, partition);
	this->exponent = tmpExp;
	return partition;
}

} /* namespace NetworKit */
