/*
 * SelectiveDissimilarityMeasure.cpp
 *
 *  Created on: 18.06.2013
 *      Author: cls
 */

#include "SelectiveDissimilarityMeasure.h"

namespace NetworKit {

SelectiveDissimilarityMeasure::SelectiveDissimilarityMeasure() {
	// TODO Auto-generated constructor stub

}

SelectiveDissimilarityMeasure::~SelectiveDissimilarityMeasure() {
	// TODO Auto-generated destructor stub
}


JaccardIndex::JaccardIndex() {
}

JaccardIndex::~JaccardIndex() {
}

double JaccardIndex::localDissimilarity(const node seedNode,
		const std::unordered_set<node>& community,
		const Clustering& groundTruth) {

	int intersection = 0;
	int clusterSize = 0;
	cluster cluster = groundTruth.clusterOf(seedNode);

	if (cluster == -1) {
		return 0;
	}
	if (community.find(seedNode) == community.end()) {
		return 0;
	}
	for(node u = 0; u < groundTruth.numberOfEntries() ; ++u) {
		if(groundTruth.clusterOf(u) == cluster) {
			clusterSize++;
		}
	}
	for(node u : community) {
		if (groundTruth.clusterOf(u) == cluster) {
			intersection++;
		}
	}
	return ((double)(intersection)) / ((double)(community.size() + clusterSize - intersection));
}

Precision::Precision() {
}

Precision::~Precision() {
}

double Precision::localDissimilarity(const node seedNode,
		const std::unordered_set<node>& community,
		const Clustering& groundTruth) {

	int counter = 0;
	cluster cluster = groundTruth.clusterOf(seedNode);
	if (community.find(seedNode) == community.end()) {
		return 0;
	}
	if (cluster == -1) {
		return 0;
	}
	for(node u : community) {
		if (groundTruth.clusterOf(u) == cluster) {
			counter++;
		}
	}
	return ((double)counter) / ((double)community.size());
}

Recall::Recall() {
}

Recall::~Recall() {
}

double Recall::localDissimilarity(const node seedNode,
		const std::unordered_set<node>& community,
		const Clustering& groundTruth) {

	int counter = 0;
	int clusterSize = 0;
	cluster cluster = groundTruth.clusterOf(seedNode);

	if (community.find(seedNode) == community.end()) {
		return 0;
	}
	if (cluster == -1) {
		return 0;
	}
	for(node u : community) {
		if (groundTruth.clusterOf(u) == cluster) {
			counter++;
		}
	}
	for(node u = 0; u < groundTruth.numberOfEntries(); ++u) {
		if(groundTruth.clusterOf(u) == cluster) {
			clusterSize++;
		}
	}
	return ((double)counter) / ((double)clusterSize);
}

NMI::NMI() {
}

NMI::~NMI() {
}

double NMI::localDissimilarity(const node seedNode,
		const std::unordered_set<node>& community,
		const Clustering& groundTruth) {

	double counter = 0;
	double clusterSize = 0;
	cluster cluster = groundTruth.clusterOf(seedNode);

	if (community.find(seedNode) == community.end()) {
		return -1;
	}
	if (cluster == -1) {
		return -1;
	}
	for(node u : community) {
		if (groundTruth.clusterOf(u) == cluster) {
			counter++;
		}
	}
	for(node u = 0; u < groundTruth.numberOfEntries(); ++u) {
		if(groundTruth.clusterOf(u) == cluster) {
			clusterSize++;
		}
	}

return 0;
}

double JaccardIndex::getDissimilarity(
		const std::unordered_map<node, std::unordered_set<node> > communities,
		const Clustering& groundTruth) {

	int intersection = 0;
	int unionSize = 0;
	cluster cluster;
	std::unordered_set<count> relevantGroundTruthClusters;
	for (auto u : communities) {
		cluster = groundTruth.clusterOf(u.first);
		if (cluster != -1) {
			relevantGroundTruthClusters.insert(cluster);
		}
	}
	for (auto u : communities) {
		cluster = groundTruth.clusterOf(u.first);
		for(node v : u.second) {
			node tmp = groundTruth.clusterOf(v);
			if (tmp == cluster) {
				intersection++;
			} else if (tmp == -1) {
				unionSize++;
			} else if (relevantGroundTruthClusters.find(tmp) == relevantGroundTruthClusters.end()) {
				unionSize++;
			}
		}
	}

	for (node v = 0; v < groundTruth.numberOfEntries(); v += 1) {
		cluster = groundTruth.clusterOf(v) ;
		if (relevantGroundTruthClusters.find(cluster) != relevantGroundTruthClusters.end()) {
			unionSize++;
		}
	}
return ((double)intersection) / ((double)unionSize);
}


double Precision::getDissimilarity(
		const std::unordered_map<node, std::unordered_set<node> > communities,
		const Clustering& groundTruth) {

	int intersection = 0;
	int counter = 0;
	cluster cluster;
	std::unordered_set<count> relevantGroundTruthClusters;
	for (auto u : communities) {
		cluster = groundTruth.clusterOf(u.first);
		if (cluster != -1) {
			relevantGroundTruthClusters.insert(cluster);
		}
	}
	for (auto u : communities) {
		cluster = groundTruth.clusterOf(u.first);
		for(node v : u.second) {
			node tmp = groundTruth.clusterOf(v);
			if (tmp == cluster) {
				intersection++;
			}
		}
	}

	for (auto u : communities) {
		counter = counter + u.second.size();
	}
return ((double)intersection) / ((double)counter);
}

double Recall::getDissimilarity(
		const std::unordered_map<node, std::unordered_set<node> > communities,
		const Clustering& groundTruth) {

	int intersection = 0;
	int counter = 0;
	cluster cluster;
	std::unordered_set<count> relevantGroundTruthClusters;
	for (auto u : communities) {
		cluster = groundTruth.clusterOf(u.first);
		if (cluster != -1) {
			relevantGroundTruthClusters.insert(cluster);
		}
	}
	for (auto u : communities) {
		cluster = groundTruth.clusterOf(u.first);
		for(node v : u.second) {
			node tmp = groundTruth.clusterOf(v);
			if (tmp == cluster) {
				intersection++;
			}
		}
	}

	for (node v = 0; v < groundTruth.numberOfEntries(); v += 1) {
		cluster = groundTruth.clusterOf(v) ;
		if (relevantGroundTruthClusters.find(cluster) != relevantGroundTruthClusters.end()) {
			counter++;
		}
	}
	return ((double)intersection) / ((double)counter);

}

double NMI::getDissimilarity(
		const std::unordered_map<node, std::unordered_set<node> > communities,
		const Clustering& groundTruth) {
return 0;
}

} /* namespace NetworKit */
