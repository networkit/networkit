/*
 *
 */

#include "LFRGenerator.h"
#include "PowerlawDegreeSequence.h"
#include "EdgeSwitchingMarkovChainGenerator.h"
#include "PubWebGenerator.h"
#include "../auxiliary/Random.h"
#include <algorithm>
#include <random>

NetworKit::LFRGenerator::LFRGenerator(NetworKit::count n, NetworKit::count avgDegree, NetworKit::count maxDegree, double mu, double nodeDegreeExp, double communitySizeExp, NetworKit::count minCommunitySize, NetworKit::count maxCommunitySize) :
n(n), avgDegree(avgDegree), maxDegree(maxDegree), mu(mu), nodeDegreeExp(nodeDegreeExp), communitySizeExp(communitySizeExp), minCommunitySize(minCommunitySize), maxCommunitySize(maxCommunitySize), hasGraph(false), hasPartition(false) {
	if (std::ceil(mu * maxDegree) >= maxCommunitySize) {
		throw std::runtime_error("Graph not realizable, the maximum internal degree is greater than the largest possible internal degree.");
	}

	if (maxDegree >= n) {
		throw std::runtime_error("Graph not realizable, the maximum degree must not be larger or equal to n");
	}
}

void NetworKit::LFRGenerator::run() {
	std::vector<count> externalDegree;

	{
		PowerlawDegreeSequence nodeDegreeSequence(1, maxDegree, nodeDegreeExp);
		nodeDegreeSequence.setMinimumFromAverageDegree(avgDegree);
		if (std::ceil(nodeDegreeSequence.getMinimumDegree() * mu) >= minCommunitySize) {
			throw std::runtime_error("Graph not realizable, the needed minimum degree in order to realize the node degree sequence is so high that no nodes can be assigned in communities of minimum size");
		}
		nodeDegreeSequence.run();
		externalDegree = nodeDegreeSequence.getDegreeSequence(n);
	}

	std::vector<count> internalDegree(n);

	#pragma omp parallel for
	for (node u = 0; u < n; ++u) {
		if (externalDegree[u] == 0) continue;

		double intDeg = (1.0 - mu) * externalDegree[u];
		if (intDeg < 1) { // assure minimum internal degree of 1
			internalDegree[u] = 1;
		} else if (Aux::Random::probability() >= std::remainder(intDeg, 1)) {
			internalDegree[u] = count(intDeg);
		} else {
			internalDegree[u] = std::ceil(intDeg);
		}
		externalDegree[u] -= internalDegree[u];
	}

	count actualMaxInternalDegree = *std::max_element(internalDegree.begin(), internalDegree.end());
	count actualMaxCommunitySize = 0;

	std::vector<count> communitySizes;

	{ // generate community sizes according to the community size power law distribution
		PowerlawDegreeSequence communityDegreeSequenceGen(minCommunitySize, maxCommunitySize, communitySizeExp);
		communityDegreeSequenceGen.run();

		count sumCommunitySizes = 0;

		while (true) {
			count newSize = communityDegreeSequenceGen.getDegree();

			if (actualMaxCommunitySize <= actualMaxInternalDegree && sumCommunitySizes + actualMaxInternalDegree + 1 + newSize > n) {
				// this is the first time this happened,  sumCommunitySizes + actualMaxInternalDegree + 1 <= n
				// however there is also no chance that without intervention the graph will be realizable
				// so modify the degree sequence such that the graph is possibly realizable
				newSize = actualMaxInternalDegree + 1;
			}

			if (sumCommunitySizes + newSize <= n) {
				communitySizes.push_back(newSize);
				sumCommunitySizes += newSize;
				actualMaxCommunitySize = std::max(newSize, actualMaxCommunitySize);
			} else { // if the new community doesn't fit anymore, increase the smallest community to fill the gap and exit the loop
				*std::min_element(communitySizes.begin(), communitySizes.end()) += n - sumCommunitySizes;

				break;
			}
		}
	}

	// a list of nodes for each community
	std::vector<std::vector<node> > communityNodeList;

	bool assignmentSucceeded = true;

	do {
		assignmentSucceeded = true;
		communityNodeList.clear();
		communityNodeList.resize(communitySizes.size());

		std::vector<count> communitySelection;
		{ // generate a random permutation of size(c) copies of each community id c
			communitySelection.reserve(n);

			for (index i = 0; i < communitySizes.size(); ++i) {
				for (node u = 0; u < communitySizes[i]; ++u) {
					communitySelection.push_back(i);
				}
			}

			std::shuffle(communitySelection.begin(), communitySelection.end(), std::default_random_engine());
		}

		// how many nodes are still missing in the community
		std::vector<count> remainingCommunitySizes(communitySizes);

		// nodes that still need to be assigned to a community
		std::vector<node> nodesToAssign;

		// in the first round, assign each node to a randomly chosen community if possible
		for (node u = 0; u < n; ++u) {
			if (communitySizes[communitySelection[u]] > internalDegree[u]) {
				communityNodeList[communitySelection[u]].push_back(u);
				--remainingCommunitySizes[communitySelection[u]];
			} else {
				nodesToAssign.push_back(u);
			}
		}

		// now assign a random unassigned node to a random feasible community (chosen at random weighted by community size)
		// if the community is full, remove a node from it
		count attempts = 0;
		count totalAttempts = 0;
		while (!nodesToAssign.empty()) {
			++totalAttempts;
			index c = Aux::Random::choice(communitySelection);

			index i = Aux::Random::index(nodesToAssign.size());
			node u = nodesToAssign[i];

			{ // remove index u from nodesToAssign
				nodesToAssign[i] = nodesToAssign.back();
				nodesToAssign.pop_back();
			}

			// pick another community till the node's internal degree fits into the community
			while (internalDegree[u] >= communitySizes[c]) {
				// we ensured above that there is always a community in which every node will fit
				c = Aux::Random::choice(communitySelection);
			}

			communityNodeList[c].push_back(u);

			// we are lucky, the community is not full yet
			if (remainingCommunitySizes[c] > 0) {
				--remainingCommunitySizes[c];
				attempts = 0;
			} else { // remove a random node from c
				index r = Aux::Random::index(communityNodeList[c].size());
				nodesToAssign.push_back(communityNodeList[c][r]);
				communityNodeList[c][r] = communityNodeList[c].back();
				communityNodeList[c].pop_back();
				++attempts;
			}

			if (attempts > 3*n) {
				// merge two smallest communities
				WARN("Needed too long to assign nodes to communities. Either there are too many high-degree nodes or the communities are simply too small.");
				WARN("Changing the community sizes by merging the two smallest communities");
				DEBUG(attempts, " nodes were assigned to full communities, in total, ", totalAttempts, " were made");
				auto minIt = std::min_element(communitySizes.begin(), communitySizes.end());
				count minVal = *minIt;

				// delete minimum element
				*minIt = communitySizes.back();
				communitySizes.pop_back();

				minIt = std::min_element(communitySizes.begin(), communitySizes.end());
				*minIt += minVal;

				assignmentSucceeded = false;

				break;
			}
		}
	} while (!assignmentSucceeded);


	// write communityNodeList into partition instance
	zeta = Partition(n);
	zeta.setUpperBound(communityNodeList.size());

	#pragma omp parallel for
	for (index i = 0; i < communityNodeList.size(); ++i) {
		for (node u : communityNodeList[i]) {
			zeta[u] = i;
		}
	}

	// initialize the result graph
	G = Graph(n);

	// generate intra-cluster edges
	for (auto & communityNodes : communityNodeList) {
		std::vector<count> intraDeg;
		intraDeg.reserve(communityNodes.size());

		for (node u : communityNodes) {
			intraDeg.push_back(internalDegree[u]);
		}

		// check if sum of degrees is even and fix if necessary
		count intraDegSum = std::accumulate(intraDeg.begin(), intraDeg.end(), 0);

		while (intraDegSum % 2 != 0) {
			index i = Aux::Random::index(communityNodes.size());
			node u = communityNodes[i];
			if (Aux::Random::real() >= 0.5) {
				if (intraDeg[i] < communityNodes.size() - 2 && externalDegree[u] > 0) {
					TRACE("Making degree distribution even by increasing intra-degree of node ", u);
					++intraDeg[i];
					++internalDegree[u];
					++intraDegSum;
					--externalDegree[u];
				}
			} else {
				if (intraDeg[i] > 1 && externalDegree[u] < n - 2) {
					TRACE("Making degree distribution even by decreasing intra-degree of node ", u);
					--intraDeg[i];
					--internalDegree[u];
					--intraDegSum;
					++externalDegree[u];
				}
			}
		}

		EdgeSwitchingMarkovChainGenerator intraGen(intraDeg, true);
		/* even though the sum is even the degree distribution isn't necessarily realizable.
		Disabling the check means that some edges might not be created because of this but at least we will get a graph. */
		Graph intraG = intraGen.generate();

		intraG.forEdges([&](node u, node v) {
			G.addEdge(communityNodes[u], communityNodes[v]);
		});
	}

	// generate inter-cluster edges
	EdgeSwitchingMarkovChainGenerator graphGen(externalDegree, true);

	Graph interG = graphGen.generate();

	// rewire intra-cluster edges in interG as only inter-cluster edges should be generated
	std::vector<std::pair<node, node> > edgesToRemove;

	interG.forEdges([&](node u, node v) {
		if (zeta[u] == zeta[v]) {
			edgesToRemove.emplace_back(u, v);
		}
	});

	if (!edgesToRemove.empty()) {
		std::vector<node> nodeSelection;
		nodeSelection.reserve(interG.numberOfEdges() * 2);

		interG.forNodes([&](node u) {
			for (count i = 0; i < interG.degree(u); ++i) {
				nodeSelection.push_back(u);
			}
		});


		count intraRemovalAttempts = 0;
		count maxIntraRemovelAttempts = n * 10;
		while (!edgesToRemove.empty()) { // FIXME simply drop these edges at some point if rewiring is not possible
			index i = Aux::Random::index(edgesToRemove.size());
			node s1, t1;
			std::tie(s1, t1) = edgesToRemove[i];

			if (!interG.hasEdge(s1, t1)) {
				edgesToRemove[i] = edgesToRemove.back();
				edgesToRemove.pop_back();
				continue;
			}

			++intraRemovalAttempts;

			node s2 = Aux::Random::choice(nodeSelection);

			if (s2 == s1 || s2 == t1) continue;

			node t2 = interG.randomNeighbor(s2);

			if (t2 == none) continue;

			if (t1 == t2 || s1 == t2) continue;

			if (interG.hasEdge(s1, t2) || interG.hasEdge(s2, t1)) continue; // FIXME make efficient

			interG.swapEdge(s1, t1, s2, t2);

			edgesToRemove[i] = edgesToRemove.back();
			edgesToRemove.pop_back();

			if (zeta[s1] == zeta[t2]) {
				edgesToRemove.emplace_back(s1, t2);
			}

			if (zeta[s2] == zeta[t1]) {
				edgesToRemove.emplace_back(s2, t1);
			}

			// if t1-t2 is in edgesToRemove, it will be removed in a later iteration of the loop

			if (intraRemovalAttempts > maxIntraRemovelAttempts) {
				break;
			}
		}

		for (auto it = edgesToRemove.begin(); it != edgesToRemove.end(); ++it) {
			if (!interG.hasEdge(it->first, it->second)) {
				*it = edgesToRemove.back();
				edgesToRemove.pop_back();
				if (it == edgesToRemove.end()) break;
			}
		}

		if (!edgesToRemove.empty()) {
			WARN("There are ", edgesToRemove.size(), " intra-cluster edges (which might be counted twice) that should actually be rewired to be inter-cluster edges but couldn't be rewired after ",
				maxIntraRemovelAttempts, " attempts. They will be simply dropped now.");

			for (auto e : edgesToRemove) {
				if (interG.hasEdge(e.first, e.second)) { // these edges might be duplicates/in both directions
					interG.removeEdge(e.first, e.second);
				}
			}
		}
	}

	count edgesRemoved = 0;
	interG.forEdges([&](node u, node v) {
#ifndef NDEBUG // check if edge rewiring actually worked.
		if (!G.hasEdge(u, v)) {
#endif
			G.addEdge(u, v);
#ifndef NDEBUG
		} else {
			++edgesRemoved;
		}
#endif
	});

	if (edgesRemoved > 0) {
		ERROR("Removed ", edgesRemoved, " intra-community edges that already existed but should have been rewired actually!");
	}

	hasGraph = true;
	hasPartition = true;
}

NetworKit::Graph NetworKit::LFRGenerator::getGraph() const {
	if (!hasGraph) throw std::runtime_error("Run must be called first");

	return G;
}

NetworKit::Graph && NetworKit::LFRGenerator::getMoveGraph() {
	if (!hasGraph) throw std::runtime_error("Run must be called first");

	hasGraph = false;

	return std::move(G);
}

NetworKit::Partition NetworKit::LFRGenerator::getPartition() const {
	if (!hasPartition) throw std::runtime_error("Run must be called first");

	return zeta;
}

NetworKit::Partition && NetworKit::LFRGenerator::getMovePartition() {
	if (!hasPartition) throw std::runtime_error("Run must be called first");

	return std::move(zeta);
}


