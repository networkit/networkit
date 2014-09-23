/*
 * Author: Michael Hamann <michael.hamann@kit.edu>
 */

#include "CutClustering.h"
#include "../flow/EdmondsKarp.h"
#include "../properties/ConnectedComponents.h"
#include "../auxiliary/Log.h"

#include <sstream>
#include <stdexcept>
#include <limits>

NetworKit::CutClustering::CutClustering(NetworKit::edgeweight alpha) : alpha(alpha) { }

NetworKit::Partition NetworKit::CutClustering::run(const NetworKit::Graph &G) {
	Partition result(G.upperNodeIdBound());
	result.setUpperBound(G.upperNodeIdBound());
	
	Graph graph(G, true, false);
	
	node t = graph.addNode();
	
	graph.forNodes([&](node u) {
		if (u != t) {
			graph.addEdge(u, t, alpha);
		}
	});
	
	graph.indexEdges();
	
	EdmondsKarp flowAlgo;
	
	graph.forNodes([&](node u) {
		if (u != t && !result.contains(u)) {
			std::vector<node> sourceSet;
			flowAlgo.run(graph, u, t, sourceSet);
			
			for (node v : sourceSet) {
				result[v] = u;
			}
		}
	});
	
	return result;
}

std::map< double, NetworKit::Partition > NetworKit::CutClustering::getClusterHierarchy(const NetworKit::Graph &G) {
	std::map<double, Partition> result;
	double lower = 0, upper = 2;

	if (G.isWeighted()) {
		G.forEdges([&](node u, node v, edgeweight weight) {
			upper = std::max(weight, upper);
		});
		upper += 1;
	}

	ConnectedComponents connComp = ConnectedComponents(G);
	connComp.run();

	Partition lowerClusters(G.upperNodeIdBound());

	std::vector<node> componentReps(connComp.numberOfComponents(), none);

	G.parallelForNodes([&](node u) {
		componentReps[connComp.componentOfNode(u)] = u;
	});

	lowerClusters.setUpperBound(G.upperNodeIdBound());

	G.parallelForNodes([&](node u) {
		lowerClusters[u] = componentReps[connComp.componentOfNode(u)];
	});

	result.insert(std::make_pair(0, lowerClusters));

	if (connComp.numberOfComponents() > 1) {
		node rep = G.randomNode();
		Partition wholeGraph(G.upperNodeIdBound(), rep);
		wholeGraph.setUpperBound(rep + 1);

		result.insert(std::make_pair(-1, std::move(wholeGraph))); // wholeGraph won't be used anymore
	}

	Partition upperClusters(G.upperNodeIdBound());
	upperClusters.allToSingletons();

	clusterHierarchyRecursion(G, lower, std::move(lowerClusters), upper, std::move(upperClusters), result); // moved values won't be used anymore

	return result;
}

void NetworKit::CutClustering::clusterHierarchyRecursion(const NetworKit::Graph &G, double lower, NetworKit::Partition lowerClusters, double upper, NetworKit::Partition upperClusters, std::map< double, Partition > &result) {
	while (true) {
		double middle = -1;

		std::map<index, count> upperSizes, lowerSizes;
		std::map<index, edgeweight> upperCut, lowerCut;

		G.forEdges([&](node u, node v, edgeweight ew) {
			if (upperClusters[u] != upperClusters[v]) {
				upperCut[upperClusters[u]] += ew;
				upperCut[upperClusters[v]] += ew;
			}

			if (lowerClusters[u] != lowerClusters[v]) {
				lowerCut[lowerClusters[u]] += ew;
				lowerCut[lowerClusters[v]] += ew;
			}
		});

		G.forNodes([&](node u) {
			++upperSizes[upperClusters[u]];
			++lowerSizes[lowerClusters[u]];
		});

		for (auto it : lowerSizes) {
			if (it.second == upperSizes[upperClusters[it.first]])
				continue;

			count upperSize = upperSizes[upperClusters[it.first]];
			count lowerSize = it.second;

			edgeweight lowerWeight = lowerCut[it.first];
			edgeweight upperWeight = upperCut[upperClusters[it.first]];

			edgeweight possibleBreakpoint = (upperWeight - lowerWeight) / (lowerSize - upperSize);

			G.forNodes([&](node u) {
				if (lowerClusters[u] == it.first && upperClusters[u] == u && u != it.first) {
					edgeweight tmpBreakPoint = (upperCut[u] - lowerWeight) / (lowerSize - upperSizes[u]);

					if (tmpBreakPoint > possibleBreakpoint) {
						possibleBreakpoint = tmpBreakPoint;
						upperWeight = upperCut[u];
						upperSize = upperSizes[u];
					}
				}
			});

			if (possibleBreakpoint + std::numeric_limits<double>::epsilon() < upper) {
				middle = possibleBreakpoint + std::numeric_limits<double>::epsilon();
				break;
			} else if (possibleBreakpoint > upper) {
				std::stringstream stream;
				stream << "Warning: possible breakpoint " << possibleBreakpoint << " higher than the upper value: " << upper;
				WARN(stream.str());
			}
		}

		if (middle == -1) {
			if (result.count(upper) == 0)
				result.insert(std::make_pair(upper, std::move(upperClusters))); // upperClusters won't be used anymore

			break;
		}

		Partition middleClusters(CutClustering(middle).run(G));
		count numMiddleClusters = middleClusters.numberOfSubsets();

		std::stringstream stream;
		stream << "Calculated clustering for alpha value " << middle;
		INFO(stream.str());


		bool lowerIsMiddle = (lowerSizes.size() == numMiddleClusters), middleIsUpper = (numMiddleClusters == upperSizes.size());

		if (lowerIsMiddle) {
			throw std::logic_error("Error: Lower clusters are middle clusters");
		} else if (middleIsUpper) {
			upper = middle;
			result.insert(std::make_pair(upper, upperClusters));
		} else {
			clusterHierarchyRecursion(G, lower, std::move(lowerClusters) /* lowerClusters is override afterwards */, middle, middleClusters, result);
			lower = middle;
			lowerClusters = std::move(middleClusters); // middleClusters won't be used anymore
		}
	}

	if (result.count(upper) == 0)
		result.insert(std::make_pair(upper, std::move(upperClusters))); // upperClusters won't be used anymore
};


std::string NetworKit::CutClustering::toString() const {
	std::stringstream stream;
	
	stream << "CutClustering(" << alpha << ")";
	return stream.str();
}

