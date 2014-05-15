/*
 * DotPartitionWriter.cpp
 */

#include "DotPartitionWriter.h"

#include <atomic>
#include <cstdint>
#include <fstream>
#include <unordered_set>

namespace NetworKit {

std::map<index, double> DotPartitionWriter::createHueMap(Graph &graph, Partition& zeta) const {
	std::unordered_set<index> uniqueIds;

	graph.forNodes([&](node u){
		if (graph.degree(u) == 0) {
			return;
		}
		index c = zeta.subsetOf(u);
		uniqueIds.insert(c);
	});

	std::map<index, double> clusterHueMap;

	std::size_t idx = 0;
	double factor = 1.0 / uniqueIds.size();
	for(const auto& value: uniqueIds) {
		clusterHueMap.insert(std::make_pair(value, factor * idx));
		++idx;
	}

	return clusterHueMap;
}

void DotPartitionWriter::write(Graph& graph, Partition& zeta, const std::string& path) const {
	std::ofstream file{path};
	
	auto hueMap = this->createHueMap(graph, zeta);
	
	std::size_t idx = 0;
	double factor = 1.0 / countConnectedNodes(graph);
	
	file << "graph {\n";
	
	graph.forNodes([&](node u){
		if (graph.degree(u) == 0) {
			return;
		}
		index c = zeta.subsetOf(u);
		double h = hueMap.at(c);
		file << u << " [style=filled, color=\"" << h << ",0.99,0.99\", label=" << c << "];\n";
	});
	graph.forEdges([&](node u, node v){
		file << u << " -- " << v << ";\n";
	});

	file << "}\n";
}

std::size_t DotPartitionWriter::countConnectedNodes(const Graph& graph) const {
	std::atomic_size_t count{0};
	graph.parallelForNodes([&](node u){
		if (graph.degree(u) != 0) {
			++count;
		}
	});
	return count;
}

} /* namespace NetworKit */
