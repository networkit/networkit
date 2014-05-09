/*
 * DotPartitionWriter.cpp
 */

#include "DotPartitionWriter.h"

#include <map>
#include <set>

namespace NetworKit {

DotPartitionWriter::DotPartitionWriter() {
	// TODO Auto-generated constructor stub

}

DotPartitionWriter::~DotPartitionWriter() {
	// TODO Auto-generated destructor stub
}

std::map<index, double> DotPartitionWriter::createHueMap(Graph &graph, Partition& zeta) const {
    std::set<index> uniqueIds;

    graph.forNodes([&](node u){
        if (graph.degree(u) == 0) {
            return;
        }
        index c = zeta.subsetOf(u);
        uniqueIds.insert(c);
	});

    std::map<index, double> clusterHueMap;

    std::set<index>::iterator iter;
    int idx = 0;
    for(iter = uniqueIds.begin(); iter != uniqueIds.end(); ++iter) {
        clusterHueMap.insert(std::pair<index, double>(*iter, (1 / (double)uniqueIds.size()) * idx));
        idx++;
    }

    return clusterHueMap;
}

void DotPartitionWriter::write(Graph& graph, Partition& zeta, std::string path) const {
	std::ofstream file;
	file.open(path.c_str());

    std::map<index, double> hueMap = this->createHueMap(graph, zeta);

    file << "graph {" << std::endl;

    graph.forNodes([&](node u){
        if (graph.degree(u) == 0) {
            return;
        }
        index c = zeta.subsetOf(u);
        double h = hueMap.at(c);
        file << u << " [style=filled, color=\"" << h << ",0.99,0.99\", label=" << c << "];" << std::endl;
	});
    graph.forEdges([&](node u, node v){
        file << u << " -- " << v << ";" << std::endl;
    });

    file << "}" << std::endl;

	file.close();
}

} /* namespace NetworKit */
