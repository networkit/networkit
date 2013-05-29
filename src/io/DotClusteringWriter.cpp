/*
 * ClusteringWriter.cpp
 */

#include "DotClusteringWriter.h"

#include <map>
#include <set>

namespace NetworKit {

DotClusteringWriter::DotClusteringWriter() {
	// TODO Auto-generated constructor stub

}

DotClusteringWriter::~DotClusteringWriter() {
	// TODO Auto-generated destructor stub
}

std::map<cluster, double> DotClusteringWriter::createHueMap(Graph &graph, Clustering& zeta) const {
    std::set<cluster> uniqueIds;

    graph.forNodes([&](node u){
        if (graph.degree(u) == 0) {
            return;
        }
        cluster c = zeta.clusterOf(u);
        uniqueIds.insert(c);
	});

    std::map<cluster, double> clusterHueMap;

    std::set<cluster>::iterator iter;
    int index = 0;
    for(iter = uniqueIds.begin(); iter != uniqueIds.end(); ++iter) {
        clusterHueMap.insert(std::pair<cluster, double>(*iter, (1 / (double)uniqueIds.size()) * index));
        index++;
    }

    return clusterHueMap;
}

void DotClusteringWriter::write(Graph& graph, Clustering& zeta, std::string path) const {
	std::ofstream file;
	file.open(path.c_str());

    std::map<cluster, double> hueMap = this->createHueMap(graph, zeta);

    file << "graph {" << std::endl;

    graph.forNodes([&](node u){
        if (graph.degree(u) == 0) {
            return;
        }
        cluster c = zeta.clusterOf(u);
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
