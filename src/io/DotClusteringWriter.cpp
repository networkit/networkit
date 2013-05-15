/*
 * ClusteringWriter.cpp
 */

#include "DotClusteringWriter.h"

namespace NetworKit {

DotClusteringWriter::DotClusteringWriter() {
	// TODO Auto-generated constructor stub

}

DotClusteringWriter::~DotClusteringWriter() {
	// TODO Auto-generated destructor stub
}

void DotClusteringWriter::write(Graph& graph, Clustering& zeta, std::string path) const {
	std::ofstream file;
	file.open(path.c_str());

    file << "graph {" << std::endl;

    graph.forNodes([&](node u){
        cluster c = zeta.clusterOf(u);
        double h = (double)c/zeta.numberOfClusters();
        file << u << " [style=filled, color=\"" << h << ",0.99,0.25\", label=" << c << "];" << std::endl;
	});
    graph.forEdges([&](node u, node v){
        file << u << " -- " << v << ";" << std::endl;
    });

    file << "}" << std::endl;

	file.close();
}

} /* namespace NetworKit */
