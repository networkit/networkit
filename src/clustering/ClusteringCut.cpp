/*
 * ClusteringCut.cpp
 *
 *  Created on: 23.10.2013
 *      Author: Weijian Ji (weijian.ji@student.kit.edu)
 */

#include "ClusteringCut.h"

namespace NetworKit {


ClusteringCut::ClusteringCut(){
	//TODO: creator with special Parameters like ClusteringCut(const Clustering& zeta, const Graph& G)
}
ClusteringCut::~ClusteringCut(){

}

cutEdges ClusteringCut::getAllCutEdges(const Graph& G, const Clustering& zeta){
	cutEdges cuttingEdges;
	G.forEdges([&](node u, node v) {
		if (zeta.clusterOf(u) != zeta.clusterOf(v)) {
			cuttingEdges.insert(std::pair<node, node>(u, v));
		}
	});
	return cuttingEdges;
}

cuttingMap ClusteringCut::getClusterToCutMatrixMap(const Graph& G, const Clustering& zeta){

	cuttingMap cMap; // map to return

	cluster clst1, clst2; // cluster indexes for std::unordered_map key
	cutEdges cuttingEdges; // cutting edges, group of edges

	G.forEdges([&](node u, node v) {
		clst1 = zeta.clusterOf(u);
		clst2 = zeta.clusterOf(v);
		std::pair<cluster, cluster> key (clst1, clst2);
		if (zeta.clusterOf(u) != zeta.clusterOf(v)) {

			cuttingMap::const_iterator got = cMap.find(key);
			 if ( got != cMap.end() ){
				 std::cout << "yes, found";
				 ;//TODO insert to mapped value;
			 }else{
				 ;//TODO add new key to map
			 }
		}
	});
	return cMap;
}

}/*namenspace NetworKit*/
