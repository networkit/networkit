/*
 * ApproximateClusteringCoefficient.cpp
 *
 *  Created on: 16.11.2013
 *      Author: gbrueckner
 */

#include "ApproximateClusteringCoefficient_Brueckner.h"
#include <random>

namespace NetworKit {

ApproximateClusteringCoefficient_Brueckner::ApproximateClusteringCoefficient_Brueckner() {

}

ApproximateClusteringCoefficient_Brueckner::~ApproximateClusteringCoefficient_Brueckner() {

}

double ApproximateClusteringCoefficient_Brueckner::calculate(Graph& G, int k) {

    using namespace std;

    default_random_engine e;
    e.seed(random_device()());
    uniform_int_distribution<node> d(0, G.numberOfNodes() - 1);

    int l = 0;

    for (int i = 0; i < k; i++) {
        node r = d(e);
        if (G.degree(r) >= 2) {
            node u = G.randomNeighbor(r);
            node w;
            do {
                w = G.randomNeighbor(r);
            } while (u == w);
            if (G.hasEdge(u, w))
                l++;
        }
    }

    return (double) l / (double) k;
}

} /* namespace NetworKit */
