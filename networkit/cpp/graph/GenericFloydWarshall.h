/*
 * GenericFloydWarshall.h
 */

#ifndef GENERICFLOYDWARSHALL_H_
#define GENERICFLOYDWARSHALL_H_

#include "Graph.h"
#include "APAPP.h"

namespace NetworKit {

/**
 * @ingroup graph
 * Floyd-Warshall all-pair algorithm
 */
template <class T>
class GenericFloydWarshall : public APAPP<T> {

public:

    /**
     * Creates the FloydWarshall class for @a G and the source node @a
     * souce.
     *
     * @param G The graph.
     * @param source The source node.
     * @param storePaths store paths and number of paths?
     */
    GenericFloydWarshall(const Graph& G, std::vector<T> edgeWeights={}) : APAPP<T>(G, edgeWeights) {

    }

    /**
     * Performs the Floyd-Warshall all-pair algorithm on the graph given in
     * the contructor.
     */
    void run() {

        // In theory the algorithm creates a new matrix for every
        // step, where we would end up with as many matrices as there
        // are nodes. We just create two matrices and overwrite them
        // every second step.
        std::vector<std::vector<T> > distanceMatrix0(G.upperNodeIdBound(), std::vector<T>(G.upperNodeIdBound())); 
        std::vector<std::vector<T> > distanceMatrix1(G.upperNodeIdBound(), std::vector<T>(G.upperNodeIdBound())); 
        // initialize distances
        T zero = T::getZero();
        T one = T::getOne();

        // set distanceMatrix0 to all onees
        for (u_int k = 0; k < G.upperNodeIdBound(); k++) {
            for (u_int i = 0; i < G.upperNodeIdBound(); i++) {
                distanceMatrix0[k][i] = T(one);
            }
        }

        // set distanceMatrix1 to all infinity...
        for (u_int i = 0; i < G.upperNodeIdBound(); i++) {
            for (u_int j = 0; j < G.upperNodeIdBound(); j++) {
                distanceMatrix1[i][j] = T(zero);
            }
        }

        // ...except where there are edges
        for (u_int k = 0; k < G.upperNodeIdBound(); k++) {
            G.forEdges([&](node u, node v, edgeid e) {
                    distanceMatrix1[u][v] = T(edgeWeights[e]);
                    if (!G.isDirected()) {
                        distanceMatrix1[v][u] = T(edgeWeights[e]);
                    }
                    });
            distanceMatrix1[k][k] =  T(one);
        }

        // perform the main loop
        for (u_int k = 0; k < G.upperNodeIdBound(); k++) {
            for (u_int i = 0; i < G.upperNodeIdBound(); i++) {
                for (u_int j = 0; j < G.upperNodeIdBound(); j++) {
                    if (k % 2 == 0) {
                        distanceMatrix0[i][j] = distanceMatrix1[i][j] + T(distanceMatrix1[i][k] * distanceMatrix1[k][j]);
                    } else {
                        distanceMatrix1[i][j] = distanceMatrix0[i][j] + T(distanceMatrix0[i][k] * distanceMatrix0[k][j]);
                    }
                }
            }
        }

        if ((G.upperNodeIdBound() - 1) % 2 == 0) {
            distances = distanceMatrix0;
        } else {
            distances = distanceMatrix1;
        }
        for (u_int i = 0; i < G.upperNodeIdBound(); i++) {
            distances[i][i] = T::getOne();
        }
    }

    std::string toString() const {
        return "Generic all-pairs Algebraic Path Algorithm by Floyd-Warshall";
    }

private:

    using APAPP<T>::G;
    using APAPP<T>::distances;
    using APAPP<T>::edgeWeights;

};

} /* namespace NetworKit */
#endif /* GENERICFLOYDWARSHALL_H_ */
