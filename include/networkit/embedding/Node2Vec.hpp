/*
 *  Node2Vec.hpp
 *
 *
 *  Created on: 29.09.2020
 *      Author: Klaus Ahrens
 *              <ahrens@informatik.hu-berlin.de>
 *
 *  adapted and reimplemented in C++17
 *  from node2vec
 *  [https://arxiv.org/pdf/1607.00653v1.pdf]
 *
 */

#ifndef NETWORKIT_EMBEDDING_NODE2VEC_HPP_
#define NETWORKIT_EMBEDDING_NODE2VEC_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>


namespace NetworKit {

    /**
     * @ingroup embedding
     *
     * Feature extraction by node2vec algorithm
     */
    class Node2Vec final : public Algorithm {

    public:

        /**
        * Creates the Node2Vec class for @a G
        *
        * @param G     The graph.
        * @param P     Walk Return parameter (stay local).
        * @param Q     Walk In-Out parameter (drift away).
        * @param L     Walk length.
        * @param N     Walk count.
        * @param D     Dimension of embedding vectors.
        */
        Node2Vec(const Graph& G, double P = 1,
                 double Q = 1, count L = 80,
                 count N = 10, count D = 128
                );

        ~Node2Vec() override = default;

        /**
         This method computes all node embeddings.
        */
        void run() override;
        
        /**
         * Returns a string with the algorithm's name and its parameters, if there are any.
         * @return The string representation of the algorithm.
         */
        std::string toString() const override;

        /**
         * @return True if algorithm can run multi-threaded.
         */
        bool isParallel() const override;

        /**
         This method returns a vector that contains feature vectors for all nodes
        */
        std::vector<std::vector<float>> getFeatures();

    private:

        // The graph
        const Graph *G;
        // Walk parameter P
        double P;
        // Walk parameter Q
        double Q;
        // Walk length L
        count L;
        // Walk count N
        count N;
        // Feature Dimension D
        count D;
        
        std::vector<std::vector<float>> features;
    
    };

inline std::vector<std::vector<float>> Node2Vec::getFeatures() {
    assureFinished();
    return features;
}

} /* namespace NetworKit */

#endif // NETWORKIT_EMBEDDING_NODE2VEC_HPP_

