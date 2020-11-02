/*
 *  Node2Vec.cpp
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

#include <vector>
#include <iostream>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/embedding/Node2Vec.hpp>

#include "AliasSampler.hpp"
#include "BiasedRandomWalk.hpp"
#include "LearnEmbedding.hpp"

namespace NetworKit {

Node2Vec::Node2Vec
    (
        const Graph& G, double P, double Q,
        count L, count N, count D
    )
    : G(&G), P(P), Q(Q), L(L), N(N), D(D)
{}

void Node2Vec::run() {
    using namespace NetworKit::Embedding;
    
    Aux::SignalHandler handler;
    
    TRACE("preprocess transition probabilities ...");
    
    preprocessTransitionProbs(*G, P, Q);
    handler.assureRunning();

    TRACE("do biased walks ...");
    auto walks = doWalks(*G, L, N);
    handler.assureRunning();
    
    TRACE("learn embeddings ...");
    count winSize = 10;
    count iterations = 1;
    features = learnEmbeddings(walks, D, winSize, iterations);
    handler.assureRunning();

    hasRun = true;
}

std::string Node2Vec::toString() const {
    return "Node2Vec()";
}

bool Node2Vec::isParallel() const {
    return true;
}

} /* namespace NetworKit */

