/*
 * AliasSampler.cpp
 *
 *  Created on: 03.07.2020
 *      Author: Klaus Ahrens  <ahrens@informatik.hu-berlin.de>
 *
 *  adapted and reimplemented in C++17
 *  from node2vec [https://arxiv.org/pdf/1607.00653v1.pdf]
 *
 */

#include <mutex>
#include <random>
#include <vector>

#include <networkit/auxiliary/Random.hpp>
#include <networkit/Globals.hpp>

#include "AliasSampler.hpp"

namespace NetworKit {
namespace Embedding {

std::mt19937_64& getURNG() {
    static bool seedPerThread = true;
    static std::once_flag once;
    std::call_once(once, []{Aux::Random::setSeed(1, seedPerThread);});

    return Aux::Random::getURNG();
}

double uniform_real() {
    static std::uniform_real_distribution<> dist;

    return dist(getURNG());
}

void AliasSampler::unigram(std::vector<float>& probs) {
        auto N = probs.size();

        std::vector<index> underV;
        std::vector<index> overV;
        for (index i = 0; i < N; ++i) {
            U[i] = probs[i] * N;
            if ( U[i] < 1 ) {
                underV.push_back(i);
            } else {
                overV.push_back(i);
            } 
        }
        while( !underV.empty() && !overV.empty() ) {
            auto small = underV.back();
            auto large = overV.back();
            underV.pop_back();
            overV.pop_back();
            K[small] = large;
            U[large] = (U[large] + U[small]) - 1;
            if (U[large] < 1) {
                underV.push_back(large);
            } else {
                overV.push_back(large);
            }
        } 
        while( !underV.empty() ){
            auto curr = underV.back();
            underV.pop_back();
            U[curr]=1;
        }
        while( !overV.empty() ){
            auto curr = overV.back();
            overV.pop_back();
            U[curr]=1;
        }
    }

index AliasSampler::sample() {
        double rx = uniform_real();
        double ry = uniform_real();
        index x = static_cast<index>(rx*K.size());
        double y = ry;
        return y < U[x] ? x : K[x];
}
 
} // namespace Embedding
} // namespace NetworKit
