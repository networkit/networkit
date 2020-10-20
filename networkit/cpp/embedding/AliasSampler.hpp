/*
 * AliasSampler.hpp
 *
 *  Created on: 03.07.2020
 *      Author: Klaus Ahrens  <ahrens@informatik.hu-berlin.de>
 *
 *  adapted and reimplemented in C++17 
 *  from node2vec [https://arxiv.org/pdf/1607.00653v1.pdf]
 *               
 */

#ifndef ALIASSAMPLER_HPP
#define ALIASSAMPLER_HPP

#include <random>
#include <vector>

#include <networkit/Globals.hpp>

namespace NetworKit::Embedding {


std::mt19937_64& getURNG();

double uniform_real();

struct AliasSampler {
    std::vector<index> K;
    std::vector<float> U;
public:
    //compute unigram table using alias sampling method
    void unigram(std::vector<float>& probs);

    AliasSampler(count degree = 0): K(degree), U(degree) {} // allowing as mapped_type in unordered_maps
    index sample();
};

}

#endif // ALIASSAMPLER_HPP
