/*
 * LearnEmbeddings.cpp
 *
 *  Created on: 06.07.2020
 *      Author: Klaus Ahrens  <ahrens@informatik.hu-berlin.de>
 *
 *  adapted and reimplemented in C++17
 *  from node2vec [https://arxiv.org/pdf/1607.00653v1.pdf]
 *
 */

#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <vector>


#include <networkit/auxiliary/Random.hpp>

#include "AliasSampler.hpp"
#include "BiasedRandomWalk.hpp"
#include "LearnEmbedding.hpp"

//Code from https://github.com/nicholas-leonard/word2vec/blob/master/word2vec.c
//Customized for SNAP and node2vec

namespace NetworKit {
namespace Embedding {

//Number of negative samples. Value taken from original word2vec code.
const int negSamN = 5;

const int maxExp = 6;

//Learning rate for SGD. Value taken from original word2vec code.
const double startAlpha = 0.025;
const double minAlpha   = startAlpha * 0.0001;

using Vocab = std::vector<count>; 

Vocab learnVocab(AllWalks& allWalks, count nn) {
    Vocab vocab(nn);
    for (auto& walk: allWalks) {  // for every walk starting from this node
        for(auto& node: walk) {   // for every node in this walk
            ++vocab[node];
        }
    }   
    return vocab;
}

//Precompute unigram table using alias sampling method
AliasSampler vocabSampler(Vocab& vocab) {
    double trainWordsPow = 0;
    auto sz = vocab.size();
    std::vector<float> probV(sz); // NOT ... {sz}
    for (count i = 0; i < sz; ++i) {
        double e = std::pow(vocab[i], 0.75);
        probV[i] = e;
        trainWordsPow += e;
    }
    std::transform(probV.begin(), probV.end(), probV.begin(),
                    [&trainWordsPow](float f){return f/trainWordsPow;} );
    
    AliasSampler as(sz);
    as.unigram(probV);
    return as; 
}

//Initialize negative embeddings
Embeddings initNegEmb(Vocab& vocab, int dimensions) {
    return Embeddings{vocab.size(), std::vector<float>(dimensions)};
}

//Initialize positive embeddings
Embeddings initPosEmb(Vocab& vocab, int dimensions) {
    Embeddings emb{vocab.size(), std::vector<float>(dimensions)};

    for (index i = 0; i < emb.size(); ++i) {
        std::generate(emb[i].begin(), emb[i].end(),
                      [&](){return (uniform_real()-0.5)/dimensions;});
    }
    return emb;
}

void trainModel (ModelData& model, count walkNr)
{
    AllWalks& allWalks  = model.allWalks;
    count& wordCntAll   = model.wordCntAll;
    AliasSampler& as    = model.as;
    count  allWords     = model.allWords;
    count  iterations   = model.iterations;
    count  dimensions   = model.dimensions;
    count  winSize      = model.winSize;
    double alpha        = model.alpha;
    Embeddings& synPos  = model.synPos;
    Embeddings& synNeg  = model.synNeg;

    Walk& thisWalk = allWalks[walkNr]; 

    for (index wordI = 0; wordI < thisWalk.size(); ++wordI) {
        if (wordCntAll%10000 == 0) {
            alpha = startAlpha * (1 - wordCntAll / ((double)(iterations * allWords) + 1.0));
            if ( alpha < minAlpha ) { alpha = minAlpha; }
        }

        node word = thisWalk[wordI];
        std::vector<float> eV(dimensions);
        index offset = uniform_real() * winSize;

        for (index a = offset; a < winSize * 2 + 1 - offset; ++a) {
            if ( a == winSize ) { continue; }
            if ( wordI + a < winSize ) { continue; }
            if ( wordI + a >= winSize + thisWalk.size() ) { continue; }
            count currWordI = wordI + a - winSize;
            count currWord = thisWalk[currWordI];
            //negative sampling
            for (index j = 0; j < negSamN+1; ++j) {
                node target;
                int label;
                if (j == 0) {
                    target = word;
                    label = 1;
                } else {
                    target = as.sample();
                    if (target == word) { continue; }
                    label = 0;
                }             
                double product = 0;
                for (index i = 0; i < dimensions; ++i) {
                    product += synPos[currWord][i] * synNeg[target][i];
                }
                double grad;                     //Gradient multiplied by learning rate
                if (product > maxExp) { grad = (label - 1) * alpha; } 
                else if (product < -maxExp) { grad = label * alpha; }
                else { 
                    double exp = std::exp(product);
                    grad = (label - 1 + 1 / (1 + exp)) * alpha;
                }

                for (index i = 0; i < dimensions; ++i) { 
                    eV[i] += grad * synNeg[target][i];
                    synNeg[target][i] += grad * synPos[currWord][i];
                }
            }

            for (index i = 0; i < dimensions; ++i) {
                synPos[currWord][i] += eV[i];
            }
        }
        ++wordCntAll;
    }
}


Embeddings learnEmbeddings(AllWalks& allWalks, count dimensions,
                           count winSize, count iterations)
{
    count walks = allWalks.size();              // number of walks
    count nodesPerWalk = allWalks[0].size();    // number of nodes per walk

    std::unordered_map<node, node> renumbered;
    std::unordered_map<node, node> original;
    count nodes = 0;
    for (index i = 0; i < walks; ++i) {
        for (index j = 0; j < nodesPerWalk; ++j) {
            node& nd = allWalks[i][j];
            if (renumbered.count(nd)) {  
                nd = renumbered[nd];
            } else {
                renumbered[nd] = nodes;
                original[nodes] = nd; 
                nd = nodes++;
            }
        }
    }

    std::vector<count> vocab = learnVocab(allWalks, nodes);

    Embeddings embeddings(nodes, std::vector<float>(dimensions)); 

    Embeddings synNeg = initNegEmb(vocab, dimensions);
    Embeddings synPos = initPosEmb(vocab, dimensions);
    AliasSampler as = vocabSampler(vocab);

    double alpha = startAlpha;                              //learning rate
    count wordCntAll = 0;

    ModelData model(allWalks,
                    dimensions,
                    winSize, 
                    iterations,
                    as, 
                    wordCntAll,
                    alpha,
                    synNeg,
                    synPos); 

//#pragma omp parallel for schedule(dynamic) collapse(2)
    for (index iterCnt = 0; iterCnt < iterations; ++iterCnt) {
#pragma omp parallel for schedule(dynamic)
//      for (index run = 0; run < walks; ++run) { 
        for (omp_index run = 0; run < static_cast<omp_index>(walks); ++run) {
            trainModel(model, run);
        }
    }

    for (index i = 0; i < nodes; ++i) {
        embeddings[original[i]] = synPos[i];
    }

    allWalks.clear();
    return embeddings;
}

} // namespace Embedding 
} // namespace NetworKit
