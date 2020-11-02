/*
 * LearnEmbeddings.hpp
 *
 *  Created on: 03.07.2020
 *      Author: Klaus Ahrens  <ahrens@informatik.hu-berlin.de>
 *
 *  adapted and reimplemented in C++17
 *  from node2vec [https://arxiv.org/pdf/1607.00653v1.pdf]
 *
 */

// networkit-format

#ifndef LEARN_EMBEDDINGS_HPP
#define LEARN_EMBEDDINGS_HPP

#include <fstream>
#include <vector>

#include <networkit/Globals.hpp>

#include "AliasSampler.hpp"
#include "BiasedRandomWalk.hpp"

namespace NetworKit {
namespace Embedding {

using Feature = std::vector<float>;
using Embeddings = std::vector<Feature>;

struct ModelData {
    AllWalks &allWalks;
    count dimensions;
    count winSize;
    count iterations;
    count allWords;
    AliasSampler &as;
    count &wordCntAll;
    double alpha;
    Embeddings &synNeg;
    Embeddings &synPos;
    ModelData(AllWalks &allWalks, count dimensions, count winSize, count iterations,
              AliasSampler &as, count &wordCntAll, double alpha, Embeddings &synNeg,
              Embeddings &synPos)
        : allWalks(allWalks), dimensions(dimensions), winSize(winSize), iterations(iterations),
          allWords(allWalks.size() * allWalks[0].size()), as(as), wordCntAll(wordCntAll),
          alpha(alpha), synNeg(synNeg), synPos(synPos) {}
};

/// Learns embeddings using SGD, Skip-gram with negative sampling.
Embeddings learnEmbeddings(AllWalks &walks, count dimensions, count winSize, count iterations);

/*
///print embeddings
void printEmbeddings(std::ofstream& o, Embeddings& e);
*/

} // namespace Embedding
} // namespace NetworKit
#endif // LEARN_EMBEDDINGS_HPP
