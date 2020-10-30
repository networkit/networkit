/*
 * BiasedRandomWalk.hpp
 *
 *  Created on: 03.07.2020
 *      Author: Klaus Ahrens  <ahrens@informatik.hu-berlin.de>
 *
 *  adapted and reimplemented in C++17
 *  from node2vec [https://arxiv.org/pdf/1607.00653v1.pdf]
 *
 */

#ifndef BIASEDRANDOMWALK_HPP
#define BIASEDRANDOMWALK_HPP

#include <vector>
#include <iostream>

#include <networkit/graph/Graph.hpp>

namespace NetworKit {
namespace Embedding {


using Walk      = std::vector<node>;      // one walk
using AllWalks  = std::vector<Walk>;      // n walks for all nodes: 2 dimensions in one big chunk


///preprocesses transition probabilities for random walks. Has to be called once before doWalks calls
void preprocessTransitionProbs(const Graph& graph, double paramP,
                               double paramQ);

///Simulates walks from every node and writes it into walks vector
AllWalks doWalks(const Graph& graph, count walkLen, count numberOfWalks);


} // namespace Embedding
} // namespace NetworKit

#endif // BIASEDRANDOMWALK_HPP
