/*
 * Bipartit.hpp
 *
 * Created on: 18.09.2023
 *     Author: Michael Kaibel
 */

#ifndef NETWORKIT_BIPARTIT_BIPARTIT_HPP
#define NETWORKIT_BIPARTIT_BIPARTIT_HPP

#include <networkit/graph/Graph.hpp>
#include <networkit/base/Algorithm.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

class Bipartit : public Algorithm {
public:
    Bipartit(const Graph &G);

    void run() override;

    bool isBipartit();

    const Partition &getPartition();

    const std::vector<node> &getOddCircle();

protected:
    const Graph *G;

    Partition partition;

    std::vector<node> oddCircle;

    bool bipartit = false;

    void findOddCircle(std::vector<node> &parent, node v, node w);
};

} // NetworKit

#endif // NETWORKIT_BIPARTIT_BIPARTIT_HPP
