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

protected:
    const Graph *G;

    Partition partition;

    bool bipartit = false;
};

} // NetworKit

#endif // NETWORKIT_BIPARTIT_BIPARTIT_HPP
