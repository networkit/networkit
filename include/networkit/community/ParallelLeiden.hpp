#ifndef NETWORKIT_COMMUNITY_PARALLEL_LEIDEN_HPP_
#define NETWORKIT_COMMUNITY_PARALLEL_LEIDEN_HPP_

#include <atomic>
#include <condition_variable>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <omp.h>
#include <thread>
#include <tlx/unused.hpp>
#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Parallelism.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/coarsening/ParallelPartitionCoarsening.hpp>
#include <networkit/community/CommunityDetectionAlgorithm.hpp>
#include <networkit/community/Modularity.hpp>
#include <networkit/community/PLM.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

class ParallelLeiden final : public CommunityDetectionAlgorithm {
public:
    /**
     *
     * @param graph A networkit graph
     * @param iterations Number of Leiden Iterations to be run
     * @param randomize Randomize node order?
     * @param gamma Resolution parameter
     */
    explicit ParallelLeiden(const Graph &graph, int iterations = 3, bool randomize = true,
                            double gamma = 1);

    void run() override;

    int VECTOR_OVERSIZE = 10000;

private:
    inline double modularityDelta(double cutD, double degreeV, double volD) const {
        return cutD - gamma * degreeV * volD * inverseGraphVolume;
    };

    inline double modularityThreshold(double cutC, double volC, double degreeV) const {
        return cutC - gamma * (volC - degreeV) * degreeV * inverseGraphVolume;
    }

    static inline void lockLowerFirst(index a, index b, std::vector<std::mutex> &locks) {
        if (a < b) {
            locks[a].lock();
            locks[b].lock();
        } else {
            locks[b].lock();
            locks[a].lock();
        }
    }

    void flattenPartition();

    void calculateVolumes(const Graph &graph);

    void parallelMove(const Graph &graph);

    Partition parallelRefine(const Graph &graph);

    double inverseGraphVolume; // 1/vol(V)

    std::vector<double> communityVolumes;

    std::vector<std::vector<node>> mappings;

    static constexpr int WORKING_SIZE = 1000;

    double gamma; // Resolution parameter

    bool changed;

    int numberOfIterations;

    Aux::SignalHandler handler;

    bool random;
};

} // namespace NetworKit

#endif // NETWORKIT_COMMUNITY_PARALLEL_LEIDEN_HPP_
