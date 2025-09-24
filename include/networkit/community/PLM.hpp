/*
 * PLM.hpp
 *
 *  Created on: 20.11.2013
 *      Author: cls
 */

#ifndef NETWORKIT_COMMUNITY_PLM_HPP_
#define NETWORKIT_COMMUNITY_PLM_HPP_

#include <networkit/community/CommunityDetectionAlgorithm.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * Parallel Louvain Method - a multi-level modularity maximizer.
 */
class PLM final : public CommunityDetectionAlgorithm {

public:
    /**
     * @param[in] G input graph
     * @param[in] refine add a second move phase to refine the communities
     * @param[in] par parallelization strategy
     * @param[in] gammamulti-resolution modularity parameter:
     *            1.0 -> standard modularity
     *            0.0 -> one community
     *            2m -> singleton communities
     * @param[in] maxIter maximum number of iterations for move phase
     * @param[in] parallelCoarsening use parallel graph coarsening
     * @param[in] turbo faster but uses O(n) additional memory per thread
     * @param[in] recurse use recursive coarsening, see
     * http://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.049902 for some explanations
     * (default: true)
     *
     */
    PLM(const Graph &G, bool refine = false, double gamma = 1.0, std::string par = "balanced",
        count maxIter = 32, bool turbo = true, bool recurse = true);

    PLM(const Graph &G, const PLM &other);

    /**
     * @param[in] G input graph
     * @param[in] baseClustering optional; the algorithm will start from the given clustering.
     * @param[in] refine add a second move phase to refine the communities
     * @param[in] par parallelization strategy
     * @param[in] gammamulti-resolution modularity parameter:
     *            1.0 -> standard modularity
     *            0.0 -> one community
     *            2m -> singleton communities
     * @param[in] maxIter maximum number of iterations for move phase
     * @param[in] parallelCoarsening use parallel graph coarsening
     * @param[in] turbo faster but uses O(n) additional memory per thread
     * @param[in] recurse use recursive coarsening, see
     * http://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.049902 for some explanations
     * (default: true)
     *
     */
     PLM(const Graph &G, const Partition &baseClustering, bool refine = false, double gamma = 1.0, std::string par = "balanced",
        count maxIter = 32, bool turbo = true, bool recurse = true);

    /**
     * Detect communities.
     */
    void run() override;

    /**
     * Coarsens a graph based on a given partition and returns both the coarsened graph and a
     * mapping for the nodes from fine to coarse.
     *
     * @param graph The input graph
     * @param zeta Partition of the graph, which represents the desired state of the coarsened graph
     * @return pair of coarsened graph and node-mappings from fine to coarse graph
     */
    static std::pair<Graph, std::vector<node>> coarsen(const Graph &G, const Partition &zeta);

    /**
     * Calculates a partition containing the mapping of node-id from a fine graph
     * to a cluster-id from partition based on a coarse graph.
     *
     * @param Gcoarse Coarsened graph
     * @param zetaCoarse Partition, which contains information about clusters in the coarsened graph
     * @param Gfine Fine graph
     * @param nodeToMetaNode mapping for node-id from fine to coarse graph
     * @return Partition, which contains the cluster-id in the coarse graph for every node from the
     * fine graph
     */
    static Partition prolong(const Graph &Gcoarse, const Partition &zetaCoarse, const Graph &Gfine,
                             std::vector<node> nodeToMetaNode);

    /**
     * Returns fine-grained running time measurements for algorithm engineering purposes.
     */
    const std::map<std::string, std::vector<count>> &getTiming() const;

private:
    std::string parallelism;
    bool refine;
    double gamma = 1.0;
    count maxIter;
    bool turbo;
    bool recurse;
    std::map<std::string, std::vector<count>> timing; // fine-grained running time measurement
};

} /* namespace NetworKit */

#endif // NETWORKIT_COMMUNITY_PLM_HPP_
