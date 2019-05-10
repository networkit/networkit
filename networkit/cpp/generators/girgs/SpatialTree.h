/*
 * SpatialTree.h
 *
 *  Created on: 03. May 2019
 *      Author: Christopher Weyand <Christopher.Weyand@hpi.de>, Manuel Penschuck <networkit@manuel.jetzt>
 *
 * Code is adopted from https://github.com/chistopher/girgs
 */

#ifndef GENERATORS_GIRGS_SPATIAL_TREE_H_
#define GENERATORS_GIRGS_SPATIAL_TREE_H_

#include <vector>
#include <algorithm>
#include <random>
#include <limits>
#include <numeric>
#include <cassert>

#include <omp.h>
#include "../../auxiliary/Random.h"
#include "../../auxiliary/Parallel.h"
#include "../../auxiliary/SignalHandling.h"

#include "SpatialTreeCoordinateHelper.h"
#include "WeightLayer.h"

#include <tlx/math.hpp>

namespace NetworKit {
namespace girgs {

using default_random_engine = std::mt19937_64;

/**
 * @brief
 *  Internal implementation of the linear time GIRG sampling algorithm following the method object pattern.
 *  It requires a graph with sampled positions and weights.
 *
 * @tparam D
 *  Dimension of the underlying geometry.
 */
template<unsigned int D, typename EdgeCallback>
class SpatialTree
{
    using CoordinateHelper = SpatialTreeCoordinateHelper<D>;

public:
    SpatialTree(const std::vector<double>& weights, const std::vector<std::vector<double>>& positions, double alpha, EdgeCallback& edgeCallback, bool profile = false);

    /**
     * @brief
     *  Samples edges for given positions and weights.
     *
     * @param graph
     *  A graph with reasonable positions and weights.
     *  All nodes should have the same number of coordinates in [0,1).
     * @param alpha
     *  A parameter of the GIRG model. It determines the entropy of links.
     *  Infinity results in a deterministic threshold case.
     *  Zero produces a clique.
     */
    void generateEdges();

protected:

    /**
     * @brief
     *  A recursive function that samples all edges between points in cells A and B.
     *
     * @param cellA
     *  The source cell for edges.
     *  Edges sampled by this function are stored in nodes inside cellA.
     * @param cellB
     *  The target cell for edges.
     *  The reverse edges sampled by this function are not stored.
     * @param level
     *  The level from which A and B are, meaning cellA and cellB must be in the same level.
     * @param
     *  NetworKit's signal handler (thread local) to abort generation if requested
     */
    void visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level, Aux::SignalHandler& handler);

    /**
     * @brief
     *  Same as visitCellPair(unsigned int, unsigned int, unsigned int) but stops recursion before first_parallel_level.
     *  Instead, the calls that would be made in this level are saved in parallel_calls.
     *  The saved calls are grouped by their (level local) cellA parameter.
     *
     * @param cellA
     *  Same as in visitCellPair(unsigned int, unsigned int, unsigned int).
     * @param cellB
     *  Same as in visitCellPair(unsigned int, unsigned int, unsigned int).
     * @param level
     *  Same as in visitCellPair(unsigned int, unsigned int, unsigned int).
     * @param first_parallel_level
     *  The level before which we "saw off" the recursion.
     *  To get sufficient parallel cells (the outer size of parallel_calls) this should be computed as
     *  \f$ 2^{dl} \geq kt \f$ solved for l (d dimension, l first_parallel_level, t threads, k tuning parameter).
     *  We get \f$ l \geq \log_2(kt) / d \f$.
     * @param parallel_calls
     */
    void visitCellPair_sequentialStart(unsigned int cellA, unsigned int cellB, unsigned int level,
            unsigned int first_parallel_level, std::vector<std::vector<unsigned int>>& parallel_calls, Aux::SignalHandler& handler);

    /**
     * @brief
     *  Sample edges of type 1 between \f$ V_i^A V_j^B \f$.
     *  Type 1 means the cells A and B must touch or be identical.
     *
     * @param cellA
     *  Same as in visitCellPair(unsigned int, unsigned int, unsigned int).
     * @param cellB
     *  Same as in visitCellPair(unsigned int, unsigned int, unsigned int).
     * @param level
     *  Same as in visitCellPair(unsigned int, unsigned int, unsigned int).
     * @param i
     *  The weight layer for all considered nodes in cellA.
     * @param j
     *  The weight layer for all considered nodes in cellB.
     */
    void sampleTypeI(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j);

    /**
     * @brief
     *  Sample edges of type 2 between \f$ V_i^A V_j^B \f$.
     *  Type 2 means the cells A and B must not touch.
     *
     * @param cellA
     *  Same as in visitCellPair(unsigned int, unsigned int, unsigned int).
     * @param cellB
     *  Same as in visitCellPair(unsigned int, unsigned int, unsigned int).
     * @param level
     *  Same as in visitCellPair(unsigned int, unsigned int, unsigned int).
     * @param i
     *  The weight layer for all considered nodes in cellA.
     * @param j
     *  The weight layer for all considered nodes in cellB.
     */
    void sampleTypeII(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j);

    /**
     * @brief
     *  The insertion level for all nodes in specified weight layer.
     *  The minimal volume of cells in the desired level is called \f$v(i)\f$ in the paper (\f$=v(i,0)\f$).
     *
     * @param layer
     *  The index of the weight layer.
     *  Note that our layers are offset by 1, i.e. layer i in the paper is layer i+1 for us.
     * @return
     *  The level on which to insert nodes of given weight layer.
     */
    unsigned int weightLayerTargetLevel(int layer) const;

    /**
     * @brief
     *  The level on which type 1 pairs are compared in the Partitioning with volume \f$v(i,j)\f$
     *  (i.e. the partitioning corresponding to the layer-pair \f$(i,j)\f$).
     *  We call the deepest level at which cells have a volume greater than \f$v(i,j)\f$ the partitioning base level of i and j.
     *  So lets compute the partitioning base level for i and j:
     *
     *  Let \f$w_0\f$ be the minimum weight, \f$w_i\f$ the boundary for weight layer i,
     *  \f$d\f$ the dimension, \f$W\f$ the sum of weights, and \f$l\f$ the desired level.
     *  First observe that \f$w_i = 2w_{i-1} = 2^i * w_0\f$ if we split \f$w_i\f$ into a power of two and \f$w_0\f$.
     *  Now we get
     *
     *  \f$2^{-ld} = w_i w_j/W\\
     *  = 2^i w_0 \cdot 2^j w_0 / W\\
     *  = 2^i \cdot 2^j / (W/w_0^2)\\
     *  \Leftrightarrow -ld = i+j - \log_2(W/w_0^2)\\
     *  \Leftrightarrow l  = (\log_2(W/w_0^2) - (i+j)) / d\f$
     *
     *  a point with weight \f$w\f$ is inserted into layer \f$\lfloor\log_2(w/w_0)\rfloor\f$
     *  - \f$w_0\f$ is inserted in layer 0 (instead of layer 1 like in paper)
     *  - so in fact \f$w_i\f$ in our implementation equals \f$w_{i+1}\f$ in paper
     *
     *  a pair of layers is compared in level \f$\lfloor (\log_2(W/w_0^2) - (i+j+2)) / d \rfloor\f$
     *  - +2 to shift from our \f$w_i\f$ back to paper \f$w_i\f$
     *  - rounding down means a level with less depth like requested in paper
     *  - the constant \f$\log_2(W/w_0^2)\f$ is precomputed (see #m_baseLevelConstant)
     *
     * @param layer1
     *  The index of a weight layer. This is the \f$i\f$ in \f$v(i,j)\f$.
     *  Note that our layers are offset by 1, i.e. layer i in the paper is layer i-1 for us.
     * @param layer2
     *  The index of a second weight layer. This is the \f$j\f$ in \f$v(i,j)\f$.
     *  Note that our layers are offset by 1, i.e. layer i in the paper is layer i-1 for us.
     * @return
     *  The partitioning base level of the two given layers.
     */
    unsigned int partitioningBaseLevel(int layer1, int layer2) const;


    std::vector<WeightLayer<D>> buildPartition(
        const std::vector<double>& weights, const std::vector<std::vector<double>>& positions);


private:
    EdgeCallback& m_EdgeCallback; ///< called for every produced edge
    const bool m_profile;

    double m_alpha;             ///< girg model parameter, with higher alpha, long edges become less likely
    long long m_n;              ///< number of nodes in the graph

    double m_w0;                ///< minimum weight
    double m_wn;                ///< maximum weight
    double m_W;                 ///< sum of weights
    int    m_baseLevelConstant; ///< \f$\log_2(W/w_0^2)\f$ see partitioningBaseLevel(int, int) const

    unsigned int m_layers; ///< number of layers
    unsigned int m_levels; ///< number of levels

    std::vector<Node<D>>        m_nodes;            ///< nodes ordered by layer first and morton code second
    std::vector<unsigned int>   m_first_in_cell;    ///< prefix sums into nodes array
    std::vector<WeightLayer<D>> m_weight_layers;    ///< provides access to the nodes as described in paper

    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> m_layer_pairs; ///< which pairs of weight layers to check in each level

#ifndef NDEBUG
    long long m_type1_checks = 0; ///< number of node pairs that are checked via a type 1 check
    long long m_type2_checks = 0; ///< number of node pairs that are checked via a type 2 check
#endif // NDEBUG
};


/// provide automatic type deduction for constructor
template <unsigned int D, typename EdgeCallback>
SpatialTree<D,EdgeCallback> makeSpatialTree(const std::vector<double>& weights, const std::vector<std::vector<double>>& positions,
        double alpha, EdgeCallback& edgeCallback, bool profile = false) {
    return {weights, positions, alpha, edgeCallback, profile};
}


template<unsigned int D, typename EdgeCallback>
SpatialTree<D, EdgeCallback>::SpatialTree(const std::vector<double>& weights, const std::vector<std::vector<double>>& positions, double alpha, EdgeCallback& edgeCallback, bool profile)
    : m_EdgeCallback(edgeCallback)
    , m_profile(profile)
    , m_alpha(alpha)
    , m_n(weights.size())
    , m_w0(*std::min_element(weights.begin(), weights.end()))
    , m_wn(*std::max_element(weights.begin(), weights.end()))
    , m_W(std::accumulate(weights.begin(), weights.end(), 0.0))
    , m_baseLevelConstant(static_cast<int>(std::log2(m_W/m_w0/m_w0))) // log2(W/w0^2)
    , m_layers(static_cast<unsigned int>(floor(std::log2(m_wn/m_w0)))+1)
    , m_levels(partitioningBaseLevel(0,0) + 1) // (log2(W/w0^2) - 2) / d
{
    assert(weights.size() == positions.size());
    assert(positions.size() > 0 && positions.front().size() == D);


    // determine which layer pairs to sample in which level
    m_layer_pairs.resize(m_levels);
    for (auto i = 0u; i < m_layers; ++i)
        for (auto j = 0u; j < m_layers; ++j)
            m_layer_pairs[partitioningBaseLevel(i, j)].emplace_back(i,j);

    // sort weights into exponentially growing layers
    m_weight_layers = buildPartition(weights, positions);
}


template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::generateEdges() {
    // one random generator and distribution for each thread
    const auto num_threads = omp_get_max_threads();

#ifndef NDEBUG
    // ensure that all node pairs are compared either type 1 or type 2
    m_type1_checks = 0;
    m_type2_checks = 0;
#endif // NDEBUG

    // sample all edges
    if (num_threads == 1) {
        // sequential
        Aux::SignalHandler handler;
        visitCellPair(0, 0, 0, handler);
        assert(m_type1_checks + m_type2_checks == m_n*(m_n - 1ll));
        return;
    }

    // parallel see docs for visitCellPair_sequentialStart
    const auto first_parallel_level = static_cast<unsigned int>(std::ceil(std::log2(4.0*num_threads) / D));
    const auto parallel_cells = SpatialTreeCoordinateHelper<D>::numCellsInLevel(first_parallel_level);
    const auto first_parallel_cell = SpatialTreeCoordinateHelper<D>::firstCellOfLevel(first_parallel_level);

    // saw off recursion before "first_parallel_level" and save all calls that would be made
    auto parallel_calls = std::vector<std::vector<unsigned int>>(parallel_cells);
    {
        Aux::SignalHandler handler;
        visitCellPair_sequentialStart(0, 0, 0, first_parallel_level, parallel_calls, handler);
    }

    // do the collected calls in parallel
    int abort_caught = 0;
    #pragma omp parallel for schedule(static), num_threads(num_threads), reduction(+:abort_caught) // dynamic scheduling would be better but not reproducible
    for (int i = 0; i < parallel_cells; ++i) {
        try {
            Aux::SignalHandler handler;
            auto current_cell = first_parallel_cell + i;

            for (auto each : parallel_calls[i])
                visitCellPair(current_cell, each, first_parallel_level, handler);
        } catch (Aux::SignalHandling::InterruptException) {
            abort_caught++;
        }
    }

    if (abort_caught)
        throw Aux::SignalHandling::InterruptException();

    assert(m_type1_checks + m_type2_checks == m_n*(m_n - 1ll));
}


template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level, Aux::SignalHandler& handler) {

    handler.assureRunning();

    if(!CoordinateHelper::touching(cellA, cellB, level)) { // not touching
        // sample all type 2 occurrences with this cell pair
        #ifdef NDEBUG
        if (m_alpha == std::numeric_limits<double>::infinity()) return; // dont trust compilter optimization
        #endif // NDEBUG
        for(auto l=level; l<m_levels; ++l)
            for(auto& layer_pair : m_layer_pairs[l])
                sampleTypeII(cellA, cellB, level, layer_pair.first, layer_pair.second);
        return;
    }

    // sample all type 1 occurrences with this cell pair
    for(auto& layer_pair : m_layer_pairs[level]){
        if(cellA != cellB || layer_pair.first <= layer_pair.second)
            sampleTypeI(cellA, cellB, level, layer_pair.first, layer_pair.second);
    }

    // break if last level reached
    if(level == m_levels-1) // if we are at the last level we don't need recursive calls
        return;

    // recursive call for all children pairs (a,b) where a in A and b in B
    // these will be type 1 if a and b touch or type 2 if they don't
    for(auto a = CoordinateHelper::firstChild(cellA); a<=CoordinateHelper::lastChild(cellA); ++a)
        for(auto b = cellA == cellB ? a : CoordinateHelper::firstChild(cellB); b<=CoordinateHelper::lastChild(cellB); ++b)
            visitCellPair(a, b, level+1, handler);
}



template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::visitCellPair_sequentialStart(unsigned int cellA, unsigned int cellB, unsigned int level,
                                                                 unsigned int first_parallel_level,
                                                                 std::vector<std::vector<unsigned int>> &parallel_calls,
                                                                 Aux::SignalHandler& handler) {

    handler.assureRunning();

    if(!CoordinateHelper::touching(cellA, cellB, level)) { // not touching
        // sample all type 2 occurrences with this cell pair
        #ifdef NDEBUG
        if (m_alpha == std::numeric_limits<double>::infinity()) return; // dont trust compilter optimization
        #endif // NDEBUG
        for(auto l=level; l<m_levels; ++l)
            for(auto& layer_pair : m_layer_pairs[l])
                sampleTypeII(cellA, cellB, level, layer_pair.first, layer_pair.second);
        return;
    }

    // sample all type 1 occurrences with this cell pair
    for(auto& layer_pair : m_layer_pairs[level]){
        if(cellA != cellB || layer_pair.first <= layer_pair.second)
            sampleTypeI(cellA, cellB, level, layer_pair.first, layer_pair.second);
    }

    // break if last level reached
    if (level == m_levels - 1) // if we are at the last level we don't need recursive calls
        return;

    // recursive call for all children pairs (a,b) where a in A and b in B
    // these will be type 1 if a and b touch or type 2 if they don't
    for(auto a = CoordinateHelper::firstChild(cellA); a<=CoordinateHelper::lastChild(cellA); ++a)
        for(auto b = cellA == cellB ? a : CoordinateHelper::firstChild(cellB); b<=CoordinateHelper::lastChild(cellB); ++b){
            if(level+1 == first_parallel_level)
                parallel_calls[a-CoordinateHelper::firstCellOfLevel(first_parallel_level)].push_back(b);
            else
                visitCellPair_sequentialStart(a, b, level+1, first_parallel_level, parallel_calls, handler);
        }
}



template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::sampleTypeI(
    unsigned int cellA, unsigned int cellB, unsigned int level,
    unsigned int i, unsigned int j)
{
    assert(partitioningBaseLevel(i, j) == level || !CoordinateHelper::touching(cellA, cellB, level)); // in this case we were redirected from typeII with maxProb==1.0

    auto rangeA = m_weight_layers[i].cellIterators(cellA, level);
    auto rangeB = m_weight_layers[j].cellIterators(cellB, level);

    if (rangeA.first == rangeA.second || rangeB.first == rangeB.second)
        return;


#ifndef NDEBUG
    {
        const auto sizeV_i_A = std::distance(rangeA.first, rangeA.second);
        const auto sizeV_j_B = std::distance(rangeB.first, rangeB.second);

        #pragma omp atomic
        m_type1_checks += (cellA == cellB && i == j) ? sizeV_i_A * (sizeV_i_A - 1)  // all pairs in AxA without {v,v}
                                                     : sizeV_i_A * sizeV_j_B * 2; // all pairs in AxB and BxA
    }
#endif // NDEBUG

    std::uniform_real_distribution<> dist;
    const auto threadId = omp_get_thread_num();
    auto& gen = Aux::Random::getURNG();

    const auto inThresholdMode = m_alpha == std::numeric_limits<double>::infinity();

    int kA = 0;
    for(auto pointerA = rangeA.first; pointerA != rangeA.second; ++kA, ++pointerA) {
        auto offset = (cellA == cellB && i==j) ? kA+1 : 0;
        for (auto pointerB = rangeB.first + offset; pointerB != rangeB.second; ++pointerB) {

            const auto& nodeInA = *pointerA;
            const auto& nodeInB = *pointerB;

            // pointer magic gives same results
            assert(nodeInA.index == m_weight_layers[i].kthPoint(cellA, level, kA).index);
            assert(nodeInB.index == m_weight_layers[j].kthPoint(cellB, level, std::distance(rangeB.first, pointerB)).index);

            // points are in correct cells
            assert(cellA - CoordinateHelper::firstCellOfLevel(level) == CoordinateHelper::cellForPoint(nodeInA.coord, level));
            assert(cellB - CoordinateHelper::firstCellOfLevel(level) == CoordinateHelper::cellForPoint(nodeInB.coord, level));

            // points are in correct weight layer
            assert(i == static_cast<unsigned int>(std::log2(nodeInA.weight/m_w0)));
            assert(j == static_cast<unsigned int>(std::log2(nodeInB.weight/m_w0)));

            assert(nodeInA.index != nodeInB.index);
            const auto distance = nodeInA.distance(nodeInB);
            const auto w_term = nodeInA.weight*nodeInB.weight/m_W;
            const auto d_term = tlx::power_to_the<D>(distance);

            if(inThresholdMode) {
                if(d_term < w_term)
                    m_EdgeCallback(nodeInA.index, nodeInB.index, threadId);
            } else {
                auto edge_prob = std::pow(w_term/d_term, m_alpha); // we don't need min with 1.0 here
                if(dist(gen) < edge_prob)
                    m_EdgeCallback(nodeInA.index, nodeInB.index, threadId);
            }
        }
    }
}


template<unsigned int D, typename EdgeCallback>
void SpatialTree<D, EdgeCallback>::sampleTypeII(
    unsigned int cellA, unsigned int cellB, unsigned int level,
    unsigned int i, unsigned int j)
{
    assert(partitioningBaseLevel(i, j) >= level);

    auto rangeA = m_weight_layers[i].cellIterators(cellA, level);
    auto rangeB = m_weight_layers[j].cellIterators(cellB, level);

    if (rangeA.first == rangeA.second || rangeB.first == rangeB.second)
        return;

    const auto sizeV_i_A = std::distance(rangeA.first, rangeA.second);
    const auto sizeV_j_B = std::distance(rangeB.first, rangeB.second);

    // get upper bound for probability
    const auto w_upper_bound = m_w0*(1<<(i+1)) * m_w0*(1<<(j+1)) / m_W;
    const auto cell_distance = CoordinateHelper::dist(cellA, cellB, level);
    const auto dist_lower_bound = tlx::power_to_the<D>(cell_distance);
    const auto max_connection_prob = std::min(std::pow(w_upper_bound/dist_lower_bound, m_alpha), 1.0);
    assert(dist_lower_bound > w_upper_bound); // in threshold model we would not sample anything
    const auto num_pairs = sizeV_i_A * sizeV_j_B;
    const auto expected_samples = num_pairs * max_connection_prob;

    // skipping over points is actually quite expensive as it messes up
    // branch predictions and prefetching. Hence low expected skip distances
    // it's cheapter to throw a coin each time!
    if (max_connection_prob > 0.2) {
        return sampleTypeI(cellA, cellB, level, i, j);
    }

#ifndef NDEBUG
#pragma omp atomic
    m_type2_checks += 2 * sizeV_i_A * sizeV_j_B;
#endif // NDEBUG

    if(expected_samples < 1e-6)
        return;

    // init geometric distribution
    auto threadId = omp_get_thread_num();
    auto& gen = Aux::Random::getURNG();
    auto geo = std::geometric_distribution<unsigned long long>(max_connection_prob);
    auto dist = std::uniform_real_distribution<>(0, max_connection_prob);

    for (auto r = geo(gen); r < num_pairs; r += 1 + geo(gen)) {
        // determine the r-th pair
        const Node<D>& nodeInA = rangeA.first[r%sizeV_i_A];
        const Node<D>& nodeInB = rangeB.first[r/sizeV_i_A];

        nodeInB.prefetch();
        nodeInA.prefetch();

        // points are in correct weight layer
        assert(i == static_cast<unsigned int>(std::log2(nodeInA.weight/m_w0)));
        assert(j == static_cast<unsigned int>(std::log2(nodeInB.weight/m_w0)));

        const auto rnd = dist(gen);

        // get actual connection probability
        const auto distance = nodeInA.distance(nodeInB);
        const auto w_term = nodeInA.weight*nodeInB.weight/m_W;
        const auto d_term = tlx::power_to_the<D>(distance);
        const auto connection_prob = std::pow(w_term/d_term, m_alpha); // we don't need min with 1.0 here
        assert(w_term < w_upper_bound);
        assert(d_term >= dist_lower_bound);

        if(rnd < connection_prob) {
            m_EdgeCallback(nodeInA.index, nodeInB.index, threadId);
        }
    }
}


template<unsigned int D, typename EdgeCallback>
unsigned int SpatialTree<D, EdgeCallback>::weightLayerTargetLevel(int layer) const {
    // -1 coz w0 is the upper bound for layer 0 in paper and our layers are shifted by -1
    auto result = std::max((m_baseLevelConstant - layer - 1) / (int)D, 0);
#ifndef NDEBUG
    {   // a lot of assertions that we have the correct insertion level
        assert(0 <= layer && layer < m_layers);
        assert(0 <= result && result <= m_levels); // note the result may be one larger than the deepest level (hence the <= m_levels)
        auto volume_requested  = m_w0*m_w0*std::pow(2,layer+1)/m_W; // v(i) = w0*wi/W
        auto volume_current    = std::pow(2.0, -(result+0.0)*D); // in paper \mu with v <= \mu < O(v)
        auto volume_one_deeper = std::pow(2.0, -(result+1.0)*D);
        assert(volume_requested <= volume_current || volume_requested >= 1.0); // current level has more volume than requested
        assert(volume_requested >  volume_one_deeper);    // but is the smallest such level
    }
#endif // NDEBUG
    return static_cast<unsigned int>(result);
}


template<unsigned int D, typename EdgeCallback>
unsigned int SpatialTree<D, EdgeCallback>::partitioningBaseLevel(int layer1, int layer2) const {

    // we do the computation on signed ints but cast back after the max with 0
    // m_baseLevelConstant is just log(W/w0^2)
    auto result = std::max((m_baseLevelConstant - layer1 - layer2 - 2) / (int)D, 0);
#ifndef NDEBUG
    {   // a lot of assertions that we have the correct comparison level
        assert(0 <= layer1 && layer1 < m_layers);
        assert(0 <= layer2 && layer2 < m_layers);
        auto volume_requested  = m_w0*std::pow(2,layer1+1) * m_w0*std::pow(2,layer2+1) / m_W; // v(i,j) = wi*wj/W
        auto volume_current    = std::pow(2.0, -(result+0.0)*D); // in paper \mu with v <= \mu < O(v)
        auto volume_one_deeper = std::pow(2.0, -(result+1.0)*D);
        assert(volume_requested <= volume_current || volume_requested >= 1.0); // current level has more volume than requested
        assert(volume_requested >  volume_one_deeper);    // but is the smallest such level
    }
#endif // NDEBUG
    return static_cast<unsigned int>(result);
}

template<unsigned int D, typename EdgeCallback>
std::vector<WeightLayer<D>> SpatialTree<D, EdgeCallback>::buildPartition(const std::vector<double>& weights, const std::vector<std::vector<double>>& positions) {

    const auto n = weights.size();
    assert(positions.size() == n);

    auto weight_to_layer = [=] (double weight) {
        return std::log2(weight / m_w0);
    };

    const auto first_cell_of_layer = [&] {
        std::vector<unsigned int> first_cell_of_layer(m_layers + 1);
        unsigned int sum = 0;
        for (auto l = 0; l < m_layers; ++l) {
            first_cell_of_layer[l] = sum;
            sum += CoordinateHelper::numCellsInLevel(weightLayerTargetLevel(l));
        }
        first_cell_of_layer.back() = sum;
        return first_cell_of_layer;
    }();
    const auto max_cell_id = first_cell_of_layer.back();

    // Node<D> should incur no init overhead; checked on godbolt
    m_nodes = std::vector<Node<D>>(n);

    // compute the cell a point belongs to
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        const auto layer = weight_to_layer(weights[i]);
        const auto level = weightLayerTargetLevel(layer);
        m_nodes[i] = Node<D>(positions[i], weights[i], i);
        m_nodes[i].cell_id = first_cell_of_layer[layer] + CoordinateHelper::cellForPoint(m_nodes[i].coord, level);
        assert(m_nodes[i].cell_id < max_cell_id);
    }

    // Sort points by cell-ids
        Aux::Parallel::sort(m_nodes.begin(), m_nodes.end(),
            [](const Node<D> &a, const Node<D> &b) { return a.cell_id < b.cell_id; });

    // compute pointers into points
    constexpr auto gap_cell_indicator = std::numeric_limits<unsigned int>::max();
    m_first_in_cell = std::vector<unsigned int>(max_cell_id + 1, gap_cell_indicator);

    {
        m_first_in_cell[max_cell_id] = n;

        // First, we mark the begin of cells that actually contain points
        // and repair the gaps (i.e., empty cells) later. In the mean time,
        // the values of those gaps will remain at gap_cell_indicator.
        m_first_in_cell[m_nodes[0].cell_id] = 0;
        #pragma omp parallel for
        for (int i = 1; i < n; ++i) {
            if (m_nodes[i - 1].cell_id != m_nodes[i].cell_id) {
                m_first_in_cell[m_nodes[i].cell_id] = i;
            }
        }

        // Now repair gaps: since first_point_in_cell shell contain
        // a prefix sum, we simply replace any "gap_cell_indicator"
        // with its nearest non-gap successor. In the main loop,
        // this is always the direct successors since we're iterating
        // from right to left.
        #pragma omp parallel
        {
            const auto threads = omp_get_num_threads();
            const auto rank = omp_get_thread_num();
            const auto chunk_size = (max_cell_id + threads - 1) / threads; // = ceil(max_cell_id / threads)

            // Fix right-most element (if gap) of each thread's elements by looking into chunk of next thread.
            // We do not need an end of array check, since it's guaranteed that the last element is n.
            // We're using on this very short code block to avoid UB even if we're only performing word-wise updates.
            #pragma omp single
            {
                for (int r = 0; r < threads; r++) {
                    const auto end = std::min(max_cell_id, chunk_size * (r + 1));
                    int first_non_invalid = end - 1;
                    while (m_first_in_cell[first_non_invalid] == gap_cell_indicator)
                        first_non_invalid++;
                    m_first_in_cell[end - 1] = m_first_in_cell[first_non_invalid];
                }
            }

            const auto begin = std::min(max_cell_id, chunk_size * rank);

            auto i = std::min(max_cell_id, begin + chunk_size);
            while (i-- > begin) {
                m_first_in_cell[i] = std::min(
                    m_first_in_cell[i],
                    m_first_in_cell[i + 1]);
            }
        }

        #ifndef NDEBUG
        {
            assert(m_nodes[n-1].cell_id < max_cell_id);

            // assert that we have a prefix sum starting at 0 and ending in n
            assert(m_first_in_cell[0] == 0);
            assert(m_first_in_cell[max_cell_id] == n);
            assert(std::is_sorted(m_first_in_cell.cbegin(), m_first_in_cell.cend()));

            // check that each point is in its right cell (and that the cell boundaries are correct)
            for (auto cid = 0u; cid != max_cell_id; ++cid) {
                const auto begin = m_nodes.begin() + m_first_in_cell[cid];
                const auto end = m_nodes.begin() + m_first_in_cell[cid + 1];
                for (auto it = begin; it != end; ++it)
                    assert(it->cell_id == cid);
            }
        }
        #endif
    }

    // build spatial structure and find insertion level for each layer based on lower bound on radius for current and smallest layer
    std::vector<WeightLayer<D>> weight_layers;
    weight_layers.reserve(m_layers);
    for (auto layer = 0u; layer < m_layers; ++layer) {
        weight_layers.emplace_back(weightLayerTargetLevel(layer), m_nodes.data(), m_first_in_cell.data() + first_cell_of_layer[layer]);
    }

    return weight_layers;
}

} // namespace girgs
} // namespace NetworKit

#endif // GENERATORS_GIRGS_SPATIAL_TREE_H_
