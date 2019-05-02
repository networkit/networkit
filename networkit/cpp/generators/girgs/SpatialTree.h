
#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <limits>
#include <numeric>
#include <cassert>

#include <omp.h>

#include <girgs/SpatialTreeCoordinateHelper.h>
#include <girgs/WeightLayer.h>


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
     * @param seed
     *  The seed for the edge sampling.
     *  If OpenMP is given more than one thread, thread i uses seed+i.
     *  This means that results are only reproducible for a combination of seed and thread number.
     */
    void generateEdges(int seed);

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
     */
    void visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level);

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
            unsigned int first_parallel_level, std::vector<std::vector<unsigned int>>& parallel_calls);

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

    std::vector<WeightLayer<D>> m_weight_layers;    ///< stores all nodes of one weight layer and provides the data structure described in paper
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> m_layer_pairs; ///< which pairs of weight layers to check in each level


    std::vector<default_random_engine> m_gens; ///< random generators for each thread

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


} // namespace girgs

#include <girgs/SpatialTree.inl>
