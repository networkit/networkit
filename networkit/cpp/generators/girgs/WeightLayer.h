#pragma once

#include <memory>

#include <girgs/Node.h>
#include <girgs/SpatialTreeCoordinateHelper.h>


namespace girgs {


/**
 * @brief
 *  This class implements the data structure to manage point access described in the paper (Lemma 4.1).
 *  The partitioning of the ground space (Lemma 4.2) implicitly results from the implementation of SpatialTree.
 *
 * @tparam D
 *  the dimension of the geometry
 */
template<unsigned int D>
class WeightLayer {
    using Helper = SpatialTreeCoordinateHelper<D>;

public:
    WeightLayer() = delete;

    WeightLayer(const WeightLayer&) = delete;
    WeightLayer& operator=(const WeightLayer&) = delete;

    WeightLayer(WeightLayer&&) = default;
    WeightLayer& operator=(WeightLayer&&) = default;

    WeightLayer(unsigned int targetLevel,
                std::shared_ptr<Node<D>[]>& base,
                std::shared_ptr<unsigned int[]>& prefix_sum_ownership,
                const unsigned int* prefix_sum)
        : m_target_level{targetLevel},
          m_base{base}, m_prefix_sum_ownership(prefix_sum_ownership), m_prefix_sums{prefix_sum}
    {}

    /**
     * @brief
     *  Returns the number of points of this weight layer in a cell.
     *  In the notation of the paper, this function returns \f$ |V_i^{cell}| \f$, where i is the index of this weight layer.
     *
     * @param cell
     *  The cell that contains the points.
     * @param level
     *  The level of the given cell. This should be less or equal to the target level of this weight layer.
     * @return
     *  Returns how many points there are in cells {begin..end} using prefix sums. Begin and end are the first/last descendants of cell in target level.
     */
    int pointsInCell(unsigned int cell, unsigned int level) const {
        auto cellBoundaries = levelledCell(cell, level);
        assert(cellBoundaries.first  + Helper::firstCellOfLevel(level) < Helper::firstCellOfLevel(m_target_level+1));
        assert(cellBoundaries.second + Helper::firstCellOfLevel(level) < Helper::firstCellOfLevel(m_target_level+1));

        return m_prefix_sums[cellBoundaries.second+1] - m_prefix_sums[cellBoundaries.first];
    }


    /**
     * @brief
     *  Implements the second operation required for the data structure in the paper (Lemma 4.1).
     *  The method finds the first descendant of the given cell in the target level.
     *  Then #m_prefix_sums and #m_A are used to find the requested point.
     *
     * @param cell
     *  The cell that contains the points.
     * @param level
     *  The level of the given cell.
     *  This should be less or equal to the target level of this weight layer.
     * @param k
     *  The point we want to access.
     *  This should be less than the number of points of this weight layer in the given cell
     *  (i.e. less than what was returned by pointsInCell(unsigned int, unsigned int) const ).
     * @return
     *  Returns the requested node.
     */
    const Node<D>& kthPoint(unsigned int cell, unsigned int level, int k) const {
        auto cellBoundaries = levelledCell(cell, level);
        return m_base.get()[m_prefix_sums[cellBoundaries.first] + k];
    }


    /**
     * @brief
     *  Return points to the first element of the cell ("begin") and the first element after the cell ("end")
     *
     * @param cell
     *  The cell that contains the points.
     * @param level
     *  The level of the given cell.
     *  This should be less or equal to the target level of this weight layer.
     * @return
     *  {begin, end}
     */
    std::pair<const Node<D>*, const Node<D>*> cellIterators(unsigned int cell, unsigned int level) const {
        auto cellBoundaries = levelledCell(cell, level);
        const auto begin_end = std::make_pair(m_base.get() + m_prefix_sums[cellBoundaries.first],
                                              m_base.get() + m_prefix_sums[cellBoundaries.second+1]);
        assert(begin_end.first <= begin_end.second);
        return begin_end;
    }

protected:
    const unsigned int m_target_level;      ///< the insertion level for the current weight layer (v(i) = wiw0/W)
    std::shared_ptr<Node<D>[]> m_base;                        ///< Pointer to the first point stored in this layer
    std::shared_ptr<unsigned int[]> m_prefix_sum_ownership; ///< not used directly, simply keep memory pointed into by m_prefix_sums alive
    const unsigned int* m_prefix_sums;                      ///< for each cell c in target level: the sum of points of this layer in all cells <c

    std::pair<unsigned int, unsigned int> levelledCell(unsigned int cell, unsigned int level) const {
        assert(level <= m_target_level);
        assert(Helper::firstCellOfLevel(level) <= cell && cell < Helper::firstCellOfLevel(level + 1)); // cell is from fromLevel

        // we want the begin-th and end-th cell in level targetLevel to be the first and last descendant of cell in this level
        // we could apply the firstChild function to find the first descendant but this is in O(1)
        auto descendants = Helper::numCellsInLevel(m_target_level - level);
        auto localIndexCell = cell - Helper::firstCellOfLevel(level);
        auto localIndexDescendant = localIndexCell * descendants; // each cell before the parent splits in 2^D cells in the next layer that are all before our descendant
        auto begin = localIndexDescendant;
        auto end = begin + descendants - 1;

        assert(begin <= end);

        return {begin, end};
    }
};

} // namespace girgs
