/*
 * SpatialTreeCoordinateHelper.h
 *
 *  Created on: 03. May 2019
 *      Author: Christopher Weyand <Christopher.Weyand@hpi.de>, Manuel Penschuck <networkit@manuel.jetzt>
 *
 * Code is adopted from https://github.com/chistopher/girgs
 */

#ifndef GENERATORS_GIRGS_SPATIAL_TREE_COORDINATE_HELPER_H_
#define GENERATORS_GIRGS_SPATIAL_TREE_COORDINATE_HELPER_H_

#include <array>
#include <cassert>

#include <tlx/define/likely.hpp>
#include <tlx/math/clz.hpp>

#include "BitManipulation.hpp"

namespace NetworKit {
namespace girgs {

template<unsigned int D>
class SpatialTreeCoordinateHelper
{
public:
    static constexpr unsigned int numChildren() noexcept {
        return 1u<<D;
    }

    static constexpr unsigned int numCellsInLevel(unsigned int level) noexcept {
        return 1u<<(D*level);
    }

    static constexpr unsigned int firstCellOfLevel(unsigned int level) noexcept {
        return ((1u<<(D*level))-1)/(numChildren()-1);
    }

    static constexpr unsigned int parent(unsigned int cell) noexcept {
        return (cell-1) / numChildren();
    }

    static constexpr unsigned int firstChild(unsigned int cell) noexcept {
        return numChildren() * cell + 1;
    }

    static constexpr unsigned int lastChild(unsigned int cell) noexcept {
        return firstChild(cell) + numChildren() - 1;
    }


    static unsigned int cellOfLevel(unsigned cell) noexcept;

    static unsigned int cellForPoint(const std::array<double, D>& position, unsigned int targetLevel) noexcept;

    static std::array<std::pair<double,double>, D> bounds(unsigned int cell, unsigned int level) noexcept;

    static bool touching(unsigned int cellA, unsigned int cellB, unsigned int level) noexcept;

    static double dist(unsigned int cellA, unsigned int cellB, unsigned int level) noexcept;

    SpatialTreeCoordinateHelper() = delete; // we want to support static accesses only
};


template<unsigned int D>
unsigned int SpatialTreeCoordinateHelper<D>::cellOfLevel(unsigned cell) noexcept {
    // sets all bits below the most significant bit set in x
    auto assertLower = [] (uint32_t x) {
        if (TLX_UNLIKELY(x == 0))
            return 0u; // __builtin_clz(x) is undefined for x == 0

        return static_cast<uint32_t>(
            (1llu << (32 - tlx::clz(x))) - 1);
    };

    constexpr auto mask = BitPattern<D>::kEveryDthBit;

    auto firstCellInLayer = mask & assertLower(cell);

    if (cell < firstCellInLayer)
        firstCellInLayer >>= D;

    return cell - firstCellInLayer;
}

template<unsigned int D>
std::array<std::pair<double, double>, D> SpatialTreeCoordinateHelper<D>::bounds(unsigned int cell, unsigned int level) noexcept {
    const auto diameter = 1.0 / (1<<level);
    const auto coord = BitManipulation<D>::extract(cellOfLevel(cell));

    auto result = std::array<std::pair<double, double>, D>();
    for(auto d=0u; d<D; ++d)
        result[d]= {coord[d]*diameter, (coord[d]+1)*diameter };

    return result;
}

template<unsigned int D>
unsigned int SpatialTreeCoordinateHelper<D>::cellForPoint(const std::array<double, D>& position, unsigned int targetLevel) noexcept {
    const auto diameter = static_cast<double>(1 << targetLevel);

    std::array<uint32_t, D> coords;
    for (auto d = 0u; d < D; ++d)
        coords[d] = static_cast<uint32_t>(position[d] * diameter);

    return BitManipulation<D>::deposit(coords);
}

template<unsigned int D>
bool SpatialTreeCoordinateHelper<D>::touching(unsigned int cellA, unsigned int cellB, unsigned int level) noexcept  {
    const auto coordA = BitManipulation<D>::extract(cellOfLevel(cellA));
    const auto coordB = BitManipulation<D>::extract(cellOfLevel(cellB));

    auto touching = true;
    for(auto d=0u; d<D; ++d){
        auto dist = std::abs(static_cast<int>(coordA[d]) - static_cast<int>(coordB[d]));
        dist = std::min(dist, (1<<level) - dist);
        touching &= (dist <= 1);
    }

    return touching;
}

template<unsigned int D>
double SpatialTreeCoordinateHelper<D>::dist(unsigned int cellA, unsigned int cellB, unsigned int level) noexcept  {
    // first work with integer d dimensional index
    const auto coordA = BitManipulation<D>::extract(cellOfLevel(cellA));
    const auto coordB = BitManipulation<D>::extract(cellOfLevel(cellB));

    auto result = 0;
    for(auto d=0u; d<D; ++d){
        auto dist = std::abs(static_cast<int>(coordA[d]) - static_cast<int>(coordB[d]));
        dist = std::min(dist, (1<<level) - dist);
        result = std::max(result, dist);
    }

    // then apply the diameter
    auto diameter = 1.0 / (1<<level);
    return std::max(0.0, (result-1) * diameter); // TODO if cellA and cellB are not touching, this max is irrelevant
}


} // namespace girgs
} // namespace NetworKit

#endif // GENERATORS_GIRGS_SPATIAL_TREE_COORDINATE_HELPER_H_
