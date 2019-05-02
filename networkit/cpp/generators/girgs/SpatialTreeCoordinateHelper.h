#pragma once

#include <array>

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


} // namespace girgs

#include <girgs/SpatialTreeCoordinateHelper.inl>
