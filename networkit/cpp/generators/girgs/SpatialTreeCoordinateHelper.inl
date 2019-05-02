#include <girgs/BitManipulation.h>
#include <cassert>

namespace girgs {

template<unsigned int D>
unsigned int SpatialTreeCoordinateHelper<D>::cellOfLevel(unsigned cell) noexcept {
    // sets all bits below the most significant bit set in x
    auto assertLower = [] (uint32_t x) {
#if defined(__GNUC__) || defined(__clang__)
        if (__builtin_expect(!x, 0))
            return 0u; // __builtin_clz(x) is undefined for x == 0

        return static_cast<uint32_t>(
            (1llu << (32 - __builtin_clz(x))) - 1);
#else
        x |= x >> 1;
        x |= x >> 2;
        x |= x >> 4;
        x |= x >> 8;
        x |= x >> 16;
        return x;
#endif
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
