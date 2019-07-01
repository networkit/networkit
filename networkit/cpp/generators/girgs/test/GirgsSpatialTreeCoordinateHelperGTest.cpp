#include <array>
#include <iostream>
#include <sstream>
#include <string>
#include <random>
#include <utility>

#include <gtest/gtest.h>

#include "../SpatialTreeCoordinateHelper.hpp"

class GirgsSpatialTreeCoordinateHelperGTest: public testing::Test {};

template <typename T, unsigned D>
std::string stringify(const std::array<T, D>& arr) {
    std::stringstream ss;
    ss << "[" << arr[0];
    for(int i=1; i < D; ++i)
        ss << " " << arr[i];
    ss << "]";
    return ss.str();
}

template <typename T, unsigned D>
std::string stringify(const std::array<std::pair<T, T>, D>& arr) {
    std::stringstream ss;
    ss << "[";
    for(int i=0; i < D; ++i) {
        if (!i) ss << " ";
        ss << "{" << arr[i].first << ", " << arr[i].second << "}";
    }
    ss << "]";
    return ss.str();
}


template<unsigned int D>
void testTreeStructure(const unsigned max_level) {
    using Tree = NetworKit::girgs::SpatialTreeCoordinateHelper<D>;

    // check the number of children on all layers and helper functions
    auto cells = Tree::firstCellOfLevel(max_level+1);
    auto sumCells = 1; // root in level 0
    auto numChildren = std::vector<int>(cells, 0);
    for(auto l=1u; l<=max_level; ++l) {
        EXPECT_EQ(sumCells, Tree::firstCellOfLevel(l));
        EXPECT_EQ(sumCells, Tree::parent(Tree::firstCellOfLevel(l+1)));
        EXPECT_EQ(sumCells, Tree::firstChild(Tree::firstCellOfLevel(l-1)));
        EXPECT_EQ(sumCells, Tree::firstCellOfLevel(l+1) - Tree::numCellsInLevel(l));
        sumCells += Tree::numCellsInLevel(l);
        for(auto i = Tree::firstCellOfLevel(l); i<Tree::firstCellOfLevel(l+1); ++i) {
            numChildren.at(Tree::parent(i)) += 1; // use at to get bounds check
        }
    }
    EXPECT_EQ(cells, sumCells);

    // check that all but the last level have correct number of children
    for(auto i=0; i<Tree::firstCellOfLevel(max_level); ++i)
        EXPECT_EQ(numChildren[i], Tree::numChildren());

    for(auto i=Tree::firstCellOfLevel(max_level); i<cells; ++i)
        EXPECT_EQ(numChildren[i], 0);
}

template<unsigned int D>
void testCoordMapping(const unsigned max_level) {
    using Tree = NetworKit::girgs::SpatialTreeCoordinateHelper<D>;

    // generate some points and check their cells on all levels
    std::mt19937 gen(1337);
    std::uniform_real_distribution<> dist(0.0, 1.0);
    for(auto i=0; i<100; ++i) {

        // generate random point in [0..1)^D
        auto point = std::array<double, D>();
        for(auto d=0u; d<D; ++d)
            point[d] = dist(gen);

        // compute containing cell in all levels and check if point is in their bounds
        auto containingCells = std::vector<unsigned int>(max_level);
        containingCells[0] = 0;
        for(auto l=1u; l < max_level; ++l){
            const auto firstPoint = Tree::firstCellOfLevel(l);
            const auto cellPoint = Tree::cellForPoint(point, l);
            containingCells[l] = firstPoint + cellPoint;

            auto cell = containingCells[l];
            auto bounds = Tree::bounds(cell, l);
            for(auto d=0u; d<D; ++d){
                // point is in cell bounds
                EXPECT_LE(bounds[d].first, point[d])
                    << " l=" << l
                    << " point=" << stringify<double, D>(point)
                    << " cell=" << cell
                    << " bounds=" << stringify<double, D>(bounds)
                    << " firstPoint=" << firstPoint
                    << " cellPoint=" << cellPoint;

                EXPECT_LT(point[d], bounds[d].second)
                    << " l=" << l
                    << " point=" << stringify<double, D>(point)
                    << " cell=" << cell
                    << " bounds=" << stringify<double, D>(bounds)
                    << " firstPoint=" << firstPoint
                    << " cellPoint=" << cellPoint;
            }
        }

        // check that all containing cells have the same parents
        for(auto l=1; l < max_level; ++l)
            EXPECT_EQ(containingCells[l-1], Tree::parent(containingCells[l]));
    }
}


TEST_F(GirgsSpatialTreeCoordinateHelperGTest, testTreeStructure)
{
    testTreeStructure<1>(12);
    testTreeStructure<2>( 6);
    testTreeStructure<3>( 4);
    testTreeStructure<4>( 3);
    testTreeStructure<5>( 2);
}


TEST_F(GirgsSpatialTreeCoordinateHelperGTest, testCoordMapping)
{
    testCoordMapping<1>(12);
    testCoordMapping<2>( 6);
    testCoordMapping<3>( 4);
    testCoordMapping<4>( 3);
    testCoordMapping<5>( 2);
}


TEST_F(GirgsSpatialTreeCoordinateHelperGTest, testTouching)
{
    // TODO write better test for touching
    using Tree = NetworKit::girgs::SpatialTreeCoordinateHelper<1>;

    auto a = Tree::firstCellOfLevel(5);
    auto b = Tree::firstCellOfLevel(5) + Tree::numCellsInLevel(5) - 1;
    EXPECT_TRUE(Tree::touching(a,b,5));
}


TEST_F(GirgsSpatialTreeCoordinateHelperGTest, testDistance)
{
    // TODO write better tests for dist of cells
    {
        using Tree = NetworKit::girgs::SpatialTreeCoordinateHelper<1>;
        EXPECT_EQ(Tree::dist(10, 11, 3), 0);
        EXPECT_EQ(Tree::dist(10, 12, 3), (1.0 / 8) * 1);
        EXPECT_EQ(Tree::dist(10, 13, 3), (1.0 / 8) * 2);
        EXPECT_EQ(Tree::dist(10, 14, 3), (1.0 / 8) * 3);
    }

    {
        using Tree = NetworKit::girgs::SpatialTreeCoordinateHelper<2>;
        EXPECT_EQ(Tree::dist(10, 8, 2), (1.0 / 4) * 1);
        EXPECT_EQ(Tree::dist(10, 9, 2), 0);
        EXPECT_EQ(Tree::dist(10, 14, 2), (1.0 / 4) * 1);
        EXPECT_EQ(Tree::dist(10, 17, 2), (1.0 / 4) * 1);
    }
}
