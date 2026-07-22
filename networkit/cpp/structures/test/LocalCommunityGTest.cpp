#include <set>
#include <unordered_map>
#include <utility>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/LocalCommunity.hpp>

namespace NetworKit {
namespace {

using ::testing::Pair;
using ::testing::UnorderedElementsAre;

template <bool ShellMaintainsExtDeg, bool MaintainBoundary, bool AllowRemoval>
struct LocalCommunityConfig {
    static constexpr bool shellMaintainsExtDeg = ShellMaintainsExtDeg;
    static constexpr bool maintainBoundary = MaintainBoundary;
    static constexpr bool allowRemoval = AllowRemoval;

    using Community = LocalCommunity<ShellMaintainsExtDeg, MaintainBoundary, AllowRemoval>;
};

struct DegSnapshot {
    double intDeg;
    double extDeg;
};

bool operator==(const DegSnapshot &lhs, const DegSnapshot &rhs) {
    return lhs.intDeg == rhs.intDeg && lhs.extDeg == rhs.extDeg;
}

template <typename LocalCommunity, typename F>
auto collectShell(LocalCommunity &community, F getValue) {
    using Value = decltype(getValue(std::declval<const typename LocalCommunity::ShellInfo &>()));
    std::unordered_map<node, Value> result;
    community.forShellNodes([&](node u, const auto &info) { result.emplace(u, getValue(info)); });
    return result;
}

template <typename LocalCommunity, typename F>
auto collectCommunity(LocalCommunity &community, F getValue) {
    using Value =
        decltype(getValue(std::declval<const typename LocalCommunity::CommunityInfo &>()));
    std::unordered_map<node, Value> result;
    community.forCommunityNodes(
        [&](node u, const auto &info) { result.emplace(u, getValue(info)); });
    return result;
}

Graph weightedGraph() {
    Graph G(6, true, false);
    G.addEdge(0, 1, 2.0);
    G.addEdge(1, 2, 3.0);
    G.addEdge(1, 3, 5.0);
    G.addEdge(2, 3, 7.0);
    G.addEdge(3, 4, 11.0);
    G.addEdge(0, 5, 13.0);
    return G;
}

template <typename Config>
class LocalCommunityGTest : public ::testing::Test {};

using LocalCommunityConfigs = ::testing::Types<
    LocalCommunityConfig<false, false, false>, LocalCommunityConfig<true, false, false>,
    LocalCommunityConfig<true, true, false>, LocalCommunityConfig<false, false, true>,
    LocalCommunityConfig<true, false, true>, LocalCommunityConfig<true, true, true>>;

TYPED_TEST_SUITE(LocalCommunityGTest, LocalCommunityConfigs);

TYPED_TEST(LocalCommunityGTest, testConstructorRejectsDirectedGraphs) {
    Graph G(2, false, true);

    EXPECT_THROW(
        { [[maybe_unused]] typename TypeParam::Community community(G); }, std::runtime_error);
}

TYPED_TEST(LocalCommunityGTest, testAddNodeMaintainsCommunityShellAndCut) {
    Graph G = weightedGraph();
    typename TypeParam::Community community(G);

    community.addNode(1);

    EXPECT_TRUE(community.contains(1));
    EXPECT_EQ(community.toSet(), std::set<node>({1}));
    EXPECT_EQ(community.internalEdgeWeight(), 0.0);
    EXPECT_EQ(community.cut(), 10.0);
    EXPECT_THAT(collectShell(community, [](const auto &info) { return *info.intDeg; }),
                UnorderedElementsAre(Pair(0, 2.0), Pair(2, 3.0), Pair(3, 5.0)));

    if constexpr (TypeParam::shellMaintainsExtDeg) {
        EXPECT_THAT(
            collectShell(community,
                         [](const auto &info) { return DegSnapshot{*info.intDeg, *info.extDeg}; }),
            UnorderedElementsAre(Pair(0, DegSnapshot{2.0, 13.0}), Pair(2, DegSnapshot{3.0, 7.0}),
                                 Pair(3, DegSnapshot{5.0, 18.0})));
    }

    if constexpr (TypeParam::maintainBoundary) {
        EXPECT_EQ(community.boundarySize(), 1u);
    }

    community.addNode(2);
    community.addNode(3);

    EXPECT_EQ(community.toSet(), std::set<node>({1, 2, 3}));
    EXPECT_EQ(community.internalEdgeWeight(), 15.0);
    EXPECT_EQ(community.cut(), 13.0);
    EXPECT_THAT(collectShell(community, [](const auto &info) { return *info.intDeg; }),
                UnorderedElementsAre(Pair(0, 2.0), Pair(4, 11.0)));

    if constexpr (TypeParam::shellMaintainsExtDeg) {
        EXPECT_THAT(
            collectShell(community,
                         [](const auto &info) { return DegSnapshot{*info.intDeg, *info.extDeg}; }),
            UnorderedElementsAre(Pair(0, DegSnapshot{2.0, 13.0}), Pair(4, DegSnapshot{11.0, 0.0})));
    }

    if constexpr (TypeParam::maintainBoundary) {
        EXPECT_EQ(community.boundarySize(), 2u);
        EXPECT_THAT(collectShell(community, [](const auto &info) { return info.boundaryChange(); }),
                    UnorderedElementsAre(Pair(0, 0), Pair(4, -1)));
    }
}

template <typename Config>
class RemovableLocalCommunityGTest : public ::testing::Test {};

using RemovableLocalCommunityConfigs = ::testing::Types<LocalCommunityConfig<false, false, true>,
                                                        LocalCommunityConfig<true, false, true>,
                                                        LocalCommunityConfig<true, true, true>>;

TYPED_TEST_SUITE(RemovableLocalCommunityGTest, RemovableLocalCommunityConfigs);

TYPED_TEST(RemovableLocalCommunityGTest, testRemoveNodeMaintainsCommunityShellAndCut) {
    Graph G = weightedGraph();
    typename TypeParam::Community community(G);

    community.addNode(1);
    community.addNode(2);
    community.addNode(3);
    community.removeNode(3);

    EXPECT_EQ(community.toSet(), std::set<node>({1, 2}));
    EXPECT_TRUE(!community.contains(3));
    EXPECT_EQ(community.internalEdgeWeight(), 3.0);
    EXPECT_EQ(community.cut(), 14.0);
    EXPECT_THAT(collectShell(community, [](const auto &info) { return *info.intDeg; }),
                UnorderedElementsAre(Pair(0, 2.0), Pair(3, 12.0)));

    if constexpr (TypeParam::shellMaintainsExtDeg) {
        EXPECT_THAT(
            collectShell(community,
                         [](const auto &info) { return DegSnapshot{*info.intDeg, *info.extDeg}; }),
            UnorderedElementsAre(Pair(0, DegSnapshot{2.0, 13.0}),
                                 Pair(3, DegSnapshot{12.0, 11.0})));
    }

    if constexpr (TypeParam::maintainBoundary) {
        EXPECT_EQ(community.boundarySize(), 2u);
    }
}

TYPED_TEST(RemovableLocalCommunityGTest, testCommunityNodeInfoSupportsRemovalScoring) {
    Graph G = weightedGraph();
    typename TypeParam::Community community(G);

    community.addNode(1);
    community.addNode(2);
    community.addNode(3);

    EXPECT_THAT(collectCommunity(community, [](const auto &info) { return *info.intDeg; }),
                UnorderedElementsAre(Pair(1, 8.0), Pair(2, 10.0), Pair(3, 12.0)));

    if constexpr (TypeParam::shellMaintainsExtDeg) {
        EXPECT_THAT(collectCommunity(
                        community,
                        [](const auto &info) { return DegSnapshot{*info.intDeg, *info.extDeg}; }),
                    UnorderedElementsAre(Pair(1, DegSnapshot{8.0, 2.0}),
                                         Pair(2, DegSnapshot{10.0, 0.0}),
                                         Pair(3, DegSnapshot{12.0, 11.0})));
    }

    if constexpr (TypeParam::maintainBoundary) {
        EXPECT_THAT(
            collectCommunity(community, [](const auto &info) { return info.boundaryChange(); }),
            UnorderedElementsAre(Pair(1, 0), Pair(2, 0), Pair(3, 0)));
    }
}

} // namespace
} // namespace NetworKit
