/*
 * StrongConnectedComponents.cpp
 *
 *  Created on: 01.06.2014
 *      Author: Klara Reichard (klara.reichard@gmail.com), Marvin Ritter (marvin.ritter@gmail.com)
 *
 *  [2016/07/14] Iterative variant avoiding stack-based recursion
 *  -- Obada Mahdi <omahdi@gmail.com>
 */

#include <stack>
#include <functional>
#include <tuple>

#include <networkit/components/StronglyConnectedComponents.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/auxiliary/Log.hpp>

namespace NetworKit {

StronglyConnectedComponents::StronglyConnectedComponents(const Graph& G, bool iterativeAlgo) : G(&G), iterativeAlgo(iterativeAlgo) {

}

void StronglyConnectedComponents::run() {
    if (iterativeAlgo) {
        runIteratively();
    } else {
        runRecursively();
    }
}

void StronglyConnectedComponents::runRecursively() {
    count z = G->upperNodeIdBound();
    component = Partition(z);

    index nextIndex = 0;
    std::vector<index> nodeIndex(z, none);
    std::vector<index> nodeLowLink(z, none);
    std::stack<node> stx;
    std::vector<bool> onStack(z, false);

    std::function<void(node)> strongConnect = [&](node v) {
        nodeIndex[v] = nextIndex++;
        nodeLowLink[v] = nodeIndex[v];
        stx.push(v);
        onStack[v] = true;

        G->forNeighborsOf(v, [&](node w) {
            if (nodeIndex[w] == none) {
                strongConnect(w);
                nodeLowLink[v] = std::min(nodeLowLink[v], nodeLowLink[w]);
            } else if (onStack[w]) {
                nodeLowLink[v] = std::min(nodeLowLink[v], nodeIndex[w]);
            }
        });

        if (nodeLowLink[v] == nodeIndex[v]) {
            component.toSingleton(v);
            while (true) {
                node w = stx.top();
                stx.pop();
                onStack[w] = false;
                if (w == v) {
                    break;
                }
                component[w] = component[v];
            }
        }
    };

    G->forNodes([&](node v) {
        if (nodeIndex[v] == none) {
            strongConnect(v);
        }
    });
}

void StronglyConnectedComponents::runIteratively() {
    count z = G->upperNodeIdBound();
    component = Partition(z);

    index nextIndex = 0;
    struct node_info {
        index i, lowLink;
    };
    std::vector<node_info> nodes(z, node_info {none, none});
    std::vector<node> stx {};
    stx.reserve(z);
    std::vector<bool> onStack(z, false);
    using state_type = typename std::tuple<node, node, int>;
    using container_type = typename std::vector<state_type>;
    //std::stack<state_type, container_type> dfss {};
    //
    // For performance tests: reserve enough memory to prevent online resizing
    // of dfss (memory-intensive operation!). Does not work easily with
    // std:stack wrapper, hence using std::vector directly.
    container_type dfss {};
    dfss.reserve(z);
    //auto max_stack_size = dfss.size();

    G->forNodes([&](node _v) {
        if (nodes[_v].i != none)
            return;
        DEBUG("DFS init: dfss.emplace(", _v, ", none, -1)");
        dfss.emplace_back(_v, none, -1L);
        while (!dfss.empty()) {
            //if (dfss.size() > max_stack_size)
            //	max_stack_size = dfss.size();
            node u, pred;
            int j;
            std::tie(u, pred, j) = dfss.back();
            DEBUG("(u=", u, ", pred=", pred, ", j=", j, ") = dfss.top()");
            dfss.pop_back();
            if (j == -1) {
                DEBUG("* new DFS branch at node u=", u, " (pred=", pred, "), index=", nextIndex);
                nodes[u].i = nextIndex;
                nodes[u].lowLink = nextIndex;
                nextIndex++;
                stx.push_back(u);
                onStack[u] = true;
                j = G->degreeOut(u);
                DEBUG("j <- ", j, ", nextIndex=", nextIndex);
            }
            if (j == 0) {
                DEBUG("< backtracking at u=", u, " (pred=", pred, ")");
                if (nodes[u].lowLink == nodes[u].i) {
                    component.toSingleton(u);
                    node w = none;
                    do {
                        w = stx.back();
                        stx.pop_back();
                        DEBUG(" - ", w, " <- stx.pop()");
                        onStack[w] = false;
                        component[w] = component[u];
                    } while (u != w);
                }
                if (pred != none) {
                    DEBUG("  ", pred, " -> ", u, ": testing pred lowlink (lowlink[pred]=", nodes[pred].lowLink, ", lowlink[u]=", nodes[u].lowLink, ")");
                    nodes[pred].lowLink = std::min(nodes[pred].lowLink, nodes[u].lowLink);
                }
            } else {
                dfss.emplace_back(u, pred, j-1);
                // getIthNeighbor() is not part of the public API; hacked in
                // manually to allow constant-time lookup of the i-th neighbor
                // (see networkit/cpp/graph/Graph.h)
                // Note: Neighbors in reverse order compared to original impl
                auto v = G->getIthNeighbor<true>(u, j-1);
                // Unoptimised version that works with the public API only:
                //auto v = G.neighbors(u)[G.degreeOut(u)-j];
                if (v == none) {
                    throw std::runtime_error("StronglyConnectedComponents::runIteratively(): unexpected value 'none' for neighbor, try the recursive implementation");
                }
                if (nodes[v].i == none) {
                    DEBUG("  ", u, " -> ", v, ": descending (current pred=", pred, ")");
                    dfss.emplace_back(v, u, -1L);
                } else if (onStack[v]) {
                    DEBUG("  ", u, " -> ", v, ": testing lowlink (lowlink=", nodes[u].lowLink, ", index at arc head is ", nodes[v].i, ")");
                    nodes[u].lowLink = std::min(nodes[u].lowLink, nodes[v].i);
                }
            }
        }
    });
    //DEBUG("max_stack_size = ", max_stack_size, ", node count = ", z);
}

Partition StronglyConnectedComponents::getPartition() {
    return this->component;
}

count StronglyConnectedComponents::numberOfComponents() {
    return this->component.numberOfSubsets();
}

count StronglyConnectedComponents::componentOfNode(node u) {
    assert (component[u] != none);
    return component[u];
}

}
