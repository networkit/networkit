/*
 * RootedTree.h
 *
 *  Created on: Apr 26, 2014
 *      Author: dhoske
 */

#ifndef ROOTED_TREE_H_
#define ROOTED_TREE_H_

#include <stack>
#include "../graph/Graph.h"

namespace NetworKit {
namespace SDD {

/**
 * Tests whether the graph @a G is connected.
 */
bool isConnected(const Graph& G);

/**
 * Tests whether the graph @a T is a tree.
 **/
bool isTree(const Graph& T);

/** @addtogroup sdd
 *  @{ */

/**
 * Rooted tree data structure with pre- and post-order traversals.
 *
 * @bug forNodesPostAbort, forNodesPost are susceptible to stack overflows
 * since they are implemented recursively.
 */
class RootedTree {
public:
  RootedTree() {
  }

  /**
   * @brief Creates a rooted tree from the undirected graph @a G rooted at @a root.
   * @param[in] G the graph representing the tree. It has to be connected and acyclic.
   * @param[in] root the node to root the tree at.
   */
  RootedTree(const Graph& G, node root);

  /** Data type for an edge from a node to its parent. */
  using Edge = std::pair<node, edgeweight>;

  /**
   * Returns the parent of @a u or -1 if @a u is the root of the tree.
   */
  node getParent(node u) const {
    return parent[u].first;
  }

  /**
   * Returns the list of children of @a u.
   */
  const std::vector<node>& getChildren(node u) const {
    return children[u];
  }

  /**
   * Returns the root of the tree.
   */
  node getRoot() const {
    return root;
  }

  /**
   * Returns whether the tree contains the edge \f$uv\f$.
   */
  bool hasEdge(node u, node v) const {
    return parent[u].first == v || parent[v].first == u;
  }

  /**
   * Returns the weight of the edge from(!) @a u to @a v.
   */
  edgeweight getWeight(node u, node v) const {
    assert(parent[v].first == u);
    return parent[v].second;
  }

  /**
   * Sets the weight of the edge from @a u to @a v to
   * @a weight.
   */
  void setWeight(node u, node v, edgeweight weight) {
    assert(parent[v].first == u);
    parent[v].second = weight;
  }

  /**
   * Returns the number of nodes of the tree.
   */
  count numberOfNodes() const {
    return children.size();
  }

  /**
   * Calls \f$f(u, v)\f$ on each edge \f$uv\f$.
   */
  template<typename F>
  void forEdges(const F& f) const {
    for (node u = 0; u < parent.size(); ++u) {
      if (u != root) {
        f(u, parent[u].first);
      }
    }
  }

  /**
   * Calls \f$f(u, v, w)\f$ for each edge \f$uv\f$
   * where \f$w\f$ is the weight of the edge.
   */
  template<typename F>
  void forWeightedEdges(const F& f) const {
    for (node u = 0; u < parent.size(); ++u) {
      if (u != root) {
        f(u, parent[u].first, parent[u].second);
      }
    }
  }

  /**
   * Recursively traverses the tree from @a u
   * with an Eulerian tour. The function @a f gets
   * the currently visited node and its depth.
   */
  template<typename F>
  void forNodesEulerian(const F& f) const {
    forNodesEulerian(f, getRoot());
  }
  template<typename F>
  void forNodesEulerian(const F& f, node u) const {
    /*
     * Recursive formulation:
     * f(u, depth);
     * for (const auto& child : children[u]) {
     *   forNodesEulerian(f, child, depth + 1);
     *   f(u, depth);
     * }
     */

    // Converted to iterative since this routine is susceptible to stack overflow.
    struct Pos {
      node u;      // node to process
      index idx;   // current index in node's children array
      count depth; // depth of the node
      Pos(node u, index idx, count depth) : u(u), idx(idx), depth(depth) {}
    };

    std::stack<Pos> st;
    st.emplace(u, 0, 0);
    while (!st.empty()) {
      Pos& p = st.top();
      f(p.u, p.depth);
      if (p.idx < children[p.u].size()) {
        st.emplace(children[p.u][p.idx++], 0, p.depth+1);
      } else {
        st.pop();
      }
    }
  }

  /**
   * Recursively traverses the tree from @a u
   * with post-order traversal and executes
   * @a f at each node. The function @a f gets
   * the currently visited node.
   */
  template<typename F>
  void forNodesPost(const F& f, node u) const {
    for (const auto& child: children[u]) {
      forNodesPost(f, child);
    }
    f(u);
  }

  template<typename F>
  void forNodesPost(const F& f) const {
    forNodesPost(f, getRoot());
  }

  /**
   * Recursively traverses the tree from @a u in pre-order traversal
   * and executes \f$f(v)\f$ for each visited node @a v.
   */
  template<typename F>
  void forNodesPre(const F& f, node u) const {
    /* Recursive formulation:
     * f(u);
     * for (const auto& child: children[u]) {
     *  forNodesPre(f, child);
     * }
     */

    std::stack<node> st;
    st.push(u);
    while (!st.empty()) {
      node cur = st.top();
      st.pop();
      f(cur);
      for (node child: children[cur]) {
        st.emplace(child);
      }
    }
  }

  /**
   * @see forNodesPre(const F&, node)
   * Pre-order traversal starting at the root.
   */
  template<typename F>
  void forNodesPre(const F& f) const {
    forNodesPre(f, getRoot());
  }

  /**
   * @see forNodesPre(const F&, node)
   * Pre-order traversal that does not expand the nodes marked
   * with true in @a is_leaf.
   */
  template<typename F>
  void forNodesPreAbort(const F& f, node u, const std::vector<bool>& is_leaf) const {
    std::stack<node> st;
    st.push(u);
    while (!st.empty()) {
      node cur = st.top();
      st.pop();
      f(cur);
      if (!is_leaf[cur]) {
        for (node child: children[cur]) {
          st.emplace(child);
        }
      }
    }
  }

  template<typename F>
  void forNodesPreAbort(const F& f, const std::vector<bool>& is_leaf) const {
    forNodesPreAbort(f, getRoot(), is_leaf);
  }

  /**
   * Traverses the tree from @a u with post-order traversal
   * and executes \f$f(v)\f$ for each visited node.
   */
  template<typename F>
  void forNodesPostAbort(const F& f, node u, const std::vector<bool>& is_leaf) const {
    if (!is_leaf[u]) {
      for (const auto& child: children[u]) {
        forNodesPostAbort(f, child, is_leaf);
      }
    }
    f(u);
  }

  template<typename F>
  void forNodesPostAbort(const F& f, const std::vector<bool>& is_leaf) const {
    forNodesPostAbort(f, getRoot(), is_leaf);
  }

  /**
   * Traverses the tree from @a u up to the root. The function
   * @a f gets every traversed edge \f$(u,v)\f$.
   *
   * @param[in] f function<void(node, node, edgeweight)> to execute for every edge
   * @param[in] u node to start at
   */
  template<typename F> void forEdgesUp(const F& f, node u) const {
    forEdgesUpUntil(f, u, getRoot());
  }
  template<typename F> void forEdgesUpUntil(const F& f, node u, node until) const {
    while (u != until) {
      const auto& edge = parent[u];
      f(u, edge.first, edge.second);
      u = edge.first;
    }
  }


  /**
   * Traverse nodes upwards from @a u and calls @a f with each node
   * until we get to @a until.
   */
  template<typename F> void forNodesUpUntil(const F& f, node u, node until) const {
    while (u != until) {
      f(u);
      u = parent[u].first;
    }
    f(until);
  }


  /**
   * Checks two rooted trees for equality. The
   * children of a node are considered to be unordered.
   */
  bool operator==(const RootedTree& that) const;

  /** Copyable and moveable */
  RootedTree& operator=(const RootedTree&) = default;
  RootedTree(const RootedTree&) = default;
  RootedTree& operator=(RootedTree&&) = default;
  RootedTree(RootedTree&&) = default;

private:
  /* Children at each node */
  std::vector<std::vector<node>> children;
  /* Single parent at each node */
  std::vector<Edge> parent;
  /* Root of the tree */
  node root;
};

/**
 * Computes the reciprocal weighted depths of the nodes in @a T.
 * By the reciprocal depth of a node @a u we mean
 * \f[\sum_{e\in P_{u}} 1/w_e \f]
 * where \f$P_{u}\f$ is the path from @a u to the root of @a T and
 * \f$w_e\f$ is the weight of the edge @a e.
 *
 * Needs once traversal of @a T.
 */
std::vector<edgeweight> computeReciprocalDepths(const RootedTree& T);

/**
 * For all edges @a e in @a edges, computes the resistance of the
 * unique path between the ends of @a e in @a T. Needs time \f$O(n\log(n)+m)\f$.
 */
std::vector<edgeweight> computeResistances(const RootedTree& T, const std::vector<Edge>& edges);

/** @} */

}
}

#endif
