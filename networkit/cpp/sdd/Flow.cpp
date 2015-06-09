/*
 * Flow.cpp
 *
 *  Created on: May 03, 2014
 *      Author: dhoske
 */

#include <stack>

#include "Config.h"
#include "Flow.h"

using namespace std;

namespace NetworKit {
namespace SDD {

RootedTree initialFlow(const RootedTree& T, const Vector& b) {
  assert(abs(vectorSum(b)) < EPSILON);
  RootedTree out = T;

  /* For each node u store how much has to flow into the
   * tree rooted at u to get equilibrium. */
  vector<edgeweight> in_flow(b.getDimension(), 0);
  T.forNodesPost([&] (node u) {
    /* How much flow do the children consume? */
    edgeweight consumed = b[u];
    for (node v: T.getChildren(u)) {
      consumed += in_flow[v];
      out.setWeight(u, v, in_flow[v]);
    }
    in_flow[u] = consumed;
  });
  /* Flow must balance at the end -> root should require zero flow at the end. */
  assert(abs(in_flow[T.getRoot()]) < EPSILON);

  return out;
}

Flow::Flow(const RootedTree& T) : Tresistances(T) {
  /* Accumulate resistances in the tree. */
  Tresistances.forWeightedEdges([&] (node u, node v, edgeweight ew) {
    Tresistances.setWeight(v, u, 1.0 / ew);
  });
}

Flow::Flow(const RootedTree& T, const Vector& b) : Flow(T) {
  assert(abs(vectorSum(b)) < EPSILON);
  assert(T.numberOfNodes() == b.getDimension());
}

Vector Flow::primalSolution() const {
  /* Generic primalSolution() just needs to call getPotential() multiple times. */
  count n = Tresistances.numberOfNodes();
  node root = Tresistances.getRoot();
  Vector out(n, 0);
  for (node u = 0; u < n; ++u) {
    out[u] = getPotential(u, root);
  }
  return out;
}

void Flow::scaleResistances(edgeweight factor) {
  Tresistances.forWeightedEdges([&] (node u, node v, edgeweight w) {
    Tresistances.setWeight(v, u, w * factor);
  });
}

void LogFlow::scaleResistances(edgeweight) {
  throw std::runtime_error("not implemented");
}

TrivialFlow::TrivialFlow(const RootedTree& T, const Vector& b)
    : Flow(T, b), flowTree(initialFlow(T, b)) {
}

Vector TrivialFlow::primalSolution() const {
  /* Improved primalSolution() directly reads the solution from the flowTree */
  Vector out(flowTree.numberOfNodes());
  flowTree.forNodesPre([&] (node u) {
    if (u != flowTree.getRoot()) {
      node v = flowTree.getParent(u);
      out[u] = out[v] + flowTree.getWeight(v, u) * Tresistances.getWeight(v, u);
    }
  });
  return out;
}

LCAFlow::LCAFlow(const RootedTree& T, const Vector& b)
    : TrivialFlow(T, b), lca(T) {
}


void LogFlow::getSizes(node root, const vector<bool>& is_leaf, vector<int>& tree_size /* out */) const {
  /* Get tree sizes in single tree traversal */
  Tresistances.forNodesPostAbort([&] (node u) {
     tree_size[u] = 1;
     if (!is_leaf[u]) {
       for (node v: Tresistances.getChildren(u)) {
         tree_size[u] += tree_size[v];
       }
     }
  }, root, is_leaf);
}

node LogFlow::findSplit(const vector<int>& tree_size, Edge e, int nnodes) const {
  node split = e.v;

  /* The lower vertex is the split vertex in small trees. */
  if (nnodes == 1) {
    return e.u;
  }
  if (nnodes == 2) {
    return e.u == e.v ? Tresistances.getChildren(e.u).front() : e.v;
  }

  /* Walk down until each subtree has size <= n/2 */
  while (true) {
    auto& children = Tresistances.getChildren(split);

    assert(children.size() > 0);
    node max_sub = *max_element(children.begin(), children.end(), [&] (node u, node v) {
      return tree_size[u] < tree_size[v];
    });
    if (tree_size[max_sub] <= nnodes/2) {
      break;
    }
    split = max_sub;
  }

  return split;
}

vector<edgeweight> LogFlow::getHeights(const vector<bool>& is_leaf, node root, node split) const {
  /* Determine nodes on root-split path */
  vector<bool> on_split(Tresistances.numberOfNodes(), false);
  Tresistances.forEdgesUp([&] (node u, node v, edgeweight) {
    on_split[u] = on_split[v] = true;
  }, split);

  /* Get height in resistances */
  vector<edgeweight> heights(Tresistances.numberOfNodes());
  Tresistances.forNodesPreAbort([&] (node u) {
    if (u != Tresistances.getRoot()) {
      node v = Tresistances.getParent(u);
      heights[u] = heights[v];
      if (on_split[u]) {
        heights[u] += Tresistances.getWeight(v, u);
      }
    }
  }, root, is_leaf);
  return move(heights);
}

/* Adds dummy elements to @a vecs so that each vector in @a vecs
 * has a length divisible by UNROLL_LOOP */
template<typename T>
static void addDummy(vector<T>& vecs) {
  for (auto& vec : vecs) {
    count n = vec.size();
    if (n % UNROLL_LOOP != 0) {
      int nadd = UNROLL_LOOP - (n % UNROLL_LOOP);
      for (int i = 0; i < nadd; ++i) {
        vec.emplace_back(0);
      }
    }
    vec.shrink_to_fit();
  }
}

LogFlow::LogFlow(const RootedTree& T)
    : Flow(T),
      update_index(T.numberOfNodes()), update_weight(T.numberOfNodes()),
      query_index(T.numberOfNodes()), query_weight(T.numberOfNodes()) {
  /* Nodes that have already been used as splitting points are considered leaves. */
  vector<bool> is_leaf(Tresistances.numberOfNodes(), false);

  /* Get the subtree sizes. */
  vector<int> tree_size(Tresistances.numberOfNodes());
  getSizes(Tresistances.getRoot(), is_leaf, tree_size);

  /* Recursive decomposition: trees stores the starting edges of the trees that are still to be processed. */
  stack<Edge> trees;
  trees.push({Tresistances.getRoot(), Tresistances.getRoot()});
  while (!trees.empty()) {
    Edge e = trees.top(); trees.pop();
    TRACE("Subtree: ", e.u, "-", e.v);

    getSizes(e.v, is_leaf, tree_size);
    int n = e.u == e.v ? tree_size[e.u] : tree_size[e.v] + 1;

    node split = findSplit(tree_size, e, n);
    TRACE("Split vertex: ", split);
    auto heights = getHeights(is_leaf, e.v, split);

    /* Subdivide further if necessary */
    if (n > 2) {
      for (node next: Tresistances.getChildren(split)) {
        trees.push({split, next});
      }
      if (split != e.u) {
        trees.push(e);
      }
    }

    /* Add correct values to state, update & query vectors */
    int idx = int(state.size()); /* Current index in state array ~ index of subtree */
    state.emplace_back(0.0);
    state.emplace_back(0.0);

    /* Determine nodes in first tree T_0.
       Temporarily make split a leaf -> somewhat ugly. */
    vector<bool> in_first(Tresistances.numberOfNodes(), false);
    bool was_leaf = is_leaf[split];
    is_leaf[split] = true;
    Tresistances.forNodesPreAbort([&] (node u) { in_first[u] = true; }, e.v, is_leaf);
    is_leaf[split] = was_leaf;

    /* Set query and update vectors appropriately (somewhat complex logic) */
    Tresistances.forNodesPreAbort([&] (node u) {
      if (u != e.u) {
        if (in_first[u] && n != 2) {
          if (abs(heights[u]) != 0) {
            query_index[u].emplace_back(idx);
            query_weight[u].emplace_back(heights[u]);
          }
        } else {
          query_index[u].emplace_back(idx + 1);
          query_weight[u].emplace_back(1);
          update_index[u].emplace_back(idx);
          update_weight[u].emplace_back(1);
        }
        if (abs(heights[u]) != 0) {
          update_index[u].emplace_back(idx + 1);
          update_weight[u].emplace_back(heights[u]);
        }
      }
    }, e.v, is_leaf);
    is_leaf[split] = true;
  }

  addDummy(update_index);
  addDummy(update_weight);
  addDummy(query_index);
  addDummy(query_weight);
}

/* Tree decomposition for log time flow. */
LogFlow::LogFlow(const RootedTree& T, const Vector& b) : LogFlow(T) {
  if (b.getDimension() != T.numberOfNodes()) {
    throw std::invalid_argument("LogFlow: b has invalid size");
  }

  // Add initial flow (somewhat expensive)
  RootedTree TFlow = initialFlow(T, b);
  TFlow.forWeightedEdges([&] (node u, node v, edgeweight flow) {
    addFlow(u, v, flow);
  });
}

count LogFlow::getTime(node u) const {
  return update_index[u].size() + query_index[u].size();
}

}
}
