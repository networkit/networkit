/*
 * SpanningTree.cpp
 *
 *  Created on: 03 May 2014
 *      Author: dhoske
 */

#include <unordered_set>
#include <algorithm>
#include <numeric>
#include <functional>

#include "Config.h"
#include "SpanningTree.h"
#include "../structures/UnionFind.h"
#include "../graph/BFS.h"

using namespace std;

namespace NetworKit {
namespace SDD {

/**
 * Warning: we must always first invert the graphs since our
 * stretch-definition is in terms of the reciprocal weight.
 */

namespace {
  // Set weights of edges in T to corresponding weights in G
  // (needed because we may calculate with inverse weights)
  static void setWeights(const Graph& G, Graph& T) {
    T.forEdges([&] (node u, node v) {
      T.setWeight(u, v, G.weight(u, v));
    });
  }
}

/**
 * Computes the shortest-distance-tree of @a G rooted at @a u.
 */
RootedTree minDistanceST(const Graph& Gin, node u) {
  assert(isConnected(Gin));
  Graph G = invertGraph(Gin);

  /* Get shortest distance tree represented by prev pointers. */
  Dijkstra dijkstra(G, u);
  dijkstra.run();
  auto prev = dijkstra.getCanonicalPredecessors();

  /* Convert to rooted tree */
  Graph T(G.numberOfNodes(), true);
  G.forNodes([&] (node u) {
    if (prev[u] != none) {
      T.addEdge(u, prev[u]);
    }
  });
  setWeights(Gin, T);
  return RootedTree(T, u);
}

/**
 * Computes the minimum spanning tree of @a G rooted at @a u using Kruskal.
 */
RootedTree minWeightST(const Graph& Gin, node u) {
  assert(isConnected(Gin));
  Graph G = invertGraph(Gin);

  // Edge list sorted by weight
  vector<WeightedEdge> edges;
  G.forEdges([&] (node u, node v, edgeweight w) {
    edges.emplace_back(u, v, w);
  });
  sort(begin(edges), end(edges));

  // Build ST with union-find data structure
  Graph T(G.numberOfNodes(), true);
  UnionFind uf(G.numberOfNodes());
  for (const auto& edge: edges) {
    if (uf.find(edge.u) != uf.find(edge.v)) {
      uf.merge(edge.u, edge.v);
      T.addEdge(edge.u, edge.v, edge.weight);
    }
  }

  setWeights(Gin, T);
  return RootedTree(T, u);
}

/*** Low-Stretch ST based on "Using Petal-decompositions to Build a Low Stretch Spanning Tree" */
namespace {

using NodeMap = unordered_map<node, node>;
using NodeSet = unordered_set<node>;
using EdgeSet = unordered_set<Edge>;

/**
 * Sets the weight of forward-edges (from candidates) to zero.
 */
void zeroForward(Graph& G, NodeSet& candidates) {
    // Determine forward edges with multi-source SSSP.
    ConeSSSP sssp(G, vector<node>(candidates.begin(), candidates.end()));
    sssp.run();
    auto dists = sssp.getDistances();

    // Set forward distances to zero
    G.forEdges([&] (node u, node v, edgeweight w) {
      if (abs(dists[u] + w - dists[v]) < EPSILON) {
        G.setWeight(u, v, 0.0);
      }
    });
}

/**
 * Grows a ball from @a center in @a G until its radius is both
 * at least @a radius_low and \f$f(cost, volume)\f$ yields true, where
 * @a cost is the cost of the outgoing edges of the current ball
 * and @a volume is its current volume.
 */
template<typename F>
vector<node> growBall(const Graph& G, const NodeSet& explored, node center, edgeweight radius_low, F f) {
  // Grow ball and remember outgoing edges as well as their cost
  vector<node> out;

  edgeweight outgoing_cost = 0.0;
  EdgeSet outgoing_edges;
  count volume = 0;

  edgeweight last_dist = 0;
  IncompleteConeSSSP sssp(&G, {center}, &explored);
  while (sssp.hasNext()) {
    node u;
    edgeweight dist;
    tie(u, dist) = sssp.next();

    // If distance > radius_low, only add the node if
    // the cost of the boundary edges of the ball is still too big
    if (abs(dist - last_dist) > EPSILON && dist > radius_low + EPSILON) {
      if (f(outgoing_cost, volume)) {
        break;
      }
    }

    out.emplace_back(u);
    G.forNeighborsOf(u, [&] (node v, edgeweight w) {
      Edge e(u, v, true);
      if (outgoing_edges.find(e) != outgoing_edges.end()) {
        // Boundary edge no longer on boundary
        outgoing_edges.erase(e);
        outgoing_cost -= 1./w;
      } else {
        // New boundary edge
        outgoing_edges.insert(e);
        outgoing_cost += 1./w;
        volume++;
      }
    });

    last_dist = dist;
  }
  return move(out);
}

/**
 * Returns a cone from @a v with parameters in (lambda_low, lambda_high).
 *
 * G is modified (tricky and not thread-safe)
 */
vector<node> coneCut(const Graph& G, node center, edgeweight radius_low, edgeweight radius_high, NodeSet& explored, NodeSet& candidates) {
  assert(radius_low < radius_high);

  // Grow the cone
  edgeweight desired_cost = -1;
  count m = G.numberOfEdges();
  return growBall(G, explored, center, radius_low, [&] (edgeweight cost, count volume) {
    if (desired_cost == -1) {
      edgeweight mu;
      if (volume == 0) {
        mu = (volume + 1)*log2(m + 1);
      } else {
        mu = volume * log2(static_cast<edgeweight>(m) / volume);
      }
      desired_cost = mu / (radius_high - radius_low);
    }
    return cost < desired_cost + EPSILON;
  });
}

/**
 * Grows a ball from @a center with initial radius \f$\rho\cdot\delta\f$ until cost
 * constraint is fulfilled.
 */
vector<node> ballCut(const Graph& G, node center, edgeweight rho, edgeweight delta) {
  assert((1 - 2*delta) > 0 && rho > 0);

  // Grow the ball
  count m = G.numberOfEdges();
  edgeweight D = log2(m + 1) / ((1 - 2*delta) * rho);
  NodeSet explored;
  return growBall(G, explored, center, rho*delta, [&] (edgeweight cost, count volume) {
     return cost < (volume + 1) * D + EPSILON;
  });
}

/**
 * Returns a (delta-epsilon) star decomposition of G.
 *
 * BUG: As noted in the paper by P´al Andr´as Papp about low-stretch spanning
 * trees, the proof that the returned tree is actually a delta-epsilon star is flawed.
 */
void starDecomposition(const Graph& G, node center, edgeweight delta, edgeweight epsilon,
                       ConeSSSP& from_center /* input: needs to have been run */,
                       vector<vector<node>>& stars, vector<pair<node, node>>& bridges) {
  // Get necessary SSSP-information (dists and prev)
  // assert: from_center has already been run
  auto dists = from_center.getDistances();
  edgeweight radius = *max_element(begin(dists), end(dists));

  // Grow inner ball
  vector<node> ball = ballCut(G, center, radius, delta);

  // Now work on a directed copy of G
  Graph Gcopy(G, true, true);

  // Remember explored nodes
  NodeSet explored;
  for (node u : ball) {
    explored.insert(u);
  }

  // Candidates are ball shell around center
  vector<node> prev(G.numberOfNodes());
  NodeSet candidates;
  for (node u : ball) {
    G.forNeighborsOf(u, [&] (node v, edgeweight w) {
      if (explored.find(v) == explored.end() && abs(dists[u] + w - dists[v]) < EPSILON) {
        prev[v] = u;
        candidates.emplace(v);
      }
    });
  }

  // Set forward-edges to zero
  zeroForward(Gcopy, candidates);

  // Now build the remaining stars
  count ncovered = ball.size();
  stars.emplace_back(move(ball));
  bridges.emplace_back(center, none);
  while (!candidates.empty()) {
    // Choose an arbitrary candidate (this is "improved" in the follow-up paper)
    node candidate = *candidates.begin();
    vector<node> star = coneCut(Gcopy, candidate, 0.0, delta*epsilon/2.0, explored, candidates);
    assert(isConnected(inducedSubgraph(G, star)));

    // Remove covered candidates and store nodes used in star
    for (node u : star) {
      candidates.erase(u);
      explored.emplace(u);
    }
    ncovered += star.size();

    // Store computed star
    stars.emplace_back(move(star));
    bridges.emplace_back(candidate, prev[candidate]);
  }

  if (ncovered != G.numberOfNodes()) {
    throw std::runtime_error("Low-Stretch-ST implementation error: ncovered is bad");
  }
}

/**
 * Builds a spanning tree of @a G by cone decomposition (Elkin et. al
 * original algorithm).
 */
void recursiveConeDecomposition(const Graph& G, Graph& T, node center, const NodeMap& cur_to_orig, edgeweight beta) {
  assert(isConnected(G));
  count n = G.numberOfNodes();

  if (n <= LOW_STRETCH_SMALL) {
    G.forEdges([&] (node u, node v, edgeweight w) {
      T.addEdge(cur_to_orig.at(u), cur_to_orig.at(v), w); // unordered_map does not have const operator[]???
    });
  } else {
    ConeSSSP sssp(G, center);
    sssp.run();

    /** @todo Add edge contraction (necessary?) */
    //auto dists = sssp.getDistances();

    /* Call star decomposition. */
    vector<vector<node>> stars;       /* Partition into stars. */
    vector<pair<node, node>> bridges; /* (x, y)-bridges towards the stars */
    starDecomposition(G, center, 1./3., beta, sssp, stars, bridges);
    assert(stars.size() == bridges.size());
    assert(stars.size() > 0);

    /* Heuristic: if cannot subdivide anymore, use Kruskal. See P´al Andr´as Papp:
     * Low-Stretch Spanning Trees for the error in the paper that make this necessary. */
    if (stars.size() == 1) {
      auto Tout = minWeightST(G, center);
      Tout.forNodesPre([&] (node u) {
        if (u != center) {
          node v = Tout.getParent(u);
          T.addEdge(cur_to_orig.at(u), cur_to_orig.at(v), Tout.getWeight(v, u));
        }
      });
      return;
    }

    /* Add bridge edges to resulting tree (the first star doesn't have a bridge) */
    for (index i = 1; i < stars.size(); ++i) {
      node x, y;
      tie(x, y) = bridges[i];
      T.addEdge(cur_to_orig.at(x), cur_to_orig.at(y), G.weight(x, y));
    }

    /* Call decomposition hierarchially on each star */
    for (index i = 0; i < stars.size(); ++i) {
      /* Create subgraph with correct names, better: node attributes
         would make this much easier */
      NodeMap cur_to_new = inducedNames(stars[i]);
      NodeMap new_to_orig;
      node new_center = cur_to_new[bridges[i].first];
      for (const auto& mapping : cur_to_new) {
        new_to_orig[mapping.second] = cur_to_orig.at(mapping.first);
      }

      Graph Gsub = inducedSubgraph(G, stars[i]);
      recursiveConeDecomposition(Gsub, T, new_center, new_to_orig, beta);
    }
  }
}

}

RootedTree lowStretchST(const Graph& Gin, node u) {
  assert(isConnected(Gin));
  assert(u <= Gin.upperNodeIdBound());
  Graph G = invertGraph(Gin);

  count n = G.numberOfNodes();
  edgeweight beta = 1.0 / (2*ceil(log(2*n+32)/log(4./3.)));
  Graph T(n, true);  /* Output tree */

  /* Identity map of node names */
  NodeMap cur_to_orig;
  G.forNodes([&] (node u) {
    cur_to_orig[u] = u;
  });

  /* Now call the recursive decomposition */
  recursiveConeDecomposition(G, T, 0, cur_to_orig, beta);

  setWeights(Gin, T);
  if (!isTree(T)) {
    throw std::runtime_error("not a tree");
  }
  return RootedTree(T, u);
}

namespace {

/* Build ST of 2D-grid recursively. Current subgrid [x1, x2] x [y1, y2] */
void gridDecomposition(const Graph& G, Graph& T, count n, index x1, index x2, index y1, index y2) {
  /** Linearised index into grid */
  auto GI = [&] (index x, index y) {
    return y*n + x;
  };

  if (x1 == x2) {
    /* Path in y-direction */
    for (index y = y1; y < y2; ++y) {
      T.addEdge(GI(x1, y), GI(x1, y + 1));
    }
  } else if (y1 == y2) {
    /* Path in x-direction */
    for (index x = x1; x < x2; ++x) {
      T.addEdge(GI(x, y1), GI(x + 1, y1));
    }
  } else {
    /* U in the center */
    index xmid = (x1 + x2) / 2;
    index ymid = (y1 + y2) / 2;
    T.addEdge(GI(xmid, ymid), GI(xmid, ymid + 1));
    T.addEdge(GI(xmid, ymid + 1), GI(xmid + 1, ymid + 1));
    T.addEdge(GI(xmid + 1, ymid + 1), GI(xmid + 1, ymid));

    /* Recursive subdivision */
    gridDecomposition(G, T, n, x1, xmid, y1, ymid);
    gridDecomposition(G, T, n, xmid+1, x2, y1, ymid);
    gridDecomposition(G, T, n, x1, xmid, ymid+1, y2);
    gridDecomposition(G, T, n, xmid+1, x2, ymid+1, y2);
  }
}

}

RootedTree specialGridST(const Graph& G, node /* ignore */) {
  count n_sqr = G.numberOfNodes();
  count n = static_cast<count>(sqrt(n_sqr));
  Graph T(n_sqr, true);
  gridDecomposition(G, T, n, 0, n-1, 0, n-1);
  setWeights(G, T);
  return RootedTree(T, (n-1)/2*n+(n-1)/2);
}

}
}
