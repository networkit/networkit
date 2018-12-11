/*
 * TopHarmonicCloseness.cpp
 *
 * Created on: 25.02.2018
 *		 Author: nemes, Eugenio Angriman
 */

#include <cmath>
#include <omp.h>

#include "../auxiliary/PrioQueue.h"
#include "../components/StronglyConnectedComponents.h"
#include "TopHarmonicCloseness.h"

namespace NetworKit {

TopHarmonicCloseness::TopHarmonicCloseness(const Graph &G, count k,
                                           bool useBFSbound)
    : G(G), k(k), useBFSbound(useBFSbound),
      allScores(G.upperNodeIdBound(), std::numeric_limits<edgeweight>::max()),
      isExact(G.upperNodeIdBound(), false),
      isValid(G.upperNodeIdBound(), false),
      cutOff(G.upperNodeIdBound(), std::numeric_limits<edgeweight>::max()),
      exactCutOff(G.upperNodeIdBound(), 0), component(G.upperNodeIdBound()),
      rOld(G.upperNodeIdBound()) {}

std::pair<edgeweight, bool>
TopHarmonicCloseness::BFScut(node v, edgeweight x, count n, count r,
                             std::vector<uint8_t> &visited,
                             std::vector<count> &distances,
                             std::vector<node> &pred, count &visitedEdges) {

  count d = 0;
  int64_t gamma = 0, nd = 0;
  distances[v] = 0;
  std::queue<node> Q;
  std::queue<node> toReset;
  Q.push(v);
  toReset.push(v);
  nd++;

  visited[v] = true;

  /* Clean up the visited vector for the next run.
   * We don't need to do this for the distances and pred vector
   * because those are only read after visited[u] is true.
   */
  auto cleanup = [&]() {
    while (!toReset.empty()) {
      node w = toReset.front();
      toReset.pop();
      visited[w] = false;
    }
  };

  edgeweight c = 0;
  edgeweight ctilde = edgeweight(n - 1);

  do {
    node u = Q.front();
    Q.pop();

    if (distances[u] > d) {
      d = d + 1;

      int64_t d2 = int64_t(r) - nd;
      ctilde = c + (edgeweight(gamma) / ((d + 2) * (d + 1))) +
               (edgeweight(d2) / (d + 2));

      if (ctilde < x) {
        exactCutOff[v] = true;
        cutOff[v] = d;
        cleanup();
        return std::make_pair(ctilde, false);
      }

      gamma = 0;
    }
    bool cont = true;
    G.forNeighborsOf(u, [&](node w) {
      if (cont) {
        visitedEdges++;
        if (!visited[w]) {
          distances[w] = distances[u] + 1;

          Q.push(w);
          toReset.push(w);
          visited[w] = true;

          c += 1.0 / distances[w];

          if (!G.isDirected()) {
            gamma += (G.degree(w) - 1); // notice: only because it's undirected
          } else {
            gamma += G.degree(w);
          }

          ++nd;
          pred[w] = u;
        } else if (distances[w] > 1 && (pred[u] != w)) {
          ctilde = ctilde - 1.0 / (d + 1) + 1.0 / (d + 2);

          if (ctilde < x) {
            cont = false;
          }
        }
      }
    });
    if (ctilde < x) {
      exactCutOff[v] = false;
      cutOff[v] = d;
      cleanup();

      return std::make_pair(ctilde, false);
    }
  } while (!Q.empty());

  cleanup();
  exactCutOff[v] = false;
  cutOff[v] = std::numeric_limits<edgeweight>::max();

  return std::make_pair(c, true);
}

void TopHarmonicCloseness::BFSbound(node source, std::vector<double> &S2,
                                    count *visEdges) {
  count r = 0;
  std::vector<std::vector<node>> levels(n);
  // nodesPerLev[i] contains the number of nodes in level i
  std::vector<count> nodesPerLev(n, 0);
  // sumLevs[i] contains the sum of the nodes in levels j <= i
  std::vector<count> sumLevs(n, 0);
  count nLevs = 0;
  levels[nLevs].clear();
  double sum_dist = 0;
  edgeweight level_bound;

  auto inverseDistance = [&](edgeweight dist) { return 1.0 / dist; };

  G.BFSfrom(source, [&](node u, count dist) {
    sum_dist += dist > 0 ? inverseDistance(dist) : 0;

    r++;
    if (dist > nLevs) {
      sumLevs[nLevs] += nodesPerLev[nLevs];
      sumLevs[nLevs + 1] = sumLevs[nLevs];
      nLevs++;
      levels[nLevs].clear();
    }
    levels[nLevs].push_back(u);
    nodesPerLev[nLevs]++;
  });
  sumLevs[nLevs] += nodesPerLev[nLevs];
  if (G.isDirected()) {
    (*visEdges) += G.numberOfEdges();
  } else {
    (*visEdges) += 2 * G.numberOfEdges();
  }

  S2[source] = sum_dist;

  // we compute the bound for the first level
  double closeNodes = 0, farNodes = 0;
  for (count j = 0; j <= nLevs; j++) {
    if (j <= 2) {
      closeNodes += nodesPerLev[j];
    } else {
      farNodes +=
          nodesPerLev[j] * inverseDistance(double(std::abs((double)j - 1.)));
    }
  }

  level_bound = inverseDistance(2.) * closeNodes + farNodes;

  for (count j = 0; j < levels[1].size(); j++) {
    node w = levels[1][j];
    // we subtract 2 not to count the node itself
    double bound = (level_bound - inverseDistance(2.) +
                    (inverseDistance(1.) - inverseDistance(2.)) * G.degree(w));

    if (bound < S2[w] &&
        (!G.isDirected() || component[w] == component[source])) {
      S2[w] = bound;
    }
  }

  // now we compute it for the other levels
  for (count i = 2; i <= nLevs; i++) {

    level_bound = 0;
    // TODO: OPTIMIZE?
    if (!G.isDirected()) {
      for (count j = 0; j <= nLevs; j++) {
        level_bound += inverseDistance(std::max(
                           2., double(std::abs((double)j - (double)i)))) *
                       nodesPerLev[j];
      }
    } else {
      for (count j = 0; j <= nLevs; j++) {
        level_bound += inverseDistance(std::max(2., (double)j - (double)i)) *
                       nodesPerLev[j];
      }
    }

    for (count j = 0; j < levels[i].size(); ++j) {
      node w = levels[i][j];
      double bound =
          (level_bound - inverseDistance(2.) +
           (inverseDistance(1.) - inverseDistance(2.)) * G.degree(w));

      if (bound < S2[w] &&
          (!G.isDirected() || component[w] == component[source])) {
        S2[w] = bound;
      }
    }
  }
}

void TopHarmonicCloseness::init() {
  n = G.upperNodeIdBound();
  assert(n >= k);

  topk.clear();
  topk.resize(k);
  topkScores.clear();
  topkScores.resize(k);
  nMinCloseness = 0;
  minCloseness = std::numeric_limits<double>::max();
  trail = 0;
}

void TopHarmonicCloseness::run() {
  init();

  std::vector<bool> toAnalyze(n, true);

  // We compute the number of reachable nodes (or an upper bound)
  // only if we use the algorithm for complex networks.
  if (!useBFSbound) {
    computeReachableNodes();
  }

  // Main priority queue with all nodes in order of decreasing degree
  Aux::PrioQueue<edgeweight, node> Q1(n);

  G.forNodes([&](node v) { Q1.insert(n + G.degree(v), v); });

  Aux::PrioQueue<edgeweight, node> top(n);

  // protects accesses to all shared variables
  omp_lock_t lock;
  omp_init_lock(&lock);

  edgeweight kth = 0;
#pragma omp parallel
  {
    std::vector<uint8_t> visited(n, false);
    std::vector<count> distances(n);
    std::vector<node> pred(n, 0);

    std::vector<edgeweight> S(n, std::numeric_limits<edgeweight>::max());
    while (Q1.size() != 0) {

      omp_set_lock(&lock);
      if (Q1.size() == 0) { // The size of Q1 might have changed.
        omp_unset_lock(&lock);
        break;
      }
      std::pair<edgeweight, node> extracted = Q1.extractMin();

      node v = extracted.second;
      toAnalyze[v] = false;

      omp_unset_lock(&lock);

      // for networks with large diameters: break if the score of the
      // current node is smaller than the k-th highest score
      if (useBFSbound && allScores[v] < kth && isValid[v]) {
        break;
      }

      if (G.degreeOut(v) == 0) {
        omp_set_lock(&lock);
        allScores[v] = 0;
        isExact[v] = 0;
        omp_unset_lock(&lock);
      } else if (useBFSbound) {
        count visEdges;

        // Perform a complete BFS from v and obtain upper bounds
        // for all nodes in the graph
        BFSbound(v, S, &visEdges);

        omp_set_lock(&lock);
        allScores[v] = S[v];
        isExact[v] = true;
        isValid[v] = true;
        omp_unset_lock(&lock);

        // Update the scores of all nodes with the bounds obtained
        // by the complete BFS
        G.forNodes([&](node u) {
          if (allScores[u] > S[u] &&
              toAnalyze[u]) { // This part must be synchronized.
            omp_set_lock(&lock);
            if (allScores[u] > S[u] &&
                toAnalyze[u]) { // Have to check again, because the variables
              // might have changed
              allScores[u] = S[u];
              isValid[u] = true;
              Q1.remove(u);
              Q1.insert(allScores[u], u);
            }
            omp_unset_lock(&lock);
          }
        });
      } else {
        count edgeCount = 0;
        // perform a pruned BFS from v in complex networks
        std::pair<edgeweight, bool> c =
            BFScut(v, kth, n, r[v], visited, distances, pred, edgeCount);

        omp_set_lock(&lock);
        allScores[v] = c.first;
        isExact[v] = c.second;
        isValid[v] = true;
        omp_unset_lock(&lock);
      }
      // Insert v into the list with the k most central nodes if
      // its score is larger than the k-th largest value
      omp_set_lock(&lock);
      if (isExact[v] && allScores[v] >= kth) {
        top.insert(allScores[v], v);
        if (top.size() > k) {
          ++trail;
          if (allScores[v] > kth) {
            if (nMinCloseness == trail) {
              while (top.size() > k) {
                top.extractMin();
              }
              trail = 0;
              nMinCloseness = 1;
              if (k > 1) {
                Aux::PrioQueue<edgeweight, node> tmp(n);
                auto last = top.extractMin();
                auto next = top.extractMin();
                minCloseness = last.first;

                if (last.first == next.first) {
                  tmp.insert(last.first, last.second);
                  while (next.first == last.first) {
                    tmp.insert(next.first, next.second);
                    ++nMinCloseness;
                    if (top.size() == 0) {
                      break;
                    }
                    next = top.extractMin();
                  }
                  if (next.first != last.first) {
                    top.insert(next.first, next.second);
                  }

                  while (tmp.size() > 0) {
                    auto elem = tmp.extractMin();
                    top.insert(elem.first, elem.second);
                  }
                } else {
                  top.insert(next.first, next.second);
                  top.insert(last.first, last.second);
                }
              }
            }
          } else {
            ++nMinCloseness;
          }
        } else if (allScores[v] < minCloseness) {
          minCloseness = allScores[v];
          nMinCloseness = 1;
        } else if (allScores[v] == minCloseness) {
          ++nMinCloseness;
        }
      }

      // Update the k-th largest value for this thread
      if (top.size() >= k) {
        std::pair<edgeweight, node> elem = top.extractMin();
        kth = elem.first;
        top.insert(elem.first, elem.second);
        if (nMinCloseness == 1) {
          minCloseness = kth;
        }
      }
      omp_unset_lock(&lock);
    }
  }

  // TODO This could go to another method
  if (trail > 0) {
    topk.resize(k + trail);
    topkScores.resize(k + trail);
  }

  // Store the nodes and their closeness centralities
  for (int64_t j = top.size() - 1; j >= 0; --j) {
    std::pair<edgeweight, node> elem = top.extractMin();
    topk[j] = elem.second;
    topkScores[j] = elem.first;
  }

  for (count i = 0; i < topk.size() - 1; ++i) {
    count toSort = 1;
    while ((i + toSort) < topk.size() &&
           topkScores[i] == topkScores[i + toSort]) {
      ++toSort;
    }
    if (toSort > 1) {
      auto begin = topk.begin() + i;
      std::sort(begin, begin + toSort);
      i += toSort - 1;
    }
  }

  hasRun = true;
}

void TopHarmonicCloseness::computeReachableNodes() {
  if (G.isDirected()) {
    computeReachableNodesDirected();
  } else {
    computeReachableNodesUndirected();
  }
}

void TopHarmonicCloseness::computeReachableNodesUndirected() {
  ConnectedComponents comps(G);

  r = std::vector<count>(n);

  comps.run();
  std::map<index, count> sizes = comps.getComponentSizes();
  G.forNodes([&](node v) {
    index cv = comps.componentOfNode(v);
    component[v] = cv;
    r[v] = sizes[cv];
  });
}

void TopHarmonicCloseness::computeReachableNodesDirected() {
  r = std::vector<count>(n);
  StronglyConnectedComponents sccs(G);
  sccs.run();

  count N = sccs.numberOfComponents();
  std::vector<count> reachU_scc(N, 0);
  std::vector<count> reachU_without_max_scc(N, 0);
  std::vector<bool> reach_from_max_scc(N, false);
  std::vector<bool> reaches_max_scc(N, false);
  std::vector<std::vector<count>> sccs_vec(N, std::vector<count>());
  Graph sccGraph(N, false, true);
  std::vector<bool> found(N, false);
  count maxSizeCC = 0;
  // We compute the vector sccs_vec, where each component contains the list of
  // its nodes
  for (count v = 0; v < n; v++) {
    component[v] = sccs.componentOfNode(v);
    sccs_vec[sccs.componentOfNode(v) - 1].push_back(v);
  }

  // We compute the SCC graph and store it in sccGraph
  for (count V = 0; V < N; V++) {
    for (node v : sccs_vec[V]) {
      G.forNeighborsOf(v, [&](node w) {
        count W = sccs.componentOfNode(w) - 1;

        if (W != V && !found[W]) {
          found[W] = true;
          sccGraph.addEdge(V, W);
        }
      });
    }
    sccGraph.forNeighborsOf(V, [&](node W) { found[W] = false; });
    if (sccGraph.degreeOut(V) > sccGraph.degreeOut(maxSizeCC)) {
      maxSizeCC = V;
    }
    // ELISABETTA: maybe the code can be made simpler by running G.forEdges to
    // scan all the edges. Would it be better to have a Graph object to store
    // the SCC graph?
  } // MICHELE: I have used a graph instead of scc_adjlist. About G.forEdges,
  // I think it is
  // a bit more complicated: I have to scan nodes, otherwise I do not know how
  // to avoid multiple edges. This scan is made using variable "found". Do you
  // have better ideas? Note that this is linear in the graph size.

  // BFS from the biggest SCC.
  std::queue<count> Q;
  Q.push(maxSizeCC);
  reach_from_max_scc[maxSizeCC] = true;
  while (!Q.empty()) {
    count V = Q.front();
    Q.pop();
    reachU_scc[maxSizeCC] += sccs_vec[V].size();
    sccGraph.forNeighborsOf(V, [&](node W) {
      if (!reach_from_max_scc[W]) {
        reach_from_max_scc[W] = true;
        Q.push(W);
      }
    });
  }
  reaches_max_scc[maxSizeCC] = true;

  // so far only the largest SCC has reach_U and reach_L > 0

  // Dynamic programming to compute number of reachable vertices
  for (count V = 0; V < N; V++) {
    if (V == maxSizeCC) {
      continue;
    }
    sccGraph.forNeighborsOf(V, [&](node W) {
      if (!reach_from_max_scc[W]) {
        reachU_without_max_scc[V] += reachU_without_max_scc[W];
      }
      reachU_scc[V] += reachU_scc[W];
      reachU_scc[V] = std::min(reachU_scc[V], n);
      reaches_max_scc[V] = reaches_max_scc[V] || reaches_max_scc[W];
    });

    if (reaches_max_scc[V]) {
      reachU_scc[V] = reachU_without_max_scc[V] + reachU_scc[V];
    }
    reachU_scc[V] += sccs_vec[V].size();
    reachU_scc[V] = std::min(reachU_scc[V], n);
  }

  for (count v = 0; v < n; v++) {
    r[v] = reachU_scc[sccs.componentOfNode(v) - 1];
  }
}
} // namespace NetworKit
