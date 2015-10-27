/*
 * CycleDistribution.cpp
 *
 *  Created on: Jun 01, 2014
 *      Author: dhoske
 */

#include <cerrno>
#include <cstring>
#include <cfenv>
#include <xmmintrin.h>
#include <sys/mman.h>

#include "Config.h"
#include "../properties/ConnectedComponents.h"
#include "../io/METISGraphReader.h"
#include "../io/DGSReader.h"
#include "../io/MatrixMarketReader.h"
#include "../io/KONECTGraphReader.h"
#include "../dynamics/GraphEventProxy.h"
#include "../auxiliary/StringTools.h"
#include "SDDSolver.h"

using namespace std;

namespace NetworKit {
namespace SDD {

double residual(const Matrix& L, const Vector& x, const Vector& b) {
  /* Currently: relative residual */
  assert(L.numberOfColumns() == x.getDimension());
  assert(L.numberOfRows() == b.getDimension());
  return (L*x - b).length() / b.length();
}

void setBenchmarkConditions() {
  // FP -> avoid denormalized
  fesetround(FE_TONEAREST);
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);

  // Lock memory pages
//  if (mlockall(MCL_FUTURE) != 0) {
//    WARNF("Locking pages failed: %s", strerror(errno));
//  }

  // From outside:
  //   - disable frequency scaling
  //   - use a single(!) NUMA node
  //   - and fix thread affinity for single-threaded testing
}

double vectorSum(const Vector& b) {
  double sum = 0.0;
  b.forElements([&] (double elem) {
    sum += elem;
  });
  return sum;
}

unordered_map<node, node> inducedNames(const vector<node>& V) {
  unordered_map<node, node> to_new;
  for (index i = 0; i < V.size(); ++i) {
    to_new[V[i]] = i;
  }
  return to_new;
}

Graph inducedSubgraph(const Graph &G, const vector<node>& V) {
  Graph Gout(V.size(), true);
  auto to_new = inducedNames(V);

  /* Add edges to new graph. */
  for (node u: V) {
    assert(0 <= u && u < G.numberOfNodes());
    G.forEdgesOf(u, [&] (node u, node v, edgeweight w) {
      if (u < v && to_new.find(v) != to_new.end()) {
        Gout.addEdge(to_new[u], to_new[v], w);
      }
    });
  }

  return Gout;
}

vector<Graph> getComponents(const Graph& G) {
  vector<Graph> out;
  ConnectedComponents con(G);
  con.run();

  // Call inducedSubgraph on each component
  auto parts = con.getPartition();
  for (const auto& part : parts.getSubsets()) {
    vector<node> nodes(part.begin(), part.end());
    out.emplace_back(inducedSubgraph(G, nodes));
  }

  return move(out);
}

Vector inducedVector(const Vector& v, const vector<node>& I) {
  auto to_new = inducedNames(I);

  Vector vout(I.size(), 0.0f);
  for (int idx: I) {
    vout[to_new[idx]] = v[idx];
  }
  return vout;
}

Vector randZeroSum(const Graph& G, size_t seed) {
  RandomEngine rand(seed);
  auto rand_value = uniform_real_distribution<double>(-1.0, 1.0);
  ConnectedComponents con(G);
  count n = G.numberOfNodes();
  con.run();
  Partition comps = con.getPartition();

  /* Fill each component randomly such that its sum is 0 */
  Vector b(n, 0.0);
  for (int id : comps.getSubsetIds()) {
    auto indexes = comps.getMembers(id);
    assert(!indexes.empty());
    double sum = 0.0;
    for (auto entry : indexes) {
      b[entry] = rand_value(rand);
      sum += b[entry];
    }
    b[*indexes.begin()] -= sum;
  }
  return b;
}

Graph invertGraph(const Graph& G) {
  Graph Gout(G, true, G.isDirected());
  Gout.forEdges([&] (node u, node v, edgeweight w) {
    Gout.setWeight(u, v, 1.0 / w);
  });
  return Gout;
}

Graph readGraph(const string& path) {
  using namespace Aux::StringTools;
  if (ends_with(path, ".graph")) {
    METISGraphReader reader;
    return reader.read(path);
  } else if (ends_with(path, ".dgs")) {
    Graph G(0, true);
    GraphEventProxy Gproxy(G);
    DGSReader reader;
    reader.read(path, Gproxy);
    return G;
  } else if (ends_with(path, ".mtx")) {
//    MatrixMarketReader reader;
//    Matrix M = reader.read(path);
//    if (!isSDD(M)) {
//      throw std::runtime_error(path + " does not contain a SDD matrix");
//    }
//    if (isLaplacian(M)) {
//      return sddToLaplacian(M);
//    } else {
//      return laplacianToGraph(M);
//    }
  } else if (ends_with(path, ".tsv")) {
    KONECTGraphReader reader;
    return reader.read(path);
  } else {
    throw std::runtime_error("unknown graph format " + path);
  }
}

Graph gen2DGrid(count n) {
  Graph G(n*n);
  for (index i = 0; i < n; ++i) {
    for (index j = 0; j < n; ++j) {
      if (i < n-1) {
        G.addEdge(i*n + j, (i+1)*n + j);
      }
      if (j < n-1) {
        G.addEdge(i*n + j, i*n + (j+1));
      }
    }
  }

  return G;
}

}
}
