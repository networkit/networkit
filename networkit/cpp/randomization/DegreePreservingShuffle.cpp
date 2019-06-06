#include <omp.h>

#include <algorithm>
#include <numeric>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/Random.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/randomization/DegreePreservingShuffle.hpp>

namespace NetworKit {

// The algorithms effectively computes for each degree (in the undirected case it's a scalar, in the directed case it's a scalar pair)
// the set of nodes which have this degree and randomly reassigns ids *within* these sets. Since directed and undirected degrees
// have a different signature, the algorithm itself is templated and compares degrees only using the < and the == operator.
// Those implementation details are hidden in the DegreePreservingShuffleDetails namespace:
namespace DegreePreservingShuffleDetails {
struct DirectedDegree {
    count in;
    count out;
    bool operator< (const DirectedDegree& o) const {return std::tie(in, out) <  std::tie(o.in, o.out);}
    bool operator==(const DirectedDegree& o) const {return std::tie(in, out) == std::tie(o.in, o.out);}
};

template <typename DegreeT>
struct NodeDegree {
    node id;
    union {
        DegreeT degree;
        node alt_id;
    };

    NodeDegree() {}
};

using UndirectedNodeDegree = NodeDegree<count>;
using DirectedNodeDegree = NodeDegree<DirectedDegree>;

template <typename DegreeT>
static std::vector<index> computePermutation(std::vector<NodeDegree<DegreeT>>& node_degrees) {
    const auto n = node_degrees.size();
    using node_degree_t = NodeDegree<DegreeT>;

    // Allocate memory for the permutation
    std::vector<index> permutation(n);
    #ifndef NDEBUG
    std::fill_n(permutation.begin(), n, std::numeric_limits<node>::max());
    #endif

    // We sort by degree and thereby move nodes of equal degrees into consecutive groups
    Aux::Parallel::sort(begin(node_degrees), end(node_degrees),
                        [] (const node_degree_t a, const node_degree_t b) {return a.degree < b.degree;});

    // Now shuffle each group
    count num_changes = 0;
    const auto requested_threads = std::max(1, std::min<int>(omp_get_max_threads(), n / 50000));
    #pragma omp parallel num_threads(requested_threads) reduction(+: num_changes)
    {
        const auto tid = omp_get_thread_num();
        const auto nthreads = omp_get_num_threads();

        const auto chunk_size = (n + nthreads - 1) / nthreads;
        auto chunk_begin = node_degrees.begin() + std::min<size_t>(chunk_size * tid, n);
        auto chunk_end   = node_degrees.begin() + std::min<size_t>(chunk_size * (tid + 1), n);

        using iter_t = decltype(chunk_begin);
        auto first_different = [] (iter_t begin, iter_t end, DegreeT degree) {
            return std::find_if_not(begin, end, [=] (node_degree_t& x) {return degree == x.degree;});
        };

        if (chunk_begin < chunk_end) {
            // find real start/end of chunk (i.e., consider groups which extends over a boundary)
            {
                // if a group started in the last chunk continues into this thread's chunk,
                // our first elements belong to the other thread and we have to skip over them.
                if (chunk_begin != node_degrees.begin()) {
                    chunk_begin = first_different(chunk_begin, chunk_end, std::prev(chunk_begin)->degree);
                }

                chunk_end = first_different(chunk_end, node_degrees.end(), std::prev(chunk_end)->degree);
            }

            // we will override the degree values in node_degrees. This barrier prevents races between
            // those writes and the previous search of chunk limits
            #pragma omp barrier

            Aux::SignalHandler handler;

            auto &prng = Aux::Random::getURNG();
            for(auto it = chunk_begin; it != chunk_end;) {
                handler.assureRunning();

                // find end of our current group
                auto group_end = first_different(it + 1, chunk_end, it->degree);

                // copy id (to make Fischer-Yates work)
                for(auto x = it; x != group_end; ++x)
                    x->alt_id = x->id;

                // go through group and randomly shuffle ids using a Fisher-Yates shuffle
                for(size_t len = std::distance(it, group_end) - 1; it != group_end; ++it, --len) {
                    const auto rand_idx = !len ? 0 : std::uniform_int_distribution<size_t>(0, len)(prng);

                    permutation[it->id] = it[rand_idx].alt_id;
                    it[rand_idx].alt_id = it->alt_id;

                    num_changes++;
                }
            }
        }
    };

    assert(num_changes == n);
    assert(!std::any_of(permutation.cbegin(), permutation.cend(),
                        [] (node i) {return i == std::numeric_limits<node>::max();}));
    assert(std::accumulate(permutation.cbegin(), permutation.cend(), 0llu) == static_cast<size_t>(n - 1) * (n) / 2);

    return permutation;
}
} // namespace DegreePreservingShuffleDetails


DegreePreservingShuffle::DegreePreservingShuffle(const Graph& G) : G{G} {}
DegreePreservingShuffle::~DegreePreservingShuffle() = default;

std::string DegreePreservingShuffle::toString() const {
    return "DegreePreservingShuffle";
}

void DegreePreservingShuffle::run() {
    const auto n = G.numberOfNodes();

    if (G.isDirected()) {
        // generate sequence of tuple (u, deg(u)) for each u
        std::vector<DegreePreservingShuffleDetails::DirectedNodeDegree> node_degrees(n);

        G.parallelForNodes([&](const node u) {
            node_degrees[u].id = u;
            node_degrees[u].degree = {G.degreeIn(u), G.degreeOut(u)};
        });

        permutation = DegreePreservingShuffleDetails::computePermutation(node_degrees);

    } else {
        // generate sequence of tuple (u, deg(u)) for each u
        std::vector<DegreePreservingShuffleDetails::UndirectedNodeDegree> node_degrees(n);

        G.parallelForNodes([&](const node u) {
            node_degrees[u].id = u;
            node_degrees[u].degree = G.degree(u);
        });

        permutation = DegreePreservingShuffleDetails::computePermutation(node_degrees);

    }
}

Graph DegreePreservingShuffle::getGraph() const {
    const auto n = G.numberOfNodes();
    assert(permutation.size() == n);

    return GraphTools::getRemappedGraph(G, n, [this] (node u) {
        return permutation[u];
    });
}

} // namespace NetworKit
