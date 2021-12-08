#ifndef NETWORKIT_GRAPHLETS_COUNTER_IMPL_HPP_
#define NETWORKIT_GRAPHLETS_COUNTER_IMPL_HPP_

#include <set>

namespace NetworKit {
namespace SubgraphsDetails {

enum class GraphletsSize {
    TWO,
    THREE,
};

/*
 * Should never be instantiated directly: always use explicit specialisations
 */
template <GraphletsSize K>
class GraphletsCounterImpl {
    static_assert(
        /* weird condition needed so that the assertion is not evaluated
         at declaration but at definition (and is therefore evaluated iff
         no explicit specialisation is found) */
        K == GraphletsSize::TWO and K != GraphletsSize::TWO,
        "Only explicit specialisations of GraphletsCounterImpl are available"
	);
};

/*
 * 2-graphlets counter. There only exists two unlabeled 2-graphs
 * (K_2 and the empty 2-graph). The number of occurrences of K_2 is simply
 * the number of edges while the number of occurrences of the empty 2-graph
 * is V(V-1)/2 - E
 */
template <>
class GraphletsCounterImpl<GraphletsSize::TWO> {
public:
    GraphletsCounterImpl(const Graph* G) noexcept: G{G} {
    }

    /*
     * Count the number of occurrences of the 2-graphlets in the input graph.
     * @return A size 2 vector containing {f(2, 1), f(2, 2)}
     */
    std::vector<count> getCounts() const;
private:
    const Graph* G;
};

std::vector<count> GraphletsCounterImpl<GraphletsSize::TWO>::getCounts() const {
    auto V{G->numberOfNodes()};
    auto E{G->numberOfEdges()};
    return {
        E,                  // number of K_2
        (V*(V-1)) / 2 - E,  // number of empty 2-graphs
    };
}

/*
 * Implementation of the Triad-Census algorithm implemented in PGD [1]
 *
 * 3-graphlets counter. There are four unlabeled 3-graphs: K_3 (E=3), the 2-star (E=2)
 * the single-edge graph (E=1) and the empty graph (E=0).
 *
 * [1] Ahmed, N. K., Neville, J., Rossi, R. A., & Duffield, N. (2015, November).
 * Efficient graphlet counting for large networks.
 * In 2015 IEEE International Conference on Data Mining (pp. 1-10). IEEE.
 */
template <>
class GraphletsCounterImpl<GraphletsSize::THREE> {
public:
    GraphletsCounterImpl(const Graph* G) noexcept: G{G} {
    }

    /*
     * Count the number of occurrences of the 3-graphlets in the input graph.
     * @return A size 4 vector containing {f(3, 1), f(3, 2), f(3, 3), f(3, 4)}
     */
    std::vector<count> getCounts() const;
private:
    const Graph* G;
};

std::vector<count> GraphletsCounterImpl<GraphletsSize::THREE>::getCounts() const {
    count nb_triangle,
          nb_2_star,
          nb_single_edge,
          nb_empty;
    auto V{G->numberOfNodes()};
    G->parallelForEdges(  // l. 3
        [this, V, &nb_triangle, &nb_2_star, &nb_single_edge](node u, node v) {
            std::vector<short> X(V, 0);
            std::set<node> star_u, star_v, triangles;
            count nb_neighbours{0};  // number of vertices w adjacent to either u or v
            G->forNeighborsOf(  // l. 5
                u,
                [v, &X, &star_u, &nb_neighbours](node w) {
                    if(v == w)  // l. 6
                        return;
                    // l. 7
                    star_u.insert(w);
                    X[w] = 1;
                    ++nb_neighbours;
                }
            );
            G->forNeighborsOf(  // l. 8
                v,
                [u, &X, &triangles, &star_u, &star_v, &nb_neighbours](node w) {
                    if(u == w)  // l. 9
                        return;
                    if(X[w] == 1) {  // l. 10
                        triangles.insert(w);
                        star_u.erase(w);
                    } else {  // l. 13
                        star_v.insert(w);
                        ++nb_neighbours;
                    }
                }
            );
            // l. 14-16
            #pragma omp atomic update
            nb_triangle += triangles.size();
            #pragma omp atomic update
            nb_2_star += star_u.size() + star_v.size();
            #pragma omp atomic update
            nb_single_edge += V - nb_neighbours;
        }
    );
    nb_empty = (V*(V-1)*(V-2)) / 6;
    nb_empty -= nb_triangle + nb_2_star + nb_single_edge;
    return {
        nb_triangle / 3,
        nb_2_star / 2,
        nb_single_edge,
        nb_empty
    };
}

}  // namespace NetworKit::SubgraphsDetails
}  // namespace NetworKit

#endif  // NETWORKIT_GRAPHLETS_COUNTER_IMPL_HPP_
