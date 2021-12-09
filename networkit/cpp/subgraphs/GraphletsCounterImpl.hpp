#ifndef NETWORKIT_GRAPHLETS_COUNTER_IMPL_HPP_
#define NETWORKIT_GRAPHLETS_COUNTER_IMPL_HPP_

#include <set>

namespace NetworKit {
namespace SubgraphsDetails {

/*
 * Compute \binom n2 = \frac {n(n-1)}2 in O(1) time
 * */
inline static constexpr count binom2(count n) {
    return (n * (n-1)) / 2;
}

/*
 * Compute \binom n3 = \frac {n(n-1)(n-2)}6 in O(1) time
 * */
inline static constexpr count binom3(count n) {
    return (n * (n-1) * (n-2)) / 6;
}

/*
 * Compute \binom n4 = \frac {n(n-1)(n-2)(n-3)}{24} in O(1) time
 */
inline static constexpr count binom4(count n) {
    return (n * (n-1) * (n-2) * (n-3)) / 24;
}

enum class GraphletsSize {
    TWO,
    THREE,
    FOUR,
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
     *
     * Requires Theta(1) in time and space.
     *
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
        E,              // number of K_2
        binom2(V) - E,  // number of empty 2-graphs
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
     *
     * Requires O(E \Delta \log \Delta) time and O(V) space (for Delta
     * the maximum degree in the graph)
     *
     * @return A size 4 vector containing {f(3, 1), f(3, 2), f(3, 3), f(3, 4)}
     */
    std::vector<count> getCounts() const;
private:
    const Graph* G;
};

std::vector<count> GraphletsCounterImpl<GraphletsSize::THREE>::getCounts() const {
    count nb_triangle{0},
          nb_2_star{0},
          nb_single_edge{0},
          nb_empty{0};
    auto V{G->numberOfNodes()};
    G->parallelForEdges(  // l. 3
        [this, V, &nb_triangle, &nb_2_star, &nb_single_edge](node u, node v) {
            std::vector<short> X(V, 0);
            std::set<node> star_u, star_v, triangles;
            // number of vertices w adjacent to either u or v
            count nb_neighbours{2};  // start at 2 because we need to include u and v
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
    nb_triangle /= 3;
    nb_2_star /= 2;
    nb_empty = binom3(V);
    nb_empty -= nb_triangle + nb_2_star + nb_single_edge;
    return {
        nb_triangle,
        nb_2_star,
        nb_single_edge,
        nb_empty
    };
}

template <>
class GraphletsCounterImpl<GraphletsSize::FOUR> {
public:
    GraphletsCounterImpl(const Graph* G) noexcept: G{G} {
    }

    /*
     * Count the number of occurrences of the 3- and 4-graphlets in the input graph.
     *
     * Requires O(E (T + S) \Delta \log \Delta) time and O(V) space (for Delta
     * the maximum degree in the graph, T the maximum number of triangles
     * incident to an edge and S the maximum number of stars incident
     * to an edge).
     *
     * @return A size 4+11 = 15 vector containing (f(3, i))_i + (f(4, i))_i
     */
    std::vector<count> getCounts() const;
private:
    const Graph* G;

    // l. 39-45
    inline count cliqueCount(std::vector<short>& X,
                             const std::set<node>& triangles) const;
    // l. 46-52
    inline count cycleCount(std::vector<short>& X,
                            const std::set<node>& stars) const;
};


std::vector<count> GraphletsCounterImpl<GraphletsSize::FOUR>::getCounts() const {
    /* Consul the paper by Ahmed et al. for a graphical representation of these graphs */
    count nb_triangle{0},           // connected 3-graphlets
          nb_2_star{0},
          nb_3_single_edge{0},        // disconnected 3-graphlets
          nb_3_empty{0},
          nb_4_clique{0},           // connected 4-graphlets
          nb_4_chordalcycle{0},
          nb_4_tailedtriangle{0},
          nb_4_cycle{0},
          nb_3_star{0},
          nb_4_path{0},
          nb_4_single_triangle{0},  // disconnected 4-graphlets
          nb_4_single_star{0},
          nb_4_two_disconnected_edges{0},
          nb_4_single_edge{0},
          nb_4_empty{0};
    /* temporary variables used to enumerate the 4-graphlets.
     * The names are not really explicit by let's use the notations
     * of the paper. */
    count N_TT{0},          // N_{T,T}
          N_SuSv{0},        // N_{S_u, S_v}
          N_TSuwedgeSv{0},  // N_{T,S_u \lor S_v}
          N_SdotSdot{0},    // N_{S_\cdot,S_\cdot}
          N_TI{0},          // N_{T,I}
          N_SuwedgeSvI{0},  // N_{S_u \lor S_v,I}
          N_II{0},          // N_{I,I}
          N_II1{0};         // N_{I,I,1}
    auto V{G->numberOfNodes()};
    auto E{G->numberOfEdges()};
    count nb_neighbours{0};
    G->parallelForEdges(  // l. 5
        [&](node u, node v) {
            std::vector<short> X(V, 0);
            std::set<node> star_u, star_v, triangles;  // l. 6
            count nb_neighbours{2};
            G->forNeighborsOf(  // l. 8
                u,
                [&](node w) {
                    if(v == w)
                        return;
                    star_u.insert(w);
                    X[w] = 1;
                    ++nb_neighbours;
                }
            );
            G->forNeighborsOf(  // l. 10
                v,
                [&](node w) {
                    if(u == w)  // l. 11
                        return;
                    if(X[w] == 1) {  // l. 12
                        triangles.insert(w);
                        X[w] = 2;
                        star_u.erase(w);
                    } else {  // l. 15
                        star_v.insert(w);
                        X[w] = 3;
                        ++nb_neighbours;
                    }
                }
            );
            auto nb_non_neighbours{V - nb_neighbours};
            // l. 16
            #pragma omp atomic update
            nb_triangle += triangles.size();
            #pragma omp atomic update
            nb_2_star += star_u.size() + star_v.size();
            #pragma omp atomic update
            nb_3_single_edge += V - nb_neighbours;
            // l. 18
            #pragma omp atomic update
            nb_4_clique += cliqueCount(X, triangles);  // count update here
            // l. 19
            #pragma omp atomic update
            nb_4_cycle += cycleCount(X, star_u);       // and here
            // l. 21-32
            #pragma omp atomic update
            N_TT += binom2(triangles.size());
            #pragma omp atomic update
            N_SuSv += star_u.size() * star_v.size();
            #pragma omp atomic update
            N_TSuwedgeSv += triangles.size() * (star_u.size() + star_v.size());
            #pragma omp atomic update
            N_SdotSdot += binom2(star_u.size()) + binom2(star_v.size());
            #pragma omp atomic update
            N_TI += triangles.size() * nb_non_neighbours;
            #pragma omp atomic update
            N_SuwedgeSvI += (star_u.size() + star_v.size()) * nb_non_neighbours;
            #pragma omp atomic update
            N_II += binom2(nb_non_neighbours);
            #pragma omp atomic update
            N_II1 += E - G->degree(u) - G->degree(v) + 1;
        }
    );
    nb_triangle /= 3;
    nb_2_star /= 2;
    nb_3_empty = binom3(V);
    nb_3_empty -= nb_triangle + nb_2_star + nb_3_single_edge;
    nb_4_clique /= 6;  // each clique is counted 6 times (edge based algorithm)
    nb_4_cycle /= 4;  // each cycle is counted 4 times
    // l. 35
    nb_4_chordalcycle = N_TT - 6*nb_4_clique;                    // Lemma 1
    nb_4_path = N_SuSv - 4*nb_4_cycle;                           // Lemma 2
    nb_4_tailedtriangle = N_TSuwedgeSv/2 - 2*nb_4_chordalcycle;  // Lemma 3
    nb_3_star = (N_SdotSdot - nb_4_tailedtriangle) / 3;          // Lemma 4
    nb_4_single_triangle = (N_TI - nb_4_tailedtriangle) / 3;     // Lemma 5
    nb_4_single_star = N_SuwedgeSvI/2 - nb_4_path;               // Lemma 6
    // l. 37
    nb_4_two_disconnected_edges = N_II1/2
                                - 3 * nb_4_clique
                                - 2 * nb_4_chordalcycle
                                - nb_4_tailedtriangle
                                - 2 * nb_4_cycle
                                - nb_4_path;                     // Eq. 13
    nb_4_single_edge = N_II - 2*nb_4_two_disconnected_edges;     // Lemma 7
    nb_4_empty = binom4(V)
               - nb_4_clique
               - nb_4_chordalcycle
               - nb_4_tailedtriangle
               - nb_4_cycle
               - nb_3_star
               - nb_4_path
               - nb_4_single_triangle
               - nb_4_single_star
               - nb_4_two_disconnected_edges
               - nb_4_single_edge;                               // Eq. 14
    return {
        nb_triangle,
        nb_2_star,
        nb_3_single_edge,
        nb_3_empty,
        nb_4_clique,
        nb_4_chordalcycle,
        nb_4_tailedtriangle,
        nb_4_cycle,
        nb_3_star,
        nb_4_path,
        nb_4_single_triangle,
        nb_4_single_star,
        nb_4_two_disconnected_edges,
        nb_4_single_edge,
        nb_4_empty
    };
}

count GraphletsCounterImpl<GraphletsSize::FOUR>::cliqueCount(std::vector<short>& X,
                                                             const std::set<node>& triangles) const {
    count ret{0};
    for(node w : triangles) {
        G->forNeighborsOf(
            w,
            [&ret, &X](node r) {
                if(X[r] == 2)
                    ++ret;
            }
        );
        X[w] = 0;
    }
    return ret;
}

count GraphletsCounterImpl<GraphletsSize::FOUR>::cycleCount(std::vector<short>& X,
                                                            const std::set<node>& stars) const {
    count ret{0};
    for(node w : stars) {
        G->forNeighborsOf(
            w,
            [&ret, &X](node r) {
                if(X[r] == 3)
                    ++ret;
            }
        );
        X[w] = 0;
    }
    return ret;
}

}  // namespace NetworKit::SubgraphsDetails
}  // namespace NetworKit

#endif  // NETWORKIT_GRAPHLETS_COUNTER_IMPL_HPP_
