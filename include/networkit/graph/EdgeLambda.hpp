#ifndef NETWORKIT_GRAPH_EDGE_LAMBDA_HPP_
#define NETWORKIT_GRAPH_EDGE_LAMBDA_HPP_

#include <cstddef>
#include <type_traits>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/FunctionTraits.hpp>

namespace NetworKit::GraphDetails {

// Return void when the callable does not have an Nth parameter.
template <class TraitsT, size_t N, bool InRange = (N < TraitsT::arity)>
struct SafeArg {
    using type = void;
};

// Enable this specialization only when N is in range.
template <class TraitsT, size_t N>
struct SafeArg<TraitsT, N, true> {
    using type = typename TraitsT::template arg<N>::type;
};

template <class F, size_t N>
using ArgType = std::decay_t<typename SafeArg<Aux::FunctionTraits<F>, N>::type>;

/* Compile-time callable dispatcher for the edge callbacks of the graph classes.
 *
 * Aux::FunctionTraits lets us inspect a callable's arity and parameter types
 * so we can invoke it with the best matching edge callback signature.
 * SafeArg prevents out-of-bounds access to Traits::arg<N> by returning void
 * when the callable has fewer than N + 1 parameters. ArgType normalizes the
 * inspected type with std::decay_t so that references and cv-qualifiers do
 * not matter. The comparison on the parameter types is needed because
 * EdgeWeightT and NodeT are implicitly convertible to each other, so checking
 * only that the call is well-formed is not enough.
 *
 * The parameters are const, so a callback that takes a non-const reference
 * does not compile instead of silently writing to a local copy.
 *
 * If no branch matches, the static_assert produces a clear compile-time error
 * instead of a long substitution failure.
 */
template <class NodeT, class EdgeWeightT, class F>
auto edgeLambda(F &f, const NodeT u, const NodeT v, const EdgeWeightT ew, const edgeid id) {
    constexpr size_t arity = Aux::FunctionTraits<F>::arity;

    // Most specific signatures first: several callable shapes overlap through
    // implicit conversions between NodeT, EdgeWeightT and edgeid. If NodeT and
    // EdgeWeightT are the same type, (u, v) takes precedence over (v, ew).
    if constexpr (arity >= 3 && std::is_same_v<ArgType<F, 2>, EdgeWeightT>
                  && std::is_same_v<ArgType<F, 3>, edgeid>) {
        return f(u, v, ew, id);
    } else if constexpr (arity >= 2 && std::is_same_v<ArgType<F, 2>, edgeid>
                         && std::is_same_v<ArgType<F, 1>, NodeT>) {
        return f(u, v, id);
    } else if constexpr (arity >= 2 && std::is_same_v<ArgType<F, 2>, EdgeWeightT>) {
        return f(u, v, ew);
    } else if constexpr (arity >= 1 && std::is_same_v<ArgType<F, 1>, NodeT>) {
        return f(u, v);
    } else if constexpr (arity >= 1 && std::is_same_v<ArgType<F, 1>, EdgeWeightT>) {
        return f(v, ew);
    } else if constexpr (std::is_same_v<ArgType<F, 0>, NodeT>) {
        return f(v);
    } else {
        static_assert(!std::is_same_v<F, F>,
                      "Your lambda does not support the required parameters or the parameters "
                      "have the wrong type.");
    }
}

} // namespace NetworKit::GraphDetails

#endif // NETWORKIT_GRAPH_EDGE_LAMBDA_HPP_
