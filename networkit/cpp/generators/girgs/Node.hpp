/*
 * Node.h
 *
 *  Created on: 03. May 2019
 *      Author: Christopher Weyand <Christopher.Weyand@hpi.de>, Manuel Penschuck <networkit@manuel.jetzt>
 *
 * Code is adopted from https://github.com/chistopher/girgs
 */

#ifndef GENERATORS_GIRGS_NODE_H_
#define GENERATORS_GIRGS_NODE_H_

#include <array>
#include <algorithm>
#include <vector>
#include <cassert>

namespace NetworKit {
namespace girgs {

template<unsigned int D>
struct Node {
    std::array<double, D>   coord;
    double                  weight;
    int                     index;
    unsigned int            cell_id;

    Node() {}; // prevent default values

    Node(const std::vector<double>& _coord, double weight, int index, int cell_id = 0)
        : weight(weight), index(index), cell_id(cell_id)
    {
        assert(_coord.size()==D);
        std::copy_n(_coord.cbegin(), D, coord.begin());
    }

    double distance(const Node& other) const {
        auto result = 0.0;
        for(auto d=0u; d<D; ++d){
            auto dist = std::abs(coord[d] - other.coord[d]);
            dist = std::min(dist, 1.0-dist);
            result = std::max(result, dist);
        }
        return result;
    }

    void prefetch() const noexcept {
#if defined(__GNUC__) || defined(__clang__)
        __builtin_prefetch(coord.data(), 0);
        __builtin_prefetch(&index, 0);
#endif
    }
};


} // namespace girgs
} // namespace NetworKit

#endif // GENERATORS_GIRGS_NODE_H_
