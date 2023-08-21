#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/StringTools.hpp>
#include <networkit/io/MTXGraphReader.hpp>
#include <networkit/io/MTXParser.hpp>

namespace NetworKit {

Graph MTXGraphReader::read(const std::string &path) {
    MTXParser parser(path);

    auto header = parser.getHeader();
    auto size = parser.getMatrixSize();
    auto weighted = true;
    auto symmetric = true;

    if (header.field == MTXParser::Field::Pattern)
        weighted = false;
    if (header.symmetry == MTXParser::Symmetry::General)
        symmetric = false;

    Graph G(std::max(size.rows, size.columns), weighted);
    std::string graphName =
        Aux::StringTools::split(Aux::StringTools::split(path, '/').back(), '.').front();

    std::optional<MTXParser::Edge> current_edge = parser.getNext(weighted);

    while (current_edge.has_value()) {
        const auto [from, to, weight] = current_edge.value();

        if (weighted && symmetric) {
            G.addPartialEdge(unsafe, from, to, *weight);
            G.addPartialEdge(unsafe, to, from, *weight);
        } else if (weighted && !symmetric) {
            G.addPartialEdge(unsafe, from, to, *weight);
        } else if (!weighted && symmetric) {
            G.addPartialEdge(unsafe, from, to);
            G.addPartialEdge(unsafe, to, from);
        } else {
            G.addPartialEdge(unsafe, from, to);
        }

        current_edge = parser.getNext(weighted);
    }
    return G;
} // namespace NetworKit

} /* namespace NetworKit */
