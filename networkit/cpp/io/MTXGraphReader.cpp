#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/StringTools.hpp>
#include <networkit/io/MTXGraphReader.hpp>
#include <networkit/io/MTXParser.hpp>

namespace NetworKit {

Graph MTXGraphReader::read(std::string_view path) {
    MTXParser parser(path.data());

    auto header = parser.getHeader();
    auto size = parser.getMatrixSize();
    auto weighted = true;
    auto symmetric = true;

    if (header.field == MTXParser::Field::Pattern)
        weighted = false;
    if (header.symmetry == MTXParser::Symmetry::General)
        symmetric = false;

    Graph G(std::max(size.rows, size.columns), weighted, !symmetric);
    std::string graphName =
        Aux::StringTools::split(Aux::StringTools::split(path, '/').back(), '.').front();

    auto current_edge = parser.getNext(weighted);

    while (current_edge.has_value()) {
        WeightedEdge e = current_edge.value();
        G.addEdge(e.u, e.v, e.weight);
        current_edge = parser.getNext(weighted);
    }
    return G;
} // namespace NetworKit

} /* namespace NetworKit */
