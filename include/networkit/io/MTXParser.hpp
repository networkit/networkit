#ifndef NETWORKIT_IO_MTX_PARSER_HPP_
#define NETWORKIT_IO_MTX_PARSER_HPP_

#include <fstream>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
/**
 * @ingroup io
 * Parser for the MTX file format.
 */
class MTXParser final {

    std::ifstream graphFile;
    std::string currentLine;

public:
    enum class Object { Matrix, Vector };
    enum class Format { Coordinate, Array };
    enum class Field { Real, Double, Complex, Integer, Pattern };
    enum class Symmetry { General, Symmetric, SkewSymmetric, Hermitian };

    std::unordered_map<std::string, Object> objectMap = {{"matrix", Object::Matrix},
                                                         {"vector", Object::Vector}};
    std::unordered_map<std::string, Format> formatMap = {{"coordinate", Format::Coordinate},
                                                         {"array", Format::Array}};
    std::unordered_map<std::string, Field> fieldMap = {{"real", Field::Real},
                                                       {"double", Field::Double},
                                                       {"integer", Field::Integer},
                                                       {"pattern", Field::Pattern}};
    std::unordered_map<std::string, Symmetry> symmetryMap = {
        {"general", Symmetry::General},
        {"symmetric", Symmetry::Symmetric},
        {"skew-symmetric", Symmetry::SkewSymmetric}};

    struct MTXHeader {
        Object object;
        Format format;
        Field field;
        Symmetry symmetry;
    };

    struct MatrixSize {
        count rows;
        count columns;
        count nonzeros;

        MatrixSize(count r, count c, count nz) : rows(r), columns(c), nonzeros(nz) {}
    };

    struct Edge {
        node from;
        node to;
        std::optional<double> weight;

        Edge(node f, node t, std::optional<double> w) : from(f), to(t), weight(w) {}
    };

    MTXParser(const std::string &path);

    /**
     * Get the MTX graph file header.
     */
    MTXHeader getHeader();

    /**
     * Get the MTX graph file matrix size header.
     */
    MatrixSize getMatrixSize();

    /**
     * Get the (weighted) edge from the next line in the MTX graph file if present.
     *
     * @param weighted
     * @return std::optional<std::pair<NetworKit::Edge, NetworKit::edgeweight>>
     */
    std::optional<NetworKit::WeightedEdge> getNext(bool weighted);
};
} /* namespace NetworKit */
#endif // NETWORKIT_IO_MTX_PARSER_HPP_
