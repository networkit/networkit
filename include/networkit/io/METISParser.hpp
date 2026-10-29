/*
 * METISParser.hpp
 *
 *  Created on: 27.11.2012
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_IO_METIS_PARSER_HPP_
#define NETWORKIT_IO_METIS_PARSER_HPP_

#include <algorithm>
#include <cctype>
#include <charconv>
#include <exception>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <networkit/Globals.hpp>
#include <networkit/auxiliary/NumberParsing.hpp>

namespace NetworKit {

namespace METISParserDetails {

inline bool isSpace(char c) {
    return std::isspace(static_cast<unsigned char>(c));
}

inline const char *skipSpaces(const char *first, const char *last) {
    return std::find_if_not(first, last, isSpace);
}

// Returns the number of whitespace-separated tokens in [first, last).
inline size_t countTokens(const char *first, const char *last) {
    size_t tokens = 0;
    for (first = skipSpaces(first, last); first != last; first = skipSpaces(first, last)) {
        ++tokens;
        first = std::find_if(first, last, isSpace);
    }
    return tokens;
}

/**
 * Parses the whitespace-separated number at the beginning of [first, last) into @a value.
 * Integral numbers must be representable by T exactly: unlike a detour through double, which
 * rounds integers above 2^53, they are parsed directly into T.
 *
 * @return The position after the number.
 * @throws std::runtime_error if the range has no number, or the number is not a valid T.
 */
template <class T>
const char *parseNumber(const char *first, const char *last, T &value) {
    first = skipSpaces(first, last);
    const char *const tokenEnd = std::find_if(first, last, isSpace);
    // A view, so that the token is only copied into a string for error messages.
    const std::string_view token(first, static_cast<size_t>(tokenEnd - first));
    if (token.empty()) {
        throw std::runtime_error("METIS line ended where a value was expected");
    }

    const char *numberEnd = first;
    if constexpr (std::is_integral_v<T>) {
        const auto [ptr, ec] = std::from_chars(first, tokenEnd, value);
        if (ec == std::errc::result_out_of_range) {
            throw std::runtime_error("METIS value " + std::string(token)
                                     + " is out of range for the requested type");
        }
        if (ec == std::errc{}) {
            numberEnd = ptr;
        }
    } else {
        // std::from_chars for floating-point types is too recent for all supported compilers.
        try {
            std::tie(value, numberEnd) =
                Aux::Parsing::strTo<T, const char *, Aux::Checkers::Enforcer>(first, tokenEnd);
        } catch (const std::exception &) {
            numberEnd = first;
        }
    }
    // E.g., the fraction of 1.5 or the sign of -1 for unsigned types are not part of the number.
    if (numberEnd != tokenEnd) {
        throw std::runtime_error("METIS value " + std::string(token)
                                 + " is not valid for the requested type");
    }
    return tokenEnd;
}

/**
 * Parses the node ids of a METIS line into @a NodeT, after skipping its first @a ignoreFirst
 * values.
 *
 * @throws std::runtime_error if the line is malformed or a node id is not a valid NodeT.
 */
template <class NodeT>
std::vector<NodeT> parseLine(const std::string &line, count ignoreFirst = 0) {
    static_assert(std::is_integral_v<NodeT>, "Integral required for NodeT.");

    const char *it = line.data();
    const char *const end = it + line.size();
    for (count i = 0; i < ignoreFirst; ++i) {
        double ignored;
        it = parseNumber(it, end, ignored);
    }

    std::vector<NodeT> adjacencies;
    adjacencies.reserve(countTokens(it, end));
    while (skipSpaces(it, end) != end) {
        NodeT v;
        it = parseNumber(it, end, v);
        adjacencies.push_back(v);
    }
    return adjacencies;
}

} // namespace METISParserDetails

class METISParser final {

    std::ifstream graphFile;

    // Returns the next line that is not a comment.
    std::string nextLine();

public:
    /**
     * Constructor for the METIS Parser.
     *
     * @param[in]  path  file path
     */
    METISParser(std::string_view path);

    /**
     * Get the METIS graph file header
     */
    std::tuple<count, count, index, count> getHeader();

    /**
     * Test if graph file has a next line.
     */
    bool hasNext();

    /**
     * Get adjacencies from the next line in the METIS graph file. The node ids are parsed
     * directly into @a NodeT, and are 1-based as in the file.
     *
     * @param[in]  ignoreFirst  number of values to ignore [in case the METIS file contains node
     * weightes]
     * @throws std::runtime_error if the line is malformed or a node id is not a valid NodeT.
     */
    template <class NodeT = node>
    std::vector<NodeT> getNext(count ignoreFirst = 0) {
        return METISParserDetails::parseLine<NodeT>(nextLine(), ignoreFirst);
    }

    /**
     * Get adjacencies with edge weights from the next line in the METIS graph file. The node ids
     * are parsed directly into @a NodeT, and are 1-based as in the file. The weights are parsed
     * directly into @a EdgeWeightT.
     *
     * @param[in]  ignoreFirst  number of values to ignore [in case the METIS file contains node
     * weightes]
     * @throws std::runtime_error if the line is malformed, a node id is not a valid NodeT, or a
     * weight is not a valid EdgeWeightT.
     */
    template <class NodeT = node, class EdgeWeightT = edgeweight>
    std::vector<std::pair<NodeT, EdgeWeightT>> getNextWithWeights(count ignoreFirst = 0);
};

template <class NodeT, class EdgeWeightT>
std::vector<std::pair<NodeT, EdgeWeightT>> METISParser::getNextWithWeights(count ignoreFirst) {
    static_assert(std::is_integral_v<NodeT>, "Integral required for NodeT.");
    using METISParserDetails::parseNumber;

    const std::string line = nextLine();
    const char *it = line.data();
    const char *const end = it + line.size();
    for (count i = 0; i < ignoreFirst; ++i) {
        double ignored;
        it = parseNumber(it, end, ignored);
    }

    std::vector<std::pair<NodeT, EdgeWeightT>> adjacencies;
    // A line with an odd number of tokens is malformed and throws below.
    adjacencies.reserve(METISParserDetails::countTokens(it, end) / 2);
    while (METISParserDetails::skipSpaces(it, end) != end) {
        NodeT v;
        EdgeWeightT weight;
        it = parseNumber(it, end, v);
        it = parseNumber(it, end, weight);
        adjacencies.emplace_back(v, weight);
    }
    return adjacencies;
}

} /* namespace NetworKit */
#endif // NETWORKIT_IO_METIS_PARSER_HPP_
