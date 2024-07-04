#include <networkit/auxiliary/Enforce.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/NumberParsing.hpp>
#include <networkit/io/MTXParser.hpp>

#include <stdexcept>

namespace NetworKit {

/**
 * Extract a (weighted) edge from a line in the file.
 */
MTXParser::Edge parseLine(const std::string &line, const bool weighted) {
    node i, j;
    std::optional<double> value;
    std::istringstream istring(line);
    istring >> i >> j;
    i--;
    j--;

    if (weighted) {
        double tmp;
        if (istring >> tmp)
            value = tmp;
    }
    return MTXParser::Edge(i, j, value);
}

MTXParser::MTXParser(const std::string &path) : graphFile(path) {
    if (!(graphFile)) {
        ERROR("invalid graph file: ", path);
        throw std::runtime_error("invalid graph file");
    }
}

MTXParser::MTXHeader MTXParser::getHeader() {

    auto tolower = [](const std::string &str) {
        std::string out;
        std::transform(str.begin(), str.end(), std::back_inserter(out), ::tolower);
        return out;
    };

    MTXParser::MTXHeader header;
    std::string identifier, object, format, field, symmetry;
    std::string line;

    Aux::enforceOpened(graphFile);
    if (std::getline(graphFile, line)) {
        std::istringstream istring(line);
        if (!(istring >> identifier >> object >> format >> field >> symmetry)) {
            throw std::runtime_error("Invalid header format.");
        }
        // the header line must start with %%MatrixMarket or %MatrixMarket
        if (identifier != "%%MatrixMarket" && identifier != "%MatrixMarket") {
            throw std::runtime_error("Invalid identifier.");
        }
        try {
            header.object = objectMap.at(tolower(object));
            header.format = formatMap.at(tolower(format));
            header.field = fieldMap.at(tolower(field));
            header.symmetry = symmetryMap.at(tolower(symmetry));
        } catch (const std::out_of_range &) {
            throw std::runtime_error("Invalid header field.");
        }
        return header;
    } else {
        throw std::runtime_error("Getting MTX header failed.");
    }
}

MTXParser::MatrixSize MTXParser::getMatrixSize() {
    count rows, columns, nonzeros;
    std::string line;

    Aux::enforceOpened(graphFile);
    if (std::getline(graphFile, line)) {
        while (line[0] == '%') {
            std::getline(graphFile, line);
        }
        std::istringstream istring(line);
        if (!(istring >> rows >> columns >> nonzeros)) {
            throw std::runtime_error("Invalid matrix size line format.");
        }
        return MTXParser::MatrixSize(rows, columns, nonzeros);
    } else {
        throw std::runtime_error("Getting MTX matrix size line failed.");
    }
}

std::optional<NetworKit::WeightedEdge> MTXParser::getNext(bool weighted) {
    std::string line;

    do {
        auto hasLine = static_cast<bool>(std::getline(graphFile, line));
        if (!hasLine)
            return std::nullopt;
        // check for comment line starting with '%'
        if (line[0] == '%') {
            throw std::runtime_error(
                "Invalid MTX file structure. No comments allowed after size line.");
        } else {
            auto e = parseLine(line, weighted);
            edgeweight w = 1.0;
            if (weighted) {
                w = e.weight.value();
            }
            return std::optional<NetworKit::WeightedEdge>({e.from, e.to, w});
        }
    } while (true);
}

} /* namespace NetworKit */
