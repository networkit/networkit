// networkit-format

#include <networkit/io/LineFileReader.hpp>

namespace NetworKit {

std::vector<std::string> read(const std::string &path) {
    std::ifstream file;
    std::string line; // the current line
    file.open(path);

    std::vector<std::string> data;

    while (file.good()) {
        std::getline(file, line);
        data.push_back(line);
    }

    file.close();

    return data;
}

} // namespace NetworKit
