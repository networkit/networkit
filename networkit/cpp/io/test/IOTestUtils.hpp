#ifndef NETWORKIT_IO_TEST_IO_TEST_UTILS_HPP_
#define NETWORKIT_IO_TEST_IO_TEST_UTILS_HPP_

#include <atomic>
#include <filesystem>
#include <fstream>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif // _WIN32

namespace NetworKit::TestUtils {

// A directory in std::filesystem::temp_directory_path() for the files of a test, removed with
// its contents when it goes out of scope. Its name is unique, so that test processes running in
// parallel (e.g., the tests with and without assertions under ctest -j) do not overwrite each
// other's files.
class TemporaryDirectory {
public:
    // Throws std::runtime_error if no directory with a unique name could be created.
    TemporaryDirectory() {
        constexpr int maxAttempts = 100;
        // The process id and the attempt counter make the name unique among running processes and
        // among the threads of this one. The random suffix skips directories left behind by a
        // crashed process whose id has been reused.
        static std::atomic<int> attempts{0};
        static thread_local std::mt19937_64 prng{std::random_device{}()};

        const std::filesystem::path tempDirectory = std::filesystem::temp_directory_path();
        const std::string prefix = "networkit_test_" + std::to_string(processId()) + "_";
        for (int i = 0; i < maxAttempts; ++i) {
            directory =
                tempDirectory
                / (prefix + std::to_string(attempts.fetch_add(1)) + "_" + std::to_string(prng()));
            // create_directory returns false if the directory already exists
            if (std::filesystem::create_directory(directory)) {
                return;
            }
        }
        throw std::runtime_error("Could not create a temporary directory in "
                                 + tempDirectory.string());
    }

    TemporaryDirectory(const TemporaryDirectory &) = delete;
    TemporaryDirectory &operator=(const TemporaryDirectory &) = delete;

    ~TemporaryDirectory() {
        std::error_code ignored; // a destructor must not throw
        std::filesystem::remove_all(directory, ignored);
    }

    // Returns the path of the file with the given name in this directory.
    std::string path(std::string_view fileName) const { return (directory / fileName).string(); }

    // Creates the file with the given name and contents in this directory and returns its path.
    std::string writeFile(std::string_view fileName, std::string_view contents) const {
        std::string filePath = path(fileName);
        std::ofstream file(filePath);
        file << contents;
        if (!file) {
            throw std::runtime_error("Could not write " + filePath);
        }
        return filePath;
    }

private:
    static int processId() {
#ifdef _WIN32
        return _getpid();
#else
        return static_cast<int>(getpid());
#endif // _WIN32
    }

    std::filesystem::path directory;
};

} // namespace NetworKit::TestUtils

#endif // NETWORKIT_IO_TEST_IO_TEST_UTILS_HPP_
