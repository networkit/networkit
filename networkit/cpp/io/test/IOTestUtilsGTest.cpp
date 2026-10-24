#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

#include <gtest/gtest.h>

#include "IOTestUtils.hpp"

namespace NetworKit {
namespace {

using TestUtils::TemporaryDirectory;

std::filesystem::path directoryOf(const TemporaryDirectory &directory) {
    return std::filesystem::path(directory.path("file")).parent_path();
}

std::string readFile(const std::string &path) {
    std::ifstream file(path);
    return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

TEST(IOTestUtilsGTest, testCreatesEmptyDirectoryInTempDir) {
    const TemporaryDirectory directory;
    const std::filesystem::path path = directoryOf(directory);
    EXPECT_TRUE(std::filesystem::is_directory(path));
    EXPECT_TRUE(std::filesystem::is_empty(path));
    EXPECT_TRUE(
        std::filesystem::equivalent(path.parent_path(), std::filesystem::temp_directory_path()));
}

TEST(IOTestUtilsGTest, testInstancesHaveDistinctDirectories) {
    const TemporaryDirectory first, second;
    EXPECT_NE(directoryOf(first), directoryOf(second));
}

TEST(IOTestUtilsGTest, testWriteFileCreatesFileAtPath) {
    const TemporaryDirectory directory;
    const std::string path = directory.writeFile("file", "first line\nsecond line\n");
    EXPECT_EQ(path, directory.path("file"));
    EXPECT_EQ(readFile(path), "first line\nsecond line\n");
}

TEST(IOTestUtilsGTest, testWriteFileOverwritesFile) {
    const TemporaryDirectory directory;
    directory.writeFile("file", "longer contents");
    EXPECT_EQ(readFile(directory.writeFile("file", "short")), "short");
}

TEST(IOTestUtilsGTest, testDestructorRemovesDirectoryWithContents) {
    std::filesystem::path path;
    {
        const TemporaryDirectory directory;
        path = directoryOf(directory);
        directory.writeFile("file", "contents");
        std::filesystem::create_directory(directory.path("subdirectory"));
        directory.writeFile("subdirectory/file", "contents");
    }
    EXPECT_FALSE(std::filesystem::exists(path));
}

} // namespace
} // namespace NetworKit
