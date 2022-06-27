/*
 * MemoryMappedFileGTest.hpp
 *
 *  Created on: 28.07.2018
 *      Author: Manuel Penschuck (networkit@manuel.jetzt)
 */

#ifndef MEMORY_MAPPED_FILE_GTEST_H_
#define MEMORY_MAPPED_FILE_GTEST_H_

#include <cstddef>
#include <gtest/gtest.h>

namespace NetworKit {

class MemoryMappedFileGTest : public ::testing::Test {};
class MemoryMappedFileIOGTest : public ::testing::TestWithParam<size_t> {};

} // namespace NetworKit

#endif // ! MEMORY_MAPPED_FILE_GTEST_H_
