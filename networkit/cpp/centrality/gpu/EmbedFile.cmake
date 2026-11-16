# Writes the bytes of the file INPUT into the header OUTPUT, as the array NetworKit::GPU::NAME.
# Usage: cmake -DINPUT=<file> -DOUTPUT=<header> -DNAME=<array> -P EmbedFile.cmake

file(READ "${INPUT}" bytes HEX)
string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," bytes "${bytes}")
# 16 bytes per line.
string(REGEX REPLACE "((0x[0-9a-f][0-9a-f],){16})" "\\1\n" bytes "${bytes}")
get_filename_component(inputName "${INPUT}" NAME)

file(WRITE "${OUTPUT}" "// Generated from ${inputName} by EmbedFile.cmake, do not edit.
#pragma once

namespace NetworKit {
namespace GPU {

constexpr unsigned char ${NAME}[] = {
${bytes}};

} // namespace GPU
} // namespace NetworKit
")
