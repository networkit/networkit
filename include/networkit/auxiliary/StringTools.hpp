/*
 * StringTools.hpp
 *
 *  Completely rewritten on: 12.05.2014
 *      Author: Florian Weber (uagws@student.kit.edu)
 */

#ifndef NETWORKIT_AUXILIARY_STRING_TOOLS_HPP_
#define NETWORKIT_AUXILIARY_STRING_TOOLS_HPP_

#include <algorithm>
#include <string>
#include <string_view>
#include <vector>

#include <tlx/define/deprecated.hpp>

namespace Aux {

/**
 * Missing string functions.
 */
namespace StringTools {

/**
 * Splits a range of characters at a delimiter into a vector of strings.
 *
 * Requirements:
 *   Character must be equality-comparable to char and it must be possible to construct a char from
 *   any Character.
 *   Iterator must be an input-iterator over Characters.
 */
template <typename Iterator, typename Character>
std::vector<std::string> split(Iterator begin, Iterator end, Character delim = Character{' '}) {

    // measurements showed that precalculating the number of tokens and
    // reserving space for them was in fact slower than just letting
    // the vector grow naturally.
    std::vector<std::string> tokens;

    auto it = begin;
    while (it != end) {
        auto tmp = std::find(it, end, delim);
        tokens.emplace_back(it, tmp);
        if (tmp == end) {
            break;
        }
        it = tmp;
        ++it;
    }
    return tokens;
}

/**
 * Split a string at delimiter and return vector of parts.
 */
inline std::vector<std::string> split(std::string_view s, char delim = ' ') {
    return split(s.begin(), s.end(), delim);
}

/**
 * DEPRECATED: use std::string_view::ends_with instead.
 * Determines whether @a str ends with @a suffix.
 */
inline bool TLX_DEPRECATED(ends_with(std::string_view str, std::string_view suffix));

/**
 * DEPRECATED: use std::string_view::starts_with instead.
 * Determines whether @a str starts with @a prefix.
 */
inline bool TLX_DEPRECATED(starts_with(std::string_view str, std::string_view prefix));

inline bool ends_with(std::string_view str, std::string_view suffix) {
    return str.ends_with(suffix);
}

inline bool starts_with(std::string_view str, std::string_view prefix) {
    return str.starts_with(prefix);
}

} /* namespace StringTools */
} /* namespace Aux */
#endif // NETWORKIT_AUXILIARY_STRING_TOOLS_HPP_
