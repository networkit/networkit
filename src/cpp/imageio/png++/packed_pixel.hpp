/*
 * Copyright (C) 2007,2008   Alex Shulgin
 *
 * This file is part of png++ the C++ wrapper for libpng.  PNG++ is free
 * software; the exact copying conditions are as follows:
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. The name of the author may not be used to endorse or promote products
 * derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
 * NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef PNGPP_PACKED_PIXEL_HPP_INCLUDED
#define PNGPP_PACKED_PIXEL_HPP_INCLUDED

#include "types.hpp"

namespace png
{

    namespace detail
    {
        template< size_t bits > class allowed_bit_depth;

        template<> class allowed_bit_depth< 1 > {};
        template<> class allowed_bit_depth< 2 > {};
        template<> class allowed_bit_depth< 4 > {};
    } // namespace detail

    /**
     * \brief The packed pixel class template.
     *
     * \see packed_gray_pixel, packed_index_pixel
     */
    template< size_t bits >
    class packed_pixel
        : detail::allowed_bit_depth< bits >
    {
    public:
        packed_pixel(byte value = 0)
            : m_value(value & get_bit_mask())
        {
        }

        operator byte() const
        {
            return m_value;
        }

        static size_t const get_bit_depth()
        {
            return bits;
        }

        static byte const get_bit_mask()
        {
            return (1 << bits) - 1;
        }

    private:
        byte m_value;
    };

} // namespace png

#endif // PNGPP_PACKED_PIXEL_HPP_INCLUDED
