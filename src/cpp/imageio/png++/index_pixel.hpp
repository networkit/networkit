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
#ifndef PNGPP_INDEX_PIXEL_HPP_INCLUDED
#define PNGPP_INDEX_PIXEL_HPP_INCLUDED

#include "types.hpp"
#include "packed_pixel.hpp"
#include "pixel_traits.hpp"

namespace png
{

    /**
     * \brief The 8-bit Indexed (colormap) pixel type.
     */
    class index_pixel
    {
    public:
        index_pixel(byte index = 0)
            : m_index(index)
        {
        }

        operator byte() const
        {
            return m_index;
        }

    private:
        byte m_index;
    };

    /**
     * \brief The packed indexed pixel class template.  The available
     * specializations are for 1-, 2- and 4-bit pixels.
     */
    template< size_t bits >
    class packed_index_pixel
        : public packed_pixel< bits >
    {
    public:
        packed_index_pixel(byte value = 0)
            : packed_pixel< bits >(value)
        {
        }
    };

    /**
     * \brief The 1-bit Indexed pixel type.
     */
    typedef packed_index_pixel< 1 > index_pixel_1;

    /**
     * \brief The 1-bit Indexed pixel type.
     */
    typedef packed_index_pixel< 2 > index_pixel_2;

    /**
     * \brief The 1-bit Indexed pixel type.
     */
    typedef packed_index_pixel< 4 > index_pixel_4;

    /**
     * \brief Pixel traits specialization for index_pixel.
     */
    template<>
    struct pixel_traits< index_pixel >
        : basic_pixel_traits< index_pixel, byte, color_type_palette >
    {
    };

    /**
     * \brief Pixel traits specialization for packed_index_pixel.
     */
    template< size_t bits >
    struct pixel_traits< packed_index_pixel< bits > >
        : basic_pixel_traits< packed_index_pixel< bits >, byte,
                              color_type_palette, /* channels = */ 1, bits >
    {
    };

} // namespace png

#endif // PNGPP_INDEX_PIXEL_HPP_INCLUDED
