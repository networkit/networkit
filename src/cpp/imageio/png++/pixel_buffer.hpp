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
#ifndef PNGPP_PIXEL_BUFFER_HPP_INCLUDED
#define PNGPP_PIXEL_BUFFER_HPP_INCLUDED

#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "packed_pixel.hpp"
#include "gray_pixel.hpp"
#include "index_pixel.hpp"

namespace png
{

    /**
     * \brief The pixel row traits class template.  Provides a common
     * way to get starting address of the row for packed and unpacked
     * row types.
     *
     * Not implemented--see specializations.
     */
    template< typename row > class row_traits;

    /**
     * \brief The basic class template to represent image pixel data.
     */
    template< typename pixel,
              typename row,
              class traits = row_traits< row > >
    class basic_pixel_buffer
    {
    public:
        /**
         * \brief A row of pixel data.
         */
        typedef row row_type;
        typedef traits row_traits;

        /**
         * \brief Constructs an empty 0x0 pixel buffer object.
         */
        basic_pixel_buffer()
            : m_width(0),
              m_height(0)
        {
        }

        /**
         * \brief Constructs an empty pixel buffer object.
         */
        basic_pixel_buffer(size_t width, size_t height)
            : m_width(0),
              m_height(0)
        {
            resize(width, height);
        }

        size_t get_width() const
        {
            return m_width;
        }

        size_t get_height() const
	    {
            return m_height;
	    }

        /**
         * \brief Resizes the pixel buffer.
         *
         * If new width or height is greater than the original,
         * expanded pixels are filled with value of \a pixel().
         */
        void resize(size_t width, size_t height)
        {
            m_width = width;
            m_height = height;
            m_rows.resize(height);
            for (typename row_vec::iterator r = m_rows.begin();
                 r != m_rows.end();
                 ++r)
            {
                r->resize(width);
            }
        }

        /**
         * \brief Returns a reference to the row of image data at
         * specified index.
         *
         * Checks the index before returning a row: an instance of
         * std::out_of_range is thrown if \c index is greater than \c
         * height.
         */
        row_type& get_row(size_t index)
        {
            return m_rows.at(index);
        }

        /**
         * \brief Returns a const reference to the row of image data at
         * specified index.
         *
         * The checking version.
         */
        row_type const& get_row(size_t index) const
        {
            return m_rows.at(index);
        }

        /**
         * \brief The non-checking version of get_row() method.
         */
        row_type& operator[](size_t index)
        {
            return m_rows[index];
        }

        /**
         * \brief The non-checking version of get_row() method.
         */
        row_type const& operator[](size_t index) const
        {
            return m_rows[index];
        }

        /**
         * \brief Replaces the row at specified index.
         */
        void put_row(size_t index, row_type const& r)
        {
            assert(r.size() == m_width);
            m_rows.at(index) = r;
        }

        /**
         * \brief Returns a pixel at (x,y) position.
         */
        pixel get_pixel(size_t x, size_t y) const
        {
            return get_row(y).at(x);
        }

        /**
         * \brief Replaces a pixel at (x,y) position.
         */
        void set_pixel(size_t x, size_t y, pixel p)
        {
            get_row(y).at(x) = p;
        }

    protected:
        size_t m_width;
        size_t m_height;
        typedef std::vector< row_type > row_vec;
        row_vec m_rows;
    };

    /**
     * \brief The row_traits specialization for unpacked pixel rows.
     */
    template< typename pixel >
    class row_traits< std::vector< pixel > >
    {
    public:
        /**
         * \brief Returns the starting address of the row.
         */
        static pixel* get_data(std::vector< pixel >& vec)
        {
            assert(vec.size());
            return & vec[0];
        }
    };

    /**
     * The pixel_buffer specialization for unpacked pixels.
     */
    template< typename pixel >
    class pixel_buffer
        : public basic_pixel_buffer< pixel, std::vector< pixel > >
    {
    public:
        pixel_buffer()
        {
        }

        pixel_buffer(size_t width, size_t height)
            : basic_pixel_buffer< pixel, std::vector< pixel > >(width, height)
        {
        }
    };

    namespace detail
    {

        template< class pixel, typename reference >
        class basic_packed_pixel_proxy
        {
        public:
            explicit basic_packed_pixel_proxy(reference ref)
                : m_ref(ref),
                  m_shift(0)
            {
            }

            basic_packed_pixel_proxy(reference ref, size_t index)
                : m_ref(ref),
                  m_shift(get_shift(index))
            {
            }

            operator pixel() const
            {
                return pixel((m_ref >> m_shift) & pixel::get_bit_mask());
            }

        protected:
            /*
             * bits: . .   .
             *    1: 7 6 5 4 3 2 1 0
             *    2:   6   4   2   0
             *    4:       4       0
             */
            static size_t get_shift(size_t index)
            {
                size_t const bits = pixel::get_bit_depth();
                return (8 - bits) - (index % get_pixels_per_byte()) * bits;
            }

            static size_t get_pixels_per_byte()
            {
                return 8 / pixel::get_bit_depth();
            }

            reference m_ref;
            size_t m_shift;
        };

        template< class pixel >
        class const_packed_pixel_proxy
            : public basic_packed_pixel_proxy< pixel, byte const& >
        {
        public:
            const_packed_pixel_proxy(byte const& ref, size_t index)
                : basic_packed_pixel_proxy< pixel, byte const& >(ref, index)
            {
            }
        };

        template< class pixel >
        class packed_pixel_proxy
            : public basic_packed_pixel_proxy< pixel, byte& >
        {
        public:
            typedef basic_packed_pixel_proxy< pixel, byte& > basic_proxy;

            packed_pixel_proxy(byte& ref, size_t index)
                : basic_proxy(ref, index)
            {
            }

            packed_pixel_proxy(packed_pixel_proxy const& other)
                : basic_proxy(other.m_ref)
            {
                this->m_shift = other.m_shift;
            }

            packed_pixel_proxy& operator=(packed_pixel_proxy const& other)
            {
                return *this = static_cast< pixel >(other);
            }

            template< typename reference >
            packed_pixel_proxy&
            operator=(basic_packed_pixel_proxy< pixel, reference > const& other)
            {
                return *this = static_cast< pixel >(other);
            }

            packed_pixel_proxy& operator=(pixel p)
            {
                this->m_ref = (this->m_ref
                    & ~(pixel::get_bit_mask() << this->m_shift))
                    | (p << this->m_shift);

                return *this;
            }
        };

    } // namespace detail

    /**
     * \brief The packed pixel row class template.
     *
     * Stores the pixel row as a std::vector of byte-s, providing
     * access to individual packed pixels via proxy objects.
     */
    template< class pixel >
    class packed_pixel_row
    {
    public:
        /**
         * \brief Constructs a pixel row object for \c size packed pixels.
         */
        explicit packed_pixel_row(size_t size = 0)
        {
            resize(size);
        }

        size_t size() const
        {
            return m_size;
        }

        /**
         * \brief Resizes the pixel row to hold up to \c size packed pixels.
         */
        void resize(size_t size)
        {
            m_vec.resize(size / get_pixels_per_byte()
                         + (size % get_pixels_per_byte() ? 1 : 0));
            m_size = size;
        }

        /**
         * \brief The immutable packed pixel proxy type.
         */
        typedef detail::const_packed_pixel_proxy< pixel > const_pixel_proxy;

        /**
         * \brief The mutable packed pixel proxy type.
         */
        typedef detail::packed_pixel_proxy< pixel > pixel_proxy;

        /**
         * \brief Returns an immutable proxy the to the pixel at \c
         * index.
         */
        const_pixel_proxy at(size_t index) const
        {
            return const_pixel_proxy(m_vec.at(index / get_pixels_per_byte()),
                                     index);
        }

        /**
         * \brief Returns a mutable proxy the to the pixel at \c
         * index.
         */
        pixel_proxy at(size_t index)
        {
            return pixel_proxy(m_vec.at(index / get_pixels_per_byte()),
                               index);
        }

        /**
         * \brief Returns an immutable proxy the to the pixel at \c
         * index.  The non-checking version.
         */
        const_pixel_proxy operator[](size_t index) const
        {
            return const_pixel_proxy(m_vec[index / get_pixels_per_byte()],
                                     index);
        }

        /**
         * \brief Returns n mutable proxy the to the pixel at \c
         * index.  The non-checking version.
         */
        pixel_proxy operator[](size_t index)
        {
            return pixel_proxy(m_vec[index / get_pixels_per_byte()],
                               index);
        }

        /**
         * \brief Returns the starting address of the row.
         */
        byte* get_data()
        {
            assert(m_vec.size());
            return & m_vec[0];
        }

    private:
        static size_t get_pixels_per_byte()
        {
            return 8 / pixel::get_bit_depth();
        }

        std::vector< byte > m_vec;
        size_t m_size;
    };

    /**
     * \brief The row_traits class template specialization for packed
     * pixel row type.
     */
    template< typename pixel >
    class row_traits< packed_pixel_row< pixel > >
    {
    public:
        /**
         * \brief Returns the starting address of the row.
         */
        static byte* get_data(packed_pixel_row< pixel >& row)
        {
            return row.get_data();
        }
    };

    /**
     * \brief The pixel buffer class template specialization for the
     * packed_gray_pixel type.
     */
    template< size_t bits >
    class pixel_buffer< packed_gray_pixel< bits > >
        : public basic_pixel_buffer< packed_gray_pixel< bits >,
                                     packed_pixel_row< packed_gray_pixel
                                                       < bits > > >
    {
    public:
        typedef packed_gray_pixel< bits > pixel_type;
        typedef packed_pixel_row< pixel_type > pixel_row_type;

        pixel_buffer()
        {
        }

        pixel_buffer(size_t width, size_t height)
            : basic_pixel_buffer< pixel_type,
                                  pixel_row_type >(width, height)
        {
        }
    };

    /**
     * \brief The pixel buffer class template specialization for the
     * packed_index_pixel type.
     */
    template< size_t bits >
    class pixel_buffer< packed_index_pixel< bits > >
        : public basic_pixel_buffer< packed_index_pixel< bits >,
                                     packed_pixel_row< packed_index_pixel
                                                       < bits > > >
    {
    public:
        typedef packed_index_pixel< bits > pixel_type;
        typedef packed_pixel_row< pixel_type > pixel_row_type;

        pixel_buffer()
        {
        }

        pixel_buffer(size_t width, size_t height)
            : basic_pixel_buffer< pixel_type,
                                  pixel_row_type >(width, height)
        {
        }
    };

} // namespace png

#endif // PNGPP_PIXEL_BUFFER_HPP_INCLUDED
