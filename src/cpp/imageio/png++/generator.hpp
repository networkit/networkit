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
#ifndef PNGPP_GENERATOR_HPP_INCLUDED
#define PNGPP_GENERATOR_HPP_INCLUDED

#include <cassert>
#include <stdexcept>
#include <iostream>
#include <ostream>

#include "config.hpp"
#include "error.hpp"
#include "streaming_base.hpp"
#include "writer.hpp"

namespace png
{

    /**
     * \brief Pixel generator class template.
     *
     * Used as a base class for custom pixel generator classes as well
     * as inside image class implementation to write pixels from the
     * pixel buffer.
     *
     * A usage example can be found in \c example/pixel_generator.cpp.
     *
     * Encapsulates PNG %image writing procedure.  In order to create
     * a custom pixel %generator use CRTP trick:
     *
     * \code
     * class pixel_generator
     *     : public png::generator< pixel, pixel_generator >
     * {
     *     ...
     * };
     * \endcode
     *
     * Your pixel %generator class should implement \c get_next_row()
     * method and \c reset() method (optional).  Their signatures are
     * as follows:
     *
     * \code
     * png::byte* get_next_row(size_t pos);
     * void reset(size_t pass);
     * \endcode
     *
     * The \c get_next_row() method is called every time a new row of
     * %image data is needed by the writer.  The position of the row
     * being written is passed as \c pos parameter.  The \c pos takes
     * values from \c 0 to \c <image_height>-1 inclusively.  The
     * method should return the starting address of a row buffer
     * storing an appropriate amount of pixels (i.e. the width of the
     * %image being written).  The address should be casted to
     * png::byte* pointer type using \c reinterpret_cast<> or a
     * C-style cast.
     *
     * The optional \c reset() method is called every time the new
     * pass of interlaced %image processing starts.  The number of
     * interlace pass is avaiable as the only parameter of the method.
     * For non-interlaced images the method is called once prior to
     * any calls to \c get_next_row().  The value of \c 0 is passed
     * for the \c pass number.  You do not have to implement this
     * method unless you are going to support interlaced %image
     * generation.
     *
     * An optional template parameter \c info_holder encapsulated
     * image_info storage policy.  Please refer to consumer class
     * documentation for the detailed description of this parameter.
     *
     * An optional \c bool template parameter \c interlacing_supported
     * specifies whether writing interlacing images is supported by
     * your %generator class.  It defaults to \c false.  An attempt to
     * write an interlaced %image will result in throwing
     * \c std::logic_error.
     *
     * In order to fully support interlacing specify \c true for \c
     * interlacing_supported parameter and implement \c reset()
     * method.  You _must_ generate the same pixels for every pass to
     * get the correct PNG %image output.
     *
     * \see image, consumer
     */
    template< typename pixel,
              class pixgen,
              class info_holder = def_image_info_holder,
              bool interlacing_supported = false >
    class generator
        : public streaming_base< pixel, info_holder >
    {
    public:
        /**
         * \brief Writes an image to the stream.
         *
         * Essentially, this method constructs a writer object and
         * instructs it to write the image to the stream.  It handles
         * writing interlaced images as long as your generator class
         * supports this.
         */
        template< typename ostream >
        void write(ostream& stream)
        {
            writer< ostream > wr(stream);
            wr.set_image_info(this->get_info());
            wr.write_info();

#if __BYTE_ORDER == __LITTLE_ENDIAN
            if (pixel_traits< pixel >::get_bit_depth() == 16)
            {
#ifdef PNG_WRITE_SWAP_SUPPORTED
                wr.set_swap();
#else
                throw error("Cannot write 16-bit image:"
                            " recompile with PNG_WRITE_SWAP_SUPPORTED.");
#endif
            }
#endif

            size_t pass_count;
            if (this->get_info().get_interlace_type() != interlace_none)
            {
#ifdef PNG_WRITE_INTERLACING_SUPPORTED
                if (interlacing_supported)
                {
                    pass_count = wr.set_interlace_handling();
                }
                else
                {
                    throw std::logic_error("Cannot write interlaced image:"
                                           " generator does not support it.");
                }
#else
                throw error("Cannot write interlaced image:"
                            " interlace handling disabled.");
#endif
            }
            else
            {
                pass_count = 1;
            }
            pixgen* pixel_gen = static_cast< pixgen* >(this);
            for (size_t pass = 0; pass < pass_count; ++pass)
            {
                pixel_gen->reset(pass);

                for (size_t pos = 0; pos < this->get_info().get_height(); ++pos)
                {
                    wr.write_row(pixel_gen->get_next_row(pos));
                }
            }

            wr.write_end_info();
        }

    protected:
        typedef streaming_base< pixel, info_holder > base;

        /**
         * \brief Constructs a generator object using passed image_info
         * object to store image information.
         */
        explicit generator(image_info& info)
            : base(info)
        {
        }

        /**
         * \brief Constructs a generator object prepared to generate
         * an image of specified width and height.
         */
        generator(size_t width, size_t height)
            : base(width, height)
        {
        }
    };

} // namespace png

#endif // PNGPP_GENERATOR_HPP_INCLUDED
