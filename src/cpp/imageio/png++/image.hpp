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
#ifndef PNGPP_IMAGE_HPP_INCLUDED
#define PNGPP_IMAGE_HPP_INCLUDED

#include <fstream>
#include "pixel_buffer.hpp"
#include "generator.hpp"
#include "consumer.hpp"
#include "convert_color_space.hpp"

namespace png
{

    /**
     * \brief Class template to represent PNG image.
     *
     * The image consists of pixel data as well as additional %image
     * %info like interlace type, compression method, palette (for
     * colormap-based images) etc.  Provides methods to read and write
     * images from/to a generic stream and to manipulate %image pixels.
     */
    template< typename pixel >
    class image
    {
    public:
        /**
         * \brief The pixel traits type for \c pixel.
         */
        typedef pixel_traits< pixel > traits;

        /**
         * \brief The pixel buffer type for \c pixel.
         */
        typedef pixel_buffer< pixel > pixbuf;

        /**
         * \brief Represents a row of image pixel data.
         */
        typedef typename pixbuf::row_type row_type;

        /**
         * \brief A transformation functor to convert any image to
         * appropriate color space.
         */
        typedef convert_color_space< pixel > transform_convert;

        /**
         * \brief The default io transformation: does nothing.
         */
        struct transform_identity
        {
            void operator()(io_base&) const {}
        };

        /**
         * \brief Constructs an empty image.
         */
        image()
            : m_info(make_image_info< pixel >())
        {
        }

        /**
         * \brief Constructs an empty image of specified width and height.
         */
        image(size_t width, size_t height)
            : m_info(make_image_info< pixel >())
        {
            resize(width, height);
        }

        /**
         * \brief Constructs an image reading data from specified file
         * using default converting transform.
         */
        explicit image(std::string const& filename)
        {
            read(filename, transform_convert());
        }

        /**
         * \brief Constructs an image reading data from specified file
         * using custom transformaton.
         */
        template< class transformation >
        image(std::string const& filename,
              transformation const& transform)
        {
            read(filename.c_str(), transform);
        }

        /**
         * \brief Constructs an image reading data from specified file
         * using default converting transform.
         */
        explicit image(char const* filename)
        {
            read(filename, transform_convert());
        }

        /**
         * \brief Constructs an image reading data from specified file
         * using custom transformaton.
         */
        template< class transformation >
        image(char const* filename, transformation const& transform)
        {
            read(filename, transform);
        }

        /**
         * \brief Constructs an image reading data from a stream using
         * default converting transform.
         */
        explicit image(std::istream& stream)
        {
            read_stream(stream, transform_convert());
        }

        /**
         * \brief Constructs an image reading data from a stream using
         * custom transformation.
         */
        template< class transformation >
        image(std::istream& stream, transformation const& transform)
        {
            read_stream(stream, transform);
        }

        /**
         * \brief Reads an image from specified file using default
         * converting transform.
         */
        void read(std::string const& filename)
        {
            read(filename, transform_convert());
        }

        /**
         * \brief Reads an image from specified file using custom
         * transformaton.
         */
        template< class transformation >
        void read(std::string const& filename, transformation const& transform)
        {
            read(filename.c_str(), transform);
        }

        /**
         * \brief Reads an image from specified file using default
         * converting transform.
         */
        void read(char const* filename)
        {
            read(filename, transform_convert());
        }

        /**
         * \brief Reads an image from specified file using custom
         * transformaton.
         */
        template< class transformation >
        void read(char const* filename, transformation const& transform)
        {
            std::ifstream stream(filename, std::ios::binary);
            if (!stream.is_open())
            {
                throw std_error(filename);
            }
            stream.exceptions(std::ios::badbit);
            read_stream(stream, transform);
        }

        /**
         * \brief Reads an image from a stream using default
         * converting transform.
         */
        void read(std::istream& stream)
        {
            read_stream(stream, transform_convert());
        }

        /**
         * \brief Reads an image from a stream using custom
         * transformation.
         */
        template< class transformation >
        void read(std::istream& stream, transformation const& transform)
        {
            read_stream(stream, transform);
        }

        /**
         * \brief Reads an image from a stream using default
         * converting transform.
         */
        template< class istream >
        void read_stream(istream& stream)
        {
            read_stream(stream, transform_convert());
        }

        /**
         * \brief Reads an image from a stream using custom
         * transformation.
         */
        template< class istream, class transformation >
        void read_stream(istream& stream, transformation const& transform)
        {
            pixel_consumer pixcon(m_info, m_pixbuf);
            pixcon.read(stream, transform);
        }

        /**
         * \brief Writes an image to specified file.
         */
        void write(std::string const& filename)
        {
            write(filename.c_str());
        }

        /**
         * \brief Writes an image to specified file.
         */
        void write(char const* filename)
        {
            std::ofstream stream(filename, std::ios::binary);
            if (!stream.is_open())
            {
                throw std_error(filename);
            }
            stream.exceptions(std::ios::badbit);
            write_stream(stream);
        }

        /**
         * \brief Writes an image to a stream.
         */
        void write_stream(std::ostream& stream)
        {
            write_stream(stream);
        }

        /**
         * \brief Writes an image to a stream.
         */
        template< class ostream >
        void write_stream(ostream& stream)
        {
            pixel_generator pixgen(m_info, m_pixbuf);
            pixgen.write(stream);
        }

        /**
         * \brief Returns a reference to image pixel buffer.
         */
        pixbuf& get_pixbuf()
        {
            return m_pixbuf;
        }
        
        /**
         * \brief Returns a const reference to image pixel buffer.
         */
        pixbuf const& get_pixbuf() const
        {
            return m_pixbuf;
        }

        /**
         * \brief Replaces the image pixel buffer.
         *
         * \param buffer  a pixel buffer object to take a copy from
         */
        void set_pixbuf(pixbuf const& buffer)
        {
            m_pixbuf = buffer;
        }

        size_t get_width() const
        {
            return m_pixbuf.get_width();
        }

        size_t get_height() const
        {
            return m_pixbuf.get_height();
        }

        /**
         * \brief Resizes the image pixel buffer.
         */
        void resize(size_t width, size_t height)
        {
            m_pixbuf.resize(width, height);
            m_info.set_width(width);
            m_info.set_height(height);
        }

        /**
         * \brief Returns a reference to the row of image data at
         * specified index.
         *
         * \see pixel_buffer::get_row()
         */
        row_type& get_row(size_t index)
        {
            return m_pixbuf.get_row(index);
        }
        
        /**
         * \brief Returns a const reference to the row of image data at
         * specified index.
         *
         * \see pixel_buffer::get_row()
         */
        row_type const& get_row(size_t index) const
        {
            return m_pixbuf.get_row(index);
        }

        /**
         * \brief The non-checking version of get_row() method.
         */
        row_type& operator[](size_t index)
        {
            return m_pixbuf[index];
        }
        
        /**
         * \brief The non-checking version of get_row() method.
         */
        row_type const& operator[](size_t index) const
        {
            return m_pixbuf[index];
        }

        /**
         * \brief Returns a pixel at (x,y) position.
         */
        pixel get_pixel(size_t x, size_t y) const
        {
            return m_pixbuf.get_pixel(x, y);
        }

        /**
         * \brief Replaces a pixel at (x,y) position.
         */
        void set_pixel(size_t x, size_t y, pixel p)
        {
            m_pixbuf.set_pixel(x, y, p);
        }

        interlace_type get_interlace_type() const
        {
            return m_info.get_interlace_type();
        }

        void set_interlace_type(interlace_type interlace)
        {
            m_info.set_interlace_type(interlace);
        }

        compression_type get_compression_type() const
        {
            return m_info.get_compression_type();
        }

        void set_compression_type(compression_type compression)
        {
            m_info.set_compression_type(compression);
        }

        filter_type get_filter_type() const
        {
            return m_info.get_filter_type();
        }

        void set_filter_type(filter_type filter)
        {
            m_info.set_filter_type(filter);
        }

        /**
         * \brief Returns a reference to the image palette.
         */
        palette& get_palette()
        {
            return m_info.get_palette();
        }

        /**
         * \brief Returns a const reference to the image palette.
         */
        palette const& get_palette() const
        {
            return m_info.get_palette();
        }

        /**
         * \brief Replaces the image palette.
         */
        void set_palette(palette const& plte)
        {
            m_info.set_palette(plte);
        }

        tRNS const& get_tRNS() const
        {
            return m_info.get_tRNS();
        }

        tRNS& get_tRNS()
        {
            return m_info.get_tRNS();
        }

        void set_tRNS(tRNS const& trns)
        {
            m_info.set_tRNS(trns);
        }

    protected:
        /**
         * \brief A common base class template for pixel_consumer and
         * pixel_generator classes.
         */
        template< typename base_impl >
        class streaming_impl
            : public base_impl
        {
        public:
            streaming_impl(image_info& info, pixbuf& pixels)
                : base_impl(info),
                  m_pixbuf(pixels)
            {
            }

            /**
             * \brief Returns the starting address of a \c pos-th row
             * in the image's pixel buffer.
             */
            byte* get_next_row(size_t pos)
            {
                typedef typename pixbuf::row_traits row_traits;
                return reinterpret_cast< byte* >
                    (row_traits::get_data(m_pixbuf.get_row(pos)));
            }

        protected:
            pixbuf& m_pixbuf;
        };

        /**
         * \brief The pixel buffer adapter for reading pixel data.
         */
        class pixel_consumer
            : public streaming_impl< consumer< pixel,
                                               pixel_consumer,
                                               image_info_ref_holder,
                                               /* interlacing = */ true > >
        {
        public:
            pixel_consumer(image_info& info, pixbuf& pixels)
                : streaming_impl< consumer< pixel,
                                            pixel_consumer,
                                            image_info_ref_holder,
                                            true > >(info, pixels)
            {
            }

            void reset(size_t pass)
            {
                if (pass == 0)
                {
                    this->m_pixbuf.resize(this->get_info().get_width(),
                                          this->get_info().get_height());
                }
            }
        };

        /**
         * \brief The pixel buffer adapter for writing pixel data.
         */
        class pixel_generator
            : public streaming_impl< generator< pixel,
                                                pixel_generator,
                                                image_info_ref_holder,
                                                /* interlacing = */ true > >
        {
        public:
            pixel_generator(image_info& info, pixbuf& pixels)
                : streaming_impl< generator< pixel,
                                             pixel_generator,
                                             image_info_ref_holder,
                                             true > >(info, pixels)
            {
            }
        };

        image_info m_info;
        pixbuf m_pixbuf;
    };

} // namespace png

#endif // PNGPP_IMAGE_HPP_INCLUDED
