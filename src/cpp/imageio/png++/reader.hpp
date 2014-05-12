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
#ifndef PNGPP_READER_HPP_INCLUDED
#define PNGPP_READER_HPP_INCLUDED

#include <cassert>
#include "io_base.hpp"

namespace png
{

    /**
     * \brief The PNG reader class template.  This is the low-level
     * reading interface--use image class or consumer class to
     * actually read images.
     *
     * The \c istream template parameter specifies the type of input
     * stream to work with.  The \c istream class should implement the
     * minimum of the following interface:
     *
     * \code
     * class my_istream
     * {
     * public:
     *     void read(char*, size_t);
     *     bool good();
     * };
     * \endcode
     *
     * With the semantics similar to the \c std::istream.  Naturally,
     * \c std::istream fits this requirement and can be used with the
     * reader class as is.
     *
     * \see image, consumer, writer, io_base
     */
    template< class istream >
    class reader
        : public io_base
    {
    public:
        /**
         * \brief Constructs a reader prepared to read PNG image from
         * a \a stream.
         */
        explicit reader(istream& stream)
            : io_base(png_create_read_struct(PNG_LIBPNG_VER_STRING,
                                             static_cast< io_base* >(this),
                                             raise_error,
                                             0))
        {
            png_set_read_fn(m_png, & stream, read_data);
        }

        ~reader()
        {
            png_destroy_read_struct(& m_png,
                                    m_info.get_png_info_ptr(),
                                    m_end_info.get_png_info_ptr());
        }

        /**
         * \brief Reads the whole PNG data stream into memory.  Not
         * particularly useful.
         */
        void read_png()
        {
            if (setjmp(png_jmpbuf(m_png)))
            {
                throw error(m_error);
            }
            png_read_png(m_png,
                         m_info.get_png_info(),
                         /* transforms = */ 0,
                         /* params = */ 0);
        }

        /**
         * \brief Reads info about PNG image.
         */
        void read_info()
        {
            if (setjmp(png_jmpbuf(m_png)))
            {
                throw error(m_error);
            }
            m_info.read();
        }

        /**
         * \brief Reads a row of image data at a time.
         */
        void read_row(byte* bytes)
        {
            if (setjmp(png_jmpbuf(m_png)))
            {
                throw error(m_error);
            }
            png_read_row(m_png, bytes, 0);
        }

        /**
         * \brief Reads ending info about PNG image.
         */
        void read_end_info()
        {
            if (setjmp(png_jmpbuf(m_png)))
            {
                throw error(m_error);
            }
            m_end_info.read();
        }

        void update_info()
        {
            m_info.update();
        }

    private:
        static void read_data(png_struct* png, byte* data, size_t length)
        {
            io_base* io = static_cast< io_base* >(png_get_error_ptr(png));
            reader* rd = static_cast< reader* >(io);
            rd->reset_error();
            istream* stream = reinterpret_cast< istream* >(png_get_io_ptr(png));
            try
            {
                stream->read(reinterpret_cast< char* >(data), length);
                if (!stream->good())
                {
                    rd->set_error("istream::read() failed");
                }
            }
            catch (std::exception const& error)
            {
                rd->set_error(error.what());
            }
            catch (...)
            {
                assert(!"read_data: caught something wrong");
                rd->set_error("read_data: caught something wrong");
            }
            if (rd->is_error())
            {
                rd->raise_error();
            }
        }
    };

} // namespace png

#endif // PNGPP_READER_HPP_INCLUDED
