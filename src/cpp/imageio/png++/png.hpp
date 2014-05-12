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
#ifndef PNGPP_PNG_HPP_INCLUDED
#define PNGPP_PNG_HPP_INCLUDED

#include <png.h>

#include "config.hpp"
#include "types.hpp"
#include "error.hpp"
#include "color.hpp"
#include "palette.hpp"
#include "tRNS.hpp"
#include "packed_pixel.hpp"
#include "rgb_pixel.hpp"
#include "rgba_pixel.hpp"
#include "gray_pixel.hpp"
#include "ga_pixel.hpp"
#include "index_pixel.hpp"
#include "info_base.hpp"
#include "info.hpp"
#include "end_info.hpp"
#include "io_base.hpp"
#include "reader.hpp"
#include "writer.hpp"
#include "generator.hpp"
#include "consumer.hpp"
#include "pixel_buffer.hpp"
#include "require_color_space.hpp"
#include "convert_color_space.hpp"
#include "image.hpp"

/**
 * \mainpage
 *
 * \section sec_intro Introduction
 *
 * This is the documentation for png++ the C++ wrapper for libpng.
 * This page documents png++ version 0.2.
 *
 * PNG++ aims to provide simple yet powerful C++ interface to libpng,
 * the PNG reference implementation library.  PNG++ is free software
 * distributed under a modified variant of BSD <a
 * href="http://www.nongnu.org/pngpp/license.html">license</a>.
 *
 * \section sec_news News
 *
 * - Added support for tRNS chunk.
 * - Added non-std IO streams support.
 * - Fixed 16-bit endianness problems.
 * - Improved test script.
 *
 * \section sec_getting_started Getting started
 *
 * The following code demonstrates how to read and write PNG images
 * using png++:
 *
 * \code
 * png::image< png::rgb_pixel > image("input.png");
 * image.write("output.png");
 * \endcode
 *
 * The code reads an image from the file named \c input.png, then
 * writes the image to a file named \c output.png.  The image class
 * template allows you to specify the desired pixel type for the image
 * data.  The available pixel types include: RGB, Grayscale and
 * Indexed pixels.  Some of the pixel types can have an alpha channel.
 *
 * The png++ naturally supports reading PNG images of any %color type
 * into RGB or Grayscale pixel buffers (with optional alpha channel).
 * This is particularly useful for an image viewer if it needs to
 * display the PNG image on, for example, RGB device regardless of the
 * image %color type.
 *
 * On the other hand one might want to read only images of particular
 * type.  With png++ you can specify it this way:
 *
 * \code
 * png::image< png::rgb_pixel > image("rgb.png", png::require_color_space< png::rgb_pixel >());
 * \endcode
 *
 * \section sec_installing Installing
 *
 * PNG++ comes as a set of header files and does not require
 * compilation to be installed.  For the same reason there are no
 * binary packages for png++.
 *
 * \subsection subs_prerequisites Prerequisites
 *
 * - png++ works with libpng-1.2.x.
 * - png++ compiles with g++-4.1 and g++-4.2.  Other version should
 *   work well too.
 * - png++ relies on GNU make for compiling tests and examples; in
 *   particular it uses "remaking makefiles" feature
 * - Documentation is produced using <a
 *   href="http://www.doxygen.org/">doxygen</a>.  See the bottom of
 *   this page for doxygen version used to compile these docs.
 *
 * \subsection subs_installing_pngpp Installing png++
 *
 * Follow these instructions in order to install png++:
 *
 * -# Unpack source package:
 * \verbatim $ tar -zxf png++-0.2.x.tar.gz -C ~/src \endverbatim
 * -# Go to your brand new png++ sources directory:
 * \verbatim $ cd ~/src/png++-0.2.x \endverbatim
 * -# Issue \c make to test how it's doing:
 * \verbatim $ make \endverbatim
 * This will compile examples in the \c example directory.  If everything
 * goes well, try \verbatim $ make test \endverbatim (or \c make \c
 * check which is the same as above) to run the test suite.  If tests
 * do not produce error messages then probably all is OK.
 * -# Now you can create documentation (optional).  Use
 * \verbatim $ make docs \endverbatim
 * to run \c doxygen in the sources directory.
 * -# Now it is time to become \c root and install png++ into your
 * system.  It's OK to issue \c make \c install under ordinary user
 * permissions if you want to install png++ into your home
 * directory.  Run the following command:
 * \verbatim $ make install PREFIX=$HOME \endverbatim
 * to copy png++ header files to \c ~/include/png++ and documentation
 * files to <tt>~/share/doc/png++-0.2.x</tt>.  Without a \c PREFIX png++
 * installs to \c /usr/local.
 *
 * \section sec_working_with_images Working with images
 *
 * In png++ you can create new images like this:
 *
 * \code
 * #include <png++/png.hpp>
 * //...
 * png::image< png::rgb_pixel > image(128, 128);
 * for (size_t y = 0; y < image.get_height(); ++y)
 * {
 *     for (size_t x = 0; x < image.get_width(); ++x)
 *     {
 *         image[y][x] = png::rgb_pixel(x, y, x + y);
 *         // non-checking equivalent of image.set_pixel(x, y, ...);
 *     }
 * }
 * image.write("rgb.png");
 * \endcode
 *
 * Optionally, you may specify
 * \code
 * image.set_interlace_type(png::interlace_adam7)
 * \endcode
 * to produce an interlaced image.
 *
 * If you are writing an indexed colors image, you should provide a
 * palette (colormap).  One of the ways to do this is the following:
 *
 * \code
 * #include <png++/png.hpp>
 * //...
 * png::image< png::index_pixel > image;
 * png::palette pal(256);
 * for (size_t i = 0; i < pal.size(); ++i)
 * {
 *     pal[i] = png::color(i, 255 - i, i);
 * }
 * image.set_palette(pal);
 * ...
 * image.write("palette.png");
 * \endcode
 *
 * It is not absolutely necessary to have the whole image data in
 * memory in order to write a PNG file.  You can use generator class
 * template to write the image row-by-row.  An example of this is
 * provided in \c example/pixel_generator.cpp bundled with the sources
 * package.
 *
 * The same holds for reading images too.  You can use consumer class
 * template in order to read the image data row-by-row.  This might
 * help in applications which have to deal with large PNG images but
 * do not want to read the entire image into memory.
 *
 * You can read or write images from/to generic IO stream, not only
 * file on disk.  Check out \c image::read(std::istream&),
 * \c image::write(std::ostream&) overloads in the reference manual.
 *
 * \section sec_compiling_user Compiling your programs
 *
 * Use the following command to compile your program:
 *
 * \verbatim $ g++ -c example.cpp `libpng-config --cflags` \endverbatim
 *
 * and the following to link it:
 *
 * \verbatim $ g++ -o example example.o `libpng-config --ldflags` \endverbatim
 *
 * When compiling you should add \c -I \c $PREFIX/include if you have
 * installed png++ to non-standard location, like your home directory.
 *
 * In your program, the line
 *
 * \code
 * #include <png++/png.hpp>
 * \endcode
 *
 * brings in all the header files in png++ which should be suitable
 * for the most of the applications.  You may include only the headers
 * you really use, for example:
 *
 * \code
 * #include <png++/image.hpp>
 * #include <png++/rgb_pixel.hpp>
 * \endcode
 *
 * If do not want to install png++ headers you still can compile your
 * programs.  Just create a subdirectory named \c png++ somewhere in
 * your project tree and copy all of the \c .hpp files in png++
 * distribution there.  Then use appropriate compiler options to add
 * this directory into the header search path.
 *
 * \section sec_further Further reading
 *
 * - To get information about specific features, use reference (can be
 *   reached from the top of this page).
 * - If you are looking for more example code, please go to the \c
 *   example/ directory of the source distribution.  You may also find
 *   sources in the \c test/ directory insightful (well, somewhat).
 * - Have a question?  Check out \ref sec_help section below.
 * - Of course, your ultimate source for learning is the source
 *   code. <tt>:-)</tt>
 *
 * \section sec_download Download
 *
 * The project is hosted at Savannah:
 * http://savannah.nongnu.org/projects/pngpp/
 *
 * Released source packages can be found here:
 * http://download.savannah.nongnu.org/releases/pngpp/
 *
 * Also, you can check out sources directly from SVN repository:
 * svn://svn.sv.nongnu.org/pngpp/trunk/ or
 * http://svn.sv.nongnu.org/pngpp/trunk/ (for people w/o outgoing svn).
 *
 * Online version of this documentation can be found here:
 * http://www.nongnu.org/pngpp/doc/0.2.3/index.html
 *
 * \section sec_bugs Bugs
 *
 * The following is a list of known bugs and limitations:
 *
 * - Lacks support for output transformations
 * - Lacks support for optional/unknown chunks in PNG data stream
 * - Documentation sucks <tt>;-)</tt>
 *
 * To report bugs, please use Savannah bug tracker:
 * http://savannah.nongnu.org/bugs/?group=pngpp&func=additem
 *
 * Do not forget to check if the bug was already filed. <tt>:-)</tt>
 *
 * \section sec_help Getting help
 *
 * There is a mailing list for developers:
 * http://lists.nongnu.org/mailman/listinfo/pngpp-devel
 *
 * You can also contact me by dropping a mail to <alex.shulgin@gmail.com>.
 *
 * Happy hacking!
 *
 * Alex Shulgin
 */

#endif // PNGPP_PNG_HPP_INCLUDED
