/*

   This file is part of the regenie software package.

   Copyright (c) 2020-2024 Joelle Mbatchou, Andrey Ziyatdinov & Jonathan Marchini

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

*/

#include <cstdio>
#include <stdexcept>
#include <string>

#include "S3_Utils.hpp"

#ifdef WITH_S3
#include "bgen/s3stream.hpp"
#endif

bool is_remote_path(std::string const& path) {
#ifdef WITH_S3
  return s3stream_is_remote(path.c_str()) != 0;
#else
  return path.compare(0, 5, "s3://") == 0 || path.compare(0, 7, "http://") == 0 ||
    path.compare(0, 8, "https://") == 0;
#endif
}

#ifdef WITH_S3

void remote_io_init() {
  if(s3stream_init() != 0)
    throw std::string("cannot initialize remote file support : ") + s3stream_last_error();
}

std::unique_ptr<std::streambuf> open_remote_streambuf(std::string const& path) {

  try {
    return std::unique_ptr<std::streambuf>( new s3stream::StreamBuf(path) );
  } catch (std::runtime_error const& e) { // rethrow in regenie's error idiom
    throw "cannot read remote file : " + path + " (" + e.what() + ")";
  }
}

#endif // WITH_S3
