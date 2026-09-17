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

#ifndef S3_UTILS_H
#define S3_UTILS_H

#include <memory>
#include <streambuf>
#include <string>

// Check if a path names a remote object ("s3://", "http://" or "https://").
// Always available, so builds without S3 can report why the path was rejected.
bool is_remote_path(std::string const& path);

#ifdef WITH_S3

// Initialize libcurl. Call once at program start, before any thread that may
// open a remote file is spawned.
void remote_io_init();

// Open a remote object for reading. Data is fetched with HTTP range requests
// as the buffer is consumed, so no temporary file is staged and opening a
// multi-gigabyte object costs one small request.
std::unique_ptr<std::streambuf> open_remote_streambuf(std::string const& path);

#endif // WITH_S3

#endif // S3_UTILS_H
