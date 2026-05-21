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

#ifdef WITH_AWS_S3

#include <string>

// Check if a path is an S3 URI (starts with s3://)
bool is_s3_path(const std::string& path);

// Parse an S3 URI into bucket and key components
void parse_s3_uri(const std::string& uri, std::string& bucket, std::string& key);

// Download an S3 object to a local temporary file, returns the local path
std::string download_s3_to_local(const std::string& s3_uri);

// Initialize the AWS SDK (call once at program start)
void aws_sdk_init();

// Shutdown the AWS SDK (call once at program end)
void aws_sdk_shutdown();

// Resolve a path: if it's an S3 URI, download it and return local path;
// otherwise return the original path unchanged
std::string resolve_s3_path(const std::string& path);

#endif // WITH_AWS_S3

#endif // S3_UTILS_H
