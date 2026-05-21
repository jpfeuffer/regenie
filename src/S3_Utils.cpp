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

#ifdef WITH_AWS_S3

#include <iostream>
#include <fstream>
#include <cstdlib>

#include <aws/core/Aws.h>
#include <aws/s3/S3Client.h>
#include <aws/s3/model/GetObjectRequest.h>
#include <aws/core/auth/AWSCredentialsProvider.h>

#include "S3_Utils.hpp"

static Aws::SDKOptions s3_sdk_options;
static bool s3_sdk_initialized = false;

bool is_s3_path(const std::string& path) {
  return path.size() > 5 && path.substr(0, 5) == "s3://";
}

void parse_s3_uri(const std::string& uri, std::string& bucket, std::string& key) {
  // uri format: s3://bucket/key
  std::string path = uri.substr(5); // remove "s3://"
  size_t slash_pos = path.find('/');
  if (slash_pos == std::string::npos) {
    bucket = path;
    key = "";
  } else {
    bucket = path.substr(0, slash_pos);
    key = path.substr(slash_pos + 1);
  }
}

void aws_sdk_init() {
  if (!s3_sdk_initialized) {
    Aws::InitAPI(s3_sdk_options);
    s3_sdk_initialized = true;
  }
}

void aws_sdk_shutdown() {
  if (s3_sdk_initialized) {
    Aws::ShutdownAPI(s3_sdk_options);
    s3_sdk_initialized = false;
  }
}

std::string download_s3_to_local(const std::string& s3_uri) {
  aws_sdk_init();

  std::string bucket, key;
  parse_s3_uri(s3_uri, bucket, key);

  if (bucket.empty() || key.empty()) {
    throw std::string("Invalid S3 URI: ") + s3_uri;
  }

  // Create a temporary file path preserving the original extension
  std::string filename = key;
  size_t last_slash = filename.rfind('/');
  if (last_slash != std::string::npos) {
    filename = filename.substr(last_slash + 1);
  }

  // Use tmpdir
  std::string tmpdir = "/tmp";
  const char* env_tmp = std::getenv("TMPDIR");
  if (env_tmp) tmpdir = env_tmp;

  std::string local_path = tmpdir + "/regenie_s3_" + filename;

  // Configure S3 client - try with default credential chain first
  Aws::S3::S3ClientConfiguration config;
  config.region = Aws::Region::US_EAST_1;

  Aws::S3::S3Client s3_client(config);

  Aws::S3::Model::GetObjectRequest request;
  request.SetBucket(bucket.c_str());
  request.SetKey(key.c_str());

  auto outcome = s3_client.GetObject(request);

  if (!outcome.IsSuccess()) {
    // If default credentials failed, try anonymous access for public buckets
    Aws::S3::S3ClientConfiguration anon_config;
    anon_config.region = Aws::Region::US_EAST_1;

    auto creds_provider = Aws::MakeShared<Aws::Auth::AnonymousAWSCredentialsProvider>("S3Alloc");
    Aws::S3::S3Client anon_client(creds_provider);

    outcome = anon_client.GetObject(request);
    if (!outcome.IsSuccess()) {
      throw std::string("Failed to download S3 object: ") + s3_uri +
        " Error: " + outcome.GetError().GetMessage().c_str();
    }
  }

  // Write to local file
  std::ofstream local_file(local_path.c_str(), std::ios::binary);
  if (!local_file.is_open()) {
    throw std::string("Cannot create temporary file: ") + local_path;
  }

  auto& body = outcome.GetResult().GetBody();
  local_file << body.rdbuf();
  local_file.close();

  return local_path;
}

std::string resolve_s3_path(const std::string& path) {
  if (is_s3_path(path)) {
    return download_s3_to_local(path);
  }
  return path;
}

#endif // WITH_AWS_S3
