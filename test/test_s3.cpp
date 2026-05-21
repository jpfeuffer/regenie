/*
  Standalone test for AWS S3 support in regenie.
  Tests downloading public S3 files using the AWS SDK.
  
  Public test files used:
  - s3://1000genomes/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
    (small text file from the 1000 Genomes project, ~3KB)
*/

#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <cstdio>

#include "S3_Utils.hpp"

// plink2 S3 functions (bundled in external_libs/pgenlib/plink2_s3.h)
namespace plink2 {
  bool IsS3Uri(const char* path);
  FILE* OpenMaybeS3(const char* path);
  void S3InitClientOnly();
  void S3ShutdownClientOnly();
}

int main(int argc, char** argv) {
  int tests_passed = 0;
  int tests_failed = 0;

  // Initialize AWS SDK once
  aws_sdk_init();
  // Initialize only the S3 client (SDK already initialized above)
  plink2::S3InitClientOnly();

  // Test 1: is_s3_path detection
  std::cout << "Test 1: is_s3_path() detection..." << std::endl;
  {
    assert(is_s3_path("s3://bucket/key") == true);
    assert(is_s3_path("s3://mybucket/path/to/file.txt") == true);
    assert(is_s3_path("/local/path/file.txt") == false);
    assert(is_s3_path("./relative/path") == false);
    assert(is_s3_path("s3:/") == false);
    assert(is_s3_path("") == false);
    std::cout << "  PASSED" << std::endl;
    tests_passed++;
  }

  // Test 2: parse_s3_uri
  std::cout << "Test 2: parse_s3_uri()..." << std::endl;
  {
    std::string bucket, key;
    
    parse_s3_uri("s3://mybucket/path/to/file.txt", bucket, key);
    assert(bucket == "mybucket");
    assert(key == "path/to/file.txt");

    parse_s3_uri("s3://1000genomes/release/20130502/file.panel", bucket, key);
    assert(bucket == "1000genomes");
    assert(key == "release/20130502/file.panel");

    parse_s3_uri("s3://bucket-only", bucket, key);
    assert(bucket == "bucket-only");
    assert(key == "");
    
    std::cout << "  PASSED" << std::endl;
    tests_passed++;
  }

  // Test 3: resolve_s3_path with local path (no-op)
  std::cout << "Test 3: resolve_s3_path() with local path..." << std::endl;
  {
    std::string local = "/tmp/some_local_file.txt";
    std::string resolved = resolve_s3_path(local);
    assert(resolved == local);
    std::cout << "  PASSED" << std::endl;
    tests_passed++;
  }

  // Test 4: Download a small public file from 1000 Genomes S3 bucket
  std::cout << "Test 4: Download public S3 file (1000genomes panel)..." << std::endl;
  {
    // This is a small (~3KB) text file listing sample IDs and populations
    std::string s3_uri = "s3://1000genomes/release/20130502/integrated_call_samples_v3.20130502.ALL.panel";
    
    try {
      std::string local_path = download_s3_to_local(s3_uri);
      std::cout << "  Downloaded to: " << local_path << std::endl;

      // Verify file exists and has content
      std::ifstream f(local_path.c_str());
      assert(f.is_open());

      std::string first_line;
      std::getline(f, first_line);
      f.close();

      std::cout << "  First line: " << first_line << std::endl;
      // The panel file has a header line with "sample", "pop", "super_pop", "gender"
      assert(first_line.find("sample") != std::string::npos);
      
      std::cout << "  PASSED" << std::endl;
      tests_passed++;
    } catch (const std::string& e) {
      std::cerr << "  FAILED: " << e << std::endl;
      tests_failed++;
    } catch (const std::exception& e) {
      std::cerr << "  FAILED: " << e.what() << std::endl;
      tests_failed++;
    }
  }

  // Test 5: Download another public file - README from 1000 Genomes
  std::cout << "Test 5: Download public S3 file (1000genomes README)..." << std::endl;
  {
    std::string s3_uri = "s3://1000genomes/README.sequence_data";

    try {
      std::string local_path = download_s3_to_local(s3_uri);
      std::cout << "  Downloaded to: " << local_path << std::endl;

      std::ifstream f(local_path.c_str());
      assert(f.is_open());

      // Check file is non-empty
      f.seekg(0, std::ios::end);
      std::streamsize size = f.tellg();
      assert(size > 0);
      std::cout << "  File size: " << size << " bytes" << std::endl;

      std::cout << "  PASSED" << std::endl;
      tests_passed++;
    } catch (const std::string& e) {
      std::cerr << "  FAILED: " << e << std::endl;
      tests_failed++;
    } catch (const std::exception& e) {
      std::cerr << "  FAILED: " << e.what() << std::endl;
      tests_failed++;
    }
  }

  // Test 6: Invalid S3 URI should throw
  std::cout << "Test 6: Invalid S3 URI error handling..." << std::endl;
  {
    bool caught = false;
    try {
      download_s3_to_local("s3://nonexistent-bucket-xyz123abc/no-such-key");
    } catch (const std::string& e) {
      caught = true;
      std::cout << "  Caught expected error: " << e.substr(0, 80) << "..." << std::endl;
    } catch (const std::exception& e) {
      caught = true;
      std::cout << "  Caught expected error: " << std::string(e.what()).substr(0, 80) << "..." << std::endl;
    }
    if (caught) {
      std::cout << "  PASSED" << std::endl;
      tests_passed++;
    } else {
      std::cerr << "  FAILED: No exception thrown for invalid S3 URI" << std::endl;
      tests_failed++;
    }
  }

  // Test 7: plink2::IsS3Uri() detection
  std::cout << "Test 7: plink2::IsS3Uri() detection..." << std::endl;
  {
    assert(plink2::IsS3Uri("s3://bucket/key") == true);
    assert(plink2::IsS3Uri("s3://mybucket/path/to/file.pgen") == true);
    assert(plink2::IsS3Uri("/local/path/file.pgen") == false);
    assert(plink2::IsS3Uri("./relative/path.pgen") == false);
    std::cout << "  PASSED" << std::endl;
    tests_passed++;
  }

  // Test 8: plink2::OpenMaybeS3() with a local file (falls back to fopen)
  std::cout << "Test 8: plink2::OpenMaybeS3() with local file..." << std::endl;
  {
    // Use an existing local file from the example data or create a temp file
    const char* tmpfile = "/tmp/regenie_pgenlib_s3_test.txt";
    FILE* fw = fopen(tmpfile, "w");
    if (fw) {
      fputs("hello pgenlib s3 test\n", fw);
      fclose(fw);
    }

    FILE* fh = plink2::OpenMaybeS3(tmpfile);
    if (fh != nullptr) {
      char buf[64];
      char* res = fgets(buf, sizeof(buf), fh);
      fclose(fh);
      if (res != nullptr && std::string(buf).find("hello") != std::string::npos) {
        std::cout << "  PASSED (read: " << std::string(buf).substr(0, 30) << ")" << std::endl;
        tests_passed++;
      } else {
        std::cerr << "  FAILED: Could not read expected content" << std::endl;
        tests_failed++;
      }
    } else {
      std::cerr << "  FAILED: OpenMaybeS3 returned nullptr for local file" << std::endl;
      tests_failed++;
    }
    remove(tmpfile);
  }

  // Test 9: plink2::OpenMaybeS3() with a real S3 pgen URI (streaming range requests)
  // Uses a small public pgen from nf-core test data (ngi-igenomes bucket, eu-west-1).
  // Requires AWS_DEFAULT_REGION=eu-west-1 (or equivalent config).
  // No credentials needed — anonymous access, equivalent to --no-sign-request.
  std::cout << "Test 9: plink2::OpenMaybeS3() streaming from S3 pgen..." << std::endl;
  {
    const char* s3_pgen = "s3://ngi-igenomes/testdata/nf-core/modules/genomics/"
                          "homo_sapiens/popgen/plink_simulated.pgen";
    FILE* fh = plink2::OpenMaybeS3(s3_pgen);
    if (fh != nullptr) {
      // pgen magic: first two bytes must be 0x6C 0x1B
      unsigned char magic[2] = {0, 0};
      size_t nr = fread(magic, 1, 2, fh);
      fclose(fh);
      if (nr == 2 && magic[0] == 0x6C && magic[1] == 0x1B) {
        std::cout << "  PASSED (pgen magic: 0x"
                  << std::hex << (int)magic[0] << " 0x" << (int)magic[1]
                  << std::dec << ")" << std::endl;
        tests_passed++;
      } else {
        std::cerr << "  FAILED: Unexpected header bytes (read " << nr
                  << " bytes, got 0x" << std::hex << (int)magic[0]
                  << " 0x" << (int)magic[1] << std::dec << ")" << std::endl;
        tests_failed++;
      }
    } else {
      std::cerr << "  FAILED: OpenMaybeS3 returned nullptr for S3 URI" << std::endl;
      tests_failed++;
    }
  }

  // Summary
  std::cout << "\n========================================" << std::endl;
  std::cout << "Results: " << tests_passed << " passed, " << tests_failed << " failed" << std::endl;
  std::cout << "========================================" << std::endl;

  // Cleanup
  plink2::S3ShutdownClientOnly();
  aws_sdk_shutdown();

  return tests_failed > 0 ? 1 : 0;
}
