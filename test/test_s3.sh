#!/usr/bin/env bash
# Test script for AWS S3 support in regenie
# Builds and runs the S3 unit test against public S3 data.
# pgenlib S3 streaming tests (tests 7-9) are always included.
#
# Public S3 test files used:
#   s3://1000genomes/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
#   s3://1000genomes/README.sequence_data
#   s3://ngi-igenomes/testdata/nf-core/modules/genomics/homo_sapiens/popgen/plink_simulated.pgen
#
# Requirements:
#   - AWS SDK for C++ installed (e.g., brew install aws-sdk-cpp)
#   - No AWS credentials needed — tests use public buckets (anonymous access)
#   - AWS_DEFAULT_REGION must be set to eu-west-1 for test 9 (ngi-igenomes bucket)
#     e.g., AWS_DEFAULT_REGION=eu-west-1 bash test_s3.sh

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BUILD_DIR="${SCRIPT_DIR}/build_s3_test"

echo "=== Regenie S3 Support Test ==="
echo ""

# Clean and create build directory
rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}"

# Configure
echo "Configuring..."
cmake -S "${SCRIPT_DIR}" -B "${BUILD_DIR}" \
  -DCMAKE_PREFIX_PATH=/opt/homebrew 2>&1 | grep -E "^--|error|S3" || true
echo ""

# Build
echo "Building..."
cmake --build "${BUILD_DIR}" 2>&1
echo ""

# Run tests
echo "Running S3 tests..."
echo ""
"${BUILD_DIR}/test_s3"

EXIT_CODE=$?
echo ""

# Cleanup temp files
rm -f /tmp/regenie_s3_* /tmp/regenie_pgenlib_s3_test.txt

if [ ${EXIT_CODE} -eq 0 ]; then
  echo "All S3 tests PASSED."
else
  echo "Some S3 tests FAILED."
fi

exit ${EXIT_CODE}
