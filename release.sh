#!/bin/bash -e

set -o pipefail

echo "Building normally .."
cargo build --release

echo "Testing release version .."
cargo test --release

echo "Now make sure git is up to date, and run cargo publish"
