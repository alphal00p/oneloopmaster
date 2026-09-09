#!/usr/bin/env bash
# Run from an isolated virtual environment with Maturin installed.
set -euo pipefail
host_source_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
maturin develop --release --manifest-path "$host_source_dir/Cargo.toml"
cd -- "$host_source_dir/../.."
ONELOOP_PYTHON_MODULE=symbolica.community.oneloop \
  python -m unittest discover -s python/tests -v
