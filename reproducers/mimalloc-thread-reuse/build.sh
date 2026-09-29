#!/usr/bin/env bash
# SOURCE is a mimalloc C source root containing include/ and src/static.c.
set -euo pipefail
if [[ $# != 2 ]]; then
    echo "usage: $0 MIMALLOC_SOURCE OUTPUT_DIRECTORY" >&2
    exit 2
fi
source_root=$(cd "$1" && pwd)
mkdir -p "$2"
output_root=$(cd "$2" && pwd)
probe_root=$(cd "$(dirname "$0")" && pwd)
compiler=${CC:-cc}
"$compiler" -shared -fPIC -g -O0 -DMI_DEBUG=0 -ftls-model=local-dynamic \
    -I"$source_root/include" "$source_root/src/static.c" "$probe_root/library.c" \
    -o "$output_root/libprobe.so" -lpthread
"$compiler" -g -O0 "$probe_root/main.c" -o "$output_root/driver" -ldl -lpthread
