#include "mimalloc.h"
#include <string.h>

int probe(void) {
    unsigned char *p = mi_malloc_aligned(64, 64);
    if (!p) return 1;
    memset(p, 42, 64);
    int result = (p[0] != 42 || p[63] != 42);
    mi_free(p);
    return result;
}
