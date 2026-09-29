#include <dlfcn.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char *library_path;
static void *library;
static int (*probe)(void);

static void *worker(void *unused) {
    (void)unused;
    if (!library) {
        library = dlopen(library_path, RTLD_NOW | RTLD_LOCAL);
        if (!library) {
            fprintf(stderr, "%s\n", dlerror());
            exit(2);
        }
        *(void **)(&probe) = dlsym(library, "probe");
        if (!probe) exit(3);
    }
    if (probe() != 0) exit(4);
    printf("thread %lu: OK\n", (unsigned long)pthread_self());
    fflush(stdout);
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 2 || argc > 3) {
        fprintf(stderr, "usage: %s LIBRARY [main]\n", argv[0]);
        return 2;
    }
    library_path = argv[1];
    if (argc == 3) {
        if (strcmp(argv[2], "main") != 0) return 2;
        worker(NULL);
    }
    for (int i = 0; i < 8; ++i) {
        pthread_t thread;
        if (pthread_create(&thread, NULL, worker, NULL) != 0) return 5;
        if (pthread_join(thread, NULL) != 0) return 6;
    }
    // Keep the library loaded through thread cleanup; no dlclose is needed.
    return 0;
}
