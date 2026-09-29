# Symbolica main regression audit

Standalone rerun of the six evaluator/startup regression files from the
historical Symbolica patch bundle, on unmodified main
`a19c760dd567c239f30d87e4e924ca2f8b8457ab` and SymJIT 2.26.0.
No OneLOop dependency or source patches. Only the original tests are copied.

```sh
RUST_MIN_STACK=134217728 cargo test --locked --release -- --test-threads=1
```

All ten tests pass (one subprocess helper is ignored by the parent harness and
run separately in three fresh processes). This covers exact/exported constants,
nested callback indices, helper registration in the Symbolica namespace,
restored JIT settings and cold transcendental initialization.
The build-metadata test is not included: it tests Symbolica's own build script
and belongs in the upstream source tree, rather than a dependent crate.
