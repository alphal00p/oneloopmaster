use std::hint::black_box;

#[global_allocator]
static ALLOCATOR: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[no_mangle]
pub extern "C" fn probe() -> i32 {
    let bytes = black_box(Box::new([42u8; 64]));
    i32::from(bytes[0] != 42 || bytes[63] != 42)
}
