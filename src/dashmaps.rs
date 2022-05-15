use std::hash::BuildHasherDefault;
use dashmap::DashMap;
use fxhash::FxHasher;

/// A custom `DashMap` w/ `FxHasher`.
///
/// ```use dashmap::DashMap;```
/// ```use fxhash::FxHasher;```
/// ```// skip```
/// ```let dashfx_hash: DashFx = DashMap::with_hasher(BuildHasherDefault::<FxHasher>::default());```
/// Useful: [Using a Custom Hash Function in Rust](https://docs.rs/hashers/1.0.1/hashers/#using-a-custom-hash-function-in-rust).
pub type DashFx = DashMap<u64, i32, BuildHasherDefault<FxHasher>>;

pub type UnpackDashFx = DashMap<Vec<u8>, i32, BuildHasherDefault<FxHasher>>;