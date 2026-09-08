//! Read-only proof for SymJIT 2.24.1 Builder::save/load serialization.
//! Decodes the documented Builder magic/fields/function-name HashSet; never
//! decodes or changes instructions, numeric constants, or other metadata.
use std::{collections::BTreeSet, fs, ops::Range};

#[derive(Debug)]
struct Table {
    builder: usize,
    entries: Range<usize>,
    names: Vec<String>,
}

fn u64_at(data: &[u8], offset: usize) -> usize {
    u64::from_le_bytes(data[offset..offset + 8].try_into().unwrap())
        .try_into()
        .unwrap()
}

fn tables(data: &[u8]) -> Vec<Table> {
    let magic = 0x12f21e25abe627bcu64.to_le_bytes();
    data.windows(8)
        .enumerate()
        .filter(|(_, bytes)| *bytes == magic)
        .map(|(builder, _)| {
            // magic, count_loops, count_stack, num_consts, f64 constants,
            // num_ft, then each UTF-8 name preceded by its one-byte length.
            let num_consts = u64_at(data, builder + 24);
            let count_offset = builder + 32 + num_consts * 8;
            let count = u64_at(data, count_offset);
            assert!(count < 100_000);
            let start = count_offset + 8;
            let mut cursor = start;
            let mut names = Vec::new();
            for _ in 0..count {
                let length = usize::from(data[cursor]);
                cursor += 1;
                let name = std::str::from_utf8(&data[cursor..cursor + length]).unwrap();
                assert!(name.is_ascii() && !name.is_empty());
                names.push(name.to_owned());
                cursor += length;
            }
            assert_eq!(names.iter().collect::<BTreeSet<_>>().len(), names.len());
            Table { builder, entries: start..cursor, names }
        })
        .collect()
}

fn canonicalize(data: &[u8], tables: &[Table]) -> Vec<u8> {
    let mut normalized = data.to_vec();
    for table in tables {
        let mut names = table.names.clone();
        names.sort();
        let mut encoded = Vec::new();
        for name in names {
            encoded.push(name.len().try_into().unwrap());
            encoded.extend_from_slice(name.as_bytes());
        }
        assert_eq!(encoded.len(), table.entries.len());
        normalized[table.entries.clone()].copy_from_slice(&encoded);
    }
    normalized
}

fn main() {
    for family in ["a0", "b0", "db0", "c0", "d0"] {
        let original = fs::read(format!(
            "/common/dev/oneloop/oneloopmaster/assets/evaluators/{family}.bin"
        )).unwrap();
        let regenerated = fs::read(format!(
            "/tmp/oneloop-assets-final.ISK1VL/assets/{family}.bin"
        )).unwrap();
        assert_eq!(original.len(), regenerated.len());
        let before = tables(&original);
        let after = tables(&regenerated);
        assert_eq!(before.len(), after.len());
        let mut permuted = Vec::new();
        for (old, new) in before.iter().zip(&after) {
            assert_eq!(old.builder, new.builder);
            assert_eq!(old.entries, new.entries);
            assert_eq!(
                old.names.iter().collect::<BTreeSet<_>>(),
                new.names.iter().collect::<BTreeSet<_>>()
            );
            if old.names != new.names {
                permuted.push((old.entries.clone(), old.names.len()));
            }
        }
        let changed = original.iter().zip(&regenerated).enumerate()
            .filter(|(_, (a, b))| a != b).map(|(i, _)| i).collect::<Vec<_>>();
        assert!(changed.iter().all(|i| permuted.iter().any(|(range, _)| range.contains(i))));
        assert_eq!(canonicalize(&original, &before), canonicalize(&regenerated, &after));
        println!(
            "{family}: {} bytes; {} Builder tables; {} differing bytes; permuted entry spans (zero-based, end-exclusive) {permuted:?}; exact name sets and all non-order bytes IDENTICAL",
            original.len(), before.len(), changed.len()
        );
    }
}
