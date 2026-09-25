//! Zero-dependency SHA-256 against FIPS 180-4 / NIST KAT vectors.

use kronrod::crypto::{sha256, to_hex, Sha256};

#[test]
fn fips_empty_message() {
    let got = to_hex(&sha256(b""));
    assert_eq!(got, "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
}

#[test]
fn fips_abc() {
    let got = to_hex(&sha256(b"abc"));
    assert_eq!(got, "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
}

#[test]
fn fips_quick_fox() {
    let got = to_hex(&sha256(b"The quick brown fox jumps over the lazy dog"));
    assert_eq!(got, "d7a8fbb307d7809469ca9abcb0082e4f8d5651e46d3cdb762d02d0bf37c9e592");
}

#[test]
fn fips_one_million_abc() {
    // NIST KAT: "abc" repeated 1,000,000 times (3 MiB).
    let mut data = Vec::with_capacity(3_000_000);
    for _ in 0..1_000_000 {
        data.extend_from_slice(b"abc");
    }
    let got = to_hex(&sha256(&data));
    assert_eq!(got, "f4096a131e7e6ebfa7a512b5c299e13b065df34d15624ee1202ab394cc4d7e90");
}

#[test]
fn incremental_matches_oneshot() {
    let data: Vec<u8> = (0..1000u16).flat_map(|i| [i as u8, (i >> 8) as u8, 0xA5]).collect();
    let oneshot = sha256(&data);

    let mut h = Sha256::new();
    // 1-byte chunks (exercises the buf fill path at every boundary).
    for b in &data {
        h.update(&[*b]);
    }
    assert_eq!(h.finalize(), oneshot);

    let mut h = Sha256::new();
    // 61-byte chunks (straddles the 64-byte block boundary mid-stream).
    for chunk in data.chunks(61) {
        h.update(chunk);
    }
    assert_eq!(h.finalize(), oneshot);
}

#[test]
fn block_boundary_lengths() {
    // Lengths around the 64-byte block and 55/56-byte padding thresholds.
    for len in [0usize, 1, 54, 55, 56, 63, 64, 65, 127, 128, 129] {
        let data: Vec<u8> = (0..len).map(|i| (i * 7 + 3) as u8).collect();
        let oneshot = sha256(&data);
        let mut h = Sha256::new();
        h.update(&data);
        assert_eq!(h.finalize(), oneshot, "len={len}");
    }
}

#[test]
fn to_hex_roundtrip() {
    let bytes: [u8; 32] = (0..32u8).map(|i| i * 3).collect::<Vec<_>>().try_into().unwrap();
    let hex = to_hex(&bytes);
    assert_eq!(hex.len(), 64);
    assert_eq!(&hex[0..8], "00030609");
    assert!(hex.chars().all(|c| c.is_ascii_hexdigit()));
}

#[test]
fn distinct_inputs_distinct_digests() {
    let a = sha256(b"kronrod 1.0.0");
    let b = sha256(b"kronrod 1.0.0 ");
    assert_ne!(a, b);
}
