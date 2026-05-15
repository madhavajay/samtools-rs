//! glibc-compatible POSIX 48-bit `srand48`/`lrand48` LCG plus HTSlib's
//! `gen_unique_id` collision-suffixing, used by `samtools merge -s SEED`
//! for `@RG`/`@PG` ID reconciliation.
//!
//! HTSlib's `hts_srand48`/`hts_lrand48` delegate to the platform
//! `srand48`/`lrand48` (glibc on Linux), the standard POSIX 48-bit linear
//! congruential generator: `Xₙ₊₁ = (A·Xₙ + C) mod 2⁴⁸`, `A=0x5DEECE66D`,
//! `C=0xB`; `srand48(s)` sets `X = (s<<16)|0x330E`; `lrand48()` returns the
//! top 31 bits (`X >> 17`).

use std::collections::HashSet;

const A: u64 = 0x5DEECE66D;
const C: u64 = 0xB;
const MASK: u64 = (1u64 << 48) - 1;

/// A deterministic `drand48`/`lrand48`-family generator.
#[derive(Clone, Debug)]
pub struct Rand48 {
    state: u64,
}

impl Rand48 {
    /// `srand48(seed)`.
    pub fn new(seed: i64) -> Self {
        Self {
            state: (((seed as u64) << 16) | 0x330E) & MASK,
        }
    }

    fn step(&mut self) -> u64 {
        self.state = (A.wrapping_mul(self.state).wrapping_add(C)) & MASK;
        self.state
    }

    /// `lrand48()` — a non-negative long in `[0, 2³¹)`.
    pub fn lrand48(&mut self) -> u32 {
        (self.step() >> 17) as u32
    }
}

/// Port of HTSlib `bam_sort.c` `gen_unique_id`: if `prefix` is unused,
/// return it as-is (no PRNG draw); otherwise loop drawing `lrand48()` and
/// trying `"{prefix}-{:08X}"` until one is unused. The chosen id is
/// inserted into `existing`.
pub fn gen_unique_id(prefix: &str, existing: &mut HashSet<String>, rng: &mut Rand48) -> String {
    if !existing.contains(prefix) {
        existing.insert(prefix.to_string());
        return prefix.to_string();
    }
    loop {
        let candidate = format!("{prefix}-{:08X}", rng.lrand48());
        if !existing.contains(&candidate) {
            existing.insert(candidate.clone());
            return candidate;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_glibc_seed_1_sequence() {
        // Verified against the upstream `merge -s 1` (`merge/2`) fixture.
        let mut r = Rand48::new(1);
        let seq: Vec<String> = (0..8).map(|_| format!("{:08X}", r.lrand48())).collect();
        assert_eq!(
            seq,
            [
                "055424A4", "3A2CCEF5", "6ADB4A65", "2B019719", "4861F4EF", "0039E5EF", "1802EEEC",
                "7EC68B3F"
            ]
        );
    }

    #[test]
    fn gen_unique_id_bare_then_suffixed() {
        let mut ids = HashSet::new();
        let mut r = Rand48::new(1);
        assert_eq!(gen_unique_id("fish", &mut ids, &mut r), "fish"); // bare, no draw
        // Collision → first draw 055424A4.
        assert_eq!(gen_unique_id("fish", &mut ids, &mut r), "fish-055424A4");
        // A fresh prefix is still bare (no draw consumed for it).
        assert_eq!(gen_unique_id("cow", &mut ids, &mut r), "cow");
        // Next collision draws the 2nd value.
        assert_eq!(gen_unique_id("fish", &mut ids, &mut r), "fish-3A2CCEF5");
    }
}
