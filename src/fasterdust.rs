use core::fmt;
use std::cmp::max;

use crate::fasta_parsing::Fasta;

#[derive(Clone)]
pub struct LCR {
    pub name: String,
    pub start: usize, // base index, inclusive
    pub end: usize,   // base index, inclusive
    pub seq: String,
}

impl LCR {
    pub fn new(name: String, start: usize, end: usize, seq: String) -> Self {
        Self { name, start, end, seq }
    }
    pub fn get_name(&self) -> &str { &self.name }
    pub fn get_start(&self) -> usize { self.start }
    pub fn get_end(&self) -> usize { self.end }
}

impl fmt::Display for LCR {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\t{}\t{}\t{}", self.name, self.start, self.end, self.seq)
    }
}

/// Fast O(n) "slowdust": reports all *good* intervals.
/// - Uses rolling 2-bit k-mer encoding
/// - Uses per-run counts in a flat array (size 4^k)
/// - Adds Δ = ln(c+1) - T for each arriving k-mer
/// - Emits the last local max when the running score falls to ≤ 0
///
/// Coordinates: [start, end] inclusive (matches your existing code)
pub fn slowdust_fast(input: &Fasta, k: usize, t: f64, output: &mut Vec<LCR>) {
    let seq = input.get_sequence().as_bytes();
    if k == 0 || seq.len() < k { return; }
    assert!(k <= 31, "k up to 31 supported (2 bits/base in u64)");

    // name: first word of FASTA header
    let name = input.get_name().split_whitespace().next().unwrap_or_default().to_owned();

    // 2-bit base map (A=0, C=1, G=2, T=3). Anything else => None (breaks run).
    #[inline]
    fn base2(b: u8) -> Option<u8> {
        match b {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' => Some(3),
            _ => None,
        }
    }

    // Rolling k-mer mask
    let mask: u64 = if 2 * k == 64 { u64::MAX } else { (1u64 << (2 * k)) - 1 };

    // Per-run counts using a flat array and an epoch trick to avoid clearing.
    let table_size = 1usize << (2 * k);
    let mut counts: Vec<u32> = vec![0; table_size];
    let mut seen:   Vec<u32> = vec![0; table_size];
    let mut epoch: u32 = 1;
    #[inline]
    fn get_and_bump(idx: usize, counts: &mut [u32], seen: &mut [u32], epoch: u32) -> u32 {
        if seen[idx] != epoch {
            seen[idx] = epoch;
            counts[idx] = 0;
        }
        let c = counts[idx];
        counts[idx] = c + 1;
        c // return previous count
    }
    #[inline]
    fn reset_counts(seen: &mut [u32], epoch: &mut u32) {
        *epoch = epoch.wrapping_add(1);
        if *epoch == 0 {
            // super rare: wrap-around, fall back to full clear
            seen.fill(0);
            *epoch = 1;
        }
    }

    // Rolling state
    let mut code: u64 = 0;
    let mut valid = 0usize; // contiguous A/C/G/T so far

    // Run state
    let mut in_run = false;
    let mut run_start_kmer_base = 0usize; // base index of the run's first k-mer
    let mut run_score = 0.0f64;

    let mut best_score = f64::NEG_INFINITY;
    let mut best_end_kmer_base: Option<usize> = None;

    // Helper to emit the last local maximum (if positive) and reset run state
    let mut close_run = |current_kmer_end_base: usize, output: &mut Vec<LCR>| {
        if let Some(end_base) = best_end_kmer_base.take() {
            if best_score > 0.0 {
                let start_base = run_start_kmer_base;
                // inclusive range [start_base ..= end_base]
                let seq_slice = unsafe {
                    // SAFE: we only slice on byte boundaries we computed from the same byte string
                    std::str::from_utf8_unchecked(&seq[start_base ..= end_base])
                }.to_owned();

                output.push(LCR {
                    name: name.clone(),
                    start: start_base,
                    end: end_base,
                    seq: seq_slice,
                });
            }
        }
        in_run = false;
        run_score = 0.0;
        best_score = f64::NEG_INFINITY;
        best_end_kmer_base = None;
        reset_counts(&mut seen, &mut epoch);
        let _ = current_kmer_end_base; // docs clarity
    };

    for (i, &b) in seq.iter().enumerate() {
        match base2(b) {
            Some(v) => {
                code = ((code << 2) | (v as u64)) & mask;
                valid += 1;

                if valid >= k {
                    let kmer_end_base = i;              // inclusive
                    let kmer_start_base = i + 1 - k;    // inclusive

                    if !in_run {
                        in_run = true;
                        run_start_kmer_base = kmer_start_base;
                        run_score = 0.0;
                        best_score = f64::NEG_INFINITY;
                        best_end_kmer_base = None;
                    }

                    // Per-k-mer increment: Δ = ln(c+1) - T, where c is previous count in THIS run
                    let idx = (code as usize) & (table_size - 1);
                    let c_prev = get_and_bump(idx, &mut counts, &mut seen, epoch);
                    let delta = ((c_prev as f64) + 1.0).ln() - t;

                    run_score += delta;

                    // Track the last local max
                    if run_score > best_score {
                        best_score = run_score;
                        best_end_kmer_base = Some(kmer_end_base);
                    }

                    // When the run falls back to ≤ 0, emit the last max and reset
                    if run_score <= 0.0 {
                        close_run(kmer_end_base, output);
                    }
                }
            }
            None => {
                // Ambiguous base breaks any run; close if needed
                if in_run {
                    // The last valid k-mer end base is i-1 (if valid >= k)
                    let last_end = if valid >= k { i - 1 } else { 0 };
                    close_run(last_end, output);
                }
                // Reset rolling state
                code = 0;
                valid = 0;
            }
        }
    }

    // End-of-sequence: close any active run
    if in_run {
        let last_end = seq.len() - 1;
        close_run(last_end, output);
    }
}
