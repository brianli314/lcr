use core::fmt;
use std::cmp::max;
use rustc_hash::FxHashMap;

use crate::fasta_parsing::Fasta;

#[derive(Clone)]
pub struct LCR {
    pub name: String,
    pub start: usize, // inclusive
    pub end: usize,   // inclusive
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

/// Same semantics as your original `slowdust`, but much faster:
/// - still loops (i, w) over all end positions and window sizes
/// - still uses `threshold` as both the per-k-mer penalty T and the min score filter
/// - still treats k-mers as exact string slices (so 'N', lowercase, etc. are allowed)
pub fn fasterdust(input: &Fasta, max_window: usize, threshold: f64, output: &mut Vec<LCR>) {
    let seq_str = input.get_sequence();
    let seq = seq_str.as_bytes();
    if seq.is_empty() { return; }

    let name = input.get_name().split_whitespace().next().unwrap_or_default().to_owned();

    // choose k to match your scoring
    let k: usize = 7;

    for i in 0..seq.len() {
        if i % 1000 == 0 {
            println!("Running on {i}th base pair for {name}");
        }

        // Per-endpoint state: we grow windows leftwards from end = i
        let mut kmer_counts: FxHashMap<&[u8], u32> = FxHashMap::default();

        // Running totals for the *current* window (which grows as w increases)
        let mut s_total: f64 = 0.0;              // S(window)
        let mut ell: usize = 0;                  // number of k-mers in window
        let mut max_pref_proper: f64 = f64::NEG_INFINITY; // max score over proper prefixes
        let mut max_suf_proper: f64 = f64::NEG_INFINITY;  // max score over proper suffixes

        // grow window size w = 1..=max, but stop at left boundary
        let max_w = max_window.min(i + 1);
        for w in 1..=max_w {
            let start = i + 1 - w;

            // If adding this base creates a *new* k-mer at the left edge, update counts/score.
            if w >= k {
                let kmer_start = start;
                let kmer_end = start + k; // exclusive
                let kmer = &seq[kmer_start..kmer_end];

                let entry = kmer_counts.entry(kmer).or_insert(0);
                let c_prev = *entry as f64;
                *entry += 1;
                ell += 1;

                // Incremental score update: Δ = ln(c_prev + 1) - threshold
                let delta = (c_prev + 1.0).ln() - threshold;

                // Before updating s_total, update max suffix based on the *previous* total.
                // The entire previous window becomes a proper suffix of the new window.
                if ell > 1 {
                    // only meaningful when there was a previous window with >=1 k-mer
                    if s_total > max_suf_proper {
                        max_suf_proper = s_total;
                    }
                }

                // Update running total
                s_total += delta;

                // Update max proper prefix: either just this delta, or delta + previous max prefix
                if ell > 1 {
                    let cand1 = delta;
                    let cand2 = delta + max_pref_proper;
                    max_pref_proper = if cand1 > cand2 { cand1 } else { cand2 };
                } else {
                    // First k-mer in the window => no proper prefix with any k-mers
                    max_pref_proper = f64::NEG_INFINITY;
                    max_suf_proper  = f64::NEG_INFINITY;
                }

                // Evaluate "good" only if it passes the min-score filter, like your code
                if s_total >= threshold {
                    // Good iff no proper prefix/suffix scores higher than the whole window
                    let good_prefix = max_pref_proper <= s_total;
                    let good_suffix = max_suf_proper  <= s_total;

                    if good_prefix && good_suffix {
                        // Push interval [start, i] inclusive, with exact sequence slice
                        let subseq = unsafe {
                            std::str::from_utf8_unchecked(&seq[start..=i])
                        }.to_owned();

                        output.push(LCR {
                            name: name.clone(),
                            start,
                            end: i,
                            seq: subseq,
                        });
                    }
                }
            } else {
                // w < k: the score remains 0 (no k-mers yet). Your original code would
                // compute 0 and then immediately filter it out if `threshold > 0`.
                // We skip "good" evaluation here for speed; behavior is identical when threshold > 0.
            }
        }
    }
}
