use rustc_hash::FxHashMap;
use crate::fasta_parsing::Fasta;
use crate::slowdust::LCR;

/// Encode A/C/G/T → 2-bit; anything else -> None
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

/// Precompute k-mer codes at each start position (length = seq.len()).
/// Positions where a full k-mer doesn't exist or includes a non-ACGT get None.
fn precompute_kmers(seq: &[u8], k: usize) -> Vec<Option<u64>> {
    let mut codes = vec![None; seq.len()];
    if k == 0 || seq.len() < k { return codes; }

    let mask: u64 = if 2*k == 64 { u64::MAX } else { (1u64 << (2*k)) - 1 };
    let mut code: u64 = 0;
    let mut valid = 0usize;

    for (i, &b) in seq.iter().enumerate() {
        match base2(b) {
            Some(v) => {
                code = ((code << 2) | (v as u64)) & mask;
                valid += 1;
                if valid >= k {
                    let start = i + 1 - k;
                    codes[start] = Some(code);
                }
            }
            None => {
                code = 0;
                valid = 0;
            }
        }
    }
    codes
}

/// EXACT algorithm per your outline (1–6).
/// k: k-mer length (set 7 to match your scoring)
/// t: threshold T (used in S_L as the per-k-mer penalty, AND as the minimum score filter)
pub fn fasterdust(
    input: &Fasta,
    k: usize,
    max_window: usize,
    t: f64,
    output: &mut Vec<LCR>,
) {
    let seq_str = input.get_sequence();
    let seq = seq_str.as_bytes();
    if seq.len() < k { return; }

    let name = input
        .get_name()
        .split_whitespace()
        .next()
        .unwrap_or_default()
        .to_owned();

    // Precompute k-mer code at each start
    let kmers = precompute_kmers(seq, k);

    // Precompute ln(n) for increments Δ = ln(c_prev+1) - t
    let max_kmers_per_window = max_window.saturating_sub(k).saturating_add(1);
    let mut ln_table = vec![0.0f64; max_kmers_per_window + 2]; // index by (c_prev+1)
    for n in 1..ln_table.len() {
        ln_table[n] = (n as f64).ln();
    }

    for end in 0..seq.len() {
        // We can only form windows with at least one k-mer if end+1 >= k
        if end + 1 < k { continue; }

        // Largest window allowed by max_window (in bases)
        let min_start_base = end.saturating_add(1).saturating_sub(max_window);

        // For this end, we expand windows leftward, adding exactly one new k-mer each step
        let first_kmer_start = end + 1 - k; // start of the rightmost k-mer in the window

        let mut win_counts: FxHashMap<u64, u32> = FxHashMap::default();
        let mut s_total = 0.0f64;  // S_L(window)
        let mut ell = 0usize;      // # of k-mers in window

        // Start from the smallest window with >=1 k-mer, and grow leftward
        // Each iteration adds the k-mer starting at `start`
        let mut start = first_kmer_start as isize;
        let stop = min_start_base as isize;

        while start >= stop {
            let s = start as usize;

            // If the k-mer at `s` is invalid (contains non-ACGT), we stop expanding this end.
            let code = match kmers.get(s).and_then(|&c| c) {
                Some(code) => code,
                None => break,
            };

            // Update window score incrementally: Δ = ln(c_prev+1) - t
            let entry = win_counts.entry(code).or_insert(0);
            let c_prev = *entry as usize;
            *entry += 1;
            ell += 1;
            if c_prev + 1 < ln_table.len() {
                s_total += ln_table[c_prev + 1] - t;
            } else {
                // Shouldn't happen with sane max_window; fallback:
                s_total += ((c_prev + 1) as f64).ln() - t;
            }

            // Only evaluate "good" if the total score passes your minimum filter
            if s_total >= t
                && is_good_window(&kmers, s, end, k, t, &ln_table, s_total) {
                    // Push [start, end] inclusive, with exact sequence slice
                    output.push(LCR {
                        name: name.clone(),
                        start: s,
                        end,
                    });
                }

            start -= 1;
        }
    }
}

/// Check "good": no proper prefix or proper suffix has higher score than S(window).
/// We recompute prefix/suffix scores **exactly** over the k-mers of this window.
/// Early-out as soon as we detect a violation.
fn is_good_window(
    kmers: &[Option<u64>],
    start_base: usize,
    end_base: usize,
    k: usize,
    t: f64,
    ln_table: &[f64],
    s_total: f64,
) -> bool {
    // The k-mers contained in [start_base ..= end_base] start at:
    //    start_k ..= last_k   where last_k = end_base + 1 - k
    let last_k = end_base + 1 - k;
    let start_k = start_base;
    debug_assert!(last_k >= start_k);

    // ---- Proper prefixes: start_k .. last_k-1 ----
    // Accumulate from left to right, stop if any prefix score exceeds s_total
    {
        let mut counts: FxHashMap<u64, u32> = FxHashMap::default();
        let mut s = 0.0f64;

        for pos in start_k..last_k { // excludes the last k-mer => proper prefix
            let code = match kmers[pos] {
                Some(c) => c,
                None => return false, // shouldn't happen if outer loop screened, but be safe
            };
            let entry = counts.entry(code).or_insert(0);
            let c_prev = *entry as usize;
            *entry += 1;
            s += ln_table[c_prev + 1] - t;

            if s > s_total {
                return false; // a proper prefix beats the window
            }
        }
    }

    // ---- Proper suffixes: start_k+1 .. last_k (right to left) ----
    {
        let mut counts: FxHashMap<u64, u32> = FxHashMap::default();
        let mut s = 0.0f64;

        // build from the rightmost k-mer backwards, but only up to a proper suffix
        // i.e., include positions last_k, last_k-1, ..., start_k+1
        for pos in (start_k + 1..=last_k).rev() {
            let code = match kmers[pos] {
                Some(c) => c,
                None => return false,
            };
            let entry = counts.entry(code).or_insert(0);
            let c_prev = *entry as usize;
            *entry += 1;
            s += ln_table[c_prev + 1] - t;

            if s > s_total {
                return false; // a proper suffix beats the window
            }
        }
    }

    true
}

