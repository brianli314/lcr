use std::collections::HashMap;
use std::fmt;

use crate::fasta_parsing::Fasta;
use crate::slowdust::LCR;

pub fn testdust(input: &Fasta, k: usize, max_window: usize, t: f64, output: &mut Vec<LCR>) {
    let seq = input.get_sequence();

    for end in k..=seq.len() {

        let mut prev_score = 0.0;

        let mut kmer_counts: HashMap<&str, f64> = HashMap::new();

        for win in k..=max_window {
            if win > end {
                break;
            }
            
            let start = end - win;
            let window = &seq[start..end];
            
            
            let new_mer = &window[..k];
            let entry = kmer_counts.entry(new_mer).or_insert(0.0);
            *entry += 1.0;

            let c_new = *entry;
            let window_score = prev_score + (c_new).ln() - t;

            prev_score = window_score;
            if window_score < t {
                continue;
            }
            

            if is_good_seq(window, window_score, k, t) {
                output.push(LCR {
                    name: input
                        .get_name()
                        .split_whitespace()
                        .next()
                        .unwrap_or_default()
                        .to_owned(),
                    start,
                    end,
                })
            }  
        }
    }
}

pub fn is_good_seq(window: &str, window_score: f64, k: usize, t: f64) -> bool {
    let mut prev_score_p = 0.0;
    let mut prev_score_s = 0.0;

    let mut kmer_counts_p: HashMap<&str, f64> = HashMap::new();
    let mut kmer_counts_s: HashMap<&str, f64> = HashMap::new();

    for i in 0..=window.len() {
        let prefix = &window[..i];

        if prefix.len() >= k{

            let new_mer_p = &prefix[prefix.len() - k..];

            let entry_p = kmer_counts_p.entry(new_mer_p).or_insert(0.0);
            *entry_p += 1.0;
            let c_new_p = *entry_p;
            let prefix_score = prev_score_p + (c_new_p).ln() - t;
            prev_score_p = prefix_score;

            if round_e12(prefix_score) > round_e12(window_score) {
                return false;
            }
        }
        let suffix = &window[window.len()-i..];

        if suffix.len() >= k{
            
            let new_mer_s = &suffix[..k];

            let entry_s = kmer_counts_s.entry(new_mer_s).or_insert(0.0);
            *entry_s += 1.0;
            let c_new_s = *entry_s;
            let suffix_score = prev_score_s + (c_new_s).ln() - t;
            prev_score_s = suffix_score;

            if round_e12(suffix_score) > round_e12(window_score) {
                return false;
            }
        }
    }
    true
}


fn round_e12(n:f64) -> f64{
    (n * 1e12).round() / 1e12
}
