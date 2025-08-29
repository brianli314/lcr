use core::fmt;
use rustc_hash::FxHashMap;
use statrs::function::factorial::ln_factorial;

use crate::fasta_parsing::Fasta;

#[derive(Clone)]
pub struct LCR {
    pub name: String,
    pub start: usize,
    pub end: usize,
}

impl LCR {
    pub fn new(name: String, start: usize, end: usize) -> Self {
        Self { name, start, end }
    }

    pub fn get_name(&self) -> &str {
        &self.name
    }
    pub fn get_start(&self) -> usize {
        self.start
    }
    pub fn get_end(&self) -> usize {
        self.end
    }
}

impl fmt::Display for LCR {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\t{}\t{}", self.name, self.start, self.end)
    }
}

pub fn slowdust(input: &Fasta, k: usize, max_window: usize, threshold: f64, output: &mut Vec<LCR>) {
    let seq = input.get_sequence();
    let name = input
        .get_name()
        .split_whitespace()
        .next()
        .unwrap_or_default();
    for i in k..=seq.len() {
        for w in k..=max_window {
            if w > i {
                break;
            }

            let window = &seq[i - w..i];

            //let window_time = Instant::now();
            let window_score = longdust_score(window, k, threshold);
            //println!("Computed window in: {:.2?}", window_time.elapsed());
            if window_score < threshold {
                continue;
            }

            let mut is_good = true;
            for j in 0..=window.len() {
                let prefix = &window[..j];
                let suffix = &window[j..];
                if longdust_score(prefix, k, threshold) > window_score
                    || longdust_score(suffix, k, threshold) > window_score
                {
                    is_good = false;
                    break;
                }
            }
            if is_good {
                output.push(LCR {
                    name: name.to_owned(),
                    start: i - w,
                    end: i,
                });
            }
        }
    }
}

pub fn longdust_score(x: &str, k: usize, threshold: f64) -> f64 {
    if x.len() < k {
        return 0.0;
    }
    let counts = count_kmers(x, k);
    let output: f64 = counts.values().map(|&c| ln_factorial(c as u64)).sum();

    output - threshold * ((x.len() - k + 1) as f64)
}

fn count_kmers(x: &str, k: usize) -> FxHashMap<&str, i32> {
    let mut counts = FxHashMap::default();
    for i in 0..x.len() - k + 1 {
        let kmer = &x[i..i + k];
        *counts.entry(kmer).or_insert(0) += 1;
    }
    counts
}

pub fn merge_intervals(mut intervals: Vec<LCR>) -> Vec<LCR> {
    if intervals.is_empty() {
        return vec![];
    }

    intervals.sort_by(|a, b| {
        a.name.cmp(&b.name)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });

    let mut merged = Vec::new();
    let mut current = intervals[0].clone();

    for next in intervals.into_iter().skip(1) {
        if next.name == current.name && next.start <= current.end {
            // same name and overlapping → merge
            current.end = current.end.max(next.end);
        } else {
            // push previous and move on
            merged.push(current);
            current = next;
        }
    }

    merged.push(current);
    merged
}

/*
pub fn slowdust2(input: Fasta, window_len: usize, threshold: f64, output: &mut Vec<LCR>){
    let seq = input.get_sequence();
    for w in 1..window_len{
        for i in 0..seq.len(){
            if i + w >= seq.len(){
                break
            }

            let window = &seq[i..i+w];
            let total_score = sdust_score(window);

            if total_score < threshold{
                continue
            }
            let is_good= is_good_seq(window, total_score, i, i+w);
            let mut start = i;
            let mut end = i+w-1;
            if is_good.0{
                output.push(LCR::new(input.get_name().to_owned(), start, end));
            } else {
                start = is_good.1;
                end = is_good.2;
                output.push(LCR::new(input.get_name().to_owned(), start, end));

            }
        }
    }
}


fn is_good_seq(window: &str, total_score: f64, window_start: usize, window_end: usize) -> (bool, usize, usize){
    let mut is_good = true;
    let mut start = window_start;
    let mut end = window_end;
    let mut best_score = total_score;
    for j in 1..window.len(){
        let prefix = &window[..j];

        let prefix_score = sdust_score(prefix);

        if prefix_score > total_score {
            is_good = false;
            end = window_start + j;
            if prefix_score > best_score{
                end = window_start + j;
                best_score = prefix_score
            }
        } else {
            let suffix = &window[j..];
            let suffix_score  = sdust_score(suffix);
            if suffix_score > total_score {
                is_good = false;
                start = window_start + j;
                if suffix_score > best_score {
                    start = window_start + j;
                    best_score = suffix_score;
                }
            }
        }
    }
    (is_good, start, end)
}

fn compute_kmer_counts(s: &str, k: usize) -> HashMap<&str, usize> {
    let mut map = HashMap::new();
    for i in 0..(s.len().saturating_sub(k)) {
        //println!("i: {}, string: {}", i, s);
        let kmer = &s[i..i + k];
        *map.entry(kmer).or_insert(0) += 1;
    }
    map
}

fn score(s: &str, f_cache: Option<&mut HashMap<u64, f64>>) -> f64 {
    let l = s.len().saturating_sub(2) + 1;
    if l == 0 { return f64::NEG_INFINITY; }

    let kmer_counts = compute_kmer_counts(s, 2);
    let log_fact_sum: f64 = kmer_counts.values().map(|&c| ln_gamma(c as f64+1.0)).sum();

    log_fact_sum
}
*/
