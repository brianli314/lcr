use std::fmt;

use crate::fasta_parsing::Fasta;

#[derive(Clone)]
pub struct LCR {
    pub name: String,
    pub start: usize, // inclusive
    pub end: usize,   // inclusive
}

impl LCR {
    pub fn new(name: String, start: usize, end: usize) -> Self {
        Self { name, start, end }
    }
    pub fn get_name(&self) -> &str { &self.name }
    pub fn get_start(&self) -> usize { self.start }
    pub fn get_end(&self) -> usize { self.end }
    pub fn add_one(self) -> Self {
        LCR { name: self.name, start: self.start, end: self.end + 1 }
    }
}

impl fmt::Display for LCR {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}\t{}\t{}", self.name, self.start, self.end + 1)
    }
}

pub fn testdust(input: &Fasta, k: usize, max_window: usize, t: f64, output: &mut Vec<LCR>){
    let seq = input.get_sequence();
    let name = input
        .get_name()
        .split_whitespace()
        .next()
        .unwrap_or_default()
        .to_owned();


    for end in k..seq.len(){
        let window_score = 0.0;

        for win in k..=max_window{
            if k > end {continue}
            
            let start = end - k;
            let window = &seq[start..=end];

            

        }
    }

}