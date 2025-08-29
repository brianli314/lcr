pub mod fasta_parsing;
pub mod slowdust;
pub mod slowdust2;
pub mod command_line;

use anyhow::{Ok, Result};
use clap::Parser;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufReader, BufWriter, Write},
    sync::{Arc, Mutex},
    time::Instant,
};
use threadpool::ThreadPool;

use crate::{
    command_line::DustArgs, fasta_parsing::{FastaIterator, BUFF_SIZE}, slowdust::{longdust_score, merge_intervals, slowdust}, slowdust2::{is_good_seq, slowdust2}
};

fn main() -> Result<()> {
    
    //print_score("ctcctctcctttcttctctccatccCCCCTCCATCCCcgtctcctttctcctctccatccccctctccatccccctctccatctccctctcctttctcctctccatccccctctcctttctccctctccatccccctctCCTTTCTTC",7, 0.6);
    //return Ok(());

    let args = DustArgs::parse();
    
    let num_threads: usize = args.threads;

    let pool = ThreadPool::new(num_threads);

    let file = File::open(args.input_file)?;
    let reader = BufReader::with_capacity(BUFF_SIZE, file);

    let output = File::create(args.output_file)?;
    let writer = Arc::new(Mutex::new(BufWriter::with_capacity(BUFF_SIZE, output)));

    {
        let mut header_guard = writer.lock().unwrap();
        let _ = writeln!(header_guard, "Name\tStart\tEnd\n");
        header_guard.flush()?;
    }

    let iterator = FastaIterator::new(reader);

    for line in iterator {
        let fasta = line?;
        let writer_clone = Arc::clone(&writer);

        pool.execute(move || {
            //let seq = fasta.get_sequence();
            let loop_now = Instant::now();
            let mut output = Vec::new();
            let name = fasta
                .get_name()
                .split_whitespace()
                .next()
                .unwrap_or_default();

            slowdust2(&fasta, 7, 5000, 0.6, &mut output);
            let merged = merge_intervals(output);
            let loop_elapsed = loop_now.elapsed();
            println!("1 Loop finished in {loop_elapsed:.2?} for {name}");
            let mut guard = writer_clone.lock().unwrap_or_else(|e| e.into_inner());
            for lcr in merged {
                let _ = writeln!(guard, "{}", lcr);
            }
            guard.flush().expect("Failed to flush writer");
        });
    }
    pool.join();
    println!("Finished running");

    Ok(())
}


//For manual sequence checking
fn print_score(seq: &str, k: usize, t: f64) {
    let score = longdust_score(seq, k, t);

    let mut kmer_counts = HashMap::new();
    let mut score_2 = 0.0;
    for i in k..=seq.len() {
        let window = &seq[..i];

        let new_mer = &window[window.len() - k..];
        let entry = kmer_counts.entry(new_mer).or_insert(0.0);
        *entry += 1.0;

        let c_new: f64 = *entry;
        let window_score = score_2 + (c_new).ln() - t;

        score_2 = window_score;
    }

    if score < t{
        println!("slowdust-below threshold: {score}")
    } else {
        println!("slowdust: {score}")
    }

    if score_2 < t {
        println!("slowdust2-below threshold: {score_2}")
    } else {
        println!("slowdust2: {score_2}")
    }

    let mut is_good = true;
    for j in 0..=seq.len() {
        let prefix = &seq[..j];
        let suffix = &seq[j..];
        if longdust_score(prefix, k, t) > score || longdust_score(suffix, k, t) > score {
            is_good = false;
            break;
        }
    }

    if is_good{
        println!("slowdust: Good")
    } else {
        println!("slowdust: Not good")
    }

    if is_good_seq(seq, score_2, k, t) {
        println!("slowdust2: Good")
    } else {
        println!("slowdust2: Not good")
    }


    if is_good_seq(seq, score, k, t) {
        println!("slowdust2 with slowdust score: Good seq")
    } else {
        println!("slowdust2 with slowdust score: Not good seq")
    }
}
