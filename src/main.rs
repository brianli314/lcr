pub mod fasta_parsing;
pub mod output;
pub mod slowdust;
pub mod fasterdust;
pub mod testdust;

use anyhow::{Ok, Result};
use std::{
    collections::HashMap, fs::File, io::{BufReader, BufWriter, Write}, sync::{Arc, Mutex}, time::Instant
};
use threadpool::ThreadPool;

use crate::{
    fasta_parsing::{FastaIterator, BUFF_SIZE},
    //slowdust::{longdust_score, merge_intervals, slowdust, LCR},
    fasterdust::fasterdust, 
    slowdust::{longdust_score, merge_intervals, slowdust},
    testdust::{is_good_seq, testdust}
};

fn main() -> Result<()> {
    const NUM_THREADS: usize = 2;

    let pool = ThreadPool::new(NUM_THREADS);

    let file = File::open("data/scoring_test.fasta")?;
    let reader = BufReader::with_capacity(BUFF_SIZE, file);

    let output = File::create("test.tsv")?;
    let writer = Arc::new(Mutex::new(BufWriter::with_capacity(BUFF_SIZE, output)));

    {
        let mut header_guard = writer.lock().unwrap();
        let _ = writeln!(header_guard, "Name\tStart\tEnd\n");
        println!("Header written");
        header_guard.flush()?;
    }

    let iterator = FastaIterator::new(reader);

    println!("Starting loop");
    for line in iterator {
        let fasta = line?;
        //let seq = fasta.get_sequence();
        let writer_clone = Arc::clone(&writer);
        
        pool.execute(move || {
            let seq = fasta.get_sequence();
            let loop_now = Instant::now();
            let mut output = Vec::new();
            let name = fasta
                .get_name()
                .split_whitespace()
                .next()
                .unwrap_or_default();

            
             
            let score = longdust_score(seq, 7, 0.6);
            println!("{}", longdust_score(seq, 7, 0.6));
            println!("{}", is_good_seq(seq, score, 7, 0.6));

            let mut kmer_counts = HashMap::new();
            let mut prev_score = 0.0;
            println!("{}", seq.len());
            for i in 0..seq.len(){
                let win = &seq[i..];
                if win.len() >= 7{
                    let new_mer = &win[..7];

                    let entry = kmer_counts.entry(new_mer).or_insert(0.0);
                    *entry += 1.0;
                    let c_new: f64 = *entry;
                    prev_score = prev_score + (c_new).ln() - 0.6;
                }
            }
            println!("{}", prev_score);
            

            testdust(&fasta, 7, 5000, 0.6, &mut output);
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
