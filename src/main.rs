pub mod fasta_parsing;
pub mod output;
pub mod slowdust;
pub mod fasterdust;

use anyhow::{Ok, Result};
use std::{
    fs::File,
    io::{BufReader, BufWriter, Write},
    sync::{Arc, Mutex},
    time::Instant,
};
use threadpool::ThreadPool;

use crate::{
    fasta_parsing::{FastaIterator, BUFF_SIZE},
    //slowdust::{longdust_score, merge_intervals, slowdust, LCR},
    fasterdust::{fasterdust, union_good_intervals}
};

fn main() -> Result<()> {
    let num_threads = 2;
    const MAX_IN_FLIGHT: usize = 4 * 2;
    
    let pool = ThreadPool::new(num_threads);

    let file = File::open("test.fasta")?;
    let reader = BufReader::with_capacity(BUFF_SIZE, file);

    let output = File::create("result_test.tsv")?;
    let writer = Arc::new(Mutex::new(BufWriter::with_capacity(BUFF_SIZE, output)));

    {
        let mut header_guard = writer.lock().unwrap();
        let _ = writeln!(header_guard, "Name\tStart\tEnd\tString\n");
        println!("Header written");
        header_guard.flush()?;
    }

    let iterator = FastaIterator::new(reader);

    println!("Starting loop");
    for line in iterator {
        let fasta = line?;
        let writer_clone = Arc::clone(&writer);
        
        while pool.queued_count() >= MAX_IN_FLIGHT {
            std::thread::sleep(std::time::Duration::from_millis(2));
        }
        
        pool.execute(move || {
            let loop_now = Instant::now();
            let mut output = Vec::new();
            let name = fasta
                .get_name()
                .split_whitespace()
                .next()
                .unwrap_or_default();
            /*
            let fasta_clone = fasta.clone();
            let seq = fasta_clone.get_sequence();
            
            slowdust(&fasta, 5000, 0.6, &mut temp);
            temp.sort_by(|lcr1, lcr2| {
                lcr1.get_start()
                    .cmp(&lcr2.get_start())
                    .then(lcr1.get_end().cmp(&lcr2.get_end()))
            });

            let merged = merge_intervals(temp, seq);

            let mut guard = writer_clone.lock().unwrap_or_else(|e| e.into_inner());
            for lcr in merged {
                let _ = writeln!(guard, "{}", lcr);
            }
            guard.flush().expect("Failed to flush writer");
            let name = fasta
                .get_name()
                .split_whitespace()
                .next()
                .unwrap_or_default();

            println!("{}", longdust_score(seq, 0.6));
             */
            fasterdust(&fasta, 7, 5000, 0.6, &mut output);
            let merged = union_good_intervals(output, true);
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
