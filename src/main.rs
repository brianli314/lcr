pub mod fasta_parsing;
pub mod slowdust;
pub mod output;

use core::num;
use std::{fs::File, io::{BufReader, BufWriter, Write}, sync::{Arc, Mutex}, time::Instant};
use anyhow::{Ok, Result};
use threadpool::ThreadPool;

use crate::{fasta_parsing::{FastaIterator, BUFF_SIZE}, slowdust::{merge_intervals, slowdust, LCR}};

fn main() -> Result<()>{

    let num_threads = 3;
    let pool = ThreadPool::new(num_threads);

    let file = File::open("test copy.fasta")?;
    let reader = BufReader::with_capacity(BUFF_SIZE, file);


    let output = File::create("result.tsv")?;
    let writer = Arc::new(Mutex::new(BufWriter::with_capacity(BUFF_SIZE, output)));

    let mut writer_guard = writer.lock().unwrap();
    let _ = writeln!(
        writer_guard,
        "Name\tStart\tEnd\tString\n"
    );
    println!("Header written");

    let iterator = FastaIterator::new(reader);
    //let mut output = Vec::new();

    println!("Starting loop");
    for line in iterator {
        let writer = Arc::clone(&writer); 
        let fasta = line?;
        
        pool.execute(move || {
            let loop_now = Instant::now();
            let mut temp = Vec::new();
            let fasta_clone = fasta.clone();
            let seq = fasta_clone.get_sequence();
            slowdust(&fasta, 5000,0.6,  &mut temp);
            
            temp.sort_by(
                |lcr1, lcr2, | 
                lcr1.get_start().cmp(&lcr2.get_start())
                .then(lcr1.get_end().cmp(&lcr2.get_end())));

            let merged = merge_intervals(temp, seq);
            let mut writer_guard = writer.lock().unwrap();
            for lcr in merged {
                let _ = writeln!(writer_guard, "{}", lcr);
            }
            let loop_elapsed = loop_now.elapsed();
            println!("1 Loop finished in: {:.2?}", loop_elapsed);
        });
        
        //output.append(&mut merge_intervals(temp, seq));
    };
    pool.join();
    println!("Finished running");

    //let mut new: Vec<LCR> = merge_intervals(output, );

    //let _ = write_lcr(&mut output, "result.tsv");
    
    Ok(())
}
