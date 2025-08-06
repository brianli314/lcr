pub mod fasta_parsing;
pub mod slowdust;
pub mod output;

use std::{fs::File, io::{BufReader, BufWriter, Write}, time::Instant};
use anyhow::{Ok, Result};

use crate::{fasta_parsing::{FastaIterator, BUFF_SIZE}, slowdust::{merge_intervals, slowdust, LCR}};

fn main() -> Result<()>{
    let file = File::open("test copy.fasta")?;
    let reader = BufReader::with_capacity(BUFF_SIZE, file);


    let output = File::create("result.tsv")?;
    let mut writer = BufWriter::with_capacity(BUFF_SIZE, output);

    let _ = writeln!(
        writer,
        "Name\tStart\tEnd\tString\n"
    );
    println!("Header written");

    let iterator = FastaIterator::new(reader);
    //let mut output = Vec::new();

    println!("Starting loop");
    for line in iterator{
        let loop_now = Instant::now();
        let mut temp = Vec::new();
        let fasta = line?;
        let fasta_clone = fasta.clone();
        let seq = fasta_clone.get_sequence();
        slowdust(fasta, 5000,0.6,  &mut temp);
        
        temp.sort_by(
            |lcr1, lcr2, | 
            lcr1.get_start().cmp(&lcr2.get_start())
            .then(lcr1.get_end().cmp(&lcr2.get_end())));

        for lcr in merge_intervals(temp, seq){
            let _ = writeln!(writer, "{}", lcr);
        }
        let loop_elapsed = loop_now.elapsed();
        println!("1 Loop finished in: {:.2?}", loop_elapsed);
        //output.append(&mut merge_intervals(temp, seq));
    }

    //let mut new: Vec<LCR> = merge_intervals(output, );

    //let _ = write_lcr(&mut output, "result.tsv");
    
    Ok(())
}
