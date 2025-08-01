pub mod fasta_parsing;
pub mod slowdust;
pub mod output;

use std::{fs::File, io::BufReader};
use anyhow::{Ok, Result};

use crate::{fasta_parsing::{FastaIterator, BUFF_SIZE}, output::write_lcr, slowdust::{merge_intervals, slowdust, LCR}};

fn main() -> Result<()>{
    let file = File::open("small.fasta")?;
    let reader = BufReader::with_capacity(BUFF_SIZE, file);

    let iterator = FastaIterator::new(reader);
    let mut output = Vec::new();

    for line in iterator{
        let mut temp = Vec::new();
        let fasta = line?;
        let fasta_clone = fasta.clone();
        let seq = fasta_clone.get_sequence();
        slowdust(fasta, 30,1.5,  &mut temp);
        temp.sort_by(
            |lcr1, lcr2, | 
            lcr1.get_start().cmp(&lcr2.get_start())
            .then(lcr1.get_end().cmp(&lcr2.get_end())));
        output.append(&mut merge_intervals(temp, seq));
    }

    //let mut new: Vec<LCR> = merge_intervals(output, );

    let _ = write_lcr(&mut output, "output1.tsv");
    
    Ok(())
}
