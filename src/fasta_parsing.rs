use anyhow::{anyhow, Ok, Result};
use std::{
    io::{BufRead, BufReader, Lines, Read},
    mem,
};

pub const BUFF_SIZE: usize = 1 << 20;

#[derive(Debug, Clone)]
pub struct Fasta {
    pub name: String,
    pub sequence: String,
}

#[derive(Debug, Clone)]
pub struct Fastq {
    pub id: String,
    pub sequence: String,
    pub quality: String,
}

impl Fastq{
    pub fn new(id: String, sequence: String, quality: String) -> Self{
        Self{id, sequence, quality}
    }
    pub fn get_sequence(&self) -> &str {
        &self.sequence
    }
    pub fn get_id(&self) -> &str {
        &self.id
    }
    pub fn get_quality(&self) -> &str {
        &self.quality
    }
    pub fn to_fasta(self) -> Fasta{
        Fasta { name: self.id, sequence: self.sequence}
    }
}

impl Fasta {
    pub fn new(name: String, sequence: String) -> Self {
        Self { name, sequence }
    }
    pub fn get_sequence(&self) -> &str {
        &self.sequence
    }
    pub fn get_name(&self) -> &str {
        &self.name
    }
    pub fn set_sequence(&mut self, seq: String){
        self.sequence = seq;
    }
}

#[derive(Debug)]
pub struct FastaIterator<T: Read> {
    lines_reader: Lines<BufReader<T>>,
    curr_name: String,
}

pub struct FastqIterator<T:Read> {
    lines_reader: Lines<BufReader<T>>,
    curr_name: String,
    counter: i32,
}

impl<T: Read> FastqIterator<T> {
    pub fn new(bufreader: BufReader<T>) -> Self {
        Self {
            lines_reader: bufreader.lines(),
            curr_name: String::new(),
            counter: 0,
        }
    }
}

impl<T: Read> Iterator for FastqIterator<T>{
    type Item = Result<Fasta>;

    fn next(&mut self) -> Option<Self::Item>{
        let mut seq: String = String::new();

        for line in self.lines_reader.by_ref(){
            let line = match line{
                Result::Ok(line) => line,
                Err(err) => return Some(Err(anyhow!("Invalid line/file format: {err}")))
            };

            if self.counter == 0 && line.starts_with('@'){
                self.counter += 1;
                let mut name: String = line.strip_prefix('@').unwrap().to_owned();
                if seq.is_empty(){
                    name.clone_into(&mut self.curr_name);
                    continue;
                }
                
                mem::swap(&mut name, &mut self.curr_name);
                return Some(Ok(Fasta {
                    name, 
                    sequence: seq
                }))
            } else if self.counter == 1{
                seq = line;
                self.counter += 1;
                continue
            } else if self.counter == 2{
                self.counter += 1;
                continue
            } else if self.counter == 3{
                self.counter = 0;
                continue
            } else {
                return Some(Err(anyhow!("Invalid fastq format")));
            }
        }
        if seq.is_empty() {
            None
        } else {
            Some(Ok(Fasta {
                name: mem::take(&mut self.curr_name),
                sequence: seq,
            }))
        }
    }
}

impl<T: Read> Iterator for FastaIterator<T> {
    type Item = Result<Fasta>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut seq: String = String::new();

        for line in self.lines_reader.by_ref() {
            let line = match line {
                Result::Ok(line) => line,
                Err(err) => return Some(Err(anyhow!("Invalid line/file format: {err}"))),
            };

            if line.starts_with('>') {
                let mut name = line.strip_prefix('>').unwrap().to_owned();
                if seq.is_empty() {
                    name.clone_into(&mut self.curr_name);
                    continue;
                }
                mem::swap(&mut name, &mut self.curr_name);
                return Some(Ok(Fasta {
                    name,
                    sequence: seq,
                }));
            //Checks for valid starting line in fasta
            } else if !self.curr_name.is_empty() {
                seq += &line;
            } else {
                return Some(Err(anyhow!("Invalid fasta format")));
            }
        }
        if seq.is_empty() {
            None
        } else {
            Some(Ok(Fasta {
                name: mem::take(&mut self.curr_name),
                sequence: seq,
            }))
        }
    }
}

impl<T: Read> FastaIterator<T> {
    pub fn new(bufreader: BufReader<T>) -> Self {
        Self {
            lines_reader: bufreader.lines(),
            curr_name: String::new()
        }
    }
}