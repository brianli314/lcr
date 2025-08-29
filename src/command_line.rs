use clap::Parser;

#[derive(Debug, Parser)]
pub struct DustArgs{
    #[arg(short, long = "input", required = true)]
    ///Input file path
    pub input_file: String,

    #[arg(short, long = "output")]
    ///Output file path.
    pub output_file: String,

    ///The file path for the list of adapter seqeuences. Must be fasta format
    #[arg(short, long, default_value_t = 1)]
    pub threads: usize,
}