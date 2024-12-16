use clap::Parser;

#[derive(Debug, Parser)]
#[clap(version)]

pub struct OverlapArgs {
    /// please provide the path to the fastq file
    pub overlap_file: String,
    /// please provide the overlap length
    pub overlap_size: usize
}
