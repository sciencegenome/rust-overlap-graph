mod args;
use args::OverlapArgs;
use clap::Parser;
use std::collections::BTreeMap;
use std::collections::HashSet;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/*
*Author Gaurav Sablok
*Universitat Potsdam
*Date 2024-12-16

rust-overlap-graphs: implementing the rust overlap graphs in rust.
It first computes the overlap based on the defined sequence lap and
also makes a offset table for the defined graph edges across all the sequences.
If the indexed table lookup graph edge is present in the node then it will
capture the offset for the same so that it can easily build a suffix array only
from those reads which can have that edge present.
* */

fn main() {
    let args = OverlapArgs::parse();
    let overlap_output = overlap_graph(&args.overlap_file, args.overlap_size).unwrap();
    println!("The bwt has been written: {:?}", overlap_output);
}

pub fn overlap_graph(path: &str, overlapsize: usize) -> Result<Vec<String>, Box<dyn Error>> {
    /* This function is searching for the overlap between the
     * adjacent pairs for the sequences stored on a BTreeMap.
     * */

    #[derive(Debug, Clone, PartialOrd, PartialEq)]
    struct ReadSeq {
        header: String,
        sequence: String,
    }

    let open_readfile = File::open(path).expect("file not present");
    let readfile = BufReader::new(open_readfile);
    let mut capture_seq: Vec<ReadSeq> = Vec::new();
    let mut header: Vec<String> = Vec::new();
    let mut seq: Vec<String> = Vec::new();
    for i in readfile.lines() {
        let line = i.expect("line not present");
        if line.starts_with("@") {
            header.push(line.clone())
        }
        if line.starts_with("A")
            || line.starts_with("T")
            || line.starts_with("G")
            || line.starts_with("C")
        {
            seq.push(line)
        }
    }

    for i in 0..header.len() {
        capture_seq.push(ReadSeq {
            header: header[i].clone(),
            sequence: seq[i].clone(),
        })
    }

    // storing the overlap and the corresponding sequence on a BTreeMap
    // so that it can used for the Binary search Tree traversal. In this
    // i made the overlap as the key in the BTreeMap, so that if the overlap
    // is there then traverse through the BTreeMap of that particular overlap.

    #[derive(Debug, Clone, PartialEq, PartialOrd)]
    struct Readtraversal {
        read1: String,
        read2: String,
        first_read_pre_suffix: String,
        second_read_pre_prefix: String,
    }

    let mut overlap_graph: BTreeMap<String, Readtraversal> = BTreeMap::new();

    for i in 0..seq.len() {
        if seq[i][seq[i].len() - overlapsize..seq[i].len()].to_string()
            == seq[i + 1][0..overlapsize].to_string()
        {
            overlap_graph.insert(
                seq[i][seq[i].len() - overlapsize..seq[i].len()].to_string(),
                Readtraversal {
                    read1: seq[i].clone(),
                    read2: seq[i + 1].clone(),
                    first_read_pre_suffix: seq[i][0..seq[i].len() - 3].to_string(),
                    second_read_pre_prefix: seq[i + 1][overlapsize..seq[i + 1].len()].to_string(),
                },
            );
        };
    }

    let mut collect_btreemap_keys: Vec<String> = Vec::new();
    for (key, _) in overlap_graph.iter() {
        collect_btreemap_keys.push(key.to_string());
    }

    // writing the overlap-graph

    let mut overlap_graph_write = File::create("overlap.txt").expect("file not present");

    for (i, j) in overlap_graph.iter() {
        write!(
            overlap_graph_write,
            "{}\t{}\t{}\t{}\t{}",
            i, j.read1, j.read2, j.first_read_pre_suffix, j.second_read_pre_prefix
        ).expect("line not present");
    }

    Ok(collect_btreemap_keys)
}

pub fn compute_frequency_table(path: &str, kmer: usize) -> HashSet<String> {
    /* computing the unique kmers from all the reads and then making a
     * frequency table. This is used later to create a graphlookup table.
     * in which each of the unique kmer obtained is compared with those
     * prefix and sufix kmer to see the location of those in the reads
     * where those kmers where not terminal and break those reads at those
     * points to create the new nodes for the suffix array.
     *
     * */

    #[derive(Debug, Clone, PartialOrd, PartialEq)]
    struct ReadSeq {
        header: String,
        sequence: String,
    }

    let open_readfile = File::open(path).expect("file not present");
    let readfile = BufReader::new(open_readfile);
    let mut header: Vec<String> = Vec::new();
    let mut seq: Vec<String> = Vec::new();
    for i in readfile.lines() {
        let line = i.expect("line not present");
        if line.starts_with("@") {
            header.push(line.clone())
        }
        if line.starts_with("A")
            || line.starts_with("T")
            || line.starts_with("G")
            || line.starts_with("C")
        {
            seq.push(line)
        }
    }

    let mut compute_frequency: Vec<Vec<String>> = Vec::new();

    for i in 0..seq.len() {
        let string_char: Vec<_> = seq[i].chars().map(|x| String::from(x)).collect::<Vec<_>>();
        let compute_freq: Vec<String> = string_char
            .windows(kmer)
            .map(|x| x.concat())
            .collect::<Vec<String>>();
        compute_frequency.push(compute_freq);
    }

    let mut uniquehash_table: HashSet<String> = HashSet::new();
    for i in compute_frequency.iter() {
        for iterable in i.iter() {
            uniquehash_table.insert(iterable.to_string());
        }
    }

    uniquehash_table
}

fn graph_lookup() {
    /*
     * The graph lookup table which will search for only those kmer that have
     * the overlap with the previous build prefix and suffix combinations so
     * that it can avoid the reads which dont have those suffix and the prefix.
     *
     *
     * */

    let args_parse = OverlapArgs::parse();

    // calling the compute_frequency_table function.

    let compute_lookup: HashSet<String> =
        compute_frequency_table(&args_parse.overlap_file, args_parse.overlap_size);

    #[derive(Debug, Clone, PartialOrd, PartialEq)]
    struct ReadSeq {
        header: String,
        sequence: String,
    }

    let open_readfile = File::open(&args_parse.overlap_file).expect("file not present");
    let readfile = BufReader::new(open_readfile);
    let mut header: Vec<String> = Vec::new();
    let mut seq: Vec<String> = Vec::new();
    for i in readfile.lines() {
        let line = i.expect("line not present");
        if line.starts_with("@") {
            header.push(line.clone())
        }
        if line.starts_with("A")
            || line.starts_with("T")
            || line.starts_with("G")
            || line.starts_with("C")
        {
            seq.push(line)
        }
    }

    // calling the overlap_graph function here to access the aligned prefix and sufix terminal
    // and find their offset in the other reads, where those are not present as prefix and sufix.

    let prefix_sufix_kmer =
        overlap_graph(&args_parse.overlap_file, args_parse.overlap_size).unwrap();

    struct Graphlookup {
        graphedge: String,
        start: usize,
        end: usize,
    }

    let mut graphlookup_table: Vec<Graphlookup> = Vec::new();

    for i in 0..seq.len() {
        for j in prefix_sufix_kmer.iter() {
            graphlookup_table.push(Graphlookup {
                graphedge: j.to_string(),
                start: seq[i].find(j).unwrap(),
                end: seq[i].find(j).unwrap(),
            })
        }
    }

    // those kmers which are unique but are not present at any offset at those prefix and the sufix
    // array.

    let mut unique_prefix_sufix_kmer: Vec<String> = Vec::new();

    for i in 0..prefix_sufix_kmer.len() {
        for j in compute_lookup.iter() {
            if prefix_sufix_kmer[i] != *j {
                unique_prefix_sufix_kmer.push(j.to_string());
            }
        }
    }
}
