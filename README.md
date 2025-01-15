# rust-overlap-graph

 - BTreeMap based approach to implement a overlap graph and then it computes the graph lookup table and search for the offset in the reads where those prefix and suffix are not present. 
 - Build a BTreeMap, for the overlap graph, 
 - Make a unique kmer has offset table.
 - search for the offset of the matched suffixes.
 - please see the last commit message and if it says compiled binary then it is completed or else still in development version.
 


```
cargo build 
```

```
➜  gauravsablok rust-overlap-graph git:(master) ✗ ./target/debug/rust-overlapgraph -h
Usage: rust-overlapgraph <OVERLAP_FILE> <OVERLAP_SIZE>

Arguments:
  <OVERLAP_FILE>  please provide the path to the fastq file
  <OVERLAP_SIZE>  please provide the overlap length

Options:
  -h, --help     Print help
  -V, --version  Print version

```

Gaurav Sablok
