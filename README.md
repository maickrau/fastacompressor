# FastaCompressor

Repeat-compressed DNA sequence index. The index can store reads composed of characters ACTG, is immutable and compressed, and has fast random access. Intended to be used for applications that need random access to repetitive read sets which don't change, eg. whole genome sequencing reads of a single sample.

### Usage

See src/example.cpp

### Advantages

Decent speed random access.

Low memory use for built index: 90Gb of sequence (30x coverage human nanopore) becomes 13Gb index (~1.3 bits per bp). 165Gb sequence becomes 23Gb (~1.1 bits per bp).

### Disadvantages

Building the index requires more memory than storing the built index: 90Gb of sequence requires peak 43Gb of memory to build the index. 165Gb sequence requires 77Gb.

No non-ATCG characters.

No base pair quality values.

Entire data structure must be kept in memory.

### Algorithm

Split the sequences into subsequences based on minimizers. Assign each distinct subsequence a number, then use a variant of byte pair encoding to compress the numbers.
