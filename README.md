# FastaCompressor

Repeat-compressed DNA sequence index. The index can store reads composed of characters ACTG, is immutable and compressed, and has fast random access. Intended to be used for applications that need random access to repetitive read sets which don't change, eg. whole genome sequencing reads of a single sample.

### Usage

See src/example.cpp

### Compression

| Dataset | ONT human 30x | ONT human 55x | HiFi human 32x |
| ----- | ----- | ----- | ----- |
| Base pairs | 90 G | 165 G | 96 G |
| Index size | 13 Gb | 23 Gb | 14 Gb |
| Bits per base pair | 1.3 | 1.1 | 1.2 |
| Index construction peak memory | 43 Gb | 77 Gb | 43 Gb |
| CPU-time (8 cores) | 3h10min | 3h30min | 3h30min |

### Limitation

No non-ATCG characters.

No base pair quality values.

Entire data structure must be kept in memory.

### Algorithm

Split the sequences into subsequences based on minimizers. Assign each distinct subsequence a number, then use a variant of byte pair encoding to compress the numbers.
