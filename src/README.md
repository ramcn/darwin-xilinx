## Intro
Darwin is a genome sequence alignment tool which implements a novel dynamic programming technique called GACT
and a novel hashing scheme called D-SOFT. More details refer the paper: https://stanford.edu/~yatisht/pubs/darwin.pdf
This repo is a clone of https://github.com/yatisht/darwin and experiments of adding few heuristics on it.
## Changes done over vanilla darwin
1. Implement edlib heuristic to optimize GACT
2. Modify the command line arguments to pass configuration parametes like bin_size, threshold, first tile size and first tile threshold
3. Add a new hash table implementation called L1 hash table as an alternative for current seed position table method.
4. L2 hash table code changes 
## Setup
Clone and build edlib library from https://github.com/Martinsos/edlib
## Getting started
Run the following commands to perform a reference-guided assembly of long reads (sampled_reads.fa) aligned to a reference (sample_ref.fa) using darwin. The sample reference is chrI of the yeast genome (sacCer3) and sample reads are generated using PBSIM (https://code.google.com/archive/p/pbsim/) for 20X coverage of the reference. Output alignments are in Multiple Alignment Format (MAF).
```
    $ make
    $ ./darwin_ref_guided data/sample_ref.fa data/sample_reads.fa dsoft_ht_enable bin_size, threshold, ft_size, ft_threshold > out.maf

    $ ./darwin_ref_guided data/sample_ref.fa data/sample_reads.fa 1 256 26 256 120 //  example1
    $ ./darwin_ref_guided data/ecoli.fa data/ecoli_reads.fa 1 256 26 256 120 //  example2
```
## Dataset

The directory data/ has 2 datasets. A simulated dataset called sample_ref and a real ecoli bacteria dataset.



