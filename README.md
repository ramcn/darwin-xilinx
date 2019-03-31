## Intro
Darwin is a genome sequence alignment tool which implements a novel dynamic programming technique called GACT
and a novel hashing scheme called D-SOFT. More details refer the paper: https://stanford.edu/~yatisht/pubs/darwin.pdf
This repo is a clone of https://github.com/yatisht/darwin and experiments of porting it to Xilinx FPGA
## Changes done over vanilla darwin
1. Implement the xilinx port of the AlignwithBT function  
## Dataset
The directory data/ has 2 datasets. A simulated dataset called sample_ref and a real ecoli bacteria dataset.



