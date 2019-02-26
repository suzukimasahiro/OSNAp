# binary_seq.py

## Introduction
Phylogeny of bacterial genomes based on whole genome sequence is difficult because of differences in their DNA sequence lengths, especially in plasmids. By replacing genome data to binary sequences consisting of 1 or 0 based on presence or absence of ORFs, differences of the lengths are resolved, and sequences become amenable to comparisons with each other. Binary sequences are generated from plasmid or chromosome sequence data by this python script program. Binary sequences can be visualized using SplitsTree4 (https://ab.inf.uni-tuebingen.de/software/splitstree4).

## Installation
Copy “binary_seq.py” to an appropriate directory, then mark the file as executable with the chmod command.

## Tested environment
- CentOS 7
- Python 2.7.15
- BLAST 2.7.1
- Prodigal 2.6.3

## Dependencies
- NCBI BLAST
- Bio python
- Prodigal 

## Quick Start
```% python binary_seq.py -g gb_file_list -o prefix_of output_file```

## Output files
| File type | Description |
|:----------|:------------|
| _all.fasta | Multi fasta containing all sequences. |
| _all.nex | Binary sequence including all sequences. |
| _binary.csv | CSV file of binary sequence without labels. |
| _dice.txt	 | Dice indices between two sequences. |
| _elementORF.fasta | Multi fasta containing elemental ORFs. |
| _no_annotation.txt | List of GenBank data without annotation. |
| _info.csv	 | Summary of input sequence files. |
| _orfinfo.txt | List of elemental ORF names. |
| _share.txt | The number of ORFs shared between two sequences. |

## Options
-g, --gb
- GeneBank file list (GB file must be designated with full PATH)  
-f, --fasta
- Fasta file list (strain<TAB>fasta_file, fasta file must be designated with full PATH).  
-o, --out
- Prefix of output file name.  
-e, --evalue
- E value (default = 0.001) setting for BLASTN'  
-c, --cover
- Cover ratio against query sequence length in BLASTN (default = 80)'  
-i, --identity
- Threshold of sequence identity in BLASTN (default = 80)'  
--predict
- Predict all ORFs using Prodigal, including gb files with annotation data.  
--num_process
- Set the number of processes (default = 4) running simultaneously.

## License
MIT
https://opensource.org/licenses/mit-license.php

## Author
- Masahiro Suzuki
- Web:http://www.fujita-hu.ac.jp/~microb/index.html
