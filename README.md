# MutView
MutView is a tool to view the reads cover a mutation in bam file, output file is html format. It's so fast and efficient, just like a lightweight IGV(Integrative Genomics Viewer).

## Requires
- python >= 3.7
- samtools

## Installation
```
pip install mutview
```
or
```
pip install git+https://github.com/TuBieJun/MutView.git
```

## Example use
```
mutview -s -n 250 -c 2 -p 25466919 -r path_of_ref_fasta -S path_of_samtools -b input_bam -o out.html
```

## Argument description
```
$ mutview
usage: mutview [-h] -p POS -c CHROM -b BAM -r REF -o OUT [-w WINDOW_SIZE] [-s]
               [-n NUM_READS] [-q BASE_Q] [-Q ALN_Q] [-f FONT] [-S SAMTOOLS]

optional arguments:
  -h, --help            show this help message and exit
  -p POS, --pos POS     the mutation position
  -c CHROM, --chrom CHROM
                        the chrom of mutation
  -b BAM, --bam BAM     the input bam file
  -r REF, --ref REF     the reference genome fasta file
  -o OUT, --out OUT     the output html file
  -w WINDOW_SIZE, --window_size WINDOW_SIZE
                        the window size of reference genome around the
                        mutation, default is 300
  -s, --sort            sort as base like igv, default is False
  -n NUM_READS, --num_reads NUM_READS
                        the max reads number to view default is 100
  -q BASE_Q, --base_q BASE_Q
                        the min base quality, reads with bq below the
                        specified threshold will be lowercases. default is 20.
  -Q ALN_Q, --aln_q ALN_Q
                        the min aln quality, bases with mq below the specified
                        threshold will be highlighted. default is 20.
  -f FONT, --font FONT  the the font of the out html file, it should be
                        monospaced font, default is consolas for windows, if
                        you view the out html file in mac, you need set the
                        value to Menlo or other monospaced fonts
  -S SAMTOOLS, --samtools SAMTOOLS
                        the samtools path, default is samtools


```

## Out html Description
- The bases that are inconsistent with the reference genome will be highlighted.
- The bases with quality below the specified threshold will be represented in lowercase letters.
- The background color of reads with alignment quality below the specified threshold will be highlighted.
- Indels are represented using the characters "I" and "_". 
- The count of each base type will be displayed at the top
- The left sidebar will display the name of the reads, pair information, and positive/negative strand information.

**Here is example out:**

![example](https://github.com/TuBieJun/MutView/raw/master/screenshots/example.png)  


