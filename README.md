# MutView
MutView is a tool to look the reads cover a mutation in Bam file, output file is html format, just like a lightweight IGV(Integrative Genomics Viewer).

## user example
python MutView.py -r refpath -p 25466919 -c 2 -b test.bam -s > test_2_25466919.html

## argument description
```
-p      the mutation postion in reference
-c      the chrom of mutation
-b      the bam file path
-r      the reference fasta path
-w      the show windows size of genome
-s      the flag of output reads  sort as base, default is Flase
-n      the max reads to show
```

## the example html look

