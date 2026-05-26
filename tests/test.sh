ref=/data/public_db/human_reference/b37/human_g1k_v37.fasta
samtools_path=samtools
mutview -c 1 -p 55299 -b test.bam -r ref.fa -s -S samtools -f Menlo -o test_out1.html 