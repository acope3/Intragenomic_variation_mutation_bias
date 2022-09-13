#!/bin/bash
INPUT="/data2/Labella2019/Genomes/cds_cleaned/$1"
OUTPUT="Results/ConstMut/$1/"
nohup Rscript --vanilla Scripts/rocAnalysis.R -i "$INPUT" -o "$OUTPUT" -d 0 -s 20000 -a 20 -t 5 -n 8 --est_csp --est_phi --est_hyp --max_num_runs 2 --codon_table 1 &> Results/ConstMut/$1.Rout &
wait
