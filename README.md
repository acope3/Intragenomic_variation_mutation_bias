# Intragenomic_variation_mutation_bias
A repository presenting the results from our analysis on the effects of intragenomic mutation bias variation within budding yeasts. 


# Processing RNA-seq data

For paired-end data, RNA-seq was processed using the following format.

```bash
fastp -i "${SPECIES_ACC}_1.fastq" -I "${SPECIES_ACC}_2.fastq" -o "${SPECIES_ACC}_1_trimmed.fastq" -O "${SPECIES_ACC}_2_trimmed.fastq" -w 8 -j "${SPECIES_ACC}_fastp.json" -h "${SPECIES_ACC}_fastp.html"
kallisto quant -i ${1}/${1}.index -o "${SPECIES_ACC}_tpm_kallisto" --bias -t 8 "${SPECIES_ACC}_1_trimmed.fastq" "${SPECIES_ACC}_2_trimmed.fastq"
```		

For single-end data, RNA-seq was processed using the following format.

```bash
fastp -i "${SPECIES_ACC}.fastq" -o "${SPECIES_ACC}_trimmed.fastq" -w 8 -j "${SPECIES_ACC}_fastp.json" -h "${SPECIES_ACC}_fastp.html"
kallisto quant -i ${1}/${1}.index -o "${SPECIES_ACC}_tpm_kallisto" --bias --single -l 200 -s 25 -t 8 "${SPECIES_ACC}_trimmed.fastq"
```

# Clustering 

Clustering was performed based on absolute codon frequencies with the `ca.R` script.
See script for command line input parameters. 

# Analysis with ROC-SEMPPR

AnaCoDa is currently available on CRAN.
However, as developers of the software project, we have implemented new functionality, performance improvements, and bug fixes that have not yet been incorporated into the version on CRAN.
We recommend using the version currently available on the `master` branch at github.com/acope3/RibmodelFramework.
We also note that a version was implemented on branch `CTG_Ser` that is able to handle the case of yeasts with CTG coding for serine as opposed to leucine.
The functionality of this branch is consistent with the `master` branch, aside from some changes to the data structures for mapping codons to amino acids.
Specifically, this code changes a static data structure (i.e. initialized at compile time) to a non-static data structure, but this seems to increase run times.
Due to this performance issue, this code has not been merged into the `master` branch, but we anticipate addressing this issue soon.
CDS files were taken from https://figshare.com/articles/dataset/Tempo_and_mode_of_genome_evolution_in_the_budding_yeast_subphylum/5854692 and processed to remove mitochondrial genes, genes including internal stop-codons, genes with a nucleotide length not divisible by 3, and genes not starting with ATG.  
Processed CDS files are available upon request, but are not included here due to their size. 


Analysis can be performed using the following commands.
The first command was used for processing CDS for yeasts with the canonical genetic code.
The second command was used for processing CDS for yeasts with CTG coding for serine by specifying `--codon_table 12`.
Note that this will only work using the branch `CTG_Ser`.

```bash
Rscript --vanilla Scripts/rocAnalysis.R -i "$INPUT" -o "$OUTPUT" -d 20 -s 20000 -a 20 -t 10 -n 8 --est_csp --est_phi --est_hyp  --max_num_runs 2
Rscript --vanilla Scripts/rocAnalysis.R -i "$INPUT" -o "$OUTPUT" -d 20 -s 20000 -a 20 -t 10 -n 8 --est_csp --est_phi --est_hyp  --max_num_runs 2 --codon_table 12
```

# Re-creating figures

Figures can be recreated using the `phi_prediction.Rmd` file. 