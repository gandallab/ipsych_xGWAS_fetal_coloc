# Colocalization between ipsych cross-disorder GWAS and developing brain xQTLs
- e/iso/sQTL, trimester-specific e/iso/sQTL
- `Snakefile` rules for each annotation
```
    - locus_egene: get a list of candidate genes for colocalization
    - write_locus_gene_list: write a list of all locus-gene pairs for parallelization
    - extract_cis_assoc_gene: for each gene, extract all cis assoc in QTL data
    - zscore: for each gene, write zscore files for QTL and GWAS for the common variants; write list of common variants for LD calculation
    - ld_qtl: write LD matrix using QTL genotype data
    - ld_gwas: subset locus LD matrix to common variants
    - ecaviar: run ecaviar
    - analyze: extract CLPP > 0.01
```
