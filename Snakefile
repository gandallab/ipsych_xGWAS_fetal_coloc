from os.path import join
import os
import numpy as np
import pandas as pd
import sys


# configfile: "config.yaml"

"""
rules:
    - locus_egene: get a list of candidate genes for colocalization
    - write_locus_gene_list: write a list of all locus-gene pairs for parallelization
    - extract_cis_assoc_gene: for each gene, extract all cis assoc in QTL data
    - zscore: for each gene, write zscore files for QTL and GWAS for the common variants; write list of common variants for LD calculation
    - ld_qtl: write LD matrix using QTL genotype data
    - ld_gwas: subset locus LD matrix to common variants
    - ecaviar: run ecaviar
    - analyze: extract CLPP > 0.01
"""

TEST_TABLE = pd.read_table("/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP.tsv")

LOCUS_GENE_TABLE = pd.read_table("locus_gene_list.txt")
LOCUS_ISO_TABLE = pd.read_table("locus_isoform_list.txt")
LOCUS_INTRON_TABLE = pd.read_table("locus_intron_list.txt")
LOCUS_GENE1_TABLE = pd.read_table("locus_gene_list_tri1.txt")
LOCUS_GENE2_TABLE = pd.read_table("locus_gene_list_tri2.txt")

LOCUS_ISO1_TABLE = pd.read_table("locus_isoform_list_tri1_35hcp.txt")
LOCUS_ISO2_TABLE = pd.read_table("locus_isoform_list_tri2_20hcp.txt")

LOCUS_INTRON1_TABLE = pd.read_table("locus_intron_list_tri1_15hcp.txt")
LOCUS_INTRON2_TABLE = pd.read_table("locus_intron_list_tri2_10hcp.txt")

# not ideal workaround to avoid wildcards ambiguity
ruleorder: analyze_ecaviar_intron_tri > analyze_ecaviar_isoform_tri > analyze_ecaviar_tri > ecaviar_intron_tri > ecaviar_isoform_tri > ecaviar_tri > ld_gwas_intron_tri > ld_gwas_isoform_tri > ld_gwas_tri > ld_qtl_intron_tri > ld_qtl_isoform_tri > ld_qtl_tri > zscore_intron_tri > zscore_isoform_tri > zscore_tri > extract_cis_assoc_intron_tri > extract_cis_assoc_isoform_tri > analyze_ecaviar_intron > ecaviar_intron > ld_gwas_intron > ld_qtl_intron > zscore_intron > extract_cis_assoc_intron > analyze_ecaviar_isoform > ecaviar_isoform > ld_gwas_isoform > ld_qtl_isoform > zscore_isoform > extract_cis_assoc_isoform > analyze_ecaviar > ecaviar > ld_gwas > ld_qtl > zscore > extract_cis_assoc_gene

rule all:
    input:
        # expand(
        #     "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}.done",
        #     zip,
        #     locus = LOCUS_GENE_TABLE.locus.values,
        #     gene = LOCUS_GENE_TABLE.gene.values,
        # ),
        # expand(
        #     "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}.done",
        #     zip,
        #     locus = LOCUS_ISO_TABLE.locus.values,
        #     gene = LOCUS_ISO_TABLE.isoform.values,
        # ),
        # expand(
        #     "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}.done",
        #     zip,
        #     locus = LOCUS_INTRON_TABLE.locus.values,
        #     gene = LOCUS_INTRON_TABLE.intron.values,
        # ),
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_tri1.done",
            zip,
            locus = LOCUS_GENE1_TABLE.locus.values,
            gene = LOCUS_GENE1_TABLE.gene.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_tri2.done",
            zip,
            locus = LOCUS_GENE2_TABLE.locus.values,
            gene = LOCUS_GENE2_TABLE.gene.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_tri1_35hcp.done",
            zip,
            locus = LOCUS_ISO1_TABLE.locus.values,
            gene = LOCUS_ISO1_TABLE.isoform.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_tri2_20hcp.done",
            zip,
            locus = LOCUS_ISO2_TABLE.locus.values,
            gene = LOCUS_ISO2_TABLE.isoform.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_tri1_15hcp.done",
            zip,
            locus = LOCUS_INTRON1_TABLE.locus.values,
            gene = LOCUS_INTRON1_TABLE.intron.values,
        ),
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_tri2_10hcp.done",
            zip,
            locus = LOCUS_INTRON2_TABLE.locus.values,
            gene = LOCUS_INTRON2_TABLE.intron.values,
        ),

############# eQTL #############
rule locus_egene:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/mixed_nominal_90hcp/significant_assoc.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_egene.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        mkdir -p /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/
        Rscript scripts/locus_egene.R \
            --locus {wildcards.locus} \
            --annot eqtl
        """

rule write_locus_gene_list:
    input:
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_egene.txt",
            locus=TEST_TABLE.locus.values,
        )
    output:
        "locus_gene_list.txt",
    shell:
        """
        cat /u/project/gandalm/cindywen/ipsych_gwas/out/locus*/locus_egene.txt > locus_gene_list.txt
        sed  -i '1i locus\tgene' locus_gene_list.txt
        """

rule extract_cis_assoc_gene:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/mixed_nominal_90hcp/gtex.allpairs.txt.gz",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_all_pairs.txt",
    shell:
        """
        awk -v a="{wildcards.gene}" '$1 == a {{print}}' <(zcat {input[0]}) > {output[0]}
        """

rule zscore:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_all_pairs.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_zscore.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_zscore.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/zscore.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene}
        """

rule ld_qtl:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}.ld",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load plink/1.90b624

        plink --bfile /u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/test \
            --r \
            --matrix \
            --extract /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/{wildcards.gene}_snps.txt \
            --out /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/{wildcards.gene}
        """

rule ld_gwas:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas.ld",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/ld_gwas.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene} \
            --annot eqtl
        """

rule ecaviar:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_zscore.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_zscore.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_ecaviar_col",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/10.2.0
        cd /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/

        /u/project/gandalm/shared/apps/caviar/CAVIAR-C++/eCAVIAR \
                -o {wildcards.gene}_ecaviar \
                -l {wildcards.gene}_gwas.ld \
                -z {wildcards.gene}_gwas_zscore.txt \
                -l {wildcards.gene}.ld \
                -z {wildcards.gene}_zscore.txt \
                -f 1 \
                -c 2 \
                -r 0.95
        """

rule analyze_ecaviar:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_ecaviar_col",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}.done"
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO

        Rscript scripts/analyze.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene} \
            --annot eqtl
        touch {output[0]}
        """


############# isoQTL #############
rule locus_isoform:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/mixed_nominal_70hcp/significant_assoc.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_isoform.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        mkdir -p /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/
        Rscript scripts/locus_egene.R \
            --locus {wildcards.locus} \
            --annot isoqtl
        """

rule write_locus_isoform_list:
    input:
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_isoform.txt",
            locus=TEST_TABLE.locus.values,
        )
    output:
        "locus_isoform_list.txt",
    shell:
        """
        cat /u/project/gandalm/cindywen/ipsych_gwas/out/locus*/locus_isoform.txt > {output[0]}
        sed  -i '1i locus\tisoform' {output[0]}
        """

rule extract_cis_assoc_isoform:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/mixed_nominal_70hcp/gtex.allpairs.uniq.txt.gz",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_all_pairs.txt",
    shell:
        """
        awk -v a="{wildcards.gene}" '$1 == a {{print}}' <(zcat {input[0]}) > {output[0]}
        """

rule zscore_isoform:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_all_pairs.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_zscore.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_zscore.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/zscore.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene}
        """

rule ld_qtl_isoform:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}.ld",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load plink/1.90b624

        plink --bfile /u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/test \
            --r \
            --matrix \
            --extract /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/{wildcards.gene}_snps.txt \
            --out /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/{wildcards.gene}
        """

rule ld_gwas_isoform:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas.ld",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/ld_gwas.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene} \
            --annot isoqtl
        """

rule ecaviar_isoform:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_zscore.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_zscore.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_ecaviar_col",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/10.2.0
        cd /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/

        /u/project/gandalm/shared/apps/caviar/CAVIAR-C++/eCAVIAR \
                -o {wildcards.gene}_ecaviar \
                -l {wildcards.gene}_gwas.ld \
                -z {wildcards.gene}_gwas_zscore.txt \
                -l {wildcards.gene}.ld \
                -z {wildcards.gene}_zscore.txt \
                -f 1 \
                -c 2 \
                -r 0.95
        """

rule analyze_ecaviar_isoform:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_ecaviar_col",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}.done"
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO

        Rscript scripts/analyze.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene} \
            --annot isoqtl
        touch {output[0]}
        """
############# sQTL #############
rule locus_intron:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/mixed_nominal_40hcp_1e6/significant_assoc.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_intron.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        mkdir -p /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/
        Rscript scripts/locus_egene.R \
            --locus {wildcards.locus} \
            --annot sqtl
        """

rule write_locus_intron_list:
    input:
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_intron.txt",
            locus=TEST_TABLE.locus.values,
        )
    output:
        "locus_intron_list.txt",
    shell:
        """
        cat /u/project/gandalm/cindywen/ipsych_gwas/out/locus*/locus_intron.txt > {output[0]}
        sed  -i '1i locus\tintron' {output[0]}
        """

rule extract_cis_assoc_intron:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/mixed_nominal_40hcp_1e6/gtex.allpairs.uniq.txt.gz"
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_all_pairs.txt",
    shell:
        """
        awk -v a="{wildcards.gene}" '$1 == a {{print}}' <(zcat {input[0]}) > {output[0]}
        """

rule zscore_intron:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_all_pairs.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_zscore.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_zscore.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/zscore.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene}
        """

rule ld_qtl_intron:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}.ld",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load plink/1.90b624

        plink --bfile /u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/test \
            --r \
            --matrix \
            --extract /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/{wildcards.gene}_snps.txt \
            --out /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/{wildcards.gene}
        """

rule ld_gwas_intron:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas.ld",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/ld_gwas.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene} \
            --annot sqtl
        """

rule ecaviar_intron:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_zscore.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_zscore.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_ecaviar_col",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/10.2.0
        cd /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/

        /u/project/gandalm/shared/apps/caviar/CAVIAR-C++/eCAVIAR \
                -o {wildcards.gene}_ecaviar \
                -l {wildcards.gene}_gwas.ld \
                -z {wildcards.gene}_gwas_zscore.txt \
                -l {wildcards.gene}.ld \
                -z {wildcards.gene}_zscore.txt \
                -f 1 \
                -c 2 \
                -r 0.95
        """

rule analyze_ecaviar_intron:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_ecaviar_col",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}.done"
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO

        Rscript scripts/analyze.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene} \
            --annot sqtl
        touch {output[0]}
        """
############# tri-eQTL #############
rule locus_egene_tri:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/eur_trimester/T{tri}-nominal_significant_assoc.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_egene_tri{tri}.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/locus_egene.R \
            --locus {wildcards.locus} \
            --annot eqtl_tri{wildcards.tri}
        """

rule write_locus_gene_list_tri:
    input:
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_egene_tri{{tri}}.txt",
            locus=TEST_TABLE.locus.values,
        )
    output:
        "locus_gene_list_tri{tri}.txt",
    shell:
        """
        cat /u/project/gandalm/cindywen/ipsych_gwas/out/locus*/locus_egene_tri{wildcards.tri}.txt > locus_gene_list_tri{wildcards.tri}.txt
        sed  -i '1i locus\tgene' locus_gene_list_tri{wildcards.tri}.txt
        """

rule extract_cis_assoc_gene_tri:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/eur_trimester/T{tri}-all.chunks.txt.gz",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_all_pairs_tri{tri}.txt",
    shell:
        """
        awk -v a="{wildcards.gene}" '$1 == a {{print}}' <(zcat {input[0]}) > {output[0]}
        """

rule zscore_tri:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_all_pairs_tri{tri}.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_zscore_tri{tri}.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_zscore_tri{tri}.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps_tri{tri}.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/zscore_tri.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene} \
            --annot tri{wildcards.tri}
        """

rule ld_qtl_tri:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps_tri{tri}.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_tri{tri}.ld",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load plink/1.90b624

        plink --bfile /u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/eur/filtered.hg19.sorted.removeRel \
            --r \
            --matrix \
            --extract {input[0]} \
            --out /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/{wildcards.gene}_tri{wildcards.tri}
        """

rule ld_gwas_tri:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps_tri{tri}.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_tri{tri}.ld",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/ld_gwas.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene} \
            --annot tri{wildcards.tri}
        """

rule ecaviar_tri:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_tri{tri}.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_tri{tri}.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_zscore_tri{tri}.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_zscore_tri{tri}.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_ecaviar_tri{tri}_col",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/10.2.0
        cd /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/

        /u/project/gandalm/shared/apps/caviar/CAVIAR-C++/eCAVIAR \
                -o {wildcards.gene}_ecaviar_tri{wildcards.tri} \
                -l {wildcards.gene}_gwas_tri{wildcards.tri}.ld \
                -z {wildcards.gene}_gwas_zscore_tri{wildcards.tri}.txt \
                -l {wildcards.gene}_tri{wildcards.tri}.ld \
                -z {wildcards.gene}_zscore_tri{wildcards.tri}.txt \
                -f 1 \
                -c 2 \
                -r 0.95
        """

rule analyze_ecaviar_tri:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_ecaviar_tri{tri}_col",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_tri{tri}.done"
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO

        Rscript scripts/analyze.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene} \
            --annot tri{wildcards.tri}
        touch {output[0]}
        """
############# tri-isoQTL #############
rule locus_isoform_tri:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/tri{tri}_nominal_{hcp}hcp/significant_assoc.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_isoform_tri{tri}_{hcp}hcp.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/locus_egene.R \
            --locus {wildcards.locus} \
            --annot isoqtl_tri{wildcards.tri}
        """

rule write_locus_isoform_list_tri:
    input:
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_isoform_tri{{tri}}_{{hcp}}hcp.txt",
            locus=TEST_TABLE.locus.values,
        )
    output:
        "locus_isoform_list_tri{tri}_{hcp}hcp.txt",
    shell:
        """
        cat /u/project/gandalm/cindywen/ipsych_gwas/out/locus*/locus_isoform_tri{wildcards.tri}_{wildcards.hcp}hcp.txt > locus_isoform_list_tri{wildcards.tri}_{wildcards.hcp}hcp.txt
        sed  -i '1i locus\tisoform' locus_isoform_list_tri{wildcards.tri}_{wildcards.hcp}hcp.txt
        """

rule extract_cis_assoc_isoform_tri:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/tri{tri}_nominal_{hcp}hcp/all.chunks.txt.gz",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_all_pairs_tri{tri}_{hcp}hcp.txt",
    shell:
        """
        awk -v a="{wildcards.gene}" '$1 == a {{print}}' <(zcat {input[0]}) > {output[0]}
        """

rule zscore_isoform_tri:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_all_pairs_tri{tri}_{hcp}hcp.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_zscore_tri{tri}_{hcp}hcp.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_zscore_tri{tri}_{hcp}hcp.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps_tri{tri}_{hcp}hcp.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/zscore_tri.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene} \
            --annot tri{wildcards.tri}_{wildcards.hcp}hcp
        """

rule ld_qtl_isoform_tri:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps_tri{tri}_{hcp}hcp.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_tri{tri}_{hcp}hcp.ld",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load plink/1.90b624

        plink --bfile /u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/eur/filtered.hg19.sorted.removeRel \
            --r \
            --matrix \
            --extract {input[0]} \
            --out /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/{wildcards.gene}_tri{wildcards.tri}_{wildcards.hcp}hcp
        """

rule ld_gwas_isoform_tri:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps_tri{tri}_{hcp}hcp.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_tri{tri}_{hcp}hcp.ld",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/ld_gwas.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene} \
            --annot tri{wildcards.tri}_{wildcards.hcp}hcp
        """

rule ecaviar_isoform_tri:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_tri{tri}_{hcp}hcp.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_tri{tri}_{hcp}hcp.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_zscore_tri{tri}_{hcp}hcp.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_zscore_tri{tri}_{hcp}hcp.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_ecaviar_tri{tri}_{hcp}hcp_col",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/10.2.0
        cd /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/

        /u/project/gandalm/shared/apps/caviar/CAVIAR-C++/eCAVIAR \
                -o {wildcards.gene}_ecaviar_tri{wildcards.tri}_{wildcards.hcp}hcp \
                -l {wildcards.gene}_gwas_tri{wildcards.tri}_{wildcards.hcp}hcp.ld \
                -z {wildcards.gene}_gwas_zscore_tri{wildcards.tri}_{wildcards.hcp}hcp.txt \
                -l {wildcards.gene}_tri{wildcards.tri}_{wildcards.hcp}hcp.ld \
                -z {wildcards.gene}_zscore_tri{wildcards.tri}_{wildcards.hcp}hcp.txt \
                -f 1 \
                -c 2 \
                -r 0.95
        """

rule analyze_ecaviar_isoform_tri:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_ecaviar_tri{tri}_{hcp}hcp_col",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_tri{tri}_{hcp}hcp.done"
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO

        Rscript scripts/analyze.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene} \
            --annot tri{wildcards.tri}_{wildcards.hcp}hcp
        touch {output[0]}
        """
############# tri-sQTL #############
rule locus_intron_tri:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/tri{tri}_nominal_{hcp}hcp/significant_assoc.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_intron_tri{tri}_{hcp}hcp.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/locus_egene.R \
            --locus {wildcards.locus} \
            --annot sqtl_tri{wildcards.tri}
        """

rule write_locus_intron_list_tri:
    input:
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_intron_tri{{tri}}_{{hcp}}hcp.txt",
            locus=TEST_TABLE.locus.values,
        )
    output:
        "locus_intron_list_tri{tri}_{hcp}hcp.txt",
    shell:
        """
        cat /u/project/gandalm/cindywen/ipsych_gwas/out/locus*/locus_intron_tri{wildcards.tri}_{wildcards.hcp}hcp.txt > locus_intron_list_tri{wildcards.tri}_{wildcards.hcp}hcp.txt
        sed  -i '1i locus\tintron' locus_intron_list_tri{wildcards.tri}_{wildcards.hcp}hcp.txt
        """

rule extract_cis_assoc_intron_tri:
    input:
        "/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/tri{tri}_nominal_{hcp}hcp/all.chunks.txt.gz",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_all_pairs_tri{tri}_{hcp}hcp.txt",
    shell:
        """
        awk -v a="{wildcards.gene}" '$1 == a {{print}}' <(zcat {input[0]}) > {output[0]}
        """

rule zscore_intron_tri:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_all_pairs_tri{tri}_{hcp}hcp.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_zscore_tri{tri}_{hcp}hcp.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_zscore_tri{tri}_{hcp}hcp.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps_tri{tri}_{hcp}hcp.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/zscore_tri.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene} \
            --annot tri{wildcards.tri}_{wildcards.hcp}hcp
        """

rule ld_qtl_intron_tri:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps_tri{tri}_{hcp}hcp.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_tri{tri}_{hcp}hcp.ld",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load plink/1.90b624

        plink --bfile /u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/eur/filtered.hg19.sorted.removeRel \
            --r \
            --matrix \
            --extract {input[0]} \
            --out /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/{wildcards.gene}_tri{wildcards.tri}_{wildcards.hcp}hcp
        """

rule ld_gwas_intron_tri:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_snps_tri{tri}_{hcp}hcp.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_tri{tri}_{hcp}hcp.ld",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO
        Rscript scripts/ld_gwas.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene} \
            --annot tri{wildcards.tri}_{wildcards.hcp}hcp
        """

rule ecaviar_intron_tri:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_tri{tri}_{hcp}hcp.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_tri{tri}_{hcp}hcp.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_zscore_tri{tri}_{hcp}hcp.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_gwas_zscore_tri{tri}_{hcp}hcp.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_ecaviar_tri{tri}_{hcp}hcp_col",
    resources:
        mem_gb=8,
        time_min=480,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/10.2.0
        cd /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/

        /u/project/gandalm/shared/apps/caviar/CAVIAR-C++/eCAVIAR \
                -o {wildcards.gene}_ecaviar_tri{wildcards.tri}_{wildcards.hcp}hcp \
                -l {wildcards.gene}_gwas_tri{wildcards.tri}_{wildcards.hcp}hcp.ld \
                -z {wildcards.gene}_gwas_zscore_tri{wildcards.tri}_{wildcards.hcp}hcp.txt \
                -l {wildcards.gene}_tri{wildcards.tri}_{wildcards.hcp}hcp.ld \
                -z {wildcards.gene}_zscore_tri{wildcards.tri}_{wildcards.hcp}hcp.txt \
                -f 1 \
                -c 2 \
                -r 0.95
        """

rule analyze_ecaviar_intron_tri:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_ecaviar_tri{tri}_{hcp}hcp_col",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{gene}_tri{tri}_{hcp}hcp.done"
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO

        Rscript scripts/analyze.R \
            --locus {wildcards.locus} \
            --gene {wildcards.gene} \
            --annot tri{wildcards.tri}_{wildcards.hcp}hcp
        touch {output[0]}
        """