from os.path import join
import os
import numpy as np
import pandas as pd
import sys

TEST_TABLE = pd.read_table("/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP.tsv")
LOCUS_PROBE_TABLE = pd.read_table("locus_probe_list.txt")

def get_1kg_eur_bfile_chr(wildcards):
    TEST_TABLE = pd.read_table("/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP.tsv").set_index(
        "locus", drop=True
    )
    chromosome = TEST_TABLE.loc[int(wildcards.locus), "CHR"]
    return (
        "/u/project/gandalm/shared/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC."
        + str(chromosome)
    )


rule all:
    input:
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/coloc_thistle.done",
            locus=TEST_TABLE.locus.values,
        ),
        "locus_probe_list.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/thistle_CLPP_sig.txt",

# Run coloc
rule coloc:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP.tsv",
        "/u/project/gandalm/shared/GenomicDatasets/THISTLE_BrainMeta/BrainMeta_cis_sqtl_summary/job.out.convert",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/coloc_thistle.done",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO

        Rscript scripts/coloc_thistle.R \
            --locus {wildcards.locus}
        touch {output[0]}
        """

# Run eCAVIAR
rule locus_probe:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP.tsv",
        "/u/project/gandalm/shared/GenomicDatasets/THISTLE_BrainMeta/BrainMeta_cis_sqtl_summary/job.out.convert",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_probe.txt",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO

        Rscript scripts/locus_probe.R \
            --locus {wildcards.locus}
        """

rule write_locus_probe_list:
    input:
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/locus_probe.txt",
            locus=TEST_TABLE.locus.values,
        ),
    output:
        "locus_probe_list.txt",
    shell:
        """
        cat /u/project/gandalm/cindywen/ipsych_gwas/out/locus*/locus_probe.txt > {output[0]}
        sed  -i '1i locus\tprobe' {output[0]}
        """

rule ld_qtl:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{probe}_thistle_snps.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{probe}_thistle.ld",
    params:
        get_1kg_eur_bfile_chr,
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load plink/1.90b624

        plink --bfile {params[0]} \
            --r \
            --matrix \
            --extract {input[0]} \
            --out /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/{wildcards.probe}_thistle
        """

rule ecaviar:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{probe}_thistle_gwas.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{probe}_thistle.ld",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{probe}_thistle_zscore.txt",
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{probe}_thistle_gwas_zscore.txt",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{probe}_thistle_ecaviar_col",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load gcc/10.2.0
        cd /u/project/gandalm/cindywen/ipsych_gwas/out/locus{wildcards.locus}/

        /u/project/gandalm/shared/apps/caviar/CAVIAR-C++/eCAVIAR \
                -o {wildcards.probe}_thistle_ecaviar \
                -l {wildcards.probe}_thistle_gwas.ld \
                -z {wildcards.probe}_thistle_gwas_zscore.txt \
                -l {wildcards.probe}_thistle.ld \
                -z {wildcards.probe}_thistle_zscore.txt \
                -f 1 \
                -c 2 \
                -r 0.95
        """

rule analyze_ecaviar:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{probe}_thistle_ecaviar_col",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{probe}_thistle.done"
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO

        Rscript scripts/analyze.R \
            --locus {wildcards.locus} \
            --gene {wildcards.probe} \
            --annot thistle
        touch {output[0]}
        """

rule concat:
    input:
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/{probe}_thistle.done",
            zip,
            locus = LOCUS_PROBE_TABLE.locus.values,
            probe = LOCUS_PROBE_TABLE.probe.values,
        ),
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/thistle_CLPP_sig.txt",
    shell:
        """
        cd /u/project/gandalm/cindywen/ipsych_gwas/out/
        awk 'FNR==1 && NR!=1{{next;}}{{print}}' locus*/*thistle_ecaviar_col_sig.txt > {output[0]}
        """