from os.path import join
import os
import numpy as np
import pandas as pd
import sys

"""
Run coloc on fetal + PsychENCODE QTLs
- fetal e/iso/sQTL, trimester e/iso/sQTL
- PEC e/isoQTL
"""

TEST_TABLE = pd.read_table("/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP.tsv")

rule all:
    input:
        expand(
            "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/coloc_{annot}.done",
            annot=["eqtl", "isoqtl", "sqtl"],
            locus=TEST_TABLE.locus.values,
        ),

rule coloc:
    input:
        "/u/project/gandalm/cindywen/ipsych_gwas/data/gwas_indexSNP.tsv",
    output:
        "/u/project/gandalm/cindywen/ipsych_gwas/out/locus{locus}/coloc_{annot}.done",
    shell:
        """
        . /u/local/Modules/default/init/modules.sh
        module load R/4.1.0-BIO

        Rscript scripts/coloc_fetal_pec.R \
            --locus {wildcards.locus} \
            --annot {wildcards.annot}
        touch {output[0]}
        """
