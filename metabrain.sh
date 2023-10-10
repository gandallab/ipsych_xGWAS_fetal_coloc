#!/bin/bash -l 
#$ -cwd
#$ -l h_data=24G,h_rt=12:00:00,highp
#$ -j y
#$ -o /u/project/gandalm/cindywen/ipsych_gwas/log/job.out.metabrain
#$ -m a

# source /u/local/apps/anaconda3/2019.03/bin/activate snakemake
module load anaconda3
source $CONDA_DIR/etc/profile.d/conda.sh
conda activate snakemake

snakemake \
    --snakefile metabrain.smk \
    --jobscript cluster_jobscript.sh \
    --cluster-sync "qsub -l h_data={resources.mem_gb}G,h_rt=00:{resources.time_min}:00,highp -pe shared {resources.num_cores} -o /u/project/gandalm/cindywen/ipsych_gwas/log/job.out.pipeline.MB" \
    --jobs 300 \
    --max-jobs-per-second 10 \
    -T 0 \
    -w 60 \
    --default-resources mem_gb=12 time_min=720 num_cores=1 

