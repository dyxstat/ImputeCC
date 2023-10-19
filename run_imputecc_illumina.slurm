#!/bin/bash
#SBATCH --job-name=hic
#SBATCH --ntasks=1
#SBATCH --mem=100000mb
#SBATCH --cpus-per-task=60
#SBATCH --time=72:00:00
#SBATCH --partition=qcb

cd /scratch1/yuxuandu/IMPUTECC/MOCK/IMPUTECC_TEST

module purge

eval "$(conda shell.bash hook)"

conda activate MetaCC_env

#python ImputeCC.py impute \
#    -t 60 \
#    ../INPUT_SHORT/final.contigs.fa \
#    ../INPUT_SHORT/contig_info.csv \
#    ../INPUT_SHORT/Normalized_contact_matrix.npz \
#    out_mock_illumina --cover -v

#python ImputeCC.py cluster \
#    -t 60 \
#    ../INPUT_SHORT/final.contigs.fa \
#    ../INPUT_SHORT/contig_info.csv \
#    ../INPUT_SHORT/Normalized_contact_matrix.npz \
#    out_mock_illumina --cover -v

python ImputeCC.py pipeline \
    -t 60 \
    ../INPUT_SHORT/final.contigs.fa \
    ../INPUT_SHORT/contig_info.csv \
    ../INPUT_SHORT/Normalized_contact_matrix.npz \
    out_mock_illumina1 --cover -v
