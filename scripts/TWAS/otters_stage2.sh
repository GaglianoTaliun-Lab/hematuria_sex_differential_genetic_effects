#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --array=6
#SBATCH --time=1:00:00
#SBATCH --job-name=female_tube_otters_stage2
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=6G

module load StdEnv/2020
module load gcc
module load python/3.9
module load tabix
module load r/4.3.1
module load plink/1.9b_6.21-x86_64 r/4

export R_LIBS=~/.local/R/$EBVERSIONR/.

project_dir="/home/fridald4/projects/def-gsarah/fridald4/hematuria_sexspecific"

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r ${project_dir}/tools/OTTERS/otters-requirements.txt

sex="FEMALE"
tissue="Tube"

chr=${SLURM_ARRAY_TASK_ID}

clump_r2=0.99
N_THREADS=1
N_TASK=10

######### Run OTTERS to prepare inputs for imputation models ############
### Set Tool Directory
OTTERS_DIR=${project_dir}/tools/OTTERS
SDPR_DIR=${project_dir}/tools/SDPR
INPUT_DATADIR=${project_dir}/otter_input

# make sure the dynamic libraries are not changed
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/MKL/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SDPR_DIR}/gsl/lib
# prevent automatically using  all available cores on a compute node
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS

### Set Required Inputs
exp_anno=${INPUT_DATADIR}/annotation_file.tsv
geno_dir=${project_dir}/reference_data/g1000_eur_chr${chr}
gwas_file=${INPUT_DATADIR}/hematuria_${sex}.tsv
out_stage1=${project_dir}/otter_output/${tissue}_${sex}/chr${chr}

# Set output directory for STAGE II
twas_dir=${project_dir}/otter_twas/${tissue}_${sex}/chr${chr}

# STAGE II
# gene-based association test using eQTL-weight trained from P+T, lassosum, SDPR and PRS-CS.

# Note from developer: 
# We have recently discovered that SDPR is not particularly sensitive to the reference panel. 
# Therefore, if you notice that many genes with significant ACAT p-values are predominantly driven by SDPR, 
# I recommend using only the p-values from lassosum.txt, P0.001.txt, P0.05.txt, PRScs.txt, 
# and conducting an ACAT analysis on them.

python3 ${OTTERS_DIR}/testing.py \
    --OTTERS_dir=${OTTERS_DIR} \
    --weight_dir=${out_stage1} \
    --models=P0.001,P0.05,lassosum,SDPR,PRScs \
    --anno_dir=${exp_anno} \
    --geno_dir=${geno_dir} \
    --out_dir=${twas_dir} \
    --gwas_file=${gwas_file} \
    --chrom=${chr} \
    --thread=${N_TASK}
