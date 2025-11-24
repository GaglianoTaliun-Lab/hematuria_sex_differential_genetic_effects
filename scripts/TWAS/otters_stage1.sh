#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --array=6
#SBATCH --time=7-00:00:00
#SBATCH --job-name=male_glom_otters_stage1
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=12G

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

sex="MALE"
tissue="Glom"

chr=${SLURM_ARRAY_TASK_ID}

clump_r2=0.99
N_THREADS=1
N_TASK=20

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
sst_file=${INPUT_DATADIR}/eQTLs_${tissue}_${sex}_MAF_filter_0.01.tsv
output_dir=${project_dir}/otter_output/${tissue}_${sex}/chr${chr}/

######### Run Imputation Models ################
python3 ${OTTERS_DIR}/training.py \
    --OTTERS_dir=${OTTERS_DIR} \
    --SDPR_dir=${SDPR_DIR} \
    --anno_dir=${exp_anno} \
    --geno_dir=${geno_dir} \
    --sst_file=${sst_file} \
    --out_dir=${output_dir} \
    --chrom=${chr} \
    --r2=${clump_r2} \
    --models=lassosum,PT,PRScs,SDPR \
    --lassosum_ld_blocks=EUR.hg19 \
    --thread=${N_TASK}
