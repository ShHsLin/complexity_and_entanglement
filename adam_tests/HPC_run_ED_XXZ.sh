#!/bin/bash
#SBATCH --partition=defq
#SBATCH --array=1-28
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=4g
#SBATCH --time=2:00:00
#below use Linux commands, which will run on compute node

echo "Running on `hostname`"
cd ${SLURM_SUBMIT_DIR}
module purge
module load anaconda-uon/3

source ~/.bashrc
conda activate complexity_env

TASK=${SLURM_ARRAY_TASK_ID}
python -u run_ED_GS_XXZ.py --filename GS_parameters_ED_XXZ.txt --option $TASK
echo "Finished job now"