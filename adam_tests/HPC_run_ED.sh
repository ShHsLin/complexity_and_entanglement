#!/bin/bash
#SBATCH --partition=defq
#SBATCH --array=1-14
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10g
#SBATCH --time=2:00:00
#below use Linux commands, which will run on compute node

echo "Running on `hostname`"
cd ${SLURM_SUBMIT_DIR}
module purge
module load anaconda-uon/3

source ~/.bashrc
conda activate complexity_env

TASK=${SLURM_ARRAY_TASK_ID}
python -u run_ED_GS.py --filename GS_parameters_ED.txt --option $TASK
echo "Finished job now"