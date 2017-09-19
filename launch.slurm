#!/bin/bash

#SBATCH -N 1			# use 1 node
#SBATCH -n 1			# run 1 MPI task per node
#SBATCH --cpus-per-task=20	# use 20 cpus per OpenMP task
#SBATCH --threads-per-core=1	# use 2 OpenMP thread per physical cpu corec (so 40 OpenMP threads in total)
#SBATCH --time=6:00:00	# use 300 hours of wall-time
#SBATCH --mail-user=francois.coppens@irsamc.ups-tlse.fr
#SBATCH --mail-type=ALL


#SBATCH -J "U*/U+"


module purge
module load intel/16.1.3
export OMP_NUM_THREADS=20
export MKL_NUM_THREADS=20
export OMP_PROC_BIND=true

export WORKDIR=$(pwd)

cd ${WORKDIR}
date
echo "Execution starting"
time srun -N 1 -n 1 -c 20 $(placement 1 20) ./loop.sh ../he-wfs/
date
echo "Execution finished"

jobinfo ${SLURM_JOBID}
