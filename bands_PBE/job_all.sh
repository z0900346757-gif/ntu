#!/bin/sh
#SBATCH --job-name=ZnO
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --nice=15
#SBATCH --output=/home/Unrequitedlove/DFT/Material/ZnO/bands_PBE/debug.o
#SBATCH --error=/home/Unrequitedlove/DFT/Material/ZnO/bands_PBE/debug.e

cd /home/Unrequitedlove/DFT/Material/ZnO/bands_PBE/

export OMP_NUM_THREADS=1
export OMPI_MCA_btl="self,vader"
export OMPI_MCA_btl_vader_single_copy_mechanism=none
export OMPI_MCA_btl_vader_eager_limit=32768


/home/Unrequitedlove/DFT//QE/mpi/bin/mpirun -n 64 /home/Unrequitedlove/DFT/QE/QE7.2/bin/pw.x < scf.in > scf.out 
sleep 10
/home/Unrequitedlove/DFT//QE/mpi/bin/mpirun -n 64 /home/Unrequitedlove/DFT/QE/QE7.2/bin/pw.x < bands.in > bands.out
sleep 10
/home/Unrequitedlove/DFT//QE/mpi/bin/mpirun -n 16 /home/Unrequitedlove/DFT/QE/QE7.2/bin/bands.x < bands_pp.in > bands_pp.out
wait

