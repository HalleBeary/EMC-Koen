#!/bin/bash
#SBATCH --job-name=2ps_1e-18_2000particles

#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8


#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=t.faber@rug.nl


module load CMake/3.20.1-GCCcore-10.3.0 
module load OpenMPI/4.1.1-GCC-10.3.0 
module load ScaLAPACK/2.1.0-gompi-2020b 


make
./EnsembleMC -fin ../examples/unitCube/2000.fma -fout 2ps-2000.fma

#cd ../EnsembleMC/Build/; make
#cd ../EnsembleMC/Build/; ./EnsembleMC -fin ../examples/unitCube/2000.fma -fout 2ps-2000.fma

