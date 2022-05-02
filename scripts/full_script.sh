#!/bin/bash

#SBATCH --job-name=
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G         # memory per cpu-core (4G per cpu-core is default
#SBATCH --time=23:59:59          # total run time limit (HH:MM:SS)


#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user= 


#Koen TO DO: CHECK OUTPUT FOLDER, CHECK INPUT FILES



#----------------------------------------- run parameters-------------------------------


jobname=
TotalTime=8001
timestep=1e-16
particles=
boxdim=1e6


#-----------------------------------------Simulation  parameters------------------------

#*** Electron-Electron interactions ***


farfield=false # If you use only near , or far field. TURN FMM OFF!!!!
nearfield=false
FMM=true
cutoffradius=0.01
order=2

# *** Electron-Phonon interaction ***

phonons=true


# *** Pulse ***

pulsesize=0.9
pulseenergy=0.1

# *** Doping ***

dopingpercentage=0
dcharge=1.0

holes=true # holes or ion doping. Do the positive charges have the same effective mass as the electrons


#----------------------------------------- Initial Conditions-------------------------------

#*** starting temperature ***
temp=300

#----------------------------------------- Material Parameters------------------------------


#*** eff. mass ***
eff_mass=10

#*** Dielectric ***
eps_infty=10
eps=30

#*** typical phonon frequency ***
wO=5e12


#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------




module purge
module load CMake/3.20.1-GCCcore-10.3.0 
module load OpenMPI/4.0.5-GCC-10.2.0 
module load ScaLAPACK/2.1.0-gompi-2020b 

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

cd ../../../data/p274598/Output/; mkdir $jobname
cd 
cd EnsembleMC/


if [ $nearfield = false ] ; then
    sed -i "s# algorithm->execute(FFmmNearField);#//  algorithm->execute(FFmmNearField);#" ./src/EnsembleMC.cpp
    
    echo "Only farfield"
fi

if [ $farfield = false ] ; then
    sed -i "s# algorithm->execute(FFmmFarField);#//  algorithm->execute(FFmmFarField);#" ./src/EnsembleMC.cpp
    echo "Only nearfield"
fi

if [ $nearfield = true ] ; then
    sed -i "s#.*//  algorithm->execute(FFmmNearField);.*#  algorithm->execute(FFmmNearField);#" ./src/EnsembleMC.cpp
    
fi


if [ $farfield = true ] ; then
    sed -i "s#.*//  algorithm->execute(FFmmFarField);.*#  algorithm->execute(FFmmFarField);#" ./src/EnsembleMC.cpp
fi



if [ $FMM = false ] ; then
    sed -i "s#algorithm->execute().*#// algorithm->execute();#" ./src/EnsembleMC.cpp
    echo "FMM is OFF"
fi

if [ $phonons = false ] ; then
    sed -i "s#.*Scatter(pa.*#        // Scatter(particles, mat_par, scat_par, idxPart);#" ./include/free_flight_scatter.hpp
    echo "No phonon scattering"
fi

if [ $FMM = true ] ; then
    sed -i "s#.*// algorithm->execute().*#    algorithm->execute();#" ./src/EnsembleMC.cpp
fi


if [ $phonons = true ] ; then
    sed -i "s#.*// Scatter(pa.*# Scatter(particles, mat_par, scat_par, idxPart);#" ./include/free_flight_scatter.hpp
fi




if [ $holes = true ] ; then
    sed -i "s#.*invmass /= 1e9; #// invmass /= 1e9; #" ./include/drift.hpp
fi



if [ $holes = false ] ; then
    sed -i "s#.*// invmass /= 1e9; #        invmass /= 1e9; #" ./include/drift.hpp
fi



sed -i "s/int totaltime.*/int totaltime = ${TotalTime};/" ./src/EnsembleMC.cpp  # set the run time

sed -i -e "/geometryClass/s/box_dim.*/box_dim{$boxdim}/" ./include/material_paramaters.hpp

sed -i "s/static double timestep.*/static double timestep = ${timestep};/" ./include/free_flight_scatter.hpp # set timestep

sed -i "s#H =.*#H = ${dopingpercentage}*N; #" ./src/EnsembleMC.cpp # set percentage of holes

sed -i "s/dopingcharge =.*/dopingcharge = ${dcharge};/" ./include/init_kspace.hpp 

sed -i "s/static constexpr unsigned ORDER.*/static constexpr unsigned ORDER = ${order};/" ./include/treeStructure.hpp

#sed -i "s/mat_par->set_eps(.*/mat_par->set_eps(${eps});/" ./include/polaroptical.hpp  #Setting of static dielectric constant

#sed -i "s/mat_par->set_eps_infty.*/mat_par->set_eps_infty(${eps_infty});/" ./include/polaroptical.hpp  #Setting of optic dielectric constant

sed -i "s/polarw0.*2*2*M_PI}/polarw0{${wO}*2*M_PI}/" ./include/material_paramaters.hpp  #Setting of typical phonon freq

sed -i -e "/mat_paramClass()/s/effmass.*}{}/effmass{${eff_mass}*9.11e-31, ${eff_mass}*9.11e-31, ${eff_mass}*9.11e-31}{}/" ./include/material_paramaters.hpp   #Setting of eff. mass

sed -i -e "/mat_paramClass()/s/eps{.* eps_/eps{$eps}, eps_/" ./include/material_paramaters.hpp   #Setting of static dielectric constant

sed -i -e "/mat_paramClass()/s/eps_inf.*, eff/eps_infty{$eps_infty}, eff/" ./include/material_paramaters.hpp  #Setting of optic dielectric constant

sed -i "s/double pulsesize = .*/double pulsesize = ${pulsesize};/" ./src/EnsembleMC.cpp

sed -i "s/pulse_energy =.*/pulse_energy = ${pulseenergy};/" ./include/init_kspace.hpp

sed -i "s#e_start =.*#e_start = - ( ${temp} * 1.38066e-23 ) / ( 1.60219e-19 * 1.5) * log(rr) ;#" ./include/init_kspace.hpp  #Setting of initial temperature. If you want to init via MB stats

#sed -i 's#ntf(file1.*tsim);#ntf(file1, "../../Output/'"${jobname}"'/%d_ps.csv", tsim);#' ./include/output.hpp

sed -i 's#nombre =.*#nombre = "../../../../data/p274598/Output/'"${jobname}"'/" + std::to_string(tsim) + ".fma";#' ./src/EnsembleMC.cpp

sed -i "s/FInterpMatrixKernelAPLUSR(const FReal inCoreWidth =.*/FInterpMatrixKernelAPLUSR(const FReal inCoreWidth = ${cutoffradius})/" ./include/scalfmm/include/Kernels/Interpolation/FInterpMatrixKernel.hpp





cd Build/; make
srun ./EnsembleMC -fin ../examples/unitCube/${particles}.fma  -fout ../Output/${jobname}.fma
#sed -i "s#../../Output/$jobname/#../../Output/#" ./include/output.hpp


