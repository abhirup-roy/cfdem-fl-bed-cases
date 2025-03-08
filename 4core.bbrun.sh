#!/bin/bash
#SBATCH --job-name=fl-bed-4core          # Job name to identify from queue
#SBATCH --mail-type=ALL              # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=16                 # Number of processes
#SBATCH --time=10-0                   # Time limit in days-hrs (can be hrs:min:sec)
#SBATCH --qos bbdefault
#SBATCH --mem=122G


# Load modules for any applications
module purge; 
module load bluebear; 
module load OpenFOAM/5.0-20180108-foss-2019a
# Change to the directory the job was submitted from
date;hostname;pwd

# Run program, using 'mpiexec' to start the job
# mpiexec automatically picks up the # of cores
# assigned to the job. No other flags are required
#  - note: don't use 'mpirun'
source $FOAM_BASH
source ~/CFDEM-JKR/CFDEM/setvars.sh
source	$CFDEM_bashrc

#- define variables
casePath=`pwd`
# check if mesh was built
if [ -f "$casePath/CFD/constant/polyMesh/points" ]; then
    echo "mesh was built before - using old mesh"
else
    echo "mesh needs to be built"
    cd $casePath/CFD
    blockMesh
fi

if [ -f "$casePath/DEM/post/restart/liggghts.restart" ];  then
    echo "LIGGGHTS init was run before - using existing restart file"
else
#--------------------------------------------------------------------------------#
#- define variables
#    casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
    logpath="$casePath"
    headerText="run_liggghts_init_DEM"
    logfileName="log_$headerText"
    solverName="in.liggghts_init"
    nrProcs="16"
    machineFileName="none"   # yourMachinefileName | none
#--------------------------------------------------------------------------------#

#- call function to run DEM case
    cd  $casePath/DEM

    mpiexec -n $nrProcs $CFDEM_LIGGGHTS_EXEC -in $solverName 
    #$CFDEM_LIGGGHTS_EXEC -in $solverName 
fi

#--------------------------------------------------------------------------------#
#- define variables
logpath=$casePath
headerText="run_parallel_cfdemSolverPisoSTM_packedBedTemp_CFDDEM"
logfileName="log_$headerText"
solverName="cfdemSolverPiso"
nrProcs=$SLURM_NTASKS
machineFileName="none"   # yourMachinefileName | none
debugMode="off"          # on | off| strict
testHarnessPath="$CFDEM_TEST_HARNESS_PATH"
runOctave="false"
postproc="false"
cleanup="false"
#--------------------------------------------------------------------------------#

#- call function to run a parallel CFD-DEM case
    rm $logpath/$logfileName
    #- change path
    cd $casePath/CFD
    #- remove old data
    rm -rf processor*
    rm -rf postProcessing

    #- decompose case
    decomposePar

    #- make proc dirs visible
    count=0
    for i in `seq $nrProcs`
    do
        let count=$i-1
        (cd $casePath/CFD/processor$count && touch file.foam)
    done

#- run applictaion
    cd $casePath/CFD
    mpiexec -n $nrProcs $solverName -parallel 
    # foamSequenceVTKFiles
    # reconstructPar 
    # rm -rf processor*
date
