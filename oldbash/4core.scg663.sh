#!/bin/bash
#SBATCH --job-name=ad3e4 #Job name
#SBATCH --mail-type=FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --nodes=1                    # Run all processes on a single node
#SBATCH --ntasks=16                   # Number of processes
#SBATCH --mem=64gb                    # Total memory limit
#SBATCH --time=100:00:00              # Time limit hrs:min:sec
#SBATCH --partition=amd_512

# Load modules for any applications
# Change to the directory the job was submitted from
source /public3/home/scg6635/OpenFOAM-5.x/etc/rebashrc
# Run program, using 'mpiexec' to start the job
# mpiexec automatically picks up the # of cores
# assigned to the job. No other flags are required
#  - note: don't use 'mpirun'
export CFDEM_VERSION=public
export CFDEM_PROJECT_DIR=$HOME/CFDEM-$CFDEM_VERSION/CFDEMcoupling-PUBLIC-$WM_PROJECT_VERSION
export CFDEM_PROJECT_USER_DIR=$HOME/CFDEM-$CFDEM_VERSION/$LOGNAME-$CFDEM_VERSION-$WM_PROJECT_VERSION
export CFDEM_bashrc=$CFDEM_PROJECT_DIR/src/lagrangian/cfdemParticle/etc/bashrc
export CFDEM_LIGGGHTS_SRC_DIR=$HOME/CFDEM-$CFDEM_VERSION/LIGGGHTS/LIGGGHTS-PUBLIC/src
export CFDEM_LIGGGHTS_MAKEFILE_NAME=auto
export CFDEM_LPP_DIR=$HOME/LIGGGHTS/lpp/src
. $CFDEM_bashrc

#source Allrun.sh >allrunoutput
#./cleanCase.sh
#- define variables
#casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
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
    #- run DEM init
    #- $casePath/parDEMrun.sh
    #- source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh

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
   # parDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName
    cd  $casePath/DEM
    
    #mpiexec -n $nrProcs $CFDEM_LIGGGHTS_EXEC -in $solverName 
    $CFDEM_LIGGGHTS_EXEC -in $solverName 
fi

#- run parallel CFD-DEM in new terminal
#. $casePath/parCFDDEMrun.sh
#--------------------------------------------------------------------------------#
#- define variables
logpath=$casePath
headerText="run_parallel_cfdemSolverPisoSTM_packedBedTemp_CFDDEM"
logfileName="log_$headerText"
#solverName="cfdemSolverIBSemi"
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
#- parCFDDEMrun $logpath $logfileName $casePath $headerText $solverName $nrProcs $machineFileName $debugMode
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

    cd $casePath/CFD
    #- run applictaion
    #mpiexec -n $nrProcs $solverName -parallel 2>&1 | tee -a $logpath/$logfileName
    mpiexec -n $nrProcs $solverName -parallel 
 
date
