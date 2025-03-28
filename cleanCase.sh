#!/bin/bash

#===================================================================#
# allrun script for testcase as part of test routine 
# run settlingTest CFD part
# Christoph Goniva - Feb. 2011
#===================================================================#

#- source CFDEM env vars
. ~/.bashrc

#--------------------------------------------------------------------------------#
#- define variables
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
#--------------------------------------------------------------------------------#



#- clean up case
echo "deleting data at: $casePath ?"
#read
rm -r $casePath/*.e*
rm -r $casePath/*.o*
rm -r $casePath/log*
rm -r $casePath/CFD/[0-9].[0-9]*
rm -r $casePath/CFD/[1-9]*
rm -r $casePath/CFD/log.*
rm -r $casePath/CFD/averageProps
rm -r $casePath/CFD/*.OpenFOAM
rm -r $casePath/CFD/*.restart*
rm -r $casePath/CFD/clockData
rm -r $casePath/CFD/*.dat
rm -r $casePath/CFD/octave/*.dat
rm -r $casePath/CFD/octave/*.eps
rm -r $casePath/CFD/octave/octave-workspace
rm -r $casePath/CFD/processor*
rm -r $casePath/CFD/VTK
rm -r $casePath/CFD/volAvU
rm -r $casePath/CFD/couplingFiles
rm -r $casePath/CFD/particles
rm -r $casePath/CFD/probe*
rm -r $casePath/CFD/postProcess*
rm -r $casePath/CFD/particleProbes
rm -r $casePath/CFD/*Restart*
rm -r $casePath/CFD/*restart*
rm -r $casePath/CFD/fluidFlowStats.res
#rm -r $casePath/CFD/constant/polyMesh/*.gz
#rm -r $casePath/CFD/constant/polyMesh/faces
#rm -r $casePath/CFD/constant/polyMesh/neighbour
#rm -r $casePath/CFD/constant/polyMesh/owner
#rm -r $casePath/CFD/constant/polyMesh/points
#rm -r $casePath/CFD/constant/polyMesh/boundary
rm -r $casePath/CFD/processor*
rm -r $casePath/DEM/postParticles/*
rm -r $casePath/DEM/postGlobal/*
rm -r $casePath/DEM/post/*.*
rm -r $casePath/DEM/post/tracer/*
#rm -r $casePath/DEM/post/restart/*.*
#rm -r $casePath/DEM/log.*
#rm -r $casePath/DEM/*Restart*
#rm -r $casePath/DEM/*restart*

rm -r $casePath/CFD/pascal/0.*
rm -r $casePath/CFD/pascal/1.*
rm -r $casePath/CFD/pascal/2.*
rm -r $casePath/CFD/codeInfo.*
rm -r $casePath/CFD/*.dat
rm -r $casePath/CFD/sequencedVTK*

rm -r $casePath/slurm*

rm -r $casePath/postprocessing/slurm*
rm -r $casePath/postprocessing/liggghts_dumps/*.liggghts_run
rm -r $casePath/postprocessing/liggghts_dumps/*.vtk
rm -r $casePath/postprocessing/liggghts_dumps/slurm*

rm -r $casePath/postprocessing/tracer_dumps/*.liggghts_run
rm -r $casePath/postprocessing/tracer_dumps/*.vtk
rm -r $casePath/postprocessing/tracer_dumps/slurm*
rm -r $casePath/DEM/post/bondnum/contactarea*

echo "done"


