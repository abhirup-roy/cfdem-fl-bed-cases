#!/bin/bash
. ~/.bashrc
source $CFDEM_SRC_DIR/lagrangian/cfdemParticle/etc/functions.sh
casePath="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
 #- keep terminal open (if started in new terminal)
    #- get VTK data from liggghts dump file
   cd $casePath/DEM/post
   python -i $CFDEM_LPP_DIR/lpp.py dump*.liggghts_run

    #- get VTK data from CFD sim
    cd $casePath/CFD
    reconstructPar
    ###foamToVTK                                                   #- serial run of foamToVTK
    
    #- start paraview
    paraview foam.foam&

    #- keep terminal open (if started in new terminal)
