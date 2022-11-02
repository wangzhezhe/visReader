#!/bin/bash
prefix=advection_streamlinesOutput.step
totalStep=59
blockNum=8
for ((i=0;i<totalStep;i++)); do
    #create the catalog file
    listFilePrefix=${prefix}${i}
    listFileName=${listFilePrefix}.visit
    touch ${listFileName}
    echo ${listFileName}
    echo "!NBLOCKS ${blockNum}" > ${listFileName}
    blockName=${prefix}${i}.rank${j}.domain0.vtk
    ls ${listFilePrefix}.*.vtk >> ${listFileName}

done