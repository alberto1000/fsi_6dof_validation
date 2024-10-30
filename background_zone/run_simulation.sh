#!/bin/bash

# Go to background_zone directory and set up the background mesh
cd ../background_zone
blockMesh
topoSet -dict system/topoSetDictR2 | tee log.topoSet
mv log.topoSet log.topoSet1
refineMesh -dict system/refineMeshDict2 -overwrite
cp -r 0.origin/* 0/

# Go to component_zone and set up the component mesh
cd ../component_zone
blockMesh
surfaceFeatureExtract 
decomposePar
mpirun -np 32 snappyHexMesh -parallel -overwrite
reconstructParMesh -constant
topoSet

# Merge component and background meshes and finalize setup
cd ../background_zone
mergeMeshes . ../component_zone -overwrite
checkMesh
topoSet -dict system/topoSetDict | tee log.topoSet
mv log.topoSet log.topoSet2
setFields

# Decompose domain and run simulation
decomposePar
mpirun -np 32 overPimpleDyMFoam -parallel | tee log.simulation

# Reconstruct final mesh and open ParaView for visualization
reconstructParMesh -constant
paraFoam
