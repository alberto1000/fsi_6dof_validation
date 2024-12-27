#!/bin/bash

# Go to background_zone directory and set up the background mesh
cd ../background_zone
./Allclean
blockMesh
#topoSet -dict system/topoSetDictR2 | tee log.topoSet
topoSet -dict system/topoSetDict | tee log.topoSet
mv log.topoSet log.topoSet1
#refineMesh -dict system/refineMeshDict2 -overwrite
cp -r 0.orig/. 0/
touch mesh_background.foam

# Go to component_zone and set up the component mesh
cd ../component_zone
./Allclean
blockMesh
surfaceFeatureExtract 
decomposePar
mpirun -np 12 snappyHexMesh -parallel -overwrite
reconstructParMesh -constant
topoSet
touch mesh_component.foam

# Merge component and background meshes and finalize setup
cd ../background_zone
mergeMeshes . ../component_zone -overwrite
checkMesh
topoSet -dict system/topoSetDict | tee log.topoSet
mv log.topoSet log.topoSet2
setFields
touch mesh_assembled.foam

# Decompose domain and run simulation
decomposePar
mpirun -np 12 overPimpleDyMFoam -parallel | tee log.simulation

# Reconstruct final mesh and open ParaView for visualization
reconstructPar
paraFoam
