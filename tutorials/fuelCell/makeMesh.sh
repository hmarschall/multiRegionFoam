#!/bin/bash

source parameter.dat

rm -r 0
cp -r 0.orig 0

mkdir constant/polyMesh
cp -r system/blockMeshDict constant/polyMesh/.

cp -rf ./system/controlDict.mesh ./system/controlDict

# =========== Create mesh ================== #
blockMesh
#construct regions and cellZones
./preprocessing.sh | tee log.pre
cp -rf ./system/controlDict.run ./system/controlDict
