#!/bin/bash

#----------------------------------------------------------------------#
# Contributor |   Shidong Zhang                                        #
# Email       |   s.zhang@fz-juelich.de                                #
#----------------------------------------------------------------------#
# Solver      |   openFuelCellFoam2                                    #
# OpenFOAM    |   OpenFOAM-v1906 or newer (ESI)                        #
#----------------------------------------------------------------------#
# Source code |   https://jugit.fz-juelich.de/s.zhang/fuelcellfoam     #
# Update from |   07.07.2021                                           #
#----------------------------------------------------------------------#

cp system/controlDict.mesh system/controlDict

rm -rf constant/*/polyMesh

rm constant/polyMesh/sets/*

setSet -batch config/make.zoneSet

rm constant/polyMesh/cellZones

rm -rf constant/polyMesh/tmp
mkdir constant/polyMesh/tmp

mv constant/polyMesh/sets/* constant/polyMesh/tmp/.
cp constant/polyMesh/tmp/air constant/polyMesh/sets/.
cp constant/polyMesh/tmp/fuel constant/polyMesh/sets/.
cp constant/polyMesh/tmp/interconnect constant/polyMesh/sets/.
cp constant/polyMesh/tmp/electrolyte constant/polyMesh/sets/.

## split

setsToZones -noFlipMap
splitMeshRegions -cellZonesOnly
cp -rf 1/. constant/.
rm -rf 1

setSet -batch config/make.setAir -region air
setSet -batch config/make.setFuel -region fuel

setsToZones -noFlipMap -region air
setsToZones -noFlipMap -region fuel

# Create faMeshes
cp -r system/air/faMesh_definition constant/air/faMesh
cp -r system/fuel/faMesh_definition constant/fuel/faMesh
cp -r system/interconnect/faMesh_definition constant/interconnect/faMesh
cp -r system/electrolyte/faMesh_definition constant/electrolyte/faMesh

makeFaMesh -region air
makeFaMesh -region fuel
makeFaMesh -region interconnect
makeFaMesh -region electrolyte


## creating protonic and electronic field
rm constant/polyMesh/cellZones

rm constant/polyMesh/sets/*

# move the other regions, keep phiEC, phiEA, phiI0, phi0, phiI0
cp constant/polyMesh/tmp/phiEA constant/polyMesh/sets/.
cp constant/polyMesh/tmp/phiEC constant/polyMesh/sets/.
cp constant/polyMesh/tmp/phi0 constant/polyMesh/sets/.
cp constant/polyMesh/tmp/phiI0 constant/polyMesh/sets/.

# clean phiEA and phiEC mesh

setsToZones -noFlipMap
splitMeshRegions -cellZonesOnly
cp -rf 1/. constant/.
rm -rf 1

setSet -batch config/make.setPhiEC -region phiEC
setsToZones -noFlipMap -region phiEC

setSet -batch config/make.setPhiEA -region phiEA
setsToZones -noFlipMap -region phiEA

rm -rf constant/phiI0
rm -rf constant/phi0

rm -rf constant/polyMesh/sets/*

# Create faMeshes
cp -r system/phiEA/faMesh_definition constant/phiEA/faMesh
cp -r system/phiEC/faMesh_definition constant/phiEC/faMesh

makeFaMesh -region phiEA
makeFaMesh -region phiEC

# creating phiI

cp constant/polyMesh/tmp/phiI constant/polyMesh/sets/.
cp constant/polyMesh/tmp/phiE0 constant/polyMesh/sets/.
cp constant/polyMesh/tmp/phiE1 constant/polyMesh/sets/.

rm constant/polyMesh/cellZones

setsToZones -noFlipMap
splitMeshRegions -cellZonesOnly
cp -rf 1/. constant/.
rm -rf 1

setSet -batch config/make.setPhiI -region phiI
setsToZones -noFlipMap -region phiI

rm -rf constant/phiE1
rm -rf constant/phiE0

cp constant/polyMesh/tmp/* constant/polyMesh/sets/.
rm -rf constant/polyMesh/tmp

setsToZones -noFlipMap

rm -rf system/phi0
rm -rf system/phiI0
rm -rf system/phiE0
rm -rf system/phiE1

# Create faMeshes
cp -r system/phiI/faMesh_definition constant/phiI/faMesh
makeFaMesh -region phiI

cp system/controlDict.run system/controlDict
