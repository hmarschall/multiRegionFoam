/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.5                                   |
|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     ;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1; // in mm

changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'use Math::Trig; print ($1)')]) dnl>

define(pi, 3.14159265)
define(phiBox, calc(0.25*pi))
define(thetaBox, calc(atan(sqrt(2))))
define(thetaArc, calc(calc(0.25*pi)))

// bubble radius
define(rIn, 1.0)
// domain length in x,y-direction as multiple of radius
define(fac, 8)
// number of cells in circumferential direction (must be divisible by 4!)
define(cellsCirc, 100)
// number of cells per bubble diameter
define(cellsDiam, 5)
// cell grading in outer domain (value bigger than 1)
define(go, 50)

define(rOut, calc(rIn*fac))

define(arcInX, calc(rIn*sin(thetaArc)))
define(arcInY, calc(rIn*cos(thetaArc)))
define(arcInZ, calc(rIn*sin(thetaArc)))

define(arcOutX, calc(rOut*sin(thetaArc)))
define(arcOutY, calc(rOut*cos(thetaArc)))
define(arcOutZ, calc(rOut*sin(thetaArc)))

define(x1, calc(rIn*sin(thetaBox)*sin(phiBox)))
define(y1, calc(rIn*cos(thetaBox)))
define(z1, calc(rIn*sin(thetaBox)*cos(phiBox)))

define(x2, calc(rOut*sin(thetaBox)*sin(phiBox)))
define(y2, calc(rOut*cos(thetaBox)))
define(z2, calc(rOut*sin(thetaBox)*cos(phiBox)))

define(c1, calc(cellsCirc/4))
define(c2, calc(fac*cellsDiam))


geometry
{
    sphereIn
    {
        type searchableSphere;
        centre (0 0 0);
        radius rIn;
    }

    sphereOut
    {
        type searchableSphere;
        centre (0 0 0);
        radius rOut;
    }
}

vertices
(
    (x1 -y1 -z1) //0
    (x2 -y2 -z2)
    (x2 y2 -z2)
    (x1 y1 -z1) //3
    (-x1 y1 -z1)
    (-x2 y2 -z2)
    (-x2 -y2 -z2)
    (-x1 -y1 -z1)  //7
    (x1 -y1 z1)
    (x2 -y2 z2)
    (x2 y2 z2)
    (x1 y1 z1)//11
    (-x1 y1 z1)
    (-x2 y2 z2)
    (-x2 -y2 z2)
    (-x1 -y1 z1)//15
);

blocks
(
    hex (3 2 5 4 11 10 13 12) (c2 c1 c1) simpleGrading (go 1 1) // upper part
    hex (4 5 6 7 12 13 14 15) (c2 c1 c1) simpleGrading (go 1 1) // equatorial part
    hex (7 6 1 0 15 14 9 8) (c2 c1 c1) simpleGrading (go 1 1) // lower part
    hex (0 1 2 3 8 9 10 11) (c2 c1 c1) simpleGrading (go 1 1) // equatorial part
    hex (0 1 6 7 3 2 5 4) (c2 c1 c1) simpleGrading (go 1 1) // equatorial part
    hex (8 9 10 11 15 14 13 12) (c2 c1 c1) simpleGrading (go 1 1) //equatorial part
);

edges
(
    arc 2 5 (0 arcOutY -arcOutZ) 
    arc 5 6 (-arcOutX 0 -arcOutZ)
    arc 5 13 (-arcOutX arcOutY 0)
    
    arc 6 1 (0 -arcOutY -arcOutZ)
    arc 6 14 (-arcOutX -arcOutY 0)

    arc 1 2 (arcOutX 0 -arcOutZ)
    arc 1 9 (arcOutX -arcOutY 0)

    arc 2 10 (arcOutX arcOutY 0)

    arc 10 13 ( 0 arcOutY arcOutZ)
    arc 13 14 (-arcOutX 0 arcOutZ)

    arc 14 9 (0 -arcOutY arcOutZ)
    arc 9 10 (arcOutX 0 arcOutZ)

    arc 3 4 (0 arcInY -arcInZ)
    arc 4 7 (-arcInX 0 -arcInZ)
    arc 7 0 (0 -arcInY -arcInZ)
    arc 0 3 (arcInX 0 -arcInZ)
    arc 11 12 ( 0 arcInY arcInZ)
    arc 12 15 (-arcInX 0 arcInZ)
    arc 15 8 (0 -arcInY arcInZ)
    arc 8 11 (arcInX 0 arcInZ)
    arc 3 11 (arcInX arcInY 0)
    arc 4 12 (-arcInX arcInY 0)
    arc 7 15 (-arcInX -arcInY 0) 
    arc 0 8  (arcInX -arcInY 0)
);

faces
(

    project (0 8 11 3) sphereIn
    project (3 11 12 4) sphereIn
    project (4 12 15 7) sphereIn
    project (7 15 8 0) sphereIn
    project (0 7 4 3) sphereIn
    project (8 11 12 15) sphereIn

    project (1 2 10 9) sphereOut 
    project (2 5 13 10) sphereOut
    project (5 6 14 13) sphereOut
    project (6 1 9 14) sphereOut
    project (1 2 5 6) sphereOut
    project (9 14 13 10) sphereOut
);

patches
(
    wall space
    (
        (1 2 10 9)
	(2 5 13 10)
	(5 6 14 13)
	(6 1 9 14)
	(1 2 5 6)
	(9 14 13 10)
    )
    wall interface
    (
	(0 8 11 3)
	(3 11 12 4)
	(4 12 15 7)
	(7 15 8 0)
	(0 7 4 3)
	(8 11 12 15)
    )
);

mergePatchPairs
(
);

// ************************************************************************* //
