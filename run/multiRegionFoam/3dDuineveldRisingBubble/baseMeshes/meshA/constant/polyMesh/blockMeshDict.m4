/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.5                                   |
|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'printf ($1)')])

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

define(arc2, calc(rOut/sqrt(2)))
define(arc1, calc(rIn/sqrt(2)))

define(i1, calc(rIn/sqrt(3))) 
define(i2, calc(rOut/sqrt(3)))

define(c1, calc(cellsCirc/4))
define(c2, calc(fac*cellsDiam))

//convertToMeters 0.001; // in mm

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
    (i1 -i1 -i1) //0
    (i2 -i2 -i2)
    (i2 i2 -i2)
    (i1 i1 -i1) //3
    (-i1 i1 -i1)
    (-i2 i2 -i2)
    (-i2 -i2 -i2)
    (-i1 -i1 -i1)  //7
    (i1 -i1 i1)
    (i2 -i2 i2)
    (i2 i2 i2)
    (i1 i1 i1)//11
    (-i1 i1 i1)
    (-i2 i2 i2)
    (-i2 -i2 i2)
    (-i1 -i1 i1)//15
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
    arc 2 5 (0 arc2 -arc2)
    arc 5 6 (-arc2 0 -arc2)
    arc 5 13 (-arc2 arc2 0)
    
    arc 6 1 (0 -arc2 -arc2)
    arc 6 14 (-arc2 -arc2 0)

    arc 1 2 (arc2 0 -arc2)
    arc 1 9 (arc2 -arc2 0)

    arc 2 10 (arc2 arc2 0)

    arc 10 13 ( 0 arc2 arc2)
    arc 13 14 (-arc2 0 arc2)

    arc 14 9 (0 -arc2 arc2)
    arc 9 10 (arc2 0 arc2)

    arc 3 4 (0 arc1 -arc1)
    arc 4 7 (-arc1 0 -arc1)
    arc 7 0 (0 -arc1 -arc1)
    arc 0 3 (arc1 0 -arc1)
    arc 11 12 ( 0 arc1 arc1)
    arc 12 15 (-arc1 0 arc1)
    arc 15 8 (0 -arc1 arc1)
    arc 8 11 (arc1 0 arc1)
    arc 3 11 (arc1 arc1 0)
    arc 4 12 (-arc1 arc1 0)
    arc 7 15 (-arc1 -arc1 0)
    arc 0 8  (arc1 -arc1 0)
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
