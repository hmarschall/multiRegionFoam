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
define(rIn, 2.0)
// domain length in x,y-direction as multiple of radius
define(fac, 8)
//pseudo height
define(z,0.005)
// number of cells in circumferential direction (must be divisible by 4!)
define(cellsCirc, 200)
// number of cells per bubble diameter (must be divisible by 4!)
define(cellsDiam, 30)
// cell grading in outer domain (value bigger than 1)
define(go, 50)

define(rOut, calc(rIn*fac))

define(i1, calc(rIn/sqrt(2))) 
define(i2, calc(rOut/sqrt(2)))

define(e1, calc(rIn*0.5)) 
define(e2, calc(rIn*0.866025404))
define(e3, calc(rOut*0.5)) 
define(e4, calc(rOut*0.866025404))

define(c1, calc(cellsCirc/4))
define(c2, calc(fac*cellsDiam/16))
define(c3, calc(cellsCirc/4))

convertToMeters 0.001; // in mm

vertices
(
    (i1 -i1 0)
    (i1 i1 0)
    (0 rIn 0)
    (0 -rIn 0)
    (i1 -i1 z)
    (i1 i1 z)
    (0 rIn z)
    (0 -rIn z)
    (i2 i2 0) //8
    (0 rOut 0)
    (0 -rOut 0)
    (i2 -i2 0)
    (i2 i2 z)
    (0 rOut z)
    (0 -rOut z)
    (i2 -i2 z)
);

blocks
(
   hex (1 8 9 2 5 12 13 6) (c1 c2 1) simpleGrading (go 1 1)
   hex (3 10 11 0 7 14 15 4) (c1 c2 1) simpleGrading (go 1 1)
   hex (0 11 8 1 4 15 12 5) (c1 c3 1) simpleGrading (go 1 1)
);

edges
(
    arc 1 2 (e1 e2 0)
    arc 3 0 (e1 -e2 0)
    arc 0 1 (rIn 0 0)
    arc 5 6 (e1 e2 z)
    arc 7 4 (e1 -e2 z)
    arc 4 5 (rIn 0 z)

// add arcs for outer domain
    arc 8 9 (e3 e4 0)
    arc 12 13 (e3 e4 z)
    arc 10 11 (e3 -e4 0)
    arc 14 15 (e3 -e4 z)
    arc 11 8 (rOut 0 0)
    arc 15 12 (rOut 0 z)
);

patches
(
    empty lowerWedge
    (
        (1 2 9 8)
        (3 0 11 10)
        (0 1 8 11)
    )
    empty upperWedge
    (
	(5 12 13 6)
        (7 14 15 4)
        (4 15 12 5)
    )
    wall space
    (
        (11 8 12 15)
        (8 9 13 12)
        (10 11 15 14)
    )
    wall interface
    (
        (0 4 5 1)
        (1 5 6 2)
        (3 7 4 0)
    )
);

mergePatchPairs
(
);

// ************************************************************************* //
