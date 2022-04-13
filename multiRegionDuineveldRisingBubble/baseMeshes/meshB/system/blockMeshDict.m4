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

// bubble radius mm
define(rIn, 2.0)
// pseudo height
define(z, 0.005) 
// domain length in x,y-direction as multiple of radius
define(fac, 8)
// number of cells in circumferential direction (must be divisible by 4!)
define(cellsCirc, 200)
// number of cells per bubble diameter
define(cellsDiam, 60)
// cell grading in inner domain (value smaller than 1)
define(gi, 0.4) //0.05

define(i1, calc(rIn/sqrt(2))) // point on bubble along the diagonal
define(i2, calc(rIn*0.5)) 
define(i3, calc(rIn*0.55)) 

define(e1, calc(rIn*0.866025404)) 
define(e2, calc(rIn*0.3)) 
define(e3, calc(rIn*0.65)) 

define(c1, calc(cellsCirc/4))
define(c2, calc(cellsDiam/4))
define(c3, calc(cellsCirc/2))

convertToMeters 0.001; // in mm

vertices
(
    (i3 -i3 0) //0
    (i1 -i1 0)
    (i1 i1 0)
    (i3 i3 0) //3
    (0 i1 0)
    (0 rIn 0)
    (0 -rIn 0)
    (0 -i1 0)  //7
    (i3 -i3 z)
    (i1 -i1 z)
    (i1 i1 z)
    (i3 i3 z)//11
    (0 i1 z)
    (0 rIn z)
    (0 -rIn z)
    (0 -i1 z)//15
);

blocks
(
    hex (0 3 4 7 8 11 12 15) (c1 c2 1) simpleGrading (1 1 1)
    hex (3 2 5 4 11 10 13 12) (c2 c2 1) simpleGrading (gi 1 1)
    hex (7 6 1 0 15 14 9 8) (c2 c2 1) simpleGrading (gi 1 1)
    hex (0 1 2 3 8 9 10 11) (c2 c1 1) simpleGrading (gi 1 1)
);

edges
(
    arc 2 5 (i2 e1 0)
    arc 6 1 (i2 -e1 0)
    arc 1 2 (rIn 0 0)
    arc 10 13 (i2 e1 z)
    arc 14 9 (i2 -e1 z)
    arc 9 10 (rIn 0 z)

    arc 3 4 (e2 e3 0)
    arc 7 0 (e2 -e3 0)
    arc 0 3 (i1 0 0)
    arc 11 12 (e2 e3 z)
    arc 15 8 (e2 -e3 z)
    arc 8 11 (i1 0 z)
);

patches
(
    empty lowerWedge
    (
        (0 7 4 3)
        (3 4 5 2)
        (0 3 2 1)
        (7 0 1 6)
    )
    empty upperWedge
    (
        (8 11 12 15)
        (11 10 13 12)
        (8 9 10 11)
        (15 14 9 8)
    )
    wall interfaceShadow
    (
        (1 2 10 9)
        (2 5 13 10)
        (6 1 9 14)
    )
);

mergePatchPairs
(
);

// ************************************************************************* //
