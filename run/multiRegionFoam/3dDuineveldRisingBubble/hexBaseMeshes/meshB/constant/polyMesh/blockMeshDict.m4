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
define(rIn, 1.0)
// domain length in x,y-direction as multiple of radius
define(fac, 8)
// number of cells in circumferential direction (must be divisible by 4!)
define(cellsCirc, 100)
// number of cells per bubble diameter
define(cellsDiam, 10)
// cell grading in inner domain (value smaller than 1)
define(gi, 0.4) //0.05


define(i1, calc(rIn/sqrt(3)/1.1)) // inner point to ensure nearly cartesian grid around origin
define(i2, calc(rIn/sqrt(3))) // point on bubble

define(arc1, calc(rIn/sqrt(2)/1.2)) // arc middle point between inner points; has to be a value between the i1-factor and (i1-factor*sqrt(3/2)) to produce curvature between 0 and circle 
define(arc2, calc(rIn/sqrt(2)))  // arc middle point between points on bubble

define(c1, calc(cellsCirc/4))
define(c2, calc(cellsDiam))

//convertToMeters 0.001; // in mm

geometry
{
    sphere
    {
        type searchableSphere;
        centre (0 0 0);
        radius rIn;
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
    hex (0 3 4 7 8 11 12 15) (c1 c1 c1) simpleGrading (1 1 1)  // inner block
    hex (3 2 5 4 11 10 13 12) (c2 c1 c1) simpleGrading (gi 1 1) // upper part
    hex (4 5 6 7 12 13 14 15) (c2 c1 c1) simpleGrading (gi 1 1) // equatorial part
    hex (7 6 1 0 15 14 9 8) (c2 c1 c1) simpleGrading (gi 1 1) // lower part
    hex (0 1 2 3 8 9 10 11) (c2 c1 c1) simpleGrading (gi 1 1) // equatorial part
    hex (0 1 6 7 3 2 5 4) (c2 c1 c1) simpleGrading (gi 1 1) // equatorial part
    hex (8 9 10 11 15 14 13 12) (c2 c1 c1) simpleGrading (gi 1 1) //equatorial part
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
    project (1 2 10 9) sphere
    project (2 5 13 10) sphere
    project (5 6 14 13) sphere
    project (6 1 9 14) sphere
    project (1 2 5 6) sphere
    project (9 14 13 10) sphere
);

patches
(
    wall interfaceShadow
    (
        (1 2 10 9)
        (2 5 13 10)
        (5 6 14 13)
        (6 1 9 14)
	(1 2 5 6)
	(9 14 13 10)
    )
);

mergePatchPairs
(
);

// ************************************************************************* //
