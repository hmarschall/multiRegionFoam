/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\    /   O peration     | Version:     3.1                                |
|   \\  /    A nd           | Web:         http://www.extend-project.de       |
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
define(calc, [esyscmd(perl -e 'print ($1)')])

//dummy z direction
define(zmin, 0) 
define(zmax, 0.00001) 

//inlet Outlet L
define(xmin, 0) 
define(xmax, 10.0) 

//d and H
define(ymin, 0.0) 
define(ymax, 10.0) 

//N of cells
define(NX, 50)
define(NY, 25)
define(NZ, 1)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

vertices
(
    (xmin ymin zmin)
    (xmax ymin zmin)
    (xmax ymax zmin)
    (xmin ymax zmin)
    (xmin ymin zmax)
    (xmax ymin zmax)
    (xmax ymax zmax)
    (xmin ymax zmax)
);

blocks
(
    hex (0 1 2 3 4 5 6 7) (NX NY NZ) simpleGrading (1 1 1)
);

edges
(
);

patches
(
    empty backFront
    (      
        (4 5 6 7)
        (0 3 2 1)      
    )
    wall inlet
    (
        (0 4 7 3)
    )
    wall outlet
    (
        (2 6 5 1)
    )
    wall Wall
    (
        (0 1 5 4)
    )
    wall interface
    (
        (3 7 6 2)
    )
);

mergePatchPairs
(
);

// ************************************************************************* //
