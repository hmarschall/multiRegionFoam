/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Global
    CourantNo

Description
    Calculates and outputs the mean and maximum Courant Numbers.

\*---------------------------------------------------------------------------*/

{
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
    scalar velMag = 0.0;

    const Time& runTime = mesh().time();

    const scalar deltaT = runTime.deltaT().value();

    const fvMesh& mesh = this->mesh();

    if (mesh.nInternalFaces())
    {
        surfaceScalarField magPhi = mag(phi_());

        surfaceScalarField SfUfbyDelta =
            mesh.surfaceInterpolation::deltaCoeffs()*magPhi;

        CoNum = max(SfUfbyDelta/mesh.magSf()).value()*deltaT;

        meanCoNum = (sum(SfUfbyDelta)/sum(mesh.magSf())).value()*deltaT;

        velMag = max(magPhi/mesh.magSf()).value();
    }

    #include "meshCourantNo.H"

    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum
        << " velocity magnitude: " << velMag
        << nl << endl;
}

// ************************************************************************* //
