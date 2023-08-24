/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 2022 AUTHOR,AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Implementation of time derivative solver for the "solidStVKPhysicsSolver" 
    momentum equation.

\*---------------------------------------------------------------------------*/

// ddt scheme check :
word ddtScheme
(
    mesh().schemesDict().ddtScheme("ddt(" + D_().name() +')')
);


// ddt scheme coefficient initialization :

dimensionedScalar Cn("Cn", dimless/dimTime, 0.0);
dimensionedScalar Co("Co", dimless/dimTime, 0.0);
dimensionedScalar Coo("Coo", dimless/dimTime, 0.0);

// Euler :
if ( ddtScheme == fv::EulerDdtScheme<vector>::typeName )
{
    Cn.value() = 1.0/mesh().time().deltaT().value();
    Co.value() = 1.0/mesh().time().deltaT().value();
}

// backward :
if (ddtScheme == fv::backwardDdtScheme<vector>::typeName)
{
    scalar deltaT = mesh().time().deltaT().value();
    scalar deltaT0 = mesh().time().deltaT0().value();

    Cn.value() = 1 + deltaT/(deltaT + deltaT0);
    Coo.value() = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    Co.value() = Cn.value() + Coo.value();

    if(mesh().time().timeIndex() == 1)
    {
        Cn.value() = 1.0;
        Co.value() = 1.0;
        Coo.value() = 0.0;
    }

    Cn.value() /= deltaT;
    Co.value() /= deltaT;
    Coo.value() /= deltaT;
}

// steadyState :
if (ddtScheme == fv::steadyStateDdtScheme<vector>::typeName)
{
    Cn.value() = 0.0;
    Co.value() = 0.0;
    Coo.value() = 0.0;
}
