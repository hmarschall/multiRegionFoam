/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    regionCoupledPressureFlux

Author
    Holger Marschall (marschall@csi.tu-darmstadt.de)

Description
    Neumann boundary condition for pressure within the ALE interface tracking.
    Normal pressure flux jump is calculated from the neighbour patch flux.

    Example usage:
        interfacePatch
            {
                type                    regionCoupledPressureFlux;
                k                       "v";
                neighbourRegionName     liquid;
                neighbourPatchName      interfacePatchShadow;
                value                   uniform 0;
            }


SourceFiles
    regionCoupledPressureFlux.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledPressureFlux_H
#define regionCoupledPressureFlux_H

#include "genericRegionCoupledFluxFvPatchField.H"
#include "capillaryInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class regionCoupledPressureFlux Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledPressureFlux
:
    public genericRegionCoupledFluxFvPatchField<scalar>
{

public:

    //- Runtime type information
    TypeName("regionCoupledPressureFlux");

    // Constructors

        //- Construct from patch and internal field
        regionCoupledPressureFlux
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupledPressureFlux
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given regionCoupledPressureFlux
        //  onto a new patch
        regionCoupledPressureFlux
        (
            const regionCoupledPressureFlux&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
            	new regionCoupledPressureFlux(*this)
            );
        }

        //- Construct as copy setting internal field reference
        regionCoupledPressureFlux
        (
            const regionCoupledPressureFlux&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new regionCoupledPressureFlux(*this, iF)
            );
        }


    // Member functions

        //- Calculate pressure flux jump
        virtual tmp<scalarField> fluxJump() const;

        const regionInterfaces::capillaryInterface& capInterface() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
