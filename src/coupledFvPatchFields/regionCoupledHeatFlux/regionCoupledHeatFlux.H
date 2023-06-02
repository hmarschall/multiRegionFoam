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
    regionCoupledHeatFlux

Author
    Constantin Habes (constantin.habes@stud.tu-darmstadt.de)

Description
    Neumann boundary condition for species within the ALE interface tracking.
    Normal species flux jump is calculated from the neighbour patch flux.

    Example usage:
        interfacePatch
            {
                type                    regionCoupledHeatFlux;
                interfaceType           heatTransferInterface;
                k                       "D";
                neighbourRegionName     liquid;
                neighbourPatchName      interfacePatchShadow;
                value                   uniform 0;
            }


SourceFiles
    regionCoupledHeatFlux.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledHeatFlux_H
#define regionCoupledHeatFlux_H

#include "genericRegionCoupledFluxFvPatchField.H"
#include "heatTransferInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class regionCoupledHeatFlux Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledHeatFlux
:
    public genericRegionCoupledFluxFvPatchField<scalar>
{

public:

    //- Runtime type information
    TypeName("regionCoupledHeatFlux");

    // Constructors

        //- Construct from patch and internal field
        regionCoupledHeatFlux
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupledHeatFlux
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given regionCoupledHeatFlux
        //  onto a new patch
        regionCoupledHeatFlux
        (
            const regionCoupledHeatFlux&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
            	new regionCoupledHeatFlux(*this)
            );
        }

        //- Construct as copy setting internal field reference
        regionCoupledHeatFlux
        (
            const regionCoupledHeatFlux&,
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
                new regionCoupledHeatFlux(*this, iF)
            );
        }


    // Member functions

        //- Calculate species flux jump
        virtual tmp<scalarField> fluxJump() const;

        const regionInterfaces::heatTransferInterface& heatTransInterface() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
