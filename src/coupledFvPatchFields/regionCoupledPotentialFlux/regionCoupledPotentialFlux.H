/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Fork:     foam-extend
    \\  /    A nd           | Version:  4.1                                 
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
    Foam::regionCoupledPotentialFlux.C

SourceFiles
    regionCoupledPotentialFlux.C

Author
    Steffen Hess(s.hess@fz-juelich.de, Affiliation A) 
    All rights reserved.

Contact
    Steffen Hess(s.hess@fz-juelich.de) 

Affiliations
    Affiliation A)
    Forschungszentrum Jülich IEK-14

Acknowledgement
    Funded by
    AIDAS – AI, Data Analytics and Scalable Simulation

Description
    
    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#ifndef regionCoupledPotentialFlux_H
#define regionCoupledPotentialFlux_H

#include "genericRegionCoupledFluxFvPatchField.H"
#include "currentTransferInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class regionCoupledPotentialFlux Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledPotentialFlux
:
    public genericRegionCoupledFluxFvPatchField<scalar>
{

public:

    //- Runtime type information
    TypeName("regionCoupledPotentialFlux");

    // Constructors

        //- Construct from patch and internal field
        regionCoupledPotentialFlux
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupledPotentialFlux
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given regionCoupledPotentialFlux
        //  onto a new patch
        regionCoupledPotentialFlux
        (
            const regionCoupledPotentialFlux&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
            	new regionCoupledPotentialFlux(*this)
            );
        }

        //- Construct as copy setting internal field reference
        regionCoupledPotentialFlux
        (
            const regionCoupledPotentialFlux&,
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
                new regionCoupledPotentialFlux(*this, iF)
            );
        }

    // Member functions

        //- Calculate species flux jump
        virtual tmp<scalarField> fluxJump() const;

        const regionInterfaces::currentTransferInterface& currentTransferInterface() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
