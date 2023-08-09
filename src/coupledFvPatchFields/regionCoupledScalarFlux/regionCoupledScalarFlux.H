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
    regionCoupledScalarFlux

SourceFiles
    regionCoupledScalarFlux.C

Author
    Constantin Habes (constantin.habes@tu-darmstadt.de, Affiliation A)
    Holger Marschall (holger.marschall@tu-darmstadt.de, Affiliation A)

Contact
    Holger Marschall (holger.marschall@tu-darmstadt.de)
    main developer and principal investigator
    TU Darmstadt

Affiliations
    Affiliation A)
    Computational Multiphase Flow
    Department of Mathematics
    Technical University of Darmstadt, Germany

Acknowledgement
    Funded by
    Hessian Ministry of Higher Education, Research, Science and the Arts,
    National High Performance Computing Center for Computational Engineering
    Science (NHR4CES)

Description
    Neumann boundary condition for species within the ALE interface tracking.
    Normal species flux jump is calculated from the neighbour patch flux.

    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#ifndef regionCoupledScalarFlux_H
#define regionCoupledScalarFlux_H

#include "genericRegionCoupledFluxFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class regionCoupledScalarFlux Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledScalarFlux
:
    public genericRegionCoupledFluxFvPatchField<scalar>
{

public:

    //- Runtime type information
    TypeName("regionCoupledScalarFlux");

    // Constructors

        //- Construct from patch and internal field
        regionCoupledScalarFlux
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupledScalarFlux
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given regionCoupledScalarFlux
        //  onto a new patch
        regionCoupledScalarFlux
        (
            const regionCoupledScalarFlux&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
            	new regionCoupledScalarFlux(*this)
            );
        }

        //- Construct as copy setting internal field reference
        regionCoupledScalarFlux
        (
            const regionCoupledScalarFlux&,
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
                new regionCoupledScalarFlux(*this, iF)
            );
        }

    // Member functions

        //- Calculate species flux jump
        virtual tmp<scalarField> fluxJump() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
