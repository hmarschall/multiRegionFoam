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
    regionCoupledTraction

Author
    Constantin Habes (constantin.habes@tu-darmstadt.de)

Description
    Neumann boundary condition for displacement within a FSI problem.
    The formulation uses the St.Venant-Kirchhoff constitutive law.
    The boundary condition computes the traction from the interacting
    fluid flow on the solid side of the FSI interface.

    Example usage:

	    interfacePatchShadow
        {
            type                    regionCoupledTraction;
            neighbourRegionName     fluid;
            neighbourPatchName      interfacePatch;
            nonOrthCorr             yes;
            secondOrder             yes;
            value                   uniform (0 0 0);
        }


SourceFiles
    regionCoupledTraction.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupledTraction_H
#define regionCoupledTraction_H

#include "genericRegionCoupledFluxFvPatchField.H"
#include "fsiInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class regionCoupledTraction Declaration
\*---------------------------------------------------------------------------*/

class regionCoupledTraction
:
    public genericRegionCoupledFluxFvPatchField<vector>
{

protected:
    //- Switch to use directly computed flux instead of flux jump
    virtual bool useDirectFlux() const
    {
        return true;
    }

public:

    //- Runtime type information
    TypeName("regionCoupledTraction");


    // Constructors

        //- Construct from patch and internal field
        regionCoupledTraction
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        regionCoupledTraction
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given regionCoupledTraction
        //  onto a new patch
        regionCoupledTraction
        (
            const regionCoupledTraction&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
            	new regionCoupledTraction(*this)
            );
        }

        //- Construct as copy setting internal field reference
        regionCoupledTraction
        (
            const regionCoupledTraction&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new regionCoupledTraction(*this, iF)
            );
        }

    // Member functions

        //- Calculate the forces
        virtual tmp<vectorField> flux() const;

        const regionInterfaces::fsiInterface& fsiInterface() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
