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
    genericInterfaceCoupledVelocityValue

Author
    Holger Marschall (marschall@csi.tu-darmstadt.de)

Description
    Dirichlet boundary condition for velocity within the ALE interface tracking.
    The no-jump condition of the velocity at the interface is taken into
    account by copying the value from the neighbour patch. Flux continuity is 
    enforced likewise. The momentum flux is calculated as the normal gradient 
    multiplied by mu.

    Example usage:
        interfacePatch
        {
            type                genericInterfaceCoupledVelocityValue;
            k                   "mu";
            relax               1.0;
            neighbourRegionName liquid;
            neighbourPatchName  interfacePatchShadow;
            nonOrthCorr         yes;
            secondOrder         yes;
            value               uniform (0 0 0);
        }


SourceFiles
    genericInterfaceCoupledVelocityValue.C

\*---------------------------------------------------------------------------*/

#ifndef genericInterfaceCoupledVelocityValue_H
#define genericInterfaceCoupledVelocityValue_H

#include "genericRegionCoupledJumpFvPatchField.H"
#include "capillaryInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class genericInterfaceCoupledVelocityValue Declaration
\*---------------------------------------------------------------------------*/

class genericInterfaceCoupledVelocityValue
:
    public genericRegionCoupledJumpFvPatchField<vector>
{

protected:

    //- Update phi 
    virtual void updatePhi();

public:

    //- Runtime type information
    TypeName("genericInterfaceCoupledVelocityValue");


    // Constructors

        //- Construct from patch and internal field
        genericInterfaceCoupledVelocityValue
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        genericInterfaceCoupledVelocityValue
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given genericInterfaceCoupledVelocityValue
        //  onto a new patch
        genericInterfaceCoupledVelocityValue
        (
            const genericInterfaceCoupledVelocityValue&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new genericInterfaceCoupledVelocityValue(*this)
            );
        }

        //- Construct as copy setting internal field reference
        genericInterfaceCoupledVelocityValue
        (
            const genericInterfaceCoupledVelocityValue&,
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
                new genericInterfaceCoupledVelocityValue(*this, iF)
            );
        }


    // Member functions

        //- Calculate velocity jump
        virtual tmp<vectorField> valueJump() const;

        const regionInterfaces::capillaryInterface& capInterface() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

