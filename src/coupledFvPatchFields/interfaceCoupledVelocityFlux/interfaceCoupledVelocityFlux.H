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
    interfaceCoupledVelocityFlux

Author
    Holger Marschall (marschall@csi.tu-darmstadt.de)

Description
	Couple boundary condition for interfacial species flux within the
	ALE interface tracking framework.
  Jump in velocity normal flux is implemented using the normal velocity
  gradient jump from tangential momentum jump condition.

    Example usage:
		interfacePatch
        {
            type                    interfaceCoupledVelocityFlux;
            neighbourPatchName      neighbourPatch;
            k                       "mu";
            value                   uniform (0 0 0);
        }

SourceFiles
    interfaceCoupledVelocityFlux.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceCoupledVelocityFlux_H
#define interfaceCoupledVelocityFlux_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedGradientCorrectedFvPatchFields.H"
#include "patchCoupleManager.H"
#include "directMappedFvPatch.H"//
#include "interfacialTransport.H"

#include "surfaceTracking.H"
#include "meshInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class interfaceCoupledVelocityFlux Declaration
\*---------------------------------------------------------------------------*/

class interfaceCoupledVelocityFlux
:
    public fixedGradientFvPatchVectorField
{
    // Private data

        //- Name of jump coefficient
		//word K_;

        //- Diffusion coefficient
		word k_;

        //- fluid properties
        dimensionedScalar muFluidA_;

        dimensionedScalar muFluidB_;

        dimensionedScalar sigma_;

        mutable areaScalarField* sigmaPtr_;

        //- Name of flux field
        word phiName_;

        //- Name of density field
        word rhoName_;

        //- NonOrthogonal correction
        Switch nonOrthCorr_;

        //- Second order accuracy correction
        Switch secondOrder_;

        //- Current time index
        label curTimeIndex_;

        // Old face centres of patch
        // for (second order) time integration
        vectorField Fc_;
        vectorField oldFc_;
        vectorField oldoldFc_;

    // Private functions

        //- Reference to meshInterface
        const meshInterface& ale() const;

        //- Reference to interfacialTransport
        const interfacialTransport& interface() const;

        //- Calc velocity jump
        tmp<vectorField> velJump() const;

        //- Couple manager
        autoPtr<patchCoupleManager> coupleManagerPtr_;

        //- Interfacial transport
        mutable autoPtr<interfacialTransport> interfacePtr_;

public:

    //- Runtime type information
    TypeName("interfaceCoupledVelocityFlux");


    // Constructors

        //- Construct from patch and internal field
        interfaceCoupledVelocityFlux
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        interfaceCoupledVelocityFlux
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given interfaceCoupledVelocityFlux
        //  onto a new patch
        interfaceCoupledVelocityFlux
        (
            const interfaceCoupledVelocityFlux&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
            	new interfaceCoupledVelocityFlux(*this)
            );
        }

        //- Construct as copy setting internal field reference
        interfaceCoupledVelocityFlux
        (
            const interfaceCoupledVelocityFlux&,
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
                new interfaceCoupledVelocityFlux(*this, iF)
            );
        }

    // Destructor

        virtual ~interfaceCoupledVelocityFlux();

    // Member functions

        //- Update the patch field coefficients
        virtual void updateCoeffs();

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::blocking
        );

        //- Return the maximum normalized coupled patch residual
        scalarField residual() const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
