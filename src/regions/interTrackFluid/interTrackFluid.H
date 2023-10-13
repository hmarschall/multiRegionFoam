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
    Foam::interTrackFluid.C

SourceFiles
    interTrackFluid.C

Author
    Constantin Habes (constantin.habes@tu-darmstadt.de, Affiliation A)

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

    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/

#ifndef interTrackFluid_H
#define interTrackFluid_H

#include "regionType.H"
#include "dimensionedScalar.H"
#include "multiRegionSystem.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"
#include "pimpleControl.H"
#include "solutionControl.H"
#include "inletOutletFvPatchFields.H"
#include "genericRegionCoupledFluxFvPatchField.H"
#include "genericRegionCoupledJumpFvPatchField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                           Class interTrackFluid Declaration
\*---------------------------------------------------------------------------*/

class interTrackFluid
:
    public regionType
{
        //- name
        word regionName_;

        //- transport properties dictionary
        IOdictionary transportProperties_;

        //- PIMPLE control
        pimpleControl pimple_;

        //- density
        autoPtr<volScalarField> rho_;

        //- viscosity
        autoPtr<volScalarField> mu_;

        //- transport
        autoPtr<volVectorField> U_;

        autoPtr<surfaceScalarField> phi_;

        autoPtr<volScalarField> p_;

        //- pcorr
        autoPtr<volScalarField> pcorr_;

        wordList pcorrTypes_;

        //- Reference pressure cell
        label pRefCell_;

        //- Reference pressure value
        scalar pRefValue_;

        //- Moving reference frame
        IOdictionary mrfProperties_;

        IOdictionary movingReferenceFrame_;

        scalar lambdaFf_;
        scalar lambdaF0_;

        word centerFromMesh_;
        vector center_;
        vector center0_;

        dimensionedVector XF_;

        dimensionedVector UF_;

        dimensionedVector aF_;

        //- closed volume?
        Switch closedVolume_;

        //- has space patch?
        Switch hasSpacePatch_;

        //- saved time index
        label myTimeIndex_;

        //- Adjustable time-step switch
        const Switch adjustTimeStep_;

        //- For adjustable time-step, this is the maximum Courant number
        const scalar maxCo_;

        //- For adjustable time-step, this is the maximum time-step
        const scalar maxDeltaT_;

        //- Cached momentum matrix
        //  shared between the momentum predictor and pressure corrector
        tmp<fvVectorMatrix> tUEqn;

        //- continuity err
        scalar cumulativeContErr_;

public:

    TypeName("interTrackFluid");

    // Constructors

        interTrackFluid
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~interTrackFluid();


    // Member Functions

        // Access
        const volVectorField& U() const
        {
            return U_();
        }

        const volScalarField& p() const
        {
            return p_();
        }

        volVectorField& U()
        {
            return U_();
        }

        volScalarField& p()
        {
            return p_();
        }

        virtual bool usesPIMPLE() const
        {
            return true;
        }

        virtual word regionTypeName() const
        {
            return interTrackFluid::typeName;
        }

        virtual void correct();

        virtual scalar getMinDeltaT();

        virtual void setCoupledEqns();

        virtual void solveRegion();

        virtual void prePredictor();

        virtual void momentumPredictor();

        virtual void pressureCorrector();

        virtual void meshMotionCorrector();

        virtual void postSolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
