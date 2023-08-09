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
    Foam::icoFluid.C

SourceFiles
    icoFluid.C

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
    
    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/

#ifndef icoFluid_H
#define icoFluid_H

#include "regionType.H"
#include "MRFZones.H"
#include "dimensionedScalar.H"
#include "multiRegionSystem.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"
#include "pimpleControl.H"
#include "solutionControl.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                           Class icoFluid Declaration
\*---------------------------------------------------------------------------*/

class icoFluid
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

        autoPtr<surfaceScalarField> phiHbyA_;

        autoPtr<volScalarField> p_;

        IOobject UfHeader_;

        surfaceVectorField Uf_;

        //- pcorr
        wordList pcorrTypes_;

        autoPtr<volScalarField> pcorr_;

        //- Under-relaxation factor
        scalar UUrf_;

        //- closed volume?
        Switch closedVolume_;

        //- has space patch?
        Switch hasSpacePatch_;

        //- Reference pressure cell
        label pRefCell_;

        //- Reference pressure value
        scalar pRefValue_;

        //- residual for inner correctors
        scalar innerResidual_;

        //- residual for pressure outer correctors
        scalar residualPressure_;

        //- Current corrector loop index
        label corr_;

        //- Current non-orthogonal corrector loop index
        label corrNonOrtho_;

        //- Moving reference frames
        MRFZones mrfZones_;
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

public:

    TypeName("icoFluid");

    // Constructors

        icoFluid
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~icoFluid();


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

        word regionTypeName() const
        {
            return icoFluid::typeName;
        }

        virtual bool usesPIMPLE() const
        {
            return true;
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
