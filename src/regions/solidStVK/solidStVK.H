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
    Foam::solidStVK.C

SourceFiles
    solidStVK.C

Author
    Holger Marschall (holger.marschall@tu-darmstadt.de, Affiliation A)
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


#ifndef solidStVK_H
#define solidStVK_H

#include "regionType.H"
#include "dimensionedScalar.H"
#include "multiRegionSystem.H"
#include "volFields.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"
#include "volPointInterpolation.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"
#include "steadyStateDdtScheme.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                           Class solidStVK Declaration
\*---------------------------------------------------------------------------*/

class solidStVK
:
    public regionType
{
        //- name
        word regionName_;

        // Interpolator for mesh deformation
        const volPointInterpolation& cpi_;

        // Properties
        IOdictionary mechanicalProperties_;
        dimensionedScalar rho_;
        dimensionedScalar E_;
        dimensionedScalar nu_;
        dimensionedScalar mu_;
        dimensionedScalar lambda_;
        dimensionedScalar threeK_;

        // Finite volume fields pointers
        autoPtr<volVectorField> D_;
        autoPtr<pointVectorField> pointD_;
        autoPtr<volVectorField> U_;
        autoPtr<volTensorField> gradD_;
        autoPtr<volSymmTensorField> sigma_;
        autoPtr<volVectorField> DPrevOuterIter_;

        //- Cached deformation matrix
        tmp<fvVectorMatrix> DEqn;

public:

    TypeName("solidStVK");

    // Constructors

        solidStVK
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~solidStVK();

    // Member Functions

        // Access
        const volVectorField& D() const
        {
            return D_();
        }

        volVectorField& D()
        {
            return D_();
        }

        const pointVectorField& pointD() const
        {
            return pointD_();
        }

        pointVectorField& pointD()
        {
            return pointD_();
        }

        const volVectorField& U() const
        {
            return U_();
        }

        volVectorField& U()
        {
            return U_();
        }

        word regionTypeName() const
        {
            return solidStVK::typeName;
        }

        virtual bool usesPIMPLE() const
        {
            return false;
        }

        virtual void correct();

        virtual scalar getMinDeltaT();

        virtual void solveRegion();

        virtual void prePredictor();

        virtual void momentumPredictor();

        virtual void pressureCorrector();

        virtual void meshMotionCorrector();

        virtual void setCoupledEqns();

        virtual void postSolve();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace regionTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
