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
    Foam::pUCoupledIcoFluid.C

SourceFiles
    pUCoupledIcoFluid.C

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

#ifndef pUCoupledIcoFluid_H 
#define pUCoupledIcoFluid_H 

#include "regionType.H"
#include "dimensionedScalar.H"
#include "multiRegionSystem.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOReferencer.H"
#include "HashPtrTable.H"
#include "MRFZones.H"
#include "fvCFD.H"
#include "fvBlockMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionTypes
{

/*---------------------------------------------------------------------------*\
                           Class pUCoupledIcoFluid Declaration
\*---------------------------------------------------------------------------*/

class pUCoupledIcoFluid
:
    public regionType
{       

        //- name
        word regionName_;      

        //- transport
        autoPtr<volVectorField> U_;

        autoPtr<surfaceScalarField> phi_;

        autoPtr<volScalarField> p_;

        //- source term of pressure poisson equation
        autoPtr<surfaceScalarField> presSource_;

        //- Block vector field for velocity (first entry) and pressure (second
        //  entry)
        autoPtr<volVector4Field> Up_;

        //- reverse coefficient field of U for pressure poisson equation
        autoPtr<volScalarField> rAU_;

        //- gradients
        autoPtr<volVectorField> gradp_;
        autoPtr<volTensorField> gradU_;

        //- pcorr
        wordList pcorrTypes_;
        autoPtr<volScalarField> pcorr_;

        //- pressure field bounds
        dimensionedScalar pMin_;
        dimensionedScalar pMax_;

        //- velocity field bounds
        dimensionedScalar UMax_;
        dimensionedScalar smallU_;

        //- transport properties dictionary 
        IOdictionary transportProperties_;

        //- density
        dimensionedScalar rhoFluid_;
        autoPtr<volScalarField> rho_;

        //- viscosity
        dimensionedScalar muFluid_;
        autoPtr<volScalarField> mu_;

        //- closed volume?
        Switch closedVolume_;

        //- has space patch?
        Switch hasSpacePatch_;

        //- Reference pressure cell
        label pRefCell_;

        //- Reference pressure value
        scalar pRefValue_;

        //- List of lable indicating if proc holds pRefCell
        labelList whichProcHasRef_;
        
        //- block solver residual
        // BlockSolverPerformance<vector4> residual_;
        // scalar maxResidual;
        // scalar convergenceCriterion;

        // //- Current corrector loop index
        // label corr_;

        //- Moving reference frames
        MRFZones mrfZones_;
        label myTimeIndex_;

        //- Adjustable time-step switch
        const Switch adjustTimeStep_;

        //- For adjustable time-step, this is the maximum Courant number
        const scalar maxCo_;

        //- For adjustable time-step, this is the maximum time-step
        const scalar maxDeltaT_;

        //- Indicator to correct phi after first solution step
        bool firstUpdate_;

public:

    TypeName("pUCoupledIcoFluid");

    // Constructors

        pUCoupledIcoFluid
        (
            const Time& runTime,
            const word& regionName
        );

    //- Destructor
    virtual ~pUCoupledIcoFluid();


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

        // Access by name
        template<class T>
        const T& getObject(const word& name)
        {
            return
            (
                mesh().thisDb().lookupObject<T>(name)
            );
        }

        template<class T>
        T& getObject(const word& name)
        {
            return
            (
                const_cast<const T&>(mesh().thisDb().lookupObject<T>(name))
            );
        }

        virtual word regionTypeName() const
        {
            return pUCoupledIcoFluid::typeName;
        }

        virtual bool usesPIMPLE() const
        {
            return false;
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
