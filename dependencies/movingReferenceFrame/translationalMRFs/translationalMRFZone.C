/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "translationalMRFZone.H"
#include "cylindricalCS.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"
#include "syncTools.H"
#include "faceSet.H"
#include "geometricOneField.H"
#include "IOstream.H"
#include "IFstream.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(translationalMRFZone, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::translationalMRFZone::centreControls,
        3
    >::names[] =
    {
        "Mesh",
        "Field",
        "Model"
    };
}

const Foam::NamedEnum
<
    Foam::translationalMRFZone::centreControls,
    3
> Foam::translationalMRFZone::centreControlNames_;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::translationalMRFZone::PIDCorr()
{
    dimensionedVector errorValueOld = errorValue_;

    errorValue_ = -(centreTarget_ - centreCurrent_);

//    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
//    {
//        errorValue_.replace(cmpt, mag(errorValue_.component(cmpt)));
//    }

    if (debug)
    {
        Info << nl << "PID error value = " << errorValue_.value() << endl;
    }
    
    dimensionedScalar deltaT = mesh().time().deltaT0();

    integralComponent_ += errorValue_*deltaT;

    dimensionedScalar KDnew = KD_*deltaT;

    PIDOut_ = 
        (KP_*errorValue_)
      + KI_*integralComponent_
      + KDnew*(errorValue_ - errorValueOld)/deltaT;

    Info<< nl
        << "PID-Controller:" << nl
        << "  P = " << (KP_*errorValue_).value() << nl
        << "  I = " << (KI_*integralComponent_).value() << nl
        << "  D = " << (KDnew*(errorValue_-errorValueOld)/deltaT).value() 
        << endl;
}

void Foam::translationalMRFZone::readMRFCentreControl()
{
    if (centreControl_ == Mesh)
    {
        dict_.lookup("centreFromMesh") >> regionName_;
    }
    else if (centreControl_ == Field)
    {
        dict_.lookup("centreFromField") >> fieldName_;
    }
    else if (centreControl_ == Model)
    {
        dict_.lookup("centreFromModel") >> modelName_;
    }
}

void Foam::translationalMRFZone::setMRFCentreControl()
{
    if (centreControl_ == Mesh)
    {   
        if (mesh().objectRegistry::parent().foundObject<fvMesh>(regionName_))
        {
            const fvMesh& centreMesh = mesh().objectRegistry::parent()
                .lookupObject<fvMesh>(regionName_);

            centreCurrent_ = centreMesh.C().weightedAverage(centreMesh.V());

//            centreCurrent_ =
//                gSum(centreMesh.C().internalField()*centreMesh.V())
//                /gSum(centreMesh.V());
        }
        else
        {
            FatalError << "Mesh for computing MRF centre " 
                << regionName_
                << " not found" 
                << endl;
            FatalErrorIn("setMRFCentreControl")
		    << abort(FatalError);
        }
    }
    else if (centreControl_ == Field)
    {
        if (mesh().foundObject<volScalarField>(fieldName_))
        {
            const volScalarField& centreField =
                mesh().lookupObject<volScalarField>(fieldName_);;

            volScalarField alpha = centreField;

            if (normalise01_)
            {
                alpha = (1.+centreField)/2.;
            }

            centreCurrent_ = mesh().C().weightedAverage(alpha*mesh().V());
        }
        else
        {
            FatalError << "Field for computing MRF centre " 
                << fieldName_
                << " not found" 
                << endl;
        }
    }
    else if (centreControl_ == Model)
    {
        const dimensionedVector centreModel =
        (
            mesh().objectRegistry::
            lookupObject<IOReferencer<dimensionedVector> >
            (modelName_)()
        );

        centreCurrent_ = centreModel;
    }

    centreCurrent_.value() =
        cmptMultiply(centreCurrent_.value(), centreMask_);

    centreTarget_.value() =
        cmptMultiply(centreTarget_.value(), centreMask_);
}

void Foam::translationalMRFZone::setPatchIDs()
{
    // TODO: list of entries so as to allow multiple patches

    // get wall patch index
    wallPatchID_ = mesh().boundaryMesh().findPatchID(wallPatchName_);

    if 
    (
        wallPatchID_ != -1
    )
    {
        if (mesh_.boundary()[wallPatchID_].size() == 0)
        {
            FatalError << "Moving wall patch " 
                << wallPatchName_
                << " of zero size" 
                << endl;
        }
    }
    else
    {
        Warning << "Moving wall patch " 
            << wallPatchName_
            << " for MRF not found" << nl
            << endl;
    }

    // get inlet patch index
    inletPatchID_ = mesh_.boundaryMesh().findPatchID(inletPatchName_);

    if 
    (
        inletPatchID_ != -1
    )
    {
        label inletPatchSize = mesh_.boundary()[inletPatchID_].size();
        reduce(inletPatchSize, sumOp<label>());

        if (inletPatchSize == 0)
        {
            FatalError << "Inlet patch of zero size" << endl;
        }
    }
    else
    {
        Warning << "Inlet patch " 
            << inletPatchName_
            << " for MRF not found" << nl
            << endl;
    }


    // get space patch index
    spacePatchID_ = mesh_.boundaryMesh().findPatchID(spacePatchName_);

    if 
    (
        spacePatchID_ != -1
    )
    {
        if (mesh_.boundary()[spacePatchID_].size() == 0)
        {
            FatalError << "Space patch of zero size" << endl;
        }
    }
    else
    {
        Warning << "Space patch " 
            << spacePatchName_
            << " for MRF not found" << nl
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::translationalMRFZone::translationalMRFZone(const fvMesh& mesh, Istream& is)
:
    name_(is),
    mesh_(mesh),
    Umrf_
    (
        IOobject
        (
            "Umrf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_, 
        dimensionedVector("Umrf", dimVelocity, vector::zero)
    ),
    dict_(is),

    KP_(dict_.lookup("KP")),
    KI_(dict_.lookup("KI")),
    KD_(dict_.lookup("KD")),
    PIDOut_("PIDOut", dimLength, vector::zero),

    wallPatchName_(dict_.lookupOrDefault<word>("wallPatchName", "")),
    inletPatchName_(dict_.lookupOrDefault<word>("inletPatchName", "")),
    spacePatchName_(dict_.lookupOrDefault<word>("spacePatchName", "")),
    wallPatchID_(-1),
    inletPatchID_(-1),
    spacePatchID_(-1),
    centreCurrent_("centreCurrent", mesh_.C().weightedAverage(mesh_.V())),
    // particle position at simulation start is set to be the target position
    centreTarget_("centreTarget", centreCurrent_),
    errorValue_("errorValue", dimLength, vector::zero),
    integralComponent_
    (
        "integralComponent",
        dimensionSet(0,1,1,0,0), 
        vector::zero
    ),
    UP_("Up", dimVelocity, vector::zero),
    UPrel_("relUp", dimVelocity, vector::zero),
    dUF_("dUF", dimVelocity, vector::zero),
    dUFMask_
    (
        vector
        (
            dict_.lookupOrDefault<vector>("dUFMask", vector::zero)
        )
    ),
    centreMask_
    (
        vector
        (
            dict_.lookupOrDefault<vector>("centreMask", vector::one)
        )
    ),
    aP_("partAcc", dimVelocity/dimTime, vector::zero),
    aPrel_
    (
        "partRelAcc", 
        dimVelocity/dimTime, 
        vector::zero
    ),
    aF_("aF", dimVelocity/dimTime, vector::zero),
    UF_("UF", dimVelocity, vector::zero),
    xF_("xF", dimLength, vector::zero), 
    Uwall_("Uwall", dimVelocity, vector::zero),
    regionName_(word::null),
    fieldName_(word::null),
    normalise01_(dict_.lookupOrDefault<Switch>("normalise01", false)),
    modelName_(word::null),
    projectUWall_(dict_.lookupOrDefault<Switch>("projectUWall", true)),
    projectUInlet_(dict_.lookupOrDefault<Switch>("projectUInlet", true)),
    centreControl_(centreControlNames_.read(dict_.lookup("centreFrom")))
{
    // set velocity mask
    if (!dict_.found("dUFMask"))
    {
        Info << "Setting velocity mask for free motion" << endl;

        const Vector<label>& dirs = mesh_.solutionD();

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            dUFMask_.replace(cmpt, pos(dirs.component(cmpt)));
        }
    }

    Info << "dUFMask = " << dUFMask_ << nl << endl;

    // read MRF centre control method
    readMRFCentreControl();

    // set MRF centre control method
    //setMRFCentreControl();

    // If restart file present, use it
    readRestart();

    // Enforce consistent target for controller
    centreTarget_ = centreCurrent_;

    // Set wall/inlet/space patch IDs
    setPatchIDs();
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::translationalMRFZone::~translationalMRFZone()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::translationalMRFZone::correctMRF()
{
    dimensionedScalar deltaT = mesh().time().deltaT0();

    // MRF parameters from old time step / iteration
    dimensionedVector centreOld = centreCurrent_;

    dimensionedVector UpOld = UP_;

    dimensionedVector relUpOld = UPrel_;

    dimensionedVector xFOld = xF_;

    // Set MRF centre control
    setMRFCentreControl();

    // Relative velocity of the particle on the mesh
    UPrel_ = (centreCurrent_ - centreOld)/deltaT;

    // Relative acceleration of the particle on the mesh
    aPrel_ = (UPrel_ - relUpOld)/deltaT;

    // Velocity increment for adjustment
    PIDCorr();

    dUF_ = PIDOut_/deltaT;

    dUF_.value() = cmptMultiply(dUF_.value(), dUFMask_);

    dimensionedVector xFAdjust = -deltaT*(0.5*dUF_ + UF_); //TODO! different from Henrik Rusche

    UF_ += dUF_;

    xF_ += xFAdjust;

    UP_ = (xF_ - xFOld)/deltaT;

    aP_ = (UP_ - UpOld)/deltaT;

//    volVectorField& U = const_cast<volVectorField&>
//        (mesh().lookupObject<volVectorField>("U"));

//    surfaceScalarField& phi = const_cast<surfaceScalarField&>
//        (mesh().lookupObject<surfaceScalarField>("phi"));

//    U += dUF_;
//    phi += (mesh().Sf() & (dUF_));

    aF_ = dUF_/deltaT;

    Info<< "Adjustment Velocity = " << dUF_.value() << nl
        << "Target Position = " << centreTarget_.value() << nl
        << "Particle Position = " << centreCurrent_.value() << nl
        << "  Velocity = " << UP_.value() << nl
        << "  Acceleration = " << aP_.value() << nl
        << "  Relative Velocity = " << UPrel_.value() << nl
        << "  Relative Acceleration = " << aPrel_.value() << nl
        << "Frame Position = " << xF_.value() << nl
        << "  Velocity = " << UF_.value() << nl
        << "  Acceleration = " << aF_.value() << nl
        << endl;
}

void Foam::translationalMRFZone::addFrameAcceleration
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho
)
{
    UEqn += rho*aF_;
}

void Foam::translationalMRFZone::addFrameAcceleration
(
    fvVectorMatrix& UEqn
)
{
    UEqn += aF_;
}


void Foam::translationalMRFZone::correctBoundaryVelocity
(
    volVectorField& U,
    surfaceScalarField& phi
)
{
    if (spacePatchID_ != -1)
    {
        if
        (
            U.boundaryField()[spacePatchID_].type()
         == inletOutletFvPatchVectorField::typeName
        )
        {
            inletOutletFvPatchVectorField& spaceU =
                refCast<inletOutletFvPatchVectorField>
                (
                    U.boundaryField()[spacePatchID_]
                );

            spaceU.refValue() = -UF_.value();

//            forAll(mesh().Cf().boundaryField()[spacePatchID_], faceI)
//            {
//                if (phi.boundaryField()[spacePatchID_][faceI] < SMALL)
//                {
//                    U.boundaryField()[spacePatchID_][faceI] = -UF_.value();

//                    phi.boundaryField()[spacePatchID_][faceI] =
//                    (
//                        mesh_.Sf().boundaryField()[spacePatchID_][faceI]
//                        &U.boundaryField()[spacePatchID_][faceI]
//                    );
//                }
//            }
        }
        else
        {
            FatalErrorIn("MRFModels::translational::correctBoundaryVelocity()")
                << "Velocity boundary condition at space patch must be "
                << inletOutletFvPatchVectorField::typeName 
                << abort(FatalError);
        }
    } 

    if (wallPatchID_ != -1)
    {
        if
        (
                U.boundaryField()[wallPatchID_].type()
             == fixedValueFvPatchVectorField::typeName
          ||
                U.boundaryField()[wallPatchID_].type()
             == noSlipWallFvPatchVectorField::typeName //recommended
        )
        {
            U.boundaryField()[wallPatchID_] == -UF_.value();

            U.correctBoundaryConditions();
        }
        else if
        (
            U.boundaryField()[wallPatchID_].type()
         == mixedFvPatchVectorField::typeName
        )
        {
            mixedFvPatchVectorField& Uwall =
                refCast<mixedFvPatchVectorField>
                (
                    U.boundaryField()[wallPatchID_]
                );

            Uwall.refValue() = -UF_.value();
            Uwall.refGrad() = vector::zero;
            Uwall.valueFraction() = 1;
        }
        else
        {
            FatalErrorIn("MRFModels::translational::correctBoundaryVelocity()")
                << "Velocity boundary condition at wall patch must be "
                << fixedValueFvPatchVectorField::typeName << " or "
                << mixedFvPatchVectorField::typeName
                << abort(FatalError);
        }
    }

    if (inletPatchID_ != -1)
    {
        if
        (
            U.boundaryField()[inletPatchID_].type()
         == fixedValueFvPatchVectorField::typeName
        )
        {
            // account for inlet boundary (if present)
            U.boundaryField()[inletPatchID_] == -UF_.value();

            phi.boundaryField()[inletPatchID_] == 
            (
                mesh_.Sf().boundaryField()[inletPatchID_]
                &U.boundaryField()[inletPatchID_]
            );
        }
        else
        {
            FatalErrorIn("MRFModels::translational::correctBoundaryVelocity()")
                << "Velocity boundary condition at inlet patch must be "
                << fixedValueFvPatchVectorField::typeName 
                << abort(FatalError);
        }
    }
}


void Foam::translationalMRFZone::writeRestart()
{
    if (mesh().time().outputTime())
    {
        const volVectorField& U = mesh().lookupObject<volVectorField>("U");

        Umrf_ = U - UF_;

        OFstream restartFile
        (
            mesh().time().path()/mesh().time().timeName()/"MRFrestart.raw"
        );

        if(restartFile.good())
        {
            restartFile << mesh().time().deltaT().value() << endl
            << centreCurrent_.value() << endl
            << centreTarget_.value() << endl
            << UF_.value() << endl
            << errorValue_.value() << endl
            << UP_.value() << endl
            << UPrel_.value() << endl
            << dUF_.value() << endl
            << integralComponent_.value() << endl
            << aP_.value() << endl
            << aPrel_.value() << endl
            << aF_.value() << endl
            << UF_.value() << endl
            << xF_.value() << endl
            << Uwall_.value() << endl;
        }
    }
}

void Foam::translationalMRFZone::readRestart()
{
    IFstream restartFile
    (
        mesh().time().path()/mesh().time().timeName()/"MRFrestart.raw"
    );

    if (restartFile)
    {
        Info << "restart from file!!! " << endl;

        if (restartFile.good())
        {
            scalar s;
            restartFile >> s;
            const_cast<Time&>(mesh().time()).setDeltaT(s, false);

            vector v;
            restartFile >> v;
            centreCurrent_.value() = v;

            restartFile >> v;
            centreTarget_.value() = v;

            restartFile >> v;
            UF_.value() = v;

            restartFile >> v;
            errorValue_.value() = v;

            restartFile >> v;
            UP_.value() = v;

            restartFile >> v;
            UPrel_.value() = v;

            restartFile >> v;
            dUF_.value() = v;

            restartFile >> v;
            integralComponent_.value() = v;

            restartFile >> v;
            aP_.value() = v;

            restartFile >> v;
            aPrel_.value() = v;

            restartFile >> v;
            aF_.value() = v;

            restartFile >> v;
            UF_.value() = v;

            restartFile >> v;
            xF_.value() = v;

            restartFile >> v;
            Uwall_.value() = v;
        }
    }
}

// ************************************************************************* //
