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

\*---------------------------------------------------------------------------*/

inline const Foam::fvPatch& Foam::regionInterface::patch() const
{
    return patch_;
};


inline const Foam::word& Foam::regionInterface::regionBName() const
{
    return regionBName_;
};

inline const Foam::word& Foam::regionInterface::patchName() const
{
    return patchName_;
};

inline const Foam::word& Foam::regionInterface::patchBName() const
{
    return patchBName_;
};


inline const Foam::word& Foam::regionInterface::BFieldName() const
{
     return BFieldName_;
};


inline const Foam::fvMesh& Foam::regionInterface::regionB() const
{
    return localRegion_.objectRegistry::parent()
        .lookupObject<fvMesh>(regionBName_);
}


inline Foam::label Foam::regionInterface::patchBID() const
{
    return regionB().boundaryMesh().findPatchID(patchBName_);
}


inline const Foam::fvPatch& Foam::regionInterface::patchB() const
{
    return regionB().boundary()[patchBID()];
}


template<class Type>
inline const Foam::fvPatchField<Type>&
Foam::regionInterface::patchBField() const
{
    return patchB().fvPatch::lookupPatchField
        <GeometricField<Type, fvPatchField, volMesh>, Type>
            (BFieldName_);
}


inline Foam::word& Foam::regionInterface::regionBName()
{
    return regionBName_;
};


inline Foam::word& Foam::regionInterface::patchBName()
{
    return patchBName_;
};


inline Foam::word& Foam::regionInterface::BFieldName()
{
    return BFieldName_;
};

//inline const Foam::GGIInterpolation<polyPatch, polyPatch>&
//Foam::regionInterface::ggiInterpolate
//(
//    const polyPatch& nbrPatch,
//    const polyPatch& ownPatch
//) const
//{
//    if (!ggiInterpolatePtr_)
//    {
//        ggiInterpolatePtr_ =
//            new GGIInterpolation<polyPatch, polyPatch>
//            (
//                nbrPatch,
//                ownPatch,
//                tensorField(0),
//                tensorField(0),
//                vectorField(0), // Slave-to-master separation. Bug fix
////                true,           // Patch data is complete on all processors
//                0,              // Non-overlapping face tolerances
//                0,
//                true,           // Rescale weighting factors
////                ggiInterpolation::BB_OCTREE
//                ggiInterpolation::AABB
//            );
//    }

//    return *ggiInterpolatePtr_;
//}


// ************************************************************************* //
