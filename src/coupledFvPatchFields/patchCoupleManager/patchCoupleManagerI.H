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

inline const Foam::fvPatch& Foam::patchCoupleManager::patch() const
{
    return patch_;
};


inline const Foam::word& Foam::patchCoupleManager::neighbourRegionName() const
{
    return neighbourRegionName_;
};


inline const Foam::word& Foam::patchCoupleManager::neighbourPatchName() const
{
    return neighbourPatchName_;
};


inline const Foam::word& Foam::patchCoupleManager::neighbourFieldName() const
{
     return neighbourFieldName_;
};


inline const Foam::fvMesh& Foam::patchCoupleManager::neighbourRegion() const
{
    return localRegion_.objectRegistry::parent()
        .lookupObject<fvMesh>(neighbourRegionName_);
}


inline Foam::label Foam::patchCoupleManager::neighbourPatchID() const
{
    return neighbourRegion().boundaryMesh().findPatchID(neighbourPatchName_);
}


inline const Foam::fvPatch& Foam::patchCoupleManager::neighbourPatch() const
{
    return neighbourRegion().boundary()[neighbourPatchID()];
}


template<class Type>
inline const Foam::fvPatchField<Type>&
Foam::patchCoupleManager::neighbourPatchField() const
{
    return neighbourPatch().fvPatch::lookupPatchField
        <GeometricField<Type, fvPatchField, volMesh>, Type>
            (neighbourFieldName_);
}


inline Foam::word& Foam::patchCoupleManager::neighbourRegionName()
{
    return neighbourRegionName_;
};


inline Foam::word& Foam::patchCoupleManager::neighbourPatchName()
{
    return neighbourPatchName_;
};


inline Foam::word& Foam::patchCoupleManager::neighbourFieldName()
{
    return neighbourFieldName_;
};

inline const Foam::GGIInterpolation<polyPatch, polyPatch>&
Foam::patchCoupleManager::ggiInterpolator
(
    const polyPatch& nbrPatch,
    const polyPatch& ownPatch
) const
{
    if (!ggiInterpolatorPtr_)
    {
        ggiInterpolatorPtr_ =
            new GGIInterpolation<polyPatch, polyPatch>
            (
                nbrPatch,
                ownPatch,
                tensorField(0),
                tensorField(0),
                vectorField(0), // Slave-to-master separation. Bug fix
                true,           // Patch data is complete on all processors
                SMALL,          // Non-overlapping face tolerances
                SMALL,
                true,           // Rescale weighting factors
                ggiInterpolation::BB_OCTREE
            );
//            (
//                nbrPatch,
//                ownPatch,
//                tensorField(0),
//                tensorField(0),
//                vectorField(0), // Slave-to-master separation. Bug fix
//                true,           // Patch data is complete on all processors
//                0,              // Non-overlapping face tolerances
//                0,
//                true,           // Rescale weighting factors
//                ggiInterpolation::BB_OCTREE
////                ggiInterpolation::AABB
//            );
    }

    return *ggiInterpolatorPtr_;
}


// ************************************************************************* //
