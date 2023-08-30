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
    interfaceToInterfaceCoupleManager

SourceFiles
    interfaceToInterfaceCoupleManager.C

Author
    Holger Marschall (holger.marschall@tu-darmstadt.de, Affiliation A)
    Constantin Habes (constantin.habes@tu-darmstadt.de, Affiliation A)
    Heba Alkafri (heba.alkafri@tu-darmstadt.de, Affiliation A)
    Luca Huefner (lucas.huefner@tu-darmstadt.de, Affiliation B)

Contact
    Holger Marschall (holger.marschall@tu-darmstadt.de)
    main developer and principal investigator
    TU Darmstadt

Affiliations
    Affiliation A)
    Computational Multiphase Flow
    Department of Mathematics
    Technical University of Darmstadt, Germany

    Affiliation B)
    Technische Universität Darmstadt
    AK Etzold / Technical Chemistry 1
    Department of Chemistry
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

    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


inline const Foam::fvMesh& Foam::interfaceToInterfaceCoupleManager::refMesh() const
{
    return refPatch().boundaryMesh().mesh();
};

inline const Foam::fvPatch& Foam::interfaceToInterfaceCoupleManager::refPatch() const
{
    return patch_;
};

inline Foam::label Foam::interfaceToInterfaceCoupleManager::refPatchID() const
{
    return refMesh().boundaryMesh().findPatchID(refPatch().name());
};

inline const Foam::word& Foam::interfaceToInterfaceCoupleManager::neighbourRegionName() const
{
    return neighbourRegionName_;
};

inline const Foam::word& Foam::interfaceToInterfaceCoupleManager::neighbourPatchName() const
{
    return neighbourPatchName_;
};

inline const Foam::word& Foam::interfaceToInterfaceCoupleManager::neighbourFieldName() const
{
     return neighbourFieldName_;
};

inline const Foam::fvMesh& Foam::interfaceToInterfaceCoupleManager::neighbourRegion() const
{
    return localRegion_.objectRegistry::parent()
        .lookupObject<fvMesh>(neighbourRegionName_);
}

inline Foam::label Foam::interfaceToInterfaceCoupleManager::neighbourPatchID() const
{
    return neighbourRegion().boundaryMesh().findPatchID(neighbourPatchName_);
}

inline const Foam::fvPatch& Foam::interfaceToInterfaceCoupleManager::neighbourPatch() const
{
    return neighbourRegion().boundary()[neighbourPatchID()];
}

template<class Type>
inline const Foam::fvPatchField<Type>&
Foam::interfaceToInterfaceCoupleManager::neighbourPatchField() const
{
    return neighbourPatch().fvPatch::lookupPatchField
        <GeometricField<Type, fvPatchField, volMesh>, Type>
            (neighbourFieldName_);
}

inline Foam::word& Foam::interfaceToInterfaceCoupleManager::neighbourRegionName()
{
    return neighbourRegionName_;
};

inline Foam::word& Foam::interfaceToInterfaceCoupleManager::neighbourPatchName()
{
    return neighbourPatchName_;
};

template<class Type>
Foam::tmp<Field<Type> >
Foam::interfaceToInterfaceCoupleManager::interpolateFromNbrField
(
    const Field<Type>& fromField
) const
{
    // B is neighbour
    if (refPatch().name() == rgInterface().patchA().name())
    {
        return
        (
            rgInterface().interpolateFacesFromB<Type>(fromField)
        );
    }
    // A is neighbour
    return
    (
        rgInterface().interpolateFacesFromA<Type>(fromField)
    );
}

// ************************************************************************* //
