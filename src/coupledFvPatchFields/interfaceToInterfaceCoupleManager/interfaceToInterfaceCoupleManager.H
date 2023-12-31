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
    Helper class to handle the coupling of region patches.
    Deployment with the regionInterfaceType class within the unified
    multi-region coupling framework. It can be queried to return
    the neighbour information and data.

    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/


#ifndef interfaceToInterfaceCoupleManager_H
#define interfaceToInterfaceCoupleManager_H

#include "Ostream.H"
#include "regionInterfaceType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class interfaceToInterfaceCoupleManager Declaration
\*---------------------------------------------------------------------------*/

class interfaceToInterfaceCoupleManager
{
    // Private data

        //- Reference to the local fvPatch
        const fvPatch& patch_;

        //- Region interface type name
        word typeName_;

        //- Name of neighbour region
        word neighbourRegionName_;

        //- Name of patch on the neighbour region
        word neighbourPatchName_;

        //- Name of field on the neighbour region
        word neighbourFieldName_;

        //- Reference to the local region
        const fvMesh& localRegion_;

    // Private Member Functions

        //- assemble interface name
        static word assembleName
        (
            const fvPatch& refPatch,
            const word& nbrRgName,
            const word& nbrPatchName,
            const word& typeName,
            const bool& reverseOrder=false
        );

        //- Disallow default bitwise copy construct
        // interfaceToInterfaceCoupleManager
        //     (
        //         const interfaceToInterfaceCoupleManager&
        //     );

        //- Disallow default bitwise assignment
        void operator=(const interfaceToInterfaceCoupleManager&);

protected:

        //- Reference to meshInterface
        const regionInterfaceType& rgInterface() const;

        //- Reference to neighbour mesh
        const fvMesh& nbrMesh() const;

        //- Reference to neighbour patch
        const fvPatch& nbrPatch() const;

public:

    // Constructors

        //- Construct from fvPatch
        interfaceToInterfaceCoupleManager
        (
            const fvPatch& patch,
            const word type=word("default")
        );

        //- Construct from fvPatch and dictionary
        interfaceToInterfaceCoupleManager
        (
            const fvPatch& patch,
            const dictionary& dict,
            const word type=word("default")
        );

        //- Copy constructor
        interfaceToInterfaceCoupleManager
        (
            const interfaceToInterfaceCoupleManager& pcm
        );

    // Destructor

        ~interfaceToInterfaceCoupleManager();

    // Member Functions

            //- Interpolate from neighbour patch values
            template<class Type>
            tmp<Field<Type> > interpolateFromNbrField
            (
                const Field<Type>& fromField
            ) const;

        // Access

            //- Return a reference to the local mesh
            inline const fvMesh& refMesh() const;

            //- Return a reference to the local patch
            inline const fvPatch& refPatch() const;

            //- Return a reference to the local patch
            inline label refPatchID() const;

            //- Return the name of the neighbour region
            inline const word& neighbourRegionName() const;

            //- Return the name of the patch on the neighbour region
            inline const word& neighbourPatchName() const;

            //- Return the name of the field on the neighbour region
            inline const word& neighbourFieldName() const;

            //- Return a reference to the neighbour mesh
            inline const fvMesh& neighbourRegion() const;

            //- Return the neighbour patch ID
            inline label neighbourPatchID() const;

            //- Return a reference to the neighbour patch
            inline const fvPatch& neighbourPatch() const;

            //- Return a reference to the neighbour patch field
            template<class Type>
            inline const fvPatchField<Type>& neighbourPatchField() const;

            //- Update region-interface
            void updateRegionInterface();

        // Edit

            //- Return the name of the neighbour region
            word& neighbourRegionName();

            //- Return the name of the patch on the neighbour region
            word& neighbourPatchName();

            //- Return the name of the field on the neighbour region
            word& neighbourFieldName();

        // Write

           //- Write couple to obj file
           void coupleToObj() const;

           //- Write entries for patches
           void writeEntries(Ostream& os) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interfaceToInterfaceCoupleManagerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
