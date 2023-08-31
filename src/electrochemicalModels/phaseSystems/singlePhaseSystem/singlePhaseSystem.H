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
    Foam::singlePhaseSystem.C

SourceFiles
    singlePhaseSystem.C

Author
    Steffen Hess(s.hess@fz-juelich.de, Affiliation A) 
    Constantin Habes (constantin.habes@tu-darmstadt.de, Affiliation B)
    Holger Marschall (holger.marschall@tu-darmstadt.de, Affiliation B)

Contact
    Steffen Hess(s.hess@fz-juelich.de) 
    Holger Marschall (holger.marschall@tu-darmstadt.de)
    main developer and principal investigator
    TU Darmstadt

Affiliations
    Affiliation A)
    Forschungszentrum Jülich IEK-14
    Affiliation B)
    Computational Multiphase Flow
    Department of Mathematics
    Technical University of Darmstadt, Germany

Acknowledgement
    Funded by
    AIDAS – AI, Data Analytics and Scalable Simulation and
    Hessian Ministry of Higher Education, Research, Science and the Arts,
    National High Performance Computing Center for Computational Engineering
    Science (NHR4CES) and 

Description
    
    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/

#ifndef singlePhaseSystem_H
#define singlePhaseSystem_H

#include "dictionary.H"
// #include "phaseModel.H"
#include "volFields.H"
#include "phaseSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class singlePhaseSystem Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseSystem
:
    public phaseSystem
{
    // Private data

        //- Phase model 1
        phaseModel& phase_;

        const Switch surfaceElectrochemistry_;

public:
        TypeName("singlePhaseSystem");

    // Constructors

        //- Construct from fvMesh
        singlePhaseSystem(const fvMesh&);

    //- Destructor
    virtual ~singlePhaseSystem();

        //- Solve for the two-phase-fractions
        virtual void solve();

        //- Correct two-phase properties other than turbulence
        virtual void correct();

        // Access

            //- Single phase
            virtual bool isSinglePhase() const
            {
                return true;
            }

            // //- Return relative flux magnitude
            // virtual tmp<surfaceScalarField> phirMag() const;

            // //- Return max magnitude flux
            // virtual tmp<surfaceScalarField> phiMagMax() const;

        // Surface coupling

            //- Return the TEqn from the phases
            virtual tmp<fvScalarMatrix> TEqn();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
