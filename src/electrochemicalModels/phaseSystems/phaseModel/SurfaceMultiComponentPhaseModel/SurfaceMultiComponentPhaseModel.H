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
    Foam::SurfaceMultiComponentPhaseModel.C

SourceFiles
    SurfaceMultiComponentPhaseModel.C

Author
    Steffen Hess(s.hess@fz-juelich.de, Affiliation A) 
    All rights reserved.

Contact
    Steffen Hess(s.hess@fz-juelich.de) 

Affiliations
    Affiliation A)
    Forschungszentrum Jülich IEK-14

Acknowledgement
    Funded by
    AIDAS – AI, Data Analytics and Scalable Simulation

Description
    
    This file is part of the multiRegionFoam library.

    You may refer to this software as :
    https://doi.org/10.48550/arXiv.2306.01924

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/

#ifndef SurfaceMultiComponentPhaseModel_H
#define SurfaceMultiComponentPhaseModel_H

#include "phaseModel.H"
#include "multiSpeciesTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class SurfaceMultiComponentPhaseModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseModel>
class SurfaceMultiComponentPhaseModel
:
    public BasePhaseModel
{
private:

    volScalarField rhoField_;

protected:

    // Protected data

        //- Thermo physical constitutive model
        autoPtr<multiSpeciesTransportModel> thermoPhysicalTransportModel_;

        //- Residual phase fraction
        dimensionedScalar residualAlpha_;

        //- Inert species index
        label inertIndex_;

        //- Pointer list to active species
        UPtrList<volScalarField> YActive_;

        //- mole fraction
        PtrList<volScalarField> X_;

        //- Mole fraction list
        //- Note: iDmdt is created here
        //- iDmdt is updated in electroChemicalReaction
        PtrList<volScalarField> iDmdt_;

        //- Pointer list to active species
        UPtrList<volScalarField> XActive_;


public:

    // Constructors

        SurfaceMultiComponentPhaseModel
        (
            const phaseSystem& fluid,
            const word& phaseName,
            const label index
        );


    //- Destructor
    virtual ~SurfaceMultiComponentPhaseModel();


    // Member Functions

        //- Correct the thermodynamics
        virtual void correctThermo();

        //- Correct crossover boundary condition (specie flux and velocity)
        virtual void correctBC();

        // Species

            //- Return whether the phase is pure (i.e., not multi-component)
            virtual bool pure() const;

            //- Return the species fraction equation
            virtual tmp<fvScalarMatrix> YiEqn(volScalarField& Yi);

            //- Return the source/sink term
            virtual const PtrList<volScalarField>& iDmdt() const;

            //- Return the source/sink term by name
            virtual const volScalarField& iDmdt(const word& name) const;

            //- Access the source/sink term
            virtual PtrList<volScalarField>& iDmdtRef();

            //- Return the species mass fractions
            virtual const PtrList<volScalarField>& Y() const;

            //- Return a species mass fraction by name
            virtual const volScalarField& Y(const word& name) const;

            //- Access the species mass fractions
            virtual PtrList<volScalarField>& YRef();

            //- Return the active species mass fractions
            virtual const UPtrList<volScalarField>& YActive() const;

            //- Access the active species mass fractions
            virtual UPtrList<volScalarField>& YActiveRef();

            //- Return the species mole fractions
            virtual const PtrList<volScalarField>& X() const;

            //- Return a species mole fraction by name
            virtual const volScalarField& X(const word& name) const;

            //- Access the species mole fractions
            virtual PtrList<volScalarField>& XRef();

            //- Return the active species mole fractions
            virtual const UPtrList<volScalarField>& XActive() const;

            //- Access the active species mole fractions
            virtual UPtrList<volScalarField>& XActiveRef();

            //- Return the total mass transfer
            virtual tmp<volScalarField> dmdt() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SurfaceMultiComponentPhaseModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
