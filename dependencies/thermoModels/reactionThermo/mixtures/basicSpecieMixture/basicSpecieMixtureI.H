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



inline bool Foam::basicSpecieMixture::contains
(
    const word& specieName
) const
{
    return species_.contains(specieName);
}

inline Foam::label Foam::basicSpecieMixture::defaultSpecie() const
{
    return defaultSpecieIndex_;
}


inline bool Foam::basicSpecieMixture::active(label speciei) const
{
    return active_[speciei];
}


inline const Foam::List<bool>& Foam::basicSpecieMixture::active() const
{
    return active_;
}


inline void Foam::basicSpecieMixture::setActive(label speciei)
{
    active_[speciei] = true;
}


inline void Foam::basicSpecieMixture::setInactive(label speciei)
{
    active_[speciei] = false;
}

// Return mass fractions
inline Foam::PtrList<Foam::volScalarField>&
Foam::basicSpecieMixture::Y()
{
    return Y_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::basicSpecieMixture::Y() const
{
    return Y_;
}


inline Foam::volScalarField& Foam::basicSpecieMixture::Y(const label i)
{
    return Y_[i];
}


inline const Foam::volScalarField& Foam::basicSpecieMixture::Y
(
    const label i
) const
{
    return Y_[i];
}


inline Foam::volScalarField& Foam::basicSpecieMixture::Y
(
    const word& specieName
)
{
    return Y_[species_[specieName]];
}


inline const Foam::volScalarField& Foam::basicSpecieMixture::Y
(
    const word& specieName
) const
{
    return Y_[species_[specieName]];
}

// // Return mole fractions
// inline Foam::PtrList<Foam::volScalarField>&
// Foam::basicSpecieMixture::X()
// {
//     return X_;
// }


// inline const Foam::PtrList<Foam::volScalarField>&
// Foam::basicSpecieMixture::X() const
// {
//     return X_;
// }


// inline Foam::volScalarField& Foam::basicSpecieMixture::X(const label i)
// {
//     return X_[i];
// }


// inline const Foam::volScalarField& Foam::basicSpecieMixture::X
// (
//     const label i
// ) const
// {
//     return X_[i];
// }


// inline Foam::volScalarField& Foam::basicSpecieMixture::X
// (
//     const word& specieName
// )
// {
//     return X_[species_[specieName]];
// }


// inline const Foam::volScalarField& Foam::basicSpecieMixture::X
// (
//     const word& specieName
// ) const
// {
//     return X_[species_[specieName]];
// }

// // Return molar concentration fractions
// inline Foam::PtrList<Foam::volScalarField>&
// Foam::basicSpecieMixture::C()
// {
//     return C_;
// }


// inline const Foam::PtrList<Foam::volScalarField>&
// Foam::basicSpecieMixture::C() const
// {
//     return C_;
// }


// inline Foam::volScalarField& Foam::basicSpecieMixture::C(const label i)
// {
//     return C_[i];
// }


// inline const Foam::volScalarField& Foam::basicSpecieMixture::C
// (
//     const label i
// ) const
// {
//     return C_[i];
// }


// inline Foam::volScalarField& Foam::basicSpecieMixture::C
// (
//     const word& specieName
// )
// {
//     return C_[species_[specieName]];
// }


// inline const Foam::volScalarField& Foam::basicSpecieMixture::C
// (
//     const word& specieName
// ) const
// {
//     return C_[species_[specieName]];
// }

// inline void basicSpecieMixture::updateMolarFields()
// {
//     const PtrList<volScalarField>& Y = Y();
//     const volScalarField& W = this->W();

//     forAll(Y,i)
//     {

//     }


// }


inline Foam::label Foam::basicSpecieMixture::index
(
    const volScalarField& Yi
) const
{
    // return species_[Yi.name()];
    return species_[Yi.member()];
    // Info << "Print Yi.member().substr(1) " << Yi.member().substr(1) << endl;
    // return species_[Yi.member().substr(1)];
}
