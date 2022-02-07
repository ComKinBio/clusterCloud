/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::ReactingClusterParcel<ParcelType>::constantProperties::
constantProperties()
:
    ParcelType::constantProperties(),
    TDevol_(this->dict_, 0.0),
    LDevol_(this->dict_, 0.0),
    hRetentionCoeff_(this->dict_, 0.0),
    diffCoeff0_(this->dict_, 0.0),
    diffCoeff1_(this->dict_, 0.0),
    diffCoeff2_(this->dict_, 0.0),
    constCp_(this->dict_, false)
{}


template<class ParcelType>
inline Foam::ReactingClusterParcel<ParcelType>::constantProperties::
constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    TDevol_(cp.TDevol_),
    LDevol_(cp.LDevol_),
    hRetentionCoeff_(cp.hRetentionCoeff_),
    diffCoeff0_(cp.diffCoeff0_),
    diffCoeff1_(cp.diffCoeff1_),
    diffCoeff2_(cp.diffCoeff2_),
    constCp_(cp.constCp_)
{}


template<class ParcelType>
inline Foam::ReactingClusterParcel<ParcelType>::constantProperties::
constantProperties
(
    const dictionary& parentDict
)
:
    ParcelType::constantProperties(parentDict),
    TDevol_(this->dict_, "TDevol"),
    LDevol_(this->dict_, "LDevol"),
    hRetentionCoeff_(this->dict_, "hRetentionCoeff"),
    diffCoeff0_(this->dict_, "diffCoeff0"),
    diffCoeff1_(this->dict_, "diffCoeff1", 0.0),
    diffCoeff2_(this->dict_, "diffCoeff2", 0.0),
    constCp_(this->dict_, "constantCp", false)
{}


template<class ParcelType>
inline Foam::ReactingClusterParcel<ParcelType>::ReactingClusterParcel
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti
)
:
    ParcelType(mesh, coordinates, celli, tetFacei, tetPti),
    YGas_(0),
    YLiquid_(0),
    YSolid_(0),
    canCombust_(0),
    alphaTilde_(0)
{}


template<class ParcelType>
inline Foam::ReactingClusterParcel<ParcelType>::ReactingClusterParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    ParcelType(mesh, position, celli),
    YGas_(0),
    YLiquid_(0),
    YSolid_(0),
    canCombust_(0),
    alphaTilde_(0)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar
Foam::ReactingClusterParcel<ParcelType>::constantProperties::TDevol() const
{
    return TDevol_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::ReactingClusterParcel<ParcelType>::constantProperties::LDevol() const
{
    return LDevol_.value();
}


template<class ParcelType>
inline Foam::scalar
Foam::ReactingClusterParcel<ParcelType>::constantProperties::
hRetentionCoeff() const
{
    scalar value = hRetentionCoeff_.value();

    if ((value < 0) || (value > 1))
    {
        FatalErrorInFunction
            << "hRetentionCoeff must be in the range 0 to 1" << nl
            << exit(FatalError) << endl;
    }

    return value;
}

template<class ParcelType>
inline Foam::scalar
Foam::ReactingClusterParcel<ParcelType>::constantProperties::
diffCoeff0() const
{
    scalar value = diffCoeff0_.value();

    return value;
}
template<class ParcelType>
inline Foam::scalar
Foam::ReactingClusterParcel<ParcelType>::constantProperties::
diffCoeff1() const
{
    scalar value = diffCoeff1_.value();

    return value;
}
template<class ParcelType>
inline Foam::scalar
Foam::ReactingClusterParcel<ParcelType>::constantProperties::
diffCoeff2() const
{
    scalar value = diffCoeff2_.value();

    return value;
}
template<class ParcelType>
inline bool
Foam::ReactingClusterParcel<ParcelType>::constantProperties::
constCp() const
{
    return constCp_.value();
}


// * * * * * * * * * * ThermoParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingClusterParcel<ParcelType>::
YGas() const
{
    return YGas_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingClusterParcel<ParcelType>::
YLiquid() const
{
    return YLiquid_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::ReactingClusterParcel<ParcelType>::
YSolid() const
{
    return YSolid_;
}


template<class ParcelType>
inline Foam::label
Foam::ReactingClusterParcel<ParcelType>::canCombust() const
{
    return canCombust_;
}

template<class ParcelType>
inline const Foam::scalar&
Foam::ReactingClusterParcel<ParcelType>::alphaTilde() const
{
    return alphaTilde_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::ReactingClusterParcel<ParcelType>::YGas()
{
    return YGas_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::ReactingClusterParcel<ParcelType>::YLiquid()
{
    return YLiquid_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::ReactingClusterParcel<ParcelType>::YSolid()
{
    return YSolid_;
}


template<class ParcelType>
inline Foam::label& Foam::ReactingClusterParcel<ParcelType>::canCombust()
{
    return canCombust_;
}

template<class ParcelType>
inline Foam::scalar&
Foam::ReactingClusterParcel<ParcelType>::alphaTilde()
{
    return alphaTilde_;
}


// ************************************************************************* //
