/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

template<class ParcelType>
template<class TrackCloudType>
inline Foam::ReactingClusterParcel<ParcelType>::trackingData::trackingData
(
    const TrackCloudType& cloud
)
:
    ParcelType::trackingData(cloud),
    epsilon_(cloud.epsilon()),
    k_(cloud.k()),
    epsilonInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.epsilon()
        )
    ),
    kInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.k()
        )
    ),
    epsilonc_(Zero),
    kc_(Zero)
{
}

template<class ParcelType>
inline const Foam::volScalarField&
Foam::ReactingClusterParcel<ParcelType>::trackingData::epsilon() const
{
    return epsilon_;
}


template<class ParcelType>
inline const Foam::volScalarField&
Foam::ReactingClusterParcel<ParcelType>::trackingData::k() const
{
    return k_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ReactingClusterParcel<ParcelType>::trackingData::epsilonInterp() const
{
    return epsilonInterp_();
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::ReactingClusterParcel<ParcelType>::trackingData::kInterp() const
{
    return kInterp_();
}



template<class ParcelType>
inline Foam::scalar Foam::ReactingClusterParcel<ParcelType>::trackingData::epsilonc() const
{
    return epsilonc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ReactingClusterParcel<ParcelType>::trackingData::epsilonc()
{
    return epsilonc_;
}


template<class ParcelType>
inline Foam::scalar Foam::ReactingClusterParcel<ParcelType>::trackingData::kc() const
{
    return kc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::ReactingClusterParcel<ParcelType>::trackingData::kc()
{
    return kc_;
}

// ************************************************************************* //
