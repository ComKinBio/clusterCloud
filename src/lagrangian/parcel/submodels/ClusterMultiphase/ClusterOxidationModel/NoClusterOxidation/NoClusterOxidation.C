/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "NoClusterOxidation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoClusterOxidation<CloudType>::NoClusterOxidation
(
    const dictionary&,
    CloudType& owner
)
:
    ClusterOxidationModel<CloudType>(owner)
{}


template<class CloudType>
Foam::NoClusterOxidation<CloudType>::NoClusterOxidation
(
    const NoClusterOxidation<CloudType>& srm
)
:
    ClusterOxidationModel<CloudType>(srm.owner_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoClusterOxidation<CloudType>::~NoClusterOxidation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::NoClusterOxidation<CloudType>::active() const
{
    return false;
}


template<class CloudType>
Foam::scalar Foam::NoClusterOxidation<CloudType>::calculate
(
    const scalar,
    const label,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalarField&,
    const scalarField&,
    const scalarField&,
    const scalarField&,
    const scalar,
    scalarField&,
    scalarField&,
    scalarField&,
    scalarField&
) const
{
    // Model does not stop gasification taking place
    return 0;
}


// ************************************************************************* //
