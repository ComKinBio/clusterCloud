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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::ReactingClusterCloud<CloudType>&
Foam::ReactingClusterCloud<CloudType>::cloudCopy() const
{
    return cloudCopyPtr_();
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::ReactingClusterCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::ReactingClusterCloud<CloudType>::constProps()
{
    return constProps_;
}

template<class CloudType>
inline Foam::volScalarField& 
Foam::ReactingClusterCloud<CloudType>::alphaf()
{
    return voidFraction_;
}

template<class CloudType>
inline Foam::volScalarField& 
Foam::ReactingClusterCloud<CloudType>::numberDens()
{
    return numDens_;
}

template<class CloudType>
inline Foam::volScalarField& 
Foam::ReactingClusterCloud<CloudType>::numberDensPrev()
{
    return numDensPrev_;
}

template<class CloudType>
inline Foam::volScalarField& 
Foam::ReactingClusterCloud<CloudType>::dAvg()
{
    return dAvg_;
}

template<class CloudType>
inline Foam::volScalarField& 
Foam::ReactingClusterCloud<CloudType>::dAvgPrev()
{
    return dAvgPrev_;
}

template<class CloudType>
inline const Foam::DevolatilisationModel
<
    Foam::ReactingClusterCloud<CloudType>
>&
Foam::ReactingClusterCloud<CloudType>::devolatilisation() const
{
    return devolatilisationModel_;
}


template<class CloudType>
inline Foam::DevolatilisationModel
<
    Foam::ReactingClusterCloud<CloudType>
>&
Foam::ReactingClusterCloud<CloudType>::devolatilisation()
{
    return devolatilisationModel_();
}


template<class CloudType>
inline const Foam::ClusterOxidationModel
<
    Foam::ReactingClusterCloud<CloudType>
>&
Foam::ReactingClusterCloud<CloudType>::clusterOxidation() const
{
    return clusterOxidationModel_;
}


template<class CloudType>
inline Foam::ClusterOxidationModel
<
    Foam::ReactingClusterCloud<CloudType>
>&
Foam::ReactingClusterCloud<CloudType>::clusterOxidation()
{
    return clusterOxidationModel_();
}

template<class CloudType>
inline const Foam::ClusterGasificationModel
<
    Foam::ReactingClusterCloud<CloudType>
>&
Foam::ReactingClusterCloud<CloudType>::clusterGasification() const
{
    return clusterGasificationModel_;
}


template<class CloudType>
inline Foam::ClusterGasificationModel
<
    Foam::ReactingClusterCloud<CloudType>
>&
Foam::ReactingClusterCloud<CloudType>::clusterGasification()
{
    return clusterGasificationModel_();
}

template<class CloudType>
inline const Foam::volScalarField& Foam::ReactingClusterCloud<CloudType>::epsilon() const
{
    return epsilon_;
}

template<class CloudType>
inline const Foam::volScalarField& Foam::ReactingClusterCloud<CloudType>::k() const
{
    return k_;
}



// ************************************************************************* //
