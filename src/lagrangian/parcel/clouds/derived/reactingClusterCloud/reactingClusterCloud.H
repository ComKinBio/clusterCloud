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

Class
    Foam::reactingClusterCloud

Description
    Cloud class to introduce multi-phase reacting parcels

\*---------------------------------------------------------------------------*/

#ifndef reactingClusterCloud_H
#define reactingClusterCloud_H

#include "ParcelCloudBase.H"
#include "MomentumCloud.H"
#include "ThermoCloud.H"
#include "ReactingCloud.H"
#include "ReactingClusterCloud.H"
#include "ParcelCloud.H"
#include "reactingClusterParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef ParcelCloud
    <
        ReactingClusterCloud
        <
            ReactingCloud
            <
                ThermoCloud
                <
                    MomentumCloud
                    <
                        ParcelCloudBase
                        <
                            reactingClusterParcel
                        >
                    >
                >
            >
        >
    > reactingClusterCloud;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
