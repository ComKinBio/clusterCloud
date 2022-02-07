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

#include "clusterParcelInjectionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::scalarList&
Foam::clusterParcelInjectionData::YGas() const
{
    return YGas_;
}


inline const Foam::scalarList&
Foam::clusterParcelInjectionData::YLiquid() const
{
    return YLiquid_;
}


inline const Foam::scalarList&
Foam::clusterParcelInjectionData::YSolid() const
{
    return YSolid_;
}


inline Foam::scalarList& Foam::clusterParcelInjectionData::YGas()
{
    return YGas_;
}


inline Foam::scalarList& Foam::clusterParcelInjectionData::YLiquid()
{
    return YLiquid_;
}


inline Foam::scalarList& Foam::clusterParcelInjectionData::YSolid()
{
    return YSolid_;
}


// ************************************************************************* //
