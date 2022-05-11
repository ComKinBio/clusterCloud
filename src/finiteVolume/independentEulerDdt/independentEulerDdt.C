/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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

#include "independentEulerDdtScheme.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::word Foam::fv::independentEulerDdt::rDeltaTName("rDeltaTi");
Foam::word Foam::fv::independentEulerDdt::rDeltaTfName("rDeltaTif");
Foam::word Foam::fv::independentEulerDdt::rSubDeltaTName("rSubDeltaTi");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::fv::independentEulerDdt::enabled(const fvMesh& mesh)
{
    return
        word(mesh.ddtScheme("default"))
     == fv::independentEulerDdtScheme<scalar>::typeName;
}


const Foam::volScalarField& Foam::fv::independentEulerDdt::localRDeltaT
(
    const fvMesh& mesh
)
{
    return mesh.objectRegistry::lookupObject<volScalarField>
    (
        mesh.time().subCycling() ? rSubDeltaTName : rDeltaTName
    );
}


const Foam::surfaceScalarField& Foam::fv::independentEulerDdt::localRDeltaTf
(
    const fvMesh& mesh
)
{
    return mesh.objectRegistry::lookupObject<surfaceScalarField>
    (
        rDeltaTfName
    );
}


Foam::tmp<Foam::volScalarField> Foam::fv::independentEulerDdt::localRSubDeltaT
(
    const fvMesh& mesh,
    const label nAlphaSubCycles
)
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            rSubDeltaTName,
            nAlphaSubCycles
           *mesh.objectRegistry::lookupObject<volScalarField>
            (
                rDeltaTName
            )
        )
    );
}


// ************************************************************************* //
