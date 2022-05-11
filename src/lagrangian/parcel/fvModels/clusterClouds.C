/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "clusterClouds.H"
#include "basicSpecieMixture.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(clusterClouds, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            clusterClouds,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::clusterClouds::clusterClouds
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(sourceName, modelType, dict, mesh),
    carrierThermo_
    (
        mesh.lookupObject<fluidThermo>
        (
            IOobject::groupName
            (
                basicThermo::dictName,
                dict.lookupOrDefault<word>("carrier",word::null)
            )
        )
    ),
    clusterClouds_
    (
        mesh.lookupObject<volScalarField>
        (
            carrierThermo_.phasePropertyName
            (
                word(dict.lookupOrDefault<word>("carrier",word::null))
             == word::null ?
                word("rho") : word("thermo:rho")
            )
        ),
        mesh.lookupObject<volVectorField>
        (
            carrierThermo_.phasePropertyName("U")
        ),
        mesh.lookupObject<uniformDimensionedVectorField>("g"),
        carrierThermo_
    ),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::clusterClouds::addSupFields() const
{
    word rhoName;
    if (this->mesh().foundObject<volScalarField>("rho"))
    {
        rhoName = "rho";
    }
    else
    {
        rhoName = "thermo:rho";
    }

    wordList fieldNames
    (
        {
            carrierThermo_.phasePropertyName(rhoName),
            carrierThermo_.phasePropertyName("U"),
            carrierThermo_.he().name()
        }
    );

    if (isA<basicSpecieMixture>(carrierThermo_))
    {
        const basicSpecieMixture& composition =
            refCast<const basicSpecieMixture>(carrierThermo_);

        const PtrList<volScalarField>& Y = composition.Y();

        forAll(Y, i)
        {
            if (composition.solve(i))
            {
                fieldNames.append(Y[i].name());
            }
        }
    }

    return fieldNames;
}


void Foam::fv::clusterClouds::correct()
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    clusterClouds_.evolve();

    curTimeIndex_ = mesh().time().timeIndex();
}


void Foam::fv::clusterClouds::addSup
(
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if
    (
        fieldName == carrierThermo_.phasePropertyName("thermo:rho")
     || fieldName == carrierThermo_.phasePropertyName("rho")
    )
    {
        eqn += clusterClouds_.Srho(eqn.psi());
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::clusterClouds::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if
    (
        fieldName == carrierThermo_.phasePropertyName("thermo:rho")
     || fieldName == carrierThermo_.phasePropertyName("rho")
    )
    {
        eqn += clusterClouds_.Srho(eqn.psi());
    }
    else if (fieldName == carrierThermo_.he().name())
    {
        eqn += clusterClouds_.Sh(eqn.psi());
    }
    else if
    (
        isA<basicSpecieMixture>(carrierThermo_)
     && refCast<const basicSpecieMixture>(carrierThermo_).contains
        (
            eqn.psi().name()
        )
    )
    {
        eqn += clusterClouds_.SYi
        (
            refCast<const basicSpecieMixture>(carrierThermo_).index(eqn.psi()),
            eqn.psi()
        );
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}

void Foam::fv::clusterClouds::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if
    (
        fieldName == carrierThermo_.phasePropertyName("thermo:rho")
     || fieldName == carrierThermo_.phasePropertyName("rho")
    )
    {
        eqn += clusterClouds_.Srho(eqn.psi());
    }
    else if (fieldName == carrierThermo_.he().name())
    {
        eqn += clusterClouds_.Sh(eqn.psi());
    }
    else if
    (
        isA<basicSpecieMixture>(carrierThermo_)
    && refCast<const basicSpecieMixture>(carrierThermo_).contains
        (
            eqn.psi().member()
        )
    )
    {
        volScalarField tY(eqn.psi());
        tY.rename(tY.member());
        eqn += clusterClouds_.SYi
        (
            refCast<const basicSpecieMixture>(carrierThermo_).index(tY),
            eqn.psi()
        );
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::clusterClouds::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (fieldName == carrierThermo_.phasePropertyName("U"))
    {
        eqn += clusterClouds_.SU(eqn.psi());
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}

void Foam::fv::clusterClouds::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const word& fieldName
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (fieldName == carrierThermo_.phasePropertyName("U"))
    {
        eqn += clusterClouds_.SU(eqn.psi());
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << fieldName << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::clusterClouds::preUpdateMesh()
{
    // Store the particle positions
    clusterClouds_.storeGlobalPositions();
}


// ************************************************************************* //
