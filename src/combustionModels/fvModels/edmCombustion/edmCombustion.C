/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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

#include "edmCombustion.H"
#include "fluidReactionThermo.H"
#include "fluidThermo.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(edmCombustion, 0);

    addToRunTimeSelectionTable
    (
        fvModel,
        edmCombustion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::edmCombustion::edmCombustion
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fvModel(sourceName, modelType, dict, mesh),
    phaseName_(dict.lookupOrDefault<word>("phase",word::null)),
    addSupFields_(),
    thermo_
    (
        mesh.lookupObject<fluidReactionThermo>
        (
            IOobject::groupName
            (
                basicThermo::dictName,
                phaseName_
            )
        )
    ),
    curTimeIndex_(-1)
{
    const compressibleMomentumTransportModel& turb =
        mesh.lookupObject<compressibleMomentumTransportModel>
        (
            IOobject::groupName
            (
                "momentumTransport",
                phaseName_
            )
        );

    if (thermo_.properties().found("defaultSpecie"))
    {
        word inertSpecie = thermo_.properties().lookup("defaultSpecie");

        addSupFields_.setSize(thermo_.composition().Y().size());

        label index = 0;
        forAll(thermo_.composition().Y(), specie)
        {
            if (inertSpecie != thermo_.composition().Y()[specie].member())
            {
                addSupFields_[index] = thermo_.composition().Y()[specie].name();
                index++;
            }
        }
    }
    else
    {
        addSupFields_.setSize(thermo_.composition().Y().size() + 1);
        forAll(thermo_.composition().Y(), specie)
        {
            addSupFields_[specie] = thermo_.composition().Y()[specie].name();
        }
    }

    addSupFields_[addSupFields_.size() - 1] = thermo_.he().name();

    edmCombustion_ = edmCombustionModel::New
    (
        thermo_,
        turb,
        word("edmCombustionProperties")
    ).ptr();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::edmCombustion::addSupFields() const
{
    return addSupFields_;
}


void Foam::fv::edmCombustion::correct()
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    Info<<"  Solving the fvCombustion model..."<<nl<<endl;
    edmCombustion_->correct();

    curTimeIndex_ = mesh().time().timeIndex();
}


void Foam::fv::edmCombustion::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (eqn.psi().name() == thermo_.he().name())
    {
        eqn += alpha*edmCombustion_->Qdot();
    }
    else
    {
        volScalarField& Yi = mesh().lookupObjectRef<volScalarField>
        (
            eqn.psi().name()
        );
        eqn += alpha*edmCombustion_->R(Yi);
    }
}


void Foam::fv::edmCombustion::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const word& fieldName
) const
{
    if (eqn.psi().name() == thermo_.he().name())
    {
        eqn += edmCombustion_->Qdot();
    }
    else
    {
        volScalarField& Yi = mesh().lookupObjectRef<volScalarField>
        (
            eqn.psi().name()
        );
        eqn += edmCombustion_->R(Yi);
    }
}


// ************************************************************************* //
