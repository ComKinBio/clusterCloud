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

#ifndef makeClusterDevolatilisationModel_H
#define makeClusterDevolatilisationModel_H

#include "addToRunTimeSelectionTable.H"
#include "DevolatilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeClusterDevolatilisationModel(CloudType)                               \
                                                                               \
    typedef Foam::CloudType::reactingClusterCloudType                           \
        reactingClusterCloudType;                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::DevolatilisationModel<reactingClusterCloudType>,                  \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            DevolatilisationModel<reactingClusterCloudType>,                    \
            dictionary                                                         \
        );                                                                     \
    }


#define makeClusterDevolatilisationModelType(SS, CloudType)                       \
                                                                               \
    typedef Foam::CloudType::reactingClusterCloudType                           \
        reactingClusterCloudType;                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::SS<reactingClusterCloudType>, 0);                                \
                                                                               \
    Foam::DevolatilisationModel<reactingClusterCloudType>::                     \
        adddictionaryConstructorToTable                                        \
        <Foam::SS<reactingClusterCloudType>>                                    \
        add##SS##CloudType##reactingClusterCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
