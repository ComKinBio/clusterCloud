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

#include "reactingClusterCloud.H"

#include "makeParcelCloudFunctionObjects.H"
//#include "makeAdditionalParcelCloudFunctionObjects.H"

// Momentum
#include "makeParcelForces.H"
//#include "makeAdditionalParcelForces.H"
#include "makeParcelDispersionModels.H"
#include "makeReactingMultiphaseParcelInjectionModels.H" // MP variant
#include "makeParcelPatchInteractionModels.H"
#include "makeReactingMultiphaseParcelStochasticCollisionModels.H" // MP variant
#include "makeReactingParcelSurfaceFilmModels.H" // Reacting variant

// Thermodynamic
#include "makeParcelHeatTransferModels.H"
//#include "makeAdditionalParcelHeatTransferModels.H"
//#include "makeAdditionalParcelTurbulenceDispersionModels.H"

// Reacting
#include "makeReactingMultiphaseParcelCompositionModels.H" // MP Variant
#include "makeReactingParcelPhaseChangeModels.H"

// Coal multiphase
#include "makeReactingClusterParcelDevolatilisationModels.H"
#include "makeReactingClusterParcelOxidationModels.H"
#include "makeReactingClusterParcelGasificationModels.H"

// parcelTurbulence
#include "makeParcelTurbulenceDispersionModels.H"
#include "makeThermoParcelTurbulenceForces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeParcelCloudFunctionObjects(reactingClusterCloud);
//makeAdditionalParcelCloudFunctionObjects(reactingClusterCloud);

// Momentum sub-models
makeParcelForces(reactingClusterCloud);
//makeAdditionalParcelForces(reactingClusterCloud);
makeParcelDispersionModels(reactingClusterCloud);
makeReactingMultiphaseParcelInjectionModels(reactingClusterCloud);
makeParcelPatchInteractionModels(reactingClusterCloud);
makeReactingMultiphaseParcelStochasticCollisionModels
(
    reactingClusterCloud
);
makeReactingParcelSurfaceFilmModels(reactingClusterCloud);
//makeAdditionalParcelTurbulenceDispersionModels(reactingClusterCloud);

// Thermo sub-models
makeParcelHeatTransferModels(reactingClusterCloud);
//makeAdditionalParcelHeatTransferModels(reactingClusterCloud);

// Reacting sub-models
makeReactingMultiphaseParcelCompositionModels
(
    reactingClusterCloud
);
makeReactingParcelPhaseChangeModels(reactingClusterCloud);

// Reacting multiphase sub-models
makeReactingClusterParcelDevolatilisationModels
(
    reactingClusterCloud
);
makeReactingClusterParcelOxidationModels(reactingClusterCloud);
makeReactingClusterParcelGasificationModels(reactingClusterCloud);

// parcel turbulence
makeThermoParcelTurbulenceForces(reactingClusterCloud);
makeParcelTurbulenceDispersionModels(reactingClusterCloud);

// ************************************************************************* //
