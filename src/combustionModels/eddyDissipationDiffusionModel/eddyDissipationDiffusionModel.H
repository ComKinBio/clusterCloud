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
    Foam::edmCombustionModels::eddyDissipationDiffusionModel

Description
    Eddy dissipation model based on the principle of mixed is burnt.

    The time scale is determined by the minimum of the turbulence
    and diffusion time scales.

    Difference in the rtDiff calculation compared to the Fluent 
    EddyDissipationModel/finite-rate calculation

SourceFiles
    eddyDissipationDiffusionModel.C

\*---------------------------------------------------------------------------*/

#ifndef eddyDissipationDiffusionModel_H
#define eddyDissipationDiffusionModel_H

#include "singleStepCombustion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace edmCombustionModels
{

/*---------------------------------------------------------------------------*\
                   Class eddyDissipationDiffusionModel Declaration
\*---------------------------------------------------------------------------*/

class eddyDissipationDiffusionModel
:
    public singleStepCombustion
{
    // Private Data

        //- Model constant
        scalar A_;

        //- Model constant
        scalar B_;

        //- Model constant
        scalar C_;


public:

    //- Runtime type information
    TypeName("eddyDissipationDiffusionModel");


    // Constructors

        //- Construct from components
        eddyDissipationDiffusionModel
        (
            const word& modelType,
            const fluidReactionThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        eddyDissipationDiffusionModel(const eddyDissipationDiffusionModel&);


    //- Destructor
    virtual ~eddyDissipationDiffusionModel();


    // Member Functions

        //- Return the reciprocal of the turbulent mixing time scale
        tmp<volScalarField> rtDiff() const;

        //- Return the reciprocal of the turbulent mixing time scale
        tmp<volScalarField> rtTurb() const;

        //- Correct combustion rate
        virtual void correct();

        //- Update properties
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const eddyDissipationDiffusionModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace edmCombustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
