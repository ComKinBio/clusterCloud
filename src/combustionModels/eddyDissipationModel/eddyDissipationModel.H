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
    Foam::edmCombustionModels::eddyDissipationModel

Description
    Eddy dissipation model based on the principle of mixed is burnt.

    The time scale is determined by the minimum of the turbulence
    and diffusion time scales.

    Difference in the rtDiff calculation compared to the Fluent 
    EddyDissipationModel/finite-rate calculation

SourceFiles
    eddyDissipationModel.C

\*---------------------------------------------------------------------------*/

#ifndef eddyDissipationModel_H
#define eddyDissipationModel_H

#include "singleStepCombustion.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace edmCombustionModels
{

/*---------------------------------------------------------------------------*\
                   Class eddyDissipationModel Declaration
\*---------------------------------------------------------------------------*/

class eddyDissipationModel
:
    public singleStepCombustion
{
    // Private Data

        //- Model constant
        scalar A_;

        //- Model constant
        scalar B_;

public:

    //- Runtime type information
    TypeName("eddyDissipationModel");


    // Constructors

        //- Construct from components
        eddyDissipationModel
        (
            const word& modelType,
            const fluidReactionThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        eddyDissipationModel(const eddyDissipationModel&);


    //- Destructor
    virtual ~eddyDissipationModel();


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
        void operator=(const eddyDissipationModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace edmCombustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
