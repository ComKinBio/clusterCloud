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

Class
    Foam::fv::edmCombustion

Description
    Calculates and applies the edmCombustion source to the energy equation.

Usage
    Example usage:
    \verbatim
    edmCombustionCoeffs
    {
        libs    ("libedmCombustionModels.so");
        type    edmCombustion;
    }
    \endverbatim

SourceFiles
    edmCombustion.C

\*---------------------------------------------------------------------------*/

#ifndef edmCombustion_H
#define edmCombustion_H

#include "fvModel.H"
#include "fluidReactionThermo.H"
#include "uniformDimensionedFields.H"
#include "edmCombustionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                          Class edmCombustion Declaration
\*---------------------------------------------------------------------------*/

class edmCombustion
:
    public fvModel
{
    // Private Data

        //- Phase name
        word phaseName_;

        //- The edmCombustion model pointer
        mutable autoPtr<edmCombustionModel> edmCombustion_;

        //- Stores the fvModel affected fields
        wordList addSupFields_;

        //- Reference to fluidReactionThermo
        const fluidReactionThermo& thermo_;

        //- Current time index (used for updating)
        mutable label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("edmCombustion");


    // Constructors

        //- Construct from explicit source name and mesh
        edmCombustion
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        edmCombustion(const edmCombustion&) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Correct

            //- Solve the edmCombustion chemistry and update the sources
            virtual void correct();


        // Evaluate

            //- Add explicit contribution to compressible energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;

            //- Add explicit contribution to phase momentum equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const word& fieldName
            ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const edmCombustion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
