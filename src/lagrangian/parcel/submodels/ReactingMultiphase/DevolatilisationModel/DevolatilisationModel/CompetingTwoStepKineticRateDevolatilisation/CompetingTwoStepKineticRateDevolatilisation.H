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
    Foam::CompetingTwoStepKineticRateDevolatilisation

Description
    Single kinetic rate devolatisation model.
    - acts on a per-specie basis
    - Rate given by Arrhenius eqn

        kappa = alpha1.A1.exp(- E1/R.T) + alpha2.A2.exp(- E2/R.T)

    Where:
        kappa = rate constant
        A1    = pre-exponential factor (user input)
        E1    = activation energy (user input)
        A2    = pre-exponential factor (user input)
        E2    = activation energy (user input)
        R     = universal gas constant
        T     = temperature

    Usage:

        CompetingTwoStepKineticRateDevolatilisationCoeffs
        {
            volatileData
            (
             // (name A1 E1 alpha1 A2 E2 alpha2)
                (CH4     12     0.5    0.34    20   0.8    0.42)
            );

            volatileResidualCoeff 1e-6;
        }


\*---------------------------------------------------------------------------*/

#ifndef CompetingTwoStepKineticRateDevolatilisation_H
#define CompetingTwoStepKineticRateDevolatilisation_H

#include "DevolatilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
              Class CompetingTwoStepKineticRateDevolatilisation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CompetingTwoStepKineticRateDevolatilisation
:
    public DevolatilisationModel<CloudType>
{
    // Helper class to store specie-local volatile data
    class volatileData
    {
        // Private Data

            //- Specie name
            word name_;

            //- Activation energy
            scalar A1_;

            //- Pre-exponential factor
            scalar E1_;

            //- Weighting factor
            scalar alpha1_;

            //- Activation energy
            scalar A2_;

            //- Pre-exponential factor
            scalar E2_;

            //- Weighting factor
            scalar alpha2_;


    public:

        // Constructors

            //- Null constructor
            volatileData()
            :
                name_(word::null),
                A1_(0.0),
                E1_(0.0),
                alpha1_(0.0),
                A2_(0.0),
                E2_(0.0),
                alpha2_(0.0)
            {}

            //- Construct from Istream
            volatileData(Istream& is)
            :
                name_(is),
                A1_(readScalar(is)),
                E1_(readScalar(is)),
                alpha1_(readScalar(is)),
                A2_(readScalar(is)),
                E2_(readScalar(is)),
                alpha2_(readScalar(is))
            {}


        //- Destructor
        ~volatileData()
        {}


        // Public Member Functions

            // Access

                //- Return const access to the name
                const word& name() const
                {
                    return name_;
                }

                //- Return const access to the activation energy
                scalar A1() const
                {
                    return A1_;
                }

                //- Return const access to the pre-exponential factor
                scalar E1() const
                {
                    return E1_;
                }

                //- Return const access to the weighting factor
                scalar alpha1() const
                {
                    return alpha1_;
                }

                //- Return const access to the activation energy
                scalar A2() const
                {
                    return A2_;
                }

                //- Return const access to the pre-exponential factor
                scalar E2() const
                {
                    return E2_;
                }

                //- Return const access to the weighting factor
                scalar alpha2() const
                {
                    return alpha2_;
                }


        // IOstream Operators

            //- Read from Istream
            friend Istream& operator>>(Istream& is, volatileData& vd)
            {
                is.readBeginList("volatileData");
                is  >>  vd.name_  >> vd.A1_  >> vd.E1_ >> vd.alpha1_ 
                                  >> vd.A2_  >> vd.E2_ >> vd.alpha2_;
                is.readEndList("volatileData");

                return is;
            }

            //- Write to Ostream
            friend Ostream& operator<<(Ostream& os, const volatileData& vd)
            {
                os  << token::BEGIN_LIST
                    << vd.name_ << token::SPACE
                    << vd.A1_ << token::SPACE
                    << vd.E1_ << token::SPACE
                    << vd.alpha1_ << token::SPACE
                    << vd.A2_ << token::SPACE
                    << vd.E2_ << token::SPACE
                    << vd.alpha2_ << token::SPACE
                    << token::END_LIST;

                return os;
            }
    };


    // Private Data

        // Model constants

            //- List of volatile data - (name A1 E)
            List<volatileData> volatileData_;

            //- List of initial volatile mass fractions
            List<scalar> YVolatile0_;

            //- Mapping between local and cloud gaseous species
            List<label> volatileToGasMap_;

            //- Volatile residual coefficient (0-1)
            //  When the fraction of volatiles are depleted below this
            //  threshold, combustion can occur
            const scalar residualCoeff_;


public:

    //- Runtime type information
    TypeName("competingTwoStepKineticRateDevolatilisation");


    // Constructors

        //- Construct from dictionary
        CompetingTwoStepKineticRateDevolatilisation
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct and return a clone
        virtual autoPtr<DevolatilisationModel<CloudType>> clone() const
        {
            return autoPtr<DevolatilisationModel<CloudType>>
            (
                new CompetingTwoStepKineticRateDevolatilisation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~CompetingTwoStepKineticRateDevolatilisation();


    // Member Functions

        //- Update model
        virtual void calculate
        (
            const scalar dt,
            const scalar age,
            const scalar mass0,
            const scalar mass,
            const scalar T,
            const scalarField& YGasEff,
            const scalarField& YLiquidEff,
            const scalarField& YSolidEff,
            label& canCombust,
            scalarField& dMassDV
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CompetingTwoStepKineticRateDevolatilisation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
