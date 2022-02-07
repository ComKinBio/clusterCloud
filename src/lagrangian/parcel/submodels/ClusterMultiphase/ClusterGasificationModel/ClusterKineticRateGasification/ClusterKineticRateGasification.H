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

Class
    Foam::ClusterKineticRateGasification

Description
    Kinetic rate oxidation reaction model for coal parcels.
    Reaction order is limited to unity for all educts and products.
    Limited to:

        C(s) + O2 -> CO2

    Effective reaction rate is calculated according to:

        R = eta_*A_P*p_O2*A*exp(-Ea/R/T)

    Usage:

        kineticRateGasificationCoeffs
        {
            CO2
            {
                A               497; // kg/m^2/s
                Ea        7.101e+07; // J/kmol
            }
            H2O
            {
                A               497; // kg/m^2/s
                Ea        7.101e+07; // J/kmol
            }
            H2
            {
                A               497; // kg/m^2/s
                Ea        7.101e+07; // J/kmol
            }
        }
    

\*---------------------------------------------------------------------------*/

#ifndef ClusterKineticRateGasification_H
#define ClusterKineticRateGasification_H

#include "ClusterGasificationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
template<class CloudType>
class ClusterKineticRateGasification;

/*---------------------------------------------------------------------------*\
            Class ClusterKineticRateGasification Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ClusterKineticRateGasification
:
    public ClusterGasificationModel<CloudType>
{
    // Private data

        // Model constants

            //- Kinetics limited rate pre-exponential constant, A
            scalarField A_;

            //- Kinetics limited rate activation energy
            scalarField E_;

            //- Effectiveness factor
            scalarField eta_;

            List<bool> active_;         


        // Addressing

            //- Cs positions in global/local lists
            label CsLocalId_;

            //- CH4 positions in global list
            label CH4GlobalId_;

            //- CO2 positions in global list
            label CO2GlobalId_;

            //- CO positions in global list
            label COGlobalId_;

            //- H2O position in global list
            label H2OGlobalId_;

            //- H2 positions in global list
            label H2GlobalId_;


        // Local copies of thermo properties

        // Molar masses


            //- Molecular weight of C [kg/kmol]
            scalar WC_;

            //- Molecular weight of H2 [kg/kmol]
            scalar WCH4_;

            //- Molecular weight of CO2 [kg/kmol]
            scalar WCO2_;

            //- Molecular weight of CO [kg/kmol]
            scalar WCO_;

            //- Molecular weight of H2O [kg/kmol]
            scalar WH2O_;

            //- Molecular weight of H2 [kg/kmol]
            scalar WH2_;


        // Formation enthalpies

            //- Formation enthalpy for CH4 [J/kg]
            scalar HcCH4_;

            //- Formation enthalpy for CO2 [J/kg]
            scalar HcCO2_;

            //- Formation enthalpy for CO [J/kg]
            scalar HcCO_;

            //- Formation enthalpy for H2O [J/kg]
            scalar HcH2O_;

            //- Formation enthalpy for H2 [J/kg]
            scalar HcH2_;


public:

    //- Runtime type information
    TypeName("clusterKineticRateGasification");


    // Constructors

        //- Construct from dictionary
        ClusterKineticRateGasification
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        ClusterKineticRateGasification
        (
            const ClusterKineticRateGasification<CloudType>& srm
        );

        //- Construct and return a clone
        virtual autoPtr<ClusterGasificationModel<CloudType>> clone() const
        {
            return autoPtr<ClusterGasificationModel<CloudType>>
            (
                new ClusterKineticRateGasification<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ClusterKineticRateGasification();


    // Member Functions

        //- Update surface reactions
        virtual scalar calculate
        (
            const scalar dt,
            const label celli,
            const scalar d,
            const scalar T,
            const scalar Tc,
            const scalar pc,
            const scalar rhoc,
            const scalar mass,
            const scalarField& YGas,
            const scalarField& YLiquid,
            const scalarField& YSolid,
            const scalarField& YMixture,
            const scalar N,
            scalarField& dMassGas,
            scalarField& dMassLiquid,
            scalarField& dMassSolid,
            scalarField& dMassGRCarrier
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ClusterKineticRateGasification.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
