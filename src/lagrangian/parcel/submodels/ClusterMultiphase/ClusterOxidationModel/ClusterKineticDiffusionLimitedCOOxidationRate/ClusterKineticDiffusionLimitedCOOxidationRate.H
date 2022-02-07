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
    Foam::KineticDiffusionLimitedOxidationRate

Description
    Kinetic/diffusion limited rate surface reaction model for coal parcels.
    Reaction order is limited to unity for all educts and products.
    Limited to:

        C(s) + O2 -> CO

    Effective reaction rate is calculated according to:

        R_eff = A_P*eta*p_O2*R_kin*D_O2/(R_kin + D_O2)

        R_kin = etaKin*A*exp(-Ea/R/T)

        D_O2 = alphaTilde*C1/d_P*((T_P+T_c)/2)^0.75

        alphaTilde - clustering effect - based on 

        @article{Kruger2017,
        author = {Kr{\"{u}}ger, Jonas and Haugen, Nils Erland L. and L{\o}v{\aa}s, Terese},
        doi = {10.1016/j.combustflame.2017.07.008},
        journal = {Combustion and Flame},
        keywords = {Char,Direct numerical simulation,Heterogeneous combustion,Particle clustering,Reacting multiphase flow,Turbulence},
        pages = {160--172},
        title = {{Correlation effects between turbulence and the conversion rate of pulverized char particles}},
        volume = {185},
        year = {2017}
        }

        @article{Haugen2018,
        abstract = {The effect of turbulence on the mass transfer between a fluid and embedded small heavy inertial particles that experience surface reactions is studied. For simplicity, the surface reaction, which takes place when a gas phase reactant is converted to a gas phase product at the external surface of the particles, is unimolar and isothermal. Two effects are identified. The first effect is due to the relative velocity between the fluid and the particles, and a model for the relative velocity is presented. The second effect is due to the clustering of particles, where the mass transfer rate is inhibited due to the rapid depletion of the consumed species inside the dense particle clusters. This last effect is relevant for large Damk{\"{o}}hler numbers, where the Damk{\"{o}}hler number is defined as the ratio of the turbulent and chemical time scales, and it may totally control the mass transfer rate for Damk{\"{o}}hler numbers larger than unity. A model that describes how this effect should be incorporated into existing simulation tools that utilize the Reynolds averaged Navier-Stokes approach is presented.},
        author = {Haugen, Nils Erland L. and Kr{\"{u}}ger, Jonas and Mitra, Dhrubaditya and L{\o}v{\aa}s, Terese},
        doi = {10.1017/jfm.2017.820},
        journal = {Journal of Fluid Mechanics},
        pages = {932--951},
        title = {{The effect of turbulence on mass transfer rates of small inertial particles with surface reactions}},
        volume = {836},
        year = {2018}
        }

    Usage:
        ClusterKineticDiffusionLimitedCOOxidationRateCoeffs
        {
            C1                      1; // m^2/s (?) diffusion coefficient?
            A                     497; // kg/m^2/s
            Ea              7.101e+07; // J/kmol
            eta                     1; // default 1 
            etaKin                  1; // default 1 
            beta                    0; // default 0, temperature exponent
            Hf                      0; // default 0, heat of formation for char (C)
            cluster              true; // default false, include alphaTilde or not
        }    

\*---------------------------------------------------------------------------*/

#ifndef ClusterKineticDiffusionLimitedCOOxidationRate_H
#define ClusterKineticDiffusionLimitedCOOxidationRate_H

#include "ClusterOxidationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
template<class CloudType>
class ClusterKineticDiffusionLimitedCOOxidationRate;

/*---------------------------------------------------------------------------*\
            Class ClusterKineticDiffusionLimitedCOOxidationRate Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ClusterKineticDiffusionLimitedCOOxidationRate
:
    public ClusterOxidationModel<CloudType>
{
    // Private data

        // Model constants

            //- Mass diffusion limited rate constant, C1
            const scalar C1_;

            //- Effectiveness factor
            const scalar eta_;

            //- Effectiveness factor for kinetic contribution
            const scalar etaKin_;

            //- ClusterKinetics limited rate pre-exponential constant, A
            const scalar A_;

            //- ClusterKinetics limited rate temperature dependency
            const scalar beta_;

            //- ClusterKinetics limited rate activation energy
            const scalar E_;

            //- heat of formation for char
            const scalar Hf_;

            //- activate or deactivate the cluster effects 
            const bool cluster_;

        // Addressing

            //- Cs positions in global/local lists
            label CsLocalId_;

            //- O2 position in global list
            label O2GlobalId_;

            //- CO positions in global list
            label COGlobalId_;


        // Local copies of thermo properties

            //- Molecular weight of C [kg/kmol]
            scalar WC_;

            //- Molecular weight of O2 [kg/kmol]
            scalar WO2_;

            //- Formation enthalpy for CO [J/kg]
            scalar HcCO_;


public:

    //- Runtime type information
    TypeName("clusterKineticDiffusionLimitedCOOxidationRate");


    // Constructors

        //- Construct from dictionary
        ClusterKineticDiffusionLimitedCOOxidationRate
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        ClusterKineticDiffusionLimitedCOOxidationRate
        (
            const ClusterKineticDiffusionLimitedCOOxidationRate<CloudType>& srm
        );

        //- Construct and return a clone
        virtual autoPtr<ClusterOxidationModel<CloudType>> clone() const
        {
            return autoPtr<ClusterOxidationModel<CloudType>>
            (
                new ClusterKineticDiffusionLimitedCOOxidationRate<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ClusterKineticDiffusionLimitedCOOxidationRate();


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
            const scalar alphaCluster,
            scalarField& dMassGas,
            scalarField& dMassLiquid,
            scalarField& dMassSolid,
            scalarField& dMassSRCarrier
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ClusterKineticDiffusionLimitedCOOxidationRate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
