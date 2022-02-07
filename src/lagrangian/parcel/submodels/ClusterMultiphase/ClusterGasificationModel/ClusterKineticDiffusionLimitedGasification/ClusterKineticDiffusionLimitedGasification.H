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
    Foam::ClusterKineticDiffusionLimitedGasification

Description
    ClusterKinetic/diffusion limited rate surface reaction model for coal parcels.
    Reaction order is limited to unity for all educts and products.
    Limited to:

        (1) C(s) + CO2 -> 2CO

        (2) C(s) + H2O -> CO + H2

        (3) C(s) + 2H2 -> CH4

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
        clusterKineticDiffusionLimitedGasificationCoeffs
        {
            CO2
            {
                C1           1e-5; // [m2/s]
                eta           0.7; // [] 
                etaKin          1; // [] 
                A           20230; // [m/s/K]
                E         3.304e8; // [J/mol]
                cluster      true; // switch on/off clustering effect
                active       true; // switch on/off this gasification reaction
            }
            H2O
            {
                C1           1e-5; // [m2/s]
                eta           0.7; // [] 
                etaKin          1; // [] 
                A           606.9; // [m/s/K]
                E         2.697e8; // [J/mol]
                cluster      true; // switch on/off clustering effect
                active       true; // switch on/off this gasification reaction
            }
            H2
            {
                C1           1e-5; // [m2/s]
                eta           0.7; // [] - default 1
                etaKin          1; // [] - default 1
                A               0; // [m/s/K]
                E         3.304e8; // [J/mol]
                cluster      true; // switch on/off clustering effect
                active       true; // switch on/off this gasification reaction
            }
        }    

\*---------------------------------------------------------------------------*/

#ifndef ClusterKineticDiffusionLimitedGasification_H
#define ClusterKineticDiffusionLimitedGasification_H

#include "ClusterGasificationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
template<class CloudType>
class ClusterKineticDiffusionLimitedGasification;


/*---------------------------------------------------------------------------*\
            Class ClusterKineticDiffusionLimitedGasification Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ClusterKineticDiffusionLimitedGasification
:
    public ClusterGasificationModel<CloudType>
{
    // Helper class to store specie-local gasification data
    class clusterGasificationData
    {

        // Private data

            //- Gasification agent name
            word gasifAgent_;

            //- Mass diffusion limited rate constant, C1
            scalar C1_;

            //- Effectiveness factor
            scalar eta_;

            //- Effectiveness factor only for kinetic contribution
            scalar etaKin_;

            //- Activation energy
            scalar A_;

            //- Pre-exponential factor
            scalar E_;

            bool active_;

            //- cluster switch
            bool cluster_;

    public:

        // Constructors

            //- Null constructor
            clusterGasificationData()
            :
                gasifAgent_(word::null),
                C1_(0.0),
                eta_(0.0),
                etaKin_(0.0),
                A_(0.0),
                E_(0.0),
                active_(false),
                cluster_(false)
            {}

            //- Construct from dictionary
            clusterGasificationData(const dictionary& dict)
            :
                gasifAgent_(dict.name()),
                C1_(readScalar(dict.lookup("C1"))),
                eta_(readScalar(dict.lookup("eta"))),
                etaKin_(readScalar(dict.lookup("etaKin"))),
                A_(readScalar(dict.lookup("A"))),
                E_(readScalar(dict.lookup("E"))),
                active_(readBool(dict.lookup("active"))),
                cluster_(readBool(dict.lookup("cluster")))
            {} 

            //- Construct from Istream
            clusterGasificationData(Istream& is)
            :
                gasifAgent_(is),
                C1_(readScalar(is)),
                eta_(readScalar(is)),
                etaKin_(readScalar(is)),
                A_(readScalar(is)),
                E_(readScalar(is)),
                active_(readBool(is)),
                cluster_(readBool(is))
            {}

            //- Construct as copy
            clusterGasificationData(const clusterGasificationData& gd)
            :
                gasifAgent_(gd.gasifAgent_),
                C1_(gd.C1_),
                eta_(gd.eta_),
                etaKin_(gd.etaKin_),
                A_(gd.A_),
                E_(gd.E_),
                active_(gd.active_),
                cluster_(gd.cluster_)
            {}


        //- Destructor
        ~clusterGasificationData()
        {}


        // Public Member Functions

            // Access

                //- Return const access to the gasification agent name
                word gasifAgent() const
                {
                    return gasifAgent_;
                }

                //- Return const access to the diffusion constant
                scalar C1() const
                {
                    return C1_;
                }

                //- Return const access to the effectiveness factor
                scalar eta() const
                {
                    return eta_;
                }

                //- Return const access to the effectiveness factor for kinetics
                scalar etaKin() const
                {
                    return etaKin_;
                }

                //- Return const access to the activation energy
                scalar A() const
                {
                    return A_;
                }

                //- Return const access to the pre-exponential factor
                scalar E() const
                {
                    return E_;
                }

                //- Return const access to the pre-exponential factor
                bool active() const
                {
                    return active_;
                }

                //- Return const access to the pre-exponential factor
                bool cluster() const
                {
                    return cluster_;
                }


        // IOstream Operators

            //- Read from Istream
            friend Istream& operator>>(Istream& is, clusterGasificationData& gd)
            {
                is.readBeginList("clusterGasificationData");
                is  >> gd.gasifAgent_  >> gd.C1_  >> gd.eta_  >> gd.A_  >> gd.E_;
                is.readEndList("clusterGasificationData");

                return is;
            }
            //- Write to Ostream
            friend Ostream& operator<<(Ostream& os, const clusterGasificationData& gd)
            {
                os  << token::BEGIN_LIST
                    << gd.gasifAgent_ << token::SPACE
                    << gd.C1_ << token::SPACE
                    << gd.eta_ << token::SPACE
                    << gd.etaKin_ << token::SPACE
                    << gd.A_ << token::SPACE
                    << gd.E_
                    << token::END_LIST;

                return os;
            }
    };


    // Private data

            //- List of gasification data - (C1 eta etaKin A E)
            List<clusterGasificationData> clusterGasificationData_;

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

            //- Mapping between local and cloud gaseous species
            //List<label> volatileToGasMap_;

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
    TypeName("clusterKineticDiffusionLimitedGasification");

    // Constructors

        //- Construct from dictionary
        ClusterKineticDiffusionLimitedGasification
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        ClusterKineticDiffusionLimitedGasification
        (
            const ClusterKineticDiffusionLimitedGasification<CloudType>& srm
        );

        //- Construct and return a clone
        virtual autoPtr<ClusterGasificationModel<CloudType>> clone() const
        {
            return autoPtr<ClusterGasificationModel<CloudType>>
            (
                new ClusterKineticDiffusionLimitedGasification<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ClusterKineticDiffusionLimitedGasification();


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
    #include "ClusterKineticDiffusionLimitedGasification.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
