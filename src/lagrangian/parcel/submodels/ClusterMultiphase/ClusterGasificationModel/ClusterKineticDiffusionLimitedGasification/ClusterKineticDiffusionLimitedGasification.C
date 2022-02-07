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

#include "ClusterKineticDiffusionLimitedGasification.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ClusterKineticDiffusionLimitedGasification<CloudType>::
ClusterKineticDiffusionLimitedGasification
(
    const dictionary& dict,
    CloudType& owner
)
:
    ClusterGasificationModel<CloudType>(dict, owner, typeName),
    clusterGasificationData_(3),
    CsLocalId_(-1),
    CH4GlobalId_(owner.composition().carrierId("CH4")),
    CO2GlobalId_(owner.composition().carrierId("CO2")),
    COGlobalId_(owner.composition().carrierId("CO")),
    H2OGlobalId_(owner.composition().carrierId("H2O")),
    H2GlobalId_(owner.composition().carrierId("H2")),
    WC_(0.0),
    WCH4_(0.0),
    WCO2_(0.0),
    WCO_(0.0),
    WH2O_(0.0),
    WH2_(0.0),
    HcCH4_(0.0),
    HcCO2_(0.0),
    HcCO_(0.0),
    HcH2O_(0.0),
    HcH2_(0.0)
{
    // Load clusterGasification data
    clusterGasificationData_[0]=clusterGasificationData(this->coeffDict().subDict("CO2"));
    clusterGasificationData_[1]=clusterGasificationData(this->coeffDict().subDict("H2O"));
    clusterGasificationData_[2]=clusterGasificationData(this->coeffDict().subDict("H2"));

    // Determine Cs ids
    label idSolid = owner.composition().idSolid();
    CsLocalId_ = owner.composition().localId(idSolid, "C");

    // Set local copies of thermo properties
    // Molar masses
    WCH4_ = owner.composition().carrier().Wi(CH4GlobalId_);
    WCO2_ = owner.composition().carrier().Wi(CO2GlobalId_);
    WCO_ = owner.composition().carrier().Wi(COGlobalId_);
    WH2O_ = owner.composition().carrier().Wi(H2OGlobalId_);
    WH2_ = owner.composition().carrier().Wi(H2GlobalId_);
    WC_ = WCH4_ - 2.0*WH2_;

    // Formation enthalpies
    HcCH4_ = owner.composition().carrier().Hf(CH4GlobalId_);
    HcCO2_ = owner.composition().carrier().Hf(CO2GlobalId_);
    HcCO_ = owner.composition().carrier().Hf(COGlobalId_);
    HcH2O_ = owner.composition().carrier().Hf(H2OGlobalId_);
    HcH2_ = owner.composition().carrier().Hf(H2GlobalId_);

    Info<<"    gasification rate CO2: C1="<<clusterGasificationData_[0].C1()
        <<", eta="<<clusterGasificationData_[0].eta()
        <<", etaKin="<<clusterGasificationData_[0].etaKin()
        <<", A="<<clusterGasificationData_[0].A()
        <<", E="<<clusterGasificationData_[0].E()<<nl
        <<"    gasification rate H2O: C1="<<clusterGasificationData_[1].C1()
        <<", eta="<<clusterGasificationData_[1].eta()
        <<", etaKin="<<clusterGasificationData_[1].etaKin()
        <<", A="<<clusterGasificationData_[1].A()
        <<", E="<<clusterGasificationData_[1].E()<<nl
        <<"    gasification rate H2:  C1="<<clusterGasificationData_[2].C1()
        <<", eta="<<clusterGasificationData_[2].eta()
        <<", etaKin="<<clusterGasificationData_[2].etaKin()
        <<", A="<<clusterGasificationData_[2].A()
        <<", E="<<clusterGasificationData_[2].E()<<nl
        <<endl;

    const scalar YCloc = owner.composition().Y0(idSolid)[CsLocalId_];
    const scalar YSolidTot = owner.composition().YMixture0()[idSolid];
    Info<< "    C(s): particle mass fraction = " << YCloc*YSolidTot << endl;
}


template<class CloudType>
Foam::ClusterKineticDiffusionLimitedGasification<CloudType>::
ClusterKineticDiffusionLimitedGasification
(
    const ClusterKineticDiffusionLimitedGasification<CloudType>& srm
)
:
    ClusterGasificationModel<CloudType>(srm),
    clusterGasificationData_(srm.clusterGasificationData_),
    CsLocalId_(srm.CsLocalId_),
    CH4GlobalId_(srm.CH4GlobalId_),
    CO2GlobalId_(srm.CO2GlobalId_),
    COGlobalId_(srm.COGlobalId_),
    H2OGlobalId_(srm.H2OGlobalId_),
    H2GlobalId_(srm.H2GlobalId_),
    WC_(srm.WC_),
    WCH4_(srm.WCH4_),
    WCO2_(srm.WCO2_),
    WCO_(srm.WCO_),
    WH2O_(srm.WH2O_),
    WH2_(srm.WH2_),
    HcCH4_(srm.HcCH4_),
    HcCO2_(srm.HcCO2_),
    HcCO_(srm.HcCO_),
    HcH2O_(srm.HcH2O_),
    HcH2_(srm.HcH2_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ClusterKineticDiffusionLimitedGasification<CloudType>::
~ClusterKineticDiffusionLimitedGasification()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::ClusterKineticDiffusionLimitedGasification<CloudType>::calculate
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
    scalarField& dMassGRCarrier
) const
{

    // Fraction of remaining combustible material
    const label idSolid = CloudType::parcelType::SLD;
    const scalar Ychar = YMixture[idSolid]*YSolid[CsLocalId_];

    // Surface combustion active combustible fraction is consumed
    if (Ychar < small)
    {
        return 0.0;
    }

    const parcelThermo& thermo = this->owner().thermo();
    const basicSpecieMixture& carrier = this->owner().composition().carrier();

    // Local gas mass fractions in the carrier phase
    const scalar YCO2 = carrier.Y(CO2GlobalId_)[celli];
    const scalar YH2O = carrier.Y(H2OGlobalId_)[celli];
    const scalar YH2 = carrier.Y(H2GlobalId_)[celli];

    // Particle surface area
    const scalar Ap = constant::mathematical::pi*sqr(d);

    // clusterGasification rates
    scalarField dOmega(clusterGasificationData_.size(), 0.0);

    // Molar consumption due to CO2
    if (YCO2 > rootVSmall and clusterGasificationData_[0].active())
    {
        // Diffusion rate coefficient
        scalar alphaCluster_ = 1.0;
        if(clusterGasificationData_[0].cluster())
        {
            alphaCluster_ = alphaCluster;
        }
        scalar D = alphaCluster_*clusterGasificationData_[0].C1()/d*pow(0.5*(T + Tc), 0.75);

        // ClusterKinetic rate
        scalar Rk = clusterGasificationData_[0].etaKin()
                   *clusterGasificationData_[0].A()
                   *exp(-clusterGasificationData_[0].E()/(RR*Tc));

        // Molar consumption
        // Molar consumption
        dOmega[0] = Ap*clusterGasificationData_[0].eta()
                   *D*Rk/(D+Rk)*dt
                   *rhoc*RR*Tc*YCO2/WCO2_/WC_;
    }

    if (YH2O > rootVSmall and clusterGasificationData_[1].active())
    {
        // Molar consumption due to H2O
        // Diffusion rate coefficient
        scalar alphaCluster_ = 1.0;
        if(clusterGasificationData_[1].cluster())
        {
            alphaCluster_ = alphaCluster;
        }
        scalar D = alphaCluster_*clusterGasificationData_[1].C1()/d*pow(0.5*(T + Tc), 0.75);


        // ClusterKinetic rate
        scalar Rk = clusterGasificationData_[1].etaKin()
                   *clusterGasificationData_[1].A()
                   *exp(-clusterGasificationData_[1].E()/(RR*Tc));

        // Molar consumption
        // Molar consumption
        dOmega[1] = Ap*clusterGasificationData_[1].eta()
                   *D*Rk/(D+Rk)*dt
                   *rhoc*RR*Tc*YH2O/WH2O_/WC_;
    }

    if (YH2 > rootVSmall and clusterGasificationData_[2].active())
    {
        // Molar consumption due to H2        
        // Diffusion rate coefficient
        scalar alphaCluster_ = 1.0;
        if(clusterGasificationData_[2].cluster())
        {
            alphaCluster_ = alphaCluster;
        }
        scalar D = alphaCluster_*clusterGasificationData_[2].C1()/d*pow(0.5*(T + Tc), 0.75);

        // ClusterKinetic rate
        scalar Rk = clusterGasificationData_[2].etaKin()
                   *clusterGasificationData_[2].A()
                   *exp(-clusterGasificationData_[2].E()/(RR*Tc));

        // Molar consumption
        dOmega[2] = Ap*clusterGasificationData_[2].eta()
                   *D*Rk/(D+Rk)*dt
                   *rhoc*RR*Tc*YH2/WH2_/WC_;
    }

    // Limit mass transfer by availability of C
    if (sum(dOmega)*WC_ > mass*Ychar)
    {
        dOmega *= mass*Ychar/(sum(dOmega)*WC_);
    }
    scalar dmC = sum(dOmega)*WC_;


    // Change due to CO2 clusterGasification
    // C(s) + CO2 -> 2CO
    scalar dmCO2 = -dOmega[0]*WCO2_;
    scalar dmCO  = 2.0*dOmega[0]*WCO_;

    // Change due to H2O clusterGasification
    // C(s) + H2O -> CO + H2O
    scalar dmH2O = -dOmega[1]*WH2O_;
           dmCO += dOmega[1]*WCO_;
    scalar dmH2  = dOmega[1]*WH2_;

    // Change due to H2 clusterGasification
    // C(s) + 2H2 -> CH4
           dmH2 -= 2.0*dOmega[2]*WH2_;
    scalar dmCH4 = dOmega[2]*WCH4_;

    // Update local particle C mass
    dMassSolid[CsLocalId_]       += dmC;

    // Update the carrier phase
    dMassGRCarrier[CO2GlobalId_] += dmCO2;
    dMassGRCarrier[COGlobalId_]  += dmCO;
    dMassGRCarrier[H2OGlobalId_] += dmH2O;
    dMassGRCarrier[H2GlobalId_]  += dmH2;
    dMassGRCarrier[CH4GlobalId_] += dmCH4;

    // carrier sensible enthalpy exchange handled via change in mass
    const scalar HsC = thermo.solids().properties()[CsLocalId_].Hs(T);

    // Heat of reaction [J]
    return dmC*HsC 
         - dmCO2*HcCO2_ 
         - dmCO*HcCO_
         - dmH2O*HcH2O_
         - dmH2*HcH2_
         - dmCH4*HcCH4_;
}


// ************************************************************************* //
