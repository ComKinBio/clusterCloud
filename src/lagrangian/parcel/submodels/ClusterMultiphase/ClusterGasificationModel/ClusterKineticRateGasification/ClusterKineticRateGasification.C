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

#include "ClusterKineticRateGasification.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ClusterKineticRateGasification<CloudType>::
ClusterKineticRateGasification
(
    const dictionary& dict,
    CloudType& owner
)
:
    ClusterGasificationModel<CloudType>(dict, owner, typeName),
    A_(3,0.0),
    E_(3,0.0),
    eta_(3,1.0),
    active_(3,true),
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

    // Read in kinetic parameters
    A_[0] = readScalar(this->coeffDict().subDict("CO2").lookup("A"));
    E_[0] = readScalar(this->coeffDict().subDict("CO2").lookup("E"));
    eta_[0] = readScalar(this->coeffDict().subDict("CO2").lookup("eta"));
    active_[0] = readBool(this->coeffDict().subDict("CO2").lookup("active"));

    A_[1] = readScalar(this->coeffDict().subDict("H2O").lookup("A"));
    E_[1] = readScalar(this->coeffDict().subDict("H2O").lookup("E"));
    eta_[1] = readScalar(this->coeffDict().subDict("H2O").lookup("eta"));
    active_[1] = readBool(this->coeffDict().subDict("H2O").lookup("active"));

    A_[2] = readScalar(this->coeffDict().subDict("H2").lookup("A"));
    E_[2] = readScalar(this->coeffDict().subDict("H2").lookup("E"));
    eta_[2] = readScalar(this->coeffDict().subDict("H2").lookup("eta"));
    active_[2] = readBool(this->coeffDict().subDict("H2").lookup("active"));


    Info<<"    gasification rate CO2: A="<<A_[0]<<", E="<<E_[0]<<", eta="<<eta_[0]<<nl
        <<"    gasification rate H2O: A="<<A_[1]<<", E="<<E_[1]<<", eta="<<eta_[1]<<nl
        <<"    gasification rate H2:  A="<<A_[2]<<", E="<<E_[2]<<", eta="<<eta_[2]
        <<endl;


    const scalar YCloc = owner.composition().Y0(idSolid)[CsLocalId_];
    const scalar YSolidTot = owner.composition().YMixture0()[idSolid];
    Info<< "    C(s): particle mass fraction = " << YCloc*YSolidTot << endl;
}


template<class CloudType>
Foam::ClusterKineticRateGasification<CloudType>::
ClusterKineticRateGasification
(
    const ClusterKineticRateGasification<CloudType>& srm
)
:
    ClusterGasificationModel<CloudType>(srm),
    A_(srm.A_),
    E_(srm.E_),
    eta_(srm.eta_),
    active_(srm.active_),
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
Foam::ClusterKineticRateGasification<CloudType>::
~ClusterKineticRateGasification()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::ClusterKineticRateGasification<CloudType>::calculate
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
    const scalar alphaTilde,
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

    // gasification rates
    scalarField dOmega(3, 0.0);

    // Molar consumption
    if (YCO2 > rootVSmall and active_[0])
    {
        dOmega[0] = eta_[0]*A_[0]*exp(-E_[0]/(RR*Tc))*Ap*rhoc*RR*Tc*YCO2/WCO2_*dt/WC_;
    }
    if (YH2O > rootVSmall and active_[1])
    {
        dOmega[1] = eta_[1]*A_[1]*exp(-E_[1]/(RR*Tc))*Ap*rhoc*RR*Tc*YH2O/WH2O_*dt/WC_;
    }
    if (YH2 > rootVSmall and active_[2])
    {
        dOmega[2] = eta_[2]*A_[2]*exp(-E_[2]/(RR*Tc))*Ap*rhoc*RR*Tc*YH2/WH2_*dt/WC_;
    }

    // Limit mass transfer by availability of C
    if (sum(dOmega)*WC_ > mass*Ychar)
    {
        dOmega *= mass*Ychar/(sum(dOmega)*WC_);
    }
    scalar dmC = sum(dOmega)*WC_;


    // Change due to CO2 gasification
    // C(s) + CO2 -> 2CO
    scalar dmCO2 = -dOmega[0]*WCO2_;
    scalar dmCO  = 2.0*dOmega[0]*WCO_;

    // Change due to H2O gasification
    // C(s) + H2O -> CO + H2O
    scalar dmH2O = -dOmega[1]*WH2O_;
           dmCO += dOmega[1]*WCO_;
    scalar dmH2  = dOmega[1]*WH2_;

    // Change due to H2 gasification
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
