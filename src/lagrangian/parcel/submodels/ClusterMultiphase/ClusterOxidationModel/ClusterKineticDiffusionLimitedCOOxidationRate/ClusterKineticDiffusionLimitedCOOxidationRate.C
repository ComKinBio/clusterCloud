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

#include "ClusterKineticDiffusionLimitedCOOxidationRate.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ClusterKineticDiffusionLimitedCOOxidationRate<CloudType>::
ClusterKineticDiffusionLimitedCOOxidationRate
(
    const dictionary& dict,
    CloudType& owner
)
:
    ClusterOxidationModel<CloudType>(dict, owner, typeName),
    C1_(this->coeffDict().template lookup<scalar>("C1")),
    eta_(this->coeffDict().template lookupOrDefault<scalar>("eta", 1.0)),
    etaKin_(this->coeffDict().template lookupOrDefault<scalar>("etaKin", 1.0)),
    A_(this->coeffDict().template lookup<scalar>("A")),
    beta_(this->coeffDict().template lookupOrDefault<scalar>("beta", 0.0)),
    E_(this->coeffDict().template lookup<scalar>("Ea")),
    Hf_(this->coeffDict().template lookupOrDefault<scalar>("Hf", 0.0)),
    cluster_(this->coeffDict().template lookupOrDefault<bool>("cluster", false)),
    CsLocalId_(-1),
    O2GlobalId_(owner.composition().carrierId("O2")),
    COGlobalId_(owner.composition().carrierId("CO")),
    WC_(0.0),
    WO2_(0.0),
    HcCO_(0.0)
{
    // Determine Cs ids
    label idSolid = owner.composition().idSolid();
    CsLocalId_ = owner.composition().localId(idSolid, "C");

    // Set local copies of thermo properties
    WO2_ = owner.composition().carrier().Wi(O2GlobalId_);
    const scalar WCO = owner.composition().carrier().Wi(COGlobalId_);
    WC_ = WCO - 0.5*WO2_;

    HcCO_ = owner.composition().carrier().Hf(COGlobalId_);

    const scalar YCloc = owner.composition().Y0(idSolid)[CsLocalId_];
    const scalar YSolidTot = owner.composition().YMixture0()[idSolid];
    Info<< "    C(s): particle mass fraction = " << YCloc*YSolidTot << endl;

    Info << " Heat of formation of C(s) = " << Hf_ << endl;
}


template<class CloudType>
Foam::ClusterKineticDiffusionLimitedCOOxidationRate<CloudType>::
ClusterKineticDiffusionLimitedCOOxidationRate
(
    const ClusterKineticDiffusionLimitedCOOxidationRate<CloudType>& srm
)
:
    ClusterOxidationModel<CloudType>(srm),
    C1_(srm.C1_),
    eta_(srm.eta_),
    etaKin_(srm.etaKin_),
    A_(srm.A_),
    beta_(srm.beta_),
    E_(srm.E_),
    Hf_(srm.Hf_),
    cluster_(srm.cluster_),
    CsLocalId_(srm.CsLocalId_),
    O2GlobalId_(srm.O2GlobalId_),
    COGlobalId_(srm.COGlobalId_),
    WC_(srm.WC_),
    WO2_(srm.WO2_),
    HcCO_(srm.HcCO_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ClusterKineticDiffusionLimitedCOOxidationRate<CloudType>::
~ClusterKineticDiffusionLimitedCOOxidationRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::ClusterKineticDiffusionLimitedCOOxidationRate<CloudType>::calculate
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

    // Local mass fraction of O2 in the carrier phase
    const scalar YO2 = carrier.Y(O2GlobalId_)[celli];

    // Quick exit if oxidant not present
    if (YO2 < rootVSmall)
    {
        return 0.0;
    }

    // Diffusion rate coefficient
    scalar alphaCluster_ = 1.0;
    if(cluster_)
    {
        alphaCluster_ = alphaCluster;
    }
    const scalar D0 = alphaCluster_*C1_/d*pow(0.5*(T + Tc), 0.75);

    // ClusterKinetic rate
    const scalar Rk = etaKin_*A_*pow(T, beta_)*exp(-E_/(RR*T));

    // Particle surface area
    const scalar Ap = constant::mathematical::pi*sqr(d);

    // Change in C mass [kg]
    scalar dmC = eta_*Ap*rhoc*RR*Tc*YO2/WO2_*D0*Rk/(D0 + Rk)*dt;

    // Limit mass transfer by availability of C
    dmC = min(mass*Ychar, dmC);

    // Molar consumption
    const scalar dOmega = dmC/WC_;

    // Change in O2 mass [kg]
    const scalar dmO2 = 0.5*dOmega*WO2_;

    // Mass of newly created CO [kg]
    const scalar dmCO = dOmega*(WC_ + 0.5*WO2_);

    // Update local particle C mass
    dMassSolid[CsLocalId_] += dOmega*WC_;

    // Update carrier O2 and CO mass
    dMassSRCarrier[O2GlobalId_] -= dmO2;
    dMassSRCarrier[COGlobalId_] += dmCO;

    const scalar HsC = thermo.solids().properties()[CsLocalId_].Hs(T);

    // carrier sensible enthalpy exchange handled via change in mass

    // Heat of reaction [J]
    return dmC*HsC + dmC*Hf_ - dmCO*HcCO_;
}


// ************************************************************************* //
