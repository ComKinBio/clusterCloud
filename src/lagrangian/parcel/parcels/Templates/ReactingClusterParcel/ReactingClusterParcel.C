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

#include "ReactingClusterParcel.H"
#include "CompositionModel.H"
#include "NoDevolatilisation.H"
#include "NoClusterOxidation.H"
#include "NoClusterGasification.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
const Foam::label Foam::ReactingClusterParcel<ParcelType>::GAS(0);

template<class ParcelType>
const Foam::label Foam::ReactingClusterParcel<ParcelType>::LIQ(1);

template<class ParcelType>
const Foam::label Foam::ReactingClusterParcel<ParcelType>::SLD(2);


// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::ReactingClusterParcel<ParcelType>::CpEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    // effective cp for particle can be switched to not be calculated based 
    // on gas phase thermodynamics properties but set constant
    scalar CpEff0 = cloud.constProps().Cp0(); 
    if(!cloud.constProps().constCp())
    {
        CpEff0 = this->Y_[GAS]*cloud.composition().Cp(idG, YGas_, p, T)
               + this->Y_[LIQ]*cloud.composition().Cp(idL, YLiquid_, p, T)
               + this->Y_[SLD]*cloud.composition().Cp(idS, YSolid_, p, T);
    }
    return CpEff0;        
}


template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::ReactingClusterParcel<ParcelType>::HsEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*cloud.composition().Hs(idG, YGas_, p, T)
      + this->Y_[LIQ]*cloud.composition().Hs(idL, YLiquid_, p, T)
      + this->Y_[SLD]*cloud.composition().Hs(idS, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackCloudType>
Foam::scalar Foam::ReactingClusterParcel<ParcelType>::LEff
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*cloud.composition().L(idG, YGas_, p, T)
      + this->Y_[LIQ]*cloud.composition().L(idL, YLiquid_, p, T)
      + this->Y_[SLD]*cloud.composition().L(idS, YSolid_, p, T);
}


template<class ParcelType>
Foam::scalar Foam::ReactingClusterParcel<ParcelType>::updateMassFractions
(
    const scalar mass0,
    const scalarField& dMassGas,
    const scalarField& dMassLiquid,
    const scalarField& dMassSolid
)
{
    scalarField& YMix = this->Y_;

    scalar massGas =
        this->updateMassFraction(mass0*YMix[GAS], dMassGas, YGas_);
    scalar massLiquid =
        this->updateMassFraction(mass0*YMix[LIQ], dMassLiquid, YLiquid_);
    scalar massSolid =
        this->updateMassFraction(mass0*YMix[SLD], dMassSolid, YSolid_);

    scalar massNew = max(massGas + massLiquid + massSolid, rootVSmall);

    YMix[GAS] = massGas/massNew;
    YMix[LIQ] = massLiquid/massNew;
    YMix[SLD] = 1.0 - YMix[GAS] - YMix[LIQ];

    return massNew;
}


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingClusterParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);
    tetIndices tetIs = this->currentTetIndices();
    td.epsilonc() = td.epsilonInterp().interpolate(this->coordinates(), tetIs);
    td.kc() = td.kInterp().interpolate(this->coordinates(), tetIs);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingClusterParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Re-use correction from reacting parcel
    ParcelType::cellValueSourceCorrection(cloud, td, dt);
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingClusterParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    typedef typename TrackCloudType::thermoCloudType thermoCloudType;
    const CompositionModel<thermoCloudType>& composition =
        cloud.composition();


    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const scalar np0 = this->nParticle_;
    const scalar d0 = this->d_;
    const vector& U0 = this->U_;
    const scalar T0 = this->T_;
    const scalar mass0 = this->mass();

    const scalar pc = td.pc();

    const scalarField& YMix = this->Y_;
    const label idG = composition.idGas();
    const label idL = composition.idLiquid();
    const label idS = composition.idSolid();


    // Calc surface values
    scalar Ts, rhos, mus, Prs, kappas;
    this->calcSurfaceValues(cloud, td, T0, Ts, rhos, mus, Prs, kappas);
    scalar Res = this->Re(rhos, U0, td.Uc(), d0, mus);


    // Sources
    //~~~~~~~~

    // Explicit momentum source for particle
    vector Su = Zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = Zero;

    // Explicit enthalpy source for particle
    scalar Sh = 0.0;

    // Linearised enthalpy source coefficient
    scalar Sph = 0.0;

    // Sensible enthalpy transfer from the particle to the carrier phase
    scalar dhsTrans = 0.0;


    // 1. Compute models that contribute to mass transfer - U, T held constant
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Phase change in liquid phase
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Mass transfer due to phase change
    scalarField dMassPC(YLiquid_.size(), 0.0);

    // Molar flux of species emitted from the particle (kmol/m^2/s)
    scalar Ne = 0.0;

    // Sum Ni*Cpi*Wi of emission species
    scalar NCpW = 0.0;

    // Surface concentrations of emitted species
    scalarField Cs(composition.carrier().species().size(), 0.0);

    // Calc mass and enthalpy transfer due to phase change
    this->calcPhaseChange
    (
        cloud,
        td,
        dt,
        Res,
        Prs,
        Ts,
        mus/rhos,
        d0,
        T0,
        mass0,
        idL,
        YMix[LIQ],
        YLiquid_,
        dMassPC,
        Sh,
        Ne,
        NCpW,
        Cs
    );


    // Devolatilisation
    // ~~~~~~~~~~~~~~~~

    // Mass transfer due to devolatilisation
    scalarField dMassDV(YGas_.size(), 0.0);

    // Calc mass and enthalpy transfer due to devolatilisation
    calcDevolatilisation
    (
        cloud,
        td,
        dt,
        this->age_,
        Ts,
        d0,
        T0,
        mass0,
        this->mass0_,
        YMix[GAS]*YGas_,
        YMix[LIQ]*YLiquid_,
        YMix[SLD]*YSolid_,
        canCombust_,
        dMassDV,
        Sh,
        Ne,
        NCpW,
        Cs
    );


    // Oxidation reactions
    // ~~~~~~~~~~~~~~~~~

    // Change in carrier phase composition due to oxidation reactions
    scalarField dMassORGas(YGas_.size(), 0.0);
    scalarField dMassORLiquid(YLiquid_.size(), 0.0);
    scalarField dMassORSolid(YSolid_.size(), 0.0);
    scalarField dMassORCarrier(composition.carrier().species().size(), 0.0);

    scalar ShOR = 0.0;
    scalar dhsTransOR = 0.0;

    // calculate the cluster modified reaction rate alphatilde

    //diffusion coefficient for oxygen based on given polynomials from dictionary
    //surface temperature used - based on particle and gas temperature
    const scalar DO2 = cloud.constProps().diffCoeff0() 
                      +cloud.constProps().diffCoeff1()*Ts
                      +cloud.constProps().diffCoeff2()*sqr(Ts); 
    
    // Schmidt number
    const scalar Sc = (td.muc() / td.rhoc()) / DO2; 

    // integral time scale
    const scalar tauL = 2.0/3.0*td.kc()/max(td.epsilonc(), 1e-15); 

    // particle response time
    const scalar tauP = this->rho_*sqr(d0)/(18* td.muc()); 

    // Stokes number
    const scalar St = tauP / max(tauL, 1e-10);

    // Particle Reynolds number based on Haugen's formulation
    const scalar wnumint = 11.54 * td.epsilonc() / max(pow(td.kc(), 1.5), 1e-15);
    const scalar wnumint_pow = pow(max(wnumint, 1e-15), -0.667);
    const scalar wnumkol = 6.283 * pow025(max(td.epsilonc(), 1e-15)) / pow(td.muc()/td.rhoc(), 0.75);
    const scalar wnumkol_pow = pow(max(wnumkol, 1e-15), -0.667);
    const scalar urms = 0.8165 * sqrt(max(td.kc(), 1e-15));
    const scalar limnumerator = max(1e-15, (St*wnumint_pow - wnumkol_pow));
    const scalar limdenominator = max(1e-15, (wnumint_pow - wnumkol_pow));
    const scalar urelative = min(0.41*urms, 0.41*urms*sqrt(limnumerator/limdenominator));

    const scalar partrey = urelative * d0 / (td.muc()/td.rhoc());

    //Sherwood number
    const scalar Sherwood = 2.0 + 0.69 * sqrt(partrey) * pow(Sc, 1/3);

    //particle surface area - limit dAvg to 1 mum
    const scalar Ap = pi * sqr(max(cloud.dAvgPrev()[this->cell()], 1e-6));

    scalar np(0.0);

    //particle number density
    if (cloud.solution().transient())
    {
        //this approximation is only valid for fixed value particle diameters
        np = 6 / (pi * pow3(d0)) * cloud.alphaf()[this->cell()] / (1 - cloud.alphaf()[this->cell()]); 
    }
    else
    {
        np = cloud.numberDensPrev()[this->cell()];
    }

    //alphahom = 1/tauC; Sherwood_hom = 2.0 used here
    const scalar alphahom = (np * Ap * 2.0 * DO2 / max(cloud.dAvgPrev()[this->cell()], 1e-6));

    // Damkoehler number
    const scalar Da = tauL * alphahom;

    // alpha cluster
    const scalar alphaCluster = (0.08 + St/3.0) * Sherwood / max((tauL * St), 1e-10);

    // effect of turbulence and particle clustering
    alphaTilde_ = Sherwood / 2.0 * alphaCluster * tauL / max((alphaCluster * tauL + 0.5 * Da * Sherwood), 1e-10);

    // hand alphatilde over to oxidation calculation
    // Calc mass and enthalpy transfer due to oxidation reactions
    calcClusterOxidationReactions
    (
        cloud,
        td,
        dt,
        d0,
        T0,
        mass0,
        canCombust_,
        alphaTilde_,
        YMix,
        YGas_,
        YLiquid_,
        YSolid_,
        dMassORGas,
        dMassORLiquid,
        dMassORSolid,
        dMassORCarrier,
        ShOR,
        dhsTransOR
    );

    // Gasification reactions
    // ~~~~~~~~~~~~~~~~~

    // Change in carrier phase composition due to gasification reactions
    scalarField dMassGRGas(YGas_.size(), 0.0);
    scalarField dMassGRLiquid(YLiquid_.size(), 0.0);
    scalarField dMassGRSolid(YSolid_.size(), 0.0);
    scalarField dMassGRCarrier(composition.carrier().species().size(), 0.0);
    scalar ShGR = 0.0;
    scalar dhsTransGR = 0.0;

    // same alpha Tilde used for the gasification reactions - could be improved
    // to use a different diffusion coefficient
    // Calc mass and enthalpy transfer due to gasification reactions
    calcClusterGasificationReactions
    (
        cloud,
        td,
        dt,
        d0,
        T0,
        mass0,
        canCombust_,
        alphaTilde_,
        YMix,
        YGas_,
        YLiquid_,
        YSolid_,
        dMassGRGas,
        dMassGRLiquid,
        dMassGRSolid,
        dMassGRCarrier,
        ShGR, 
        dhsTransGR 
    );

    // check total solid mass consumption and scale if necessary
    label CsLocalId = composition.localId(idS, "C");
    scalar mCmax = YMix[idS]*YSolid_[CsLocalId]*mass0;

    if 
    ( 
        dMassGRSolid[CsLocalId] 
      + dMassORSolid[CsLocalId] 
      > mCmax 
      && 
        mCmax > rootVSmall
    )
    {
        scalar scale =
            mCmax/(dMassGRSolid[CsLocalId] + dMassORSolid[CsLocalId]);
        dMassORGas *= scale;
        dMassORLiquid *= scale;
        dMassORSolid *= scale;
        dMassORCarrier *= scale;
        ShOR *= scale;
        dhsTransOR *= scale;

        dMassGRGas *= scale;
        dMassGRLiquid *= scale;
        dMassGRSolid *= scale;
        dMassGRCarrier *= scale;
        ShGR *= scale;
        dhsTransGR *= scale;
    }

    Sh += ShOR + ShGR;
    dhsTrans += dhsTransOR + dhsTransGR;


    // 2. Update the parcel properties due to change in mass
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalarField dMassGas(dMassDV + dMassORGas + dMassGRGas);
    scalarField dMassLiquid(dMassPC + dMassORLiquid + dMassGRLiquid);
    scalarField dMassSolid(dMassORSolid + dMassGRSolid);
    scalar mass1 =
        updateMassFractions(mass0, dMassGas, dMassLiquid, dMassSolid);

    this->Cp_ = CpEff(cloud, td, pc, T0, idG, idL, idS);

    // Update particle density or diameter
    if (cloud.constProps().constantVolume())
    {
        this->rho_ = mass1/this->volume();
    }
    else
    {
        this->d_ = cbrt(mass1/this->rho_*6.0/pi);
    }

    // Remove the particle when mass falls below minimum threshold
    if (np0*mass1 < cloud.constProps().minParcelMass())
    {
        td.keepParticle = false;

        if (cloud.solution().coupled())
        {
            scalar dm = np0*mass0;

            // Absorb parcel into carrier phase
            forAll(YGas_, i)
            {
                label gid = composition.localToCarrierId(GAS, i);
                cloud.rhoTrans(gid)[this->cell()] += dm*YMix[GAS]*YGas_[i];
            }
            forAll(YLiquid_, i)
            {
                label gid = composition.localToCarrierId(LIQ, i);
                cloud.rhoTrans(gid)[this->cell()] += dm*YMix[LIQ]*YLiquid_[i];
            }

            // No mapping between solid components and carrier phase
            /*
            forAll(YSolid_, i)
            {
                label gid = composition.localToCarrierId(SLD, i);
                cloud.rhoTrans(gid)[this->cell()] += dm*YMix[SLD]*YSolid_[i];
            }
            */

            cloud.UTransRef()[this->cell()] += dm*U0;

            cloud.hsTransRef()[this->cell()] +=
                dm*HsEff(cloud, td, pc, T0, idG, idL, idS);

            cloud.phaseChange().addToPhaseChangeMass(np0*mass1);
        }

        return;
    }

    // Correct surface values due to emitted species
    this->correctSurfaceValues(cloud, td, Ts, Cs, rhos, mus, Prs, kappas);
    Res = this->Re(rhos, U0, td.Uc(), this->d_, mus);


    // 3. Compute heat- and momentum transfers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Heat transfer
    // ~~~~~~~~~~~~~

    // Calculate new particle temperature
    this->T_ =
        this->calcHeatTransfer
        (
            cloud,
            td,
            dt,
            Res,
            Prs,
            kappas,
            NCpW,
            Sh,
            dhsTrans,
            Sph
        );


    this->Cp_ = CpEff(cloud, td, pc, this->T_, idG, idL, idS);


    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    this->U_ =
        this->calcVelocity(cloud, td, dt, Res, mus, mass1, Su, dUTrans, Spu);


    // 4. Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (cloud.solution().coupled())
    {
        // Transfer mass lost to carrier mass, momentum and enthalpy sources
        forAll(YGas_, i)
        {
            scalar dm = np0*dMassGas[i];
            label gid = composition.localToCarrierId(GAS, i);
            scalar hs = composition.carrier().Hs(gid, pc, T0);
            cloud.rhoTrans(gid)[this->cell()] += dm;
            cloud.UTransRef()[this->cell()] += dm*U0;
            cloud.hsTransRef()[this->cell()] += dm*hs;
        }
        forAll(YLiquid_, i)
        {
            scalar dm = np0*dMassLiquid[i];
            label gid = composition.localToCarrierId(LIQ, i);
            scalar hs = composition.carrier().Hs(gid, pc, T0);
            cloud.rhoTrans(gid)[this->cell()] += dm;
            cloud.UTransRef()[this->cell()] += dm*U0;
            cloud.hsTransRef()[this->cell()] += dm*hs;
        }

        // No mapping between solid components and carrier phase
        /*
        forAll(YSolid_, i)
        {
            scalar dm = np0*dMassSolid[i];
            label gid = composition.localToCarrierId(SLD, i);
            scalar hs = composition.carrier().Hs(gid, pc, T0);
            cloud.rhoTrans(gid)[this->cell()] += dm;
            cloud.UTransRef()[this->cell()] += dm*U0;
            cloud.hsTransRef()[this->cell()] += dm*hs;
        }
        */

        // Accumulate all carrier gas changes
        scalarField dMassCarrier(dMassORCarrier + dMassGRCarrier);

        forAll(dMassCarrier, i)
        {
            scalar dm = np0*dMassCarrier[i];
            scalar hs = composition.carrier().Hs(i, pc, T0);
            cloud.rhoTrans(i)[this->cell()] += dm;
            cloud.UTransRef()[this->cell()] += dm*U0;
            cloud.hsTransRef()[this->cell()] += dm*hs;
        }

        // Update momentum transfer
        cloud.UTransRef()[this->cell()] += np0*dUTrans;
        cloud.UCoeffRef()[this->cell()] += np0*Spu;

        // Update sensible enthalpy transfer
        cloud.hsTransRef()[this->cell()] += np0*dhsTrans;
        cloud.hsCoeffRef()[this->cell()] += np0*Sph;

        // Update radiation fields
        if (cloud.radiation())
        {
            const scalar ap = this->areaP();
            const scalar T4 = pow4(T0);
            cloud.radAreaP()[this->cell()] += dt*np0*ap;
            cloud.radT4()[this->cell()] += dt*np0*T4;
            cloud.radAreaPT4()[this->cell()] += dt*np0*ap*T4;
        }

        if (!cloud.solution().transient())
        {
            cloud.alphaf()[this->cell()] += dt*this->nParticle()*this->volume();
            cloud.numberDens()[this->cell()] += dt*this->nParticle();
            cloud.dAvg()[this->cell()] += dt*this->d_*this->nParticle();
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingClusterParcel<ParcelType>::calcDevolatilisation
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar age,
    const scalar Ts,
    const scalar d,
    const scalar T,
    const scalar mass,
    const scalar mass0,
    const scalarField& YGasEff,
    const scalarField& YLiquidEff,
    const scalarField& YSolidEff,
    label& canCombust,
    scalarField& dMassDV,
    scalar& Sh,
    scalar& N,
    scalar& NCpW,
    scalarField& Cs
) const
{
    // Check that model is active
    if
    (
        isType
        <
            NoDevolatilisation
            <
                typename TrackCloudType::reactingClusterCloudType
            >
        >
        (
            cloud.devolatilisation()
        )
    )
    {
        if (canCombust != -1)
        {
            canCombust = 1;
        }
        return;
    }

    // Initialise demand-driven constants
    (void)cloud.constProps().TDevol();
    (void)cloud.constProps().LDevol();

    // Check that the parcel temperature is within necessary limits for
    // devolatilisation to occur
    if (T < cloud.constProps().TDevol() || canCombust == -1)
    {
        return;
    }

    typedef typename TrackCloudType::thermoCloudType thermoCloudType;
    const CompositionModel<thermoCloudType>& composition =
        cloud.composition();


    // Total mass of volatiles evolved
    cloud.devolatilisation().calculate
    (
        dt,
        age,
        mass0,
        mass,
        T,
        YGasEff,
        YLiquidEff,
        YSolidEff,
        canCombust,
        dMassDV
    );

    scalar dMassTot = sum(dMassDV);

    cloud.devolatilisation().addToDevolatilisationMass
    (
        this->nParticle_*dMassTot
    );

    Sh -= dMassTot*cloud.constProps().LDevol()/dt;

    // Update molar emissions
    if (cloud.heatTransfer().BirdCorrection())
    {
        // Molar average molecular weight of carrier mix
        const scalar Wc = max(small, td.rhoc()*RR*td.Tc()/td.pc());

        // Note: hardcoded gaseous diffusivities for now
        // TODO: add to carrier thermo
        const scalar beta = sqr(cbrt(15.0) + cbrt(15.0));

        forAll(dMassDV, i)
        {
            const label id = composition.localToCarrierId(GAS, i);
            const scalar Cp = composition.carrier().Cp(id, td.pc(), Ts);
            const scalar W = composition.carrier().Wi(id);
            const scalar Ni = dMassDV[i]/(this->areaS(d)*dt*W);

            // Dab calc'd using API vapour mass diffusivity function
            const scalar Dab =
                3.6059e-3*(pow(1.8*Ts, 1.75))
               *sqrt(1.0/W + 1.0/Wc)
               /(td.pc()*beta);

            N += Ni;
            NCpW += Ni*Cp*W;
            Cs[id] += Ni*d/(2.0*Dab);
        }
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingClusterParcel<ParcelType>::calcClusterOxidationReactions
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar d,
    const scalar T,
    const scalar mass,
    const label canCombust,
    const scalar alphaTilde_,
    const scalarField& YMix,
    const scalarField& YGas,
    const scalarField& YLiquid,
    const scalarField& YSolid,
    scalarField& dMassORGas,
    scalarField& dMassORLiquid,
    scalarField& dMassORSolid,
    scalarField& dMassORCarrier,
    scalar& Sh,
    scalar& dhsTrans
) const
{
    // Check that model is active
    if
    (
        isType
        <
            NoClusterOxidation
            <
                typename TrackCloudType::reactingClusterCloudType
            >
        >
        (
            cloud.clusterOxidation()
        )
    )
    {
        return;
    }

    // Initialise demand-driven constants
    (void)cloud.constProps().hRetentionCoeff();
    (void)cloud.constProps().TMax();
    (void)cloud.constProps().diffCoeff0();
    (void)cloud.constProps().diffCoeff1();
    (void)cloud.constProps().diffCoeff2();

    // Check that model is active
    if (canCombust != 1)
    {
        return;
    }


    // Update oxidation reactions
    const scalar hReaction = cloud.clusterOxidation().calculate
    (
        dt,
        this->cell(),
        d,
        T,
        td.Tc(),
        td.pc(),
        td.rhoc(),
        mass,
        YGas,
        YLiquid,
        YSolid,
        YMix,
        alphaTilde_,
        dMassORGas,
        dMassORLiquid,
        dMassORSolid,
        dMassORCarrier
    );

    cloud.clusterOxidation().addToOxidationMass
    (
        this->nParticle_
       *(sum(dMassORGas) + sum(dMassORLiquid) + sum(dMassORSolid))
    );

    const scalar xsi = min(T/cloud.constProps().TMax(), 1.0);
    const scalar coeff =
        (1.0 - xsi*xsi)*cloud.constProps().hRetentionCoeff();

    Sh += coeff*hReaction/dt;

    dhsTrans += (1.0 - coeff)*hReaction;
}

template<class ParcelType>
template<class TrackCloudType>
void Foam::ReactingClusterParcel<ParcelType>::calcClusterGasificationReactions
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt,
    const scalar d,
    const scalar T,
    const scalar mass,
    const label canCombust,
    const scalar N,
    const scalarField& YMix,
    const scalarField& YGas,
    const scalarField& YLiquid,
    const scalarField& YSolid,
    scalarField& dMassGRGas,
    scalarField& dMassGRLiquid,
    scalarField& dMassGRSolid,
    scalarField& dMassGRCarrier,
    scalar& Sh,
    scalar& dhsTrans
) const
{
    // Check that model is active
    if
    (
        isType
        <
            NoClusterGasification
            <
                typename TrackCloudType::reactingClusterCloudType
            >
        >
        (
            cloud.clusterGasification()
        )
    )
    {
        return;
    }

    // Initialise demand-driven constants
    (void)cloud.constProps().hRetentionCoeff();
    (void)cloud.constProps().TMax();


    // Check that model is active
    if (canCombust != 1)
    {
        return;
    }

    // Update surface reactions
    const scalar hReaction = cloud.clusterGasification().calculate
    (
        dt,
        this->cell(),
        d,
        T,
        td.Tc(),
        td.pc(),
        td.rhoc(),
        mass,
        YGas,
        YLiquid,
        YSolid,
        YMix,
        N,
        dMassGRGas,
        dMassGRLiquid,
        dMassGRSolid,
        dMassGRCarrier
    );

    cloud.clusterGasification().addToGasificationMass
    (
        this->nParticle_
       *(sum(dMassGRGas) + sum(dMassGRLiquid) + sum(dMassGRSolid))
    );

    const scalar xsi = min(T/cloud.constProps().TMax(), 1.0);
    const scalar coeff =
        (1.0 - xsi*xsi)*cloud.constProps().hRetentionCoeff();

    Sh += coeff*hReaction/dt;

    dhsTrans += (1.0 - coeff)*hReaction;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ReactingClusterParcel<ParcelType>::ReactingClusterParcel
(
    const ReactingClusterParcel<ParcelType>& p
)
:
    ParcelType(p),
    YGas_(p.YGas_),
    YLiquid_(p.YLiquid_),
    YSolid_(p.YSolid_),
    canCombust_(p.canCombust_),
    alphaTilde_(p.alphaTilde_)
{}


template<class ParcelType>
Foam::ReactingClusterParcel<ParcelType>::ReactingClusterParcel
(
    const ReactingClusterParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    YGas_(p.YGas_),
    YLiquid_(p.YLiquid_),
    YSolid_(p.YSolid_),
    canCombust_(p.canCombust_),
    alphaTilde_(p.alphaTilde_)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingClusterParcelIO.C"

// ************************************************************************* //
