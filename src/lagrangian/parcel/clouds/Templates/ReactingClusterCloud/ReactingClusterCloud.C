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

#include "ReactingClusterCloud.H"

#include "DevolatilisationModel.H"
#include "ClusterOxidationModel.H"
#include "ClusterGasificationModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ReactingClusterCloud<CloudType>::setModels()
{
    devolatilisationModel_.reset
    (
        DevolatilisationModel<ReactingClusterCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    clusterOxidationModel_.reset
    (
        ClusterOxidationModel<ReactingClusterCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    clusterGasificationModel_.reset
    (
        ClusterGasificationModel<ReactingClusterCloud<CloudType>>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


template<class CloudType>
void Foam::ReactingClusterCloud<CloudType>::cloudReset
(
    ReactingClusterCloud<CloudType>& c
)
{
    CloudType::cloudReset(c);

    devolatilisationModel_.reset(c.devolatilisationModel_.ptr());
    clusterOxidationModel_.reset(c.clusterOxidationModel_.ptr());
    clusterGasificationModel_.reset(c.clusterGasificationModel_.ptr());

    dMassDevolatilisation_ = c.dMassDevolatilisation_;
    dMassOxidation_ = c.dMassOxidation_;
    dMassGasification_ = c.dMassGasification_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ReactingClusterCloud<CloudType>::ReactingClusterCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const fluidThermo& carrierThermo,
    const bool readFields
)
:
    CloudType(cloudName, rho, U, g, carrierThermo, false),
    cloudCopyPtr_(nullptr),
    voidFraction_
    (
         IOobject
         (
             this->name() + "ThetaCustom",
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::AUTO_WRITE
         ),
         this->mesh(),
         dimensionedScalar("zero", dimless, 0.0)
    ),
    numDens_
    (
         IOobject
         (
             this->name() + "numDens",
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::AUTO_WRITE
         ),
         this->mesh(),
         dimensionedScalar("zero", dimless, 0.0)
    ),
    numDensPrev_
    (
         IOobject
         (
             this->name() + "numDensPrev",
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         this->mesh(),
         dimensionedScalar("zero", dimless, 0.0)
    ),
    dAvg_
    (
         IOobject
         (
             this->name() + "dAvg",
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::AUTO_WRITE
         ),
         this->mesh(),
         dimensionedScalar("zero", dimLength, 0.0)
    ),
    dAvgPrev_
    (
         IOobject
         (
             this->name() + "dAvgPrev",
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         this->mesh(),
         dimensionedScalar("zero", dimLength, 0.0)
    ),
    constProps_(this->particleProperties()),
    devolatilisationModel_(nullptr),
    clusterOxidationModel_(nullptr),
    clusterGasificationModel_(nullptr),
    dMassDevolatilisation_(0.0),
    dMassOxidation_(0.0),
    dMassGasification_(0.0),
    epsilon_(this->mesh().objectRegistry::lookupObjectRef<volScalarField>(carrierThermo.phasePropertyName("epsilon"))),
    k_(this->mesh().objectRegistry::lookupObjectRef<volScalarField>(carrierThermo.phasePropertyName("k")))
{
    setModels();

    if (readFields)
    {
        parcelType::readFields(*this, this->composition());
        this->deleteLostParticles();
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::ReactingClusterCloud<CloudType>::ReactingClusterCloud
(
    ReactingClusterCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    cloudCopyPtr_(nullptr),
    voidFraction_(c.voidFraction_),
    numDens_(c.numDens_),
    numDensPrev_(c.numDensPrev_),
    dAvg_(c.dAvg_),
    dAvgPrev_(c.dAvgPrev_),
    constProps_(c.constProps_),
    devolatilisationModel_(c.devolatilisationModel_->clone()),
    clusterOxidationModel_(c.clusterOxidationModel_->clone()),
    clusterGasificationModel_(c.clusterGasificationModel_->clone()),
    dMassDevolatilisation_(c.dMassDevolatilisation_),
    dMassOxidation_(c.dMassOxidation_),
    dMassGasification_(c.dMassGasification_),
    epsilon_(c.epsilon_),
    k_(c.k_)
{}


template<class CloudType>
Foam::ReactingClusterCloud<CloudType>::ReactingClusterCloud
(
    const fvMesh& mesh,
    const word& name,
    const ReactingClusterCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    cloudCopyPtr_(nullptr),
    voidFraction_(c.voidFraction_),
    numDens_(c.numDens_),
    numDensPrev_(c.numDensPrev_),
    dAvg_(c.dAvg_),
    dAvgPrev_(c.dAvgPrev_),
    constProps_(c.constProps_),
    devolatilisationModel_(nullptr),
    clusterOxidationModel_(nullptr),
    clusterGasificationModel_(nullptr),
    dMassDevolatilisation_(0.0),
    dMassOxidation_(0.0),
    dMassGasification_(0.0),
    epsilon_(c.epsilon_),
    k_(c.k_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ReactingClusterCloud<CloudType>::~ReactingClusterCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ReactingClusterCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);

    label idGas = this->composition().idGas();
    label idLiquid = this->composition().idLiquid();
    label idSolid = this->composition().idSolid();

    parcel.YGas() = this->composition().Y0(idGas);
    parcel.YLiquid() = this->composition().Y0(idLiquid);
    parcel.YSolid() = this->composition().Y0(idSolid);
}


template<class CloudType>
void Foam::ReactingClusterCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);

    if (fullyDescribed)
    {
        label idGas = this->composition().idGas();
        label idLiquid = this->composition().idLiquid();
        label idSolid = this->composition().idSolid();

        this->checkSuppliedComposition
        (
            parcel.YGas(),
            this->composition().Y0(idGas),
            "YGas"
        );
        this->checkSuppliedComposition
        (
            parcel.YLiquid(),
            this->composition().Y0(idLiquid),
            "YLiquid"
        );
        this->checkSuppliedComposition
        (
            parcel.YSolid(),
            this->composition().Y0(idSolid),
            "YSolid"
        );
    }
}


template<class CloudType>
void Foam::ReactingClusterCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<ReactingClusterCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::ReactingClusterCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::ReactingClusterCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
}


template<class CloudType>
void Foam::ReactingClusterCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        if (!this->solution().transient())
        {
            voidFraction_ = Zero;
            numDens_ = Zero;
            dAvg_ = Zero;
        }
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}


template<class CloudType>
void Foam::ReactingClusterCloud<CloudType>::autoMap
(
    const mapPolyMesh& mapper
)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::ReactingClusterCloud<CloudType>::info()
{
    CloudType::info();

    if (this->solution().transient())
    {
        voidFraction_ = this->theta()();
    }
    else
    {
        voidFraction_.primitiveFieldRef() /= this->mesh().time().deltaTValue()*this->mesh().V(); 
        dAvg_.primitiveFieldRef() /= this->mesh().time().deltaTValue()
                                    *(max(numDens_, small)/this->mesh().time().deltaTValue()); 
        dAvgPrev_.primitiveFieldRef() = dAvg_;
        numDens_.primitiveFieldRef() /= this->mesh().time().deltaTValue()*this->mesh().V(); 
        numDensPrev_.primitiveFieldRef() = numDens_;
    }

    this->devolatilisation().info(Info);
    this->clusterOxidation().info(Info);
    this->clusterGasification().info(Info);
}


template<class CloudType>
void Foam::ReactingClusterCloud<CloudType>::writeFields() const
{
    if (this->compositionModel_.valid())
    {
        CloudType::particleType::writeFields(*this, this->composition());
    }
    else
    {
        CloudType::particleType::writeFields(*this);
    }
}


// ************************************************************************* //
