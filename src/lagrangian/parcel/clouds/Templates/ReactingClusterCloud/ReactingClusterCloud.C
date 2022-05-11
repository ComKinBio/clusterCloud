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
    if (sourceDistribution_)
    {
        rDeltaTi_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "rDeltaTi",
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimless/dimTime, 0.0)
            )
        );
    }
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

template<class CloudType>
template<class TrackCloudType>
void Foam::ReactingClusterCloud<CloudType>::evolveCloud
(
    TrackCloudType& cloud,
    typename parcelType::trackingData& td
)
{
    if (this->solution_.coupled())
    {
        cloud.resetSourceTerms();
    }

    if (this->solution_.transient())
    {
        label preInjectionSize = this->size();

        this->surfaceFilm().inject(cloud);

        // Update the cellOccupancy if the size of the cloud has changed
        // during the injection.
        if (preInjectionSize != this->size())
        {
            this->updateCellOccupancy();
            preInjectionSize = this->size();
        }

        // Assume that motion will update the cellOccupancy as necessary
        // before it is required.
        cloud.motion(cloud, td);

        this->injectors_.inject(cloud, td);

        this->stochasticCollision().update(td, this->solution_.trackTime());

    }
    else
    {
//        this->surfaceFilm().injectSteadyState(cloud);

        this->injectors_.injectSteadyState
        (
            cloud,
            td,
            this->solution_.trackTime()
        );

        CloudType::move(cloud, td, this->solution_.trackTime());
    }
    if(sourceDistribution_)
    {
        this->hsTransRef() = diffusion(this->hsTrans(), "heat");
        this->hsCoeffRef() = diffusion(this->hsCoeff(), "heat");

        if (this->radiation_)
        {
            this->radAreaP() = diffusion(this->radAreaP(),"radia");
            this->radT4() = diffusion(this->radT4(),"radia");
            this->radAreaPT4() = diffusion(this->radAreaPT4(),"radia");
        }

        this->UTransRef() = diffusion(this->UTrans(),"momentum");
        this->UCoeffRef() = diffusion(this->UCoeff(),"momentum");
        
        if (this->compositionModel_.valid())
        {
            label idGas = this->composition().idGas();
            label idLiquid = this->composition().idLiquid();

            forAll(this->composition().Y0(idGas), i)
            {
                label gid = this->composition().localToCarrierId(0, i);
                this->rhoTrans(gid) = diffusion(this->rhoTrans(gid), "species");
            }
            forAll(this->composition().Y0(idLiquid), i)
            {
                label gid = this->composition().localToCarrierId(1, i);
                this->rhoTrans(gid) = diffusion(this->rhoTrans(gid), "species");
            }
        }
    }
}


template<class CloudType>
void Foam::ReactingClusterCloud<CloudType>::diffusion
(
    volScalarField& s,
    word type
)
{    
    volScalarField diffWorkField
    (
        IOobject
        (
            "Smooth",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(s.dimensions(), scalar(0.0)),
        zeroGradientFvPatchScalarField::typeName
        
    );
    
    scalarField& diffWorkFieldInterFeildRef = diffWorkField.ref();
    
    scalarField& sInterFeildRef = s.ref();

    diffWorkFieldInterFeildRef = sInterFeildRef;

    scalarField& trDeltaTi = rDeltaTi_->field();
    trDeltaTi = 1.0/diffusionDeltaT(type);
    
    if(implicitFvm_)
    {
        Foam::solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT_, diffWorkField));
    }
    else
    {
        Foam::solve(fvm::ddt(diffWorkField) - fvc::laplacian(DT_, diffWorkField));
    }

    sInterFeildRef = diffWorkField;
    
    return;
}
        

template<class CloudType>
Foam::tmp<Foam::volScalarField::Internal> Foam::ReactingClusterCloud<CloudType>::diffusion
(
    const volScalarField::Internal& s,
    word type
)
{        
    tmp<volScalarField::Internal> tS
    (
        volScalarField::Internal::New
        (
            "tS",
            this->mesh(),
            dimensionedScalar(s.dimensions(), 0)
        )
    );

    scalarField& S = tS.ref();
    
    S = s;
    
    volScalarField diffWorkField
    (
        IOobject
        (
            "Smooth",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(s.dimensions(), scalar(0.0)),
        zeroGradientFvPatchScalarField::typeName
    );
    
    diffWorkField.ref() = s;
    diffWorkField.primitiveFieldRef() = s;

    scalarField& diffWorkFieldInterFeildRef = diffWorkField.ref();

    diffWorkFieldInterFeildRef = S;

    scalarField& trDeltaTi = rDeltaTi_->field();
    trDeltaTi = 1.0/diffusionDeltaT(type);

    if(implicitFvm_)
    {
        Foam::solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT_, diffWorkField));
    }
    else
    {
        Foam::solve(fvm::ddt(diffWorkField) - fvc::laplacian(DT_, diffWorkField));
    }

    S = diffWorkField.internalField();

    return tS;
}

template<class CloudType>
Foam::tmp<Foam::volVectorField::Internal> Foam::ReactingClusterCloud<CloudType>::diffusion
(
    const volVectorField::Internal& s,
    word type
)
{        
    tmp<volVectorField::Internal> tS
    (
        volVectorField::Internal::New
        (
            "tS",
            this->mesh(),
            dimensionedVector(s.dimensions(), vector::zero)
        )
    );

    vectorField& S = tS.ref();
    
    S = s;
    
    volVectorField diffWorkField
    (
        IOobject
        (
            "Smooth",
            this->mesh().time().timeName(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedVector(s.dimensions(),vector::zero),
        zeroGradientFvPatchVectorField::typeName
    );

    
    vectorField& diffWorkFieldInterFeildRef = diffWorkField.ref();
    
    diffWorkFieldInterFeildRef = S;

    scalarField& trDeltaTi = rDeltaTi_->field();
    trDeltaTi = 1.0/diffusionDeltaT(type);
    
    if(implicitFvm_)
    {
        Foam::solve(fvm::ddt(diffWorkField) - fvm::laplacian(DT_, diffWorkField));
    }
    else
    {
        Foam::solve(fvm::ddt(diffWorkField) - fvc::laplacian(DT_, diffWorkField));
    }

    S = diffWorkField.internalField();
    
    return tS;
}

template<class CloudType>
template<class TrackCloudType>
void Foam::ReactingClusterCloud<CloudType>::solve
(
    TrackCloudType& cloud,
    typename parcelType::trackingData& td
)
{
    if (this->solution_.steadyState())
    {
        cloud.storeState();

        cloud.preEvolve();

        evolveCloud(cloud, td);

        if (this->solution_.coupled())
        {
            cloud.relaxSources(cloud.cloudCopy());
        }
    }
    else
    {
        cloud.preEvolve();

        evolveCloud(cloud, td);

        if (this->solution_.coupled())
        {
            cloud.scaleSources();
        }
    }

    cloud.info();

    cloud.postEvolve();

    if (this->solution_.steadyState())
    {
        cloud.restoreState();
    }
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
    sourceDistribution_(this->solution().dict().template lookupOrDefault<bool>("sourceDistribution", false)),
    implicitFvm_(this->solution().dict().subDict("diffusion").template lookupOrDefault<bool>("useImplicitLaplacian", false)),
    diffusionBandWidth_(this->solution().dict().subDict("diffusion").template lookupOrDefault<scalar>("diffusionBandWidth", 0.024)),
    diffusionBandWidthForMassCoupling_(this->solution().dict().subDict("diffusion").template lookupOrDefault<scalar>("diffusionBandWidthMass", diffusionBandWidth_)),
    diffusionBandWidthForMomentumCoupling_(this->solution().dict().subDict("diffusion").template lookupOrDefault<scalar>("diffusionBandWidth", diffusionBandWidth_)),
    diffusionBandWidthForHeatCoupling_(this->solution().dict().subDict("diffusion").template lookupOrDefault<scalar>("diffusionBandWidthMomentum", diffusionBandWidth_)),
    diffusionBandWidthForRadiaCoupling_(this->solution().dict().subDict("diffusion").template lookupOrDefault<scalar>("diffusionBandWidthRadia", diffusionBandWidth_)),
    diffusionBandWidthForSpecieCoupling_(this->solution().dict().subDict("diffusion").template lookupOrDefault<scalar>("diffusionBandWidthSpecie", diffusionBandWidth_)),
    smoothDirection_
    (
        this->solution().dict().subDict("diffusion").template lookupOrDefault
        (
            "smoothDirection",
            tensor(1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0)
        )
    ),
    DT_("DT", dimensionSet(0, 2, -1, 0, 0), smoothDirection_),
    rDeltaTi_(nullptr),
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
    sourceDistribution_(c.sourceDistribution_),
    implicitFvm_(c.implicitFvm_),
    diffusionBandWidth_(c.diffusionBandWidth_),
    diffusionBandWidthForMassCoupling_(c.diffusionBandWidthForMassCoupling_),
    diffusionBandWidthForMomentumCoupling_(c.diffusionBandWidthForMomentumCoupling_),
    diffusionBandWidthForHeatCoupling_(c.diffusionBandWidthForHeatCoupling_),
    diffusionBandWidthForRadiaCoupling_(c.diffusionBandWidthForRadiaCoupling_),
    diffusionBandWidthForSpecieCoupling_(c.diffusionBandWidthForSpecieCoupling_),
    smoothDirection_(c.smoothDirection_),
    DT_(c.DT_),
    rDeltaTi_(nullptr),
    devolatilisationModel_(c.devolatilisationModel_->clone()),
    clusterOxidationModel_(c.clusterOxidationModel_->clone()),
    clusterGasificationModel_(c.clusterGasificationModel_->clone()),
    dMassDevolatilisation_(c.dMassDevolatilisation_),
    dMassOxidation_(c.dMassOxidation_),
    dMassGasification_(c.dMassGasification_),
    epsilon_(c.epsilon_),
    k_(c.k_)
{
    if (c.sourceDistribution_)
    {
        rDeltaTi_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "rDeltaTi",
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimless/dimTime, 0.0)
            )
        );
    }
}



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
    sourceDistribution_(c.sourceDistribution_),
    implicitFvm_(c.implicitFvm_),
    diffusionBandWidth_(c.diffusionBandWidth_),
    diffusionBandWidthForMassCoupling_(c.diffusionBandWidthForMassCoupling_),
    diffusionBandWidthForMomentumCoupling_(c.diffusionBandWidthForMomentumCoupling_),
    diffusionBandWidthForHeatCoupling_(c.diffusionBandWidthForHeatCoupling_),
    diffusionBandWidthForRadiaCoupling_(c.diffusionBandWidthForRadiaCoupling_),
    diffusionBandWidthForSpecieCoupling_(c.diffusionBandWidthForSpecieCoupling_),
    smoothDirection_(c.smoothDirection_),
    DT_(c.DT_),
    rDeltaTi_(nullptr),
    devolatilisationModel_(nullptr),
    clusterOxidationModel_(nullptr),
    clusterGasificationModel_(nullptr),
    dMassDevolatilisation_(0.0),
    dMassOxidation_(0.0),
    dMassGasification_(0.0),
    epsilon_(c.epsilon_),
    k_(c.k_)
{
    if (c.sourceDistribution_)
    {
        rDeltaTi_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "rDeltaTi",
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimless/dimTime, 0.0)
            )
        );
    }
}


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

        if (sourceDistribution_)
    {
        rDeltaTi_->field() = 0.0;
    }
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
