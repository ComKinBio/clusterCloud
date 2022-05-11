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
    Foam::ReactingClusterCloud

Description
    Templated base class for reacting cloud considering
    turbulent clustering effects

    - Adds to reacting multiphase cloud
      - turbulent clustering effect
       @article{Kruger2017,
       author = {Kr{\"{u}}ger, Jonas and Haugen, Nils Erland L. and L{\o}v{\aa}s, Terese},
       doi = {10.1016/j.combustflame.2017.07.008},
       issn = {15562921},
       journal = {Combustion and Flame},
       pages = {160--172},
       title = {{Correlation effects between turbulence and the conversion rate of pulverized char particles}},
       volume = {185},
       year = {2017}
       }

       @article{Haugen2018,
       author = {Haugen, Nils Erland L. and Kr{\"{u}}ger, Jonas and Mitra, Dhrubaditya and L{\o}v{\aa}s, Terese},
       doi = {10.1017/jfm.2017.820},
       issn = {14697645},
       journal = {Journal of Fluid Mechanics},
       pages = {932--951},
       title = {{The effect of turbulence on mass transfer rates of small inertial particles with surface reactions}},
       volume = {836},
       year = {2018}
       }


    - Additionally it adds the possibility for source term diffusion 
      defined in the cloud dictionary:

        sourceDistribution true;

        diffusion
        {
            //for the implicit solution 
            // ddt(Smooth) independentEuler;
            // needs to be specified in the fvSchemes 
            useImplicitLaplacian true;//false;
            
            //proposed to be taken as 6dp
            diffusionBandWidth  0.018;
    
            //Default to the above, if not given
            diffusionBandWidthRadia  0.018; 
            diffusionBandWidthHeat   0.018;
            diffusionBandWidthMass   0.018;
            diffusionBandWidthMomentum  0.018;
        }
        

SourceFiles
    ReactingClusterCloudI.H
    ReactingClusterCloud.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingClusterCloud_H
#define ReactingClusterCloud_H

#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "dimensionedTypes.H"
#include "fvMesh.H"
#include "fluidThermo.H"
#include "Cloud.H"

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class DevolatilisationModel;

template<class CloudType>
class ClusterOxidationModel;

template<class CloudType>
class ClusterGasificationModel;


/*---------------------------------------------------------------------------*\
                 Class ReactingClusterCloudName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(ReactingClusterCloud);


/*---------------------------------------------------------------------------*\
                   Class ReactingClusterCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ReactingClusterCloud
:
    public CloudType,
    public ReactingClusterCloudName
{
public:

    // Public Typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef ReactingClusterCloud<CloudType> reactingClusterCloudType;


private:

    // Private Data

        //- Cloud copy pointer
        autoPtr<ReactingClusterCloud<CloudType>> cloudCopyPtr_;

        //- Function object pointer void fraction
        mutable volScalarField voidFraction_;

        //- Function object pointer void fraction
        mutable volScalarField numDens_;

        //- Function object pointer void fraction from previous time step for steady state sim
        mutable volScalarField numDensPrev_;

        //- Function object pointer average diameter
        mutable volScalarField dAvg_;

        //- Function object pointer average diameter from previous time step for steady state sim
        mutable volScalarField dAvgPrev_;


protected:

    // Protected data

        //- Parcel constant properties
        typename parcelType::constantProperties constProps_;

        //- flag to indicate if source terms should be distributed / diffused
        const bool sourceDistribution_;

        //- flag to switch between explicit and implicit diffusion term
        const bool implicitFvm_;

        //- Phase diffusion band width
        const scalar diffusionBandWidth_;
        
        //- Mass diffusion band width
        const scalar diffusionBandWidthForMassCoupling_;
        
        //- Momentum diffusion band width
        const scalar diffusionBandWidthForMomentumCoupling_;
        
        //- Heat diffusion band width
        const scalar diffusionBandWidthForHeatCoupling_;
        
        //- Radiation diffusion band width
        const scalar diffusionBandWidthForRadiaCoupling_;

        //- Species diffusion band width
        const scalar diffusionBandWidthForSpecieCoupling_;
                        
        //- Diffusion tensor direction
        tensor smoothDirection_;
        
        //- Diffusion tensor
        dimensionedTensor DT_;
        autoPtr<volScalarField> rDeltaTi_;

        // References to the cloud sub-models

            //- Devolatilisation model
            autoPtr
            <
                DevolatilisationModel<ReactingClusterCloud<CloudType>>
            >
            devolatilisationModel_;

            //- Surface reaction model
            autoPtr
            <
                ClusterOxidationModel<ReactingClusterCloud<CloudType>>
            >
            clusterOxidationModel_;

            //- Surface reaction model
            autoPtr
            <
                ClusterGasificationModel<ReactingClusterCloud<CloudType>>
            >
            clusterGasificationModel_;


        // Check

            //- Total mass transferred to continuous phase via devolatilisation
            scalar dMassDevolatilisation_;

            //- Total mass transferred to continuous phase via surface
            //  reactions
            scalar dMassOxidation_;

            //- Total mass transferred to continuous phase via surface
            //  reactions
            scalar dMassGasification_;

            //- epsilonC dissipation rate from the turbulence of the continuous phase
            const volScalarField& epsilon_;

            //- kC turbulenct kinetic energy from the turbulence of the continuous phase
            const volScalarField& k_;


    // Protected Member Functions

        // Initialisation

            //- Set cloud sub-models
            void setModels();

        // Cloud evolution functions

            //- Solve the cloud - calls all evolution functions
            template<class TrackCloudType>
            void solve
            (
                TrackCloudType& cloud,
                typename parcelType::trackingData& td
            );

            //- Evolve the cloud
            template<class TrackCloudType>
            void evolveCloud
            (
                TrackCloudType& cloud,
                typename parcelType::trackingData& td
            );

            //- Return the diffused field obtained from the given volScalarField s
            void diffusion(volScalarField& s, word type);
        
            //- Return the diffused field obtained from the given volScalarField s
            tmp<volScalarField::Internal> diffusion(const volScalarField::Internal& s, word type);

            //- Return the diffused field obtained from the given volVectorField s
            tmp<volVectorField::Internal> diffusion(const volVectorField::Internal& s, word type);

        // Cloud evolution functions

            //- Reset state of cloud
            void cloudReset(ReactingClusterCloud<CloudType>& c);


public:

    // Constructors

        //- Construct given carrier fields and thermo
        ReactingClusterCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const fluidThermo& carrierThermo,
            const bool readFields = true
        );

        //- Copy constructor with new name
        ReactingClusterCloud
        (
            ReactingClusterCloud<CloudType>& c,
            const word& name
        );

        //- Copy constructor with new name - creates bare cloud
        ReactingClusterCloud
        (
            const fvMesh& mesh,
            const word& name,
            const ReactingClusterCloud<CloudType>& c
        );

        //- Disallow default bitwise copy construction
        ReactingClusterCloud(const ReactingClusterCloud&) = delete;

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new ReactingClusterCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new ReactingClusterCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~ReactingClusterCloud();


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const ReactingClusterCloud& cloudCopy() const;

            //- Return the constant properties
            inline const typename parcelType::constantProperties&
                constProps() const;

            //- Return access to the constant properties
            inline typename parcelType::constantProperties& constProps();

            //- Access to particle phase fraction
            inline volScalarField& alphaf();

            //- Access to particle number density 
            inline volScalarField& numberDens();

            //- Access to particle number density from previous timestep
            inline volScalarField& numberDensPrev();

            //- Access to particle average diameter 
            inline volScalarField& dAvg();

            //- Access to particle average diameter from previous timestep
            inline volScalarField& dAvgPrev();


            // Sub-models

                //- Return const access to devolatilisation model
                inline const DevolatilisationModel
                <
                    ReactingClusterCloud<CloudType>
                >&
                devolatilisation() const;

                //- Return reference to devolatilisation model
                inline DevolatilisationModel
                <
                    ReactingClusterCloud<CloudType>
                >&
                devolatilisation();

                //- Return const access to reacting surface reaction model
                inline const ClusterOxidationModel
                <
                    ReactingClusterCloud<CloudType>
                >&
                clusterOxidation() const;

                //- Return reference to reacting surface reaction model
                inline ClusterOxidationModel
                <
                    ReactingClusterCloud<CloudType>
                >&
                clusterOxidation();

                //- Return const access to reacting surface reaction model
                inline const ClusterGasificationModel
                <
                    ReactingClusterCloud<CloudType>
                >&
                clusterGasification() const;

                //- Return reference to reacting surface reaction model
                inline ClusterGasificationModel
                <
                    ReactingClusterCloud<CloudType>
                >&
                clusterGasification();


        // Cloud evolution functions

            //- Set parcel thermo properties
            void setParcelThermoProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt
            );

            //- Check parcel properties
            void checkParcelProperties
            (
                parcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Reset the cloud source terms
            void resetSourceTerms();

            //- Evolve the cloud
            void evolve();


        // Mapping

            //- Remap the cells of particles corresponding to the
            //  mesh topology change with a default tracking data object
            virtual void autoMap(const mapPolyMesh&);


        // I-O

            //- Print cloud information
            void info();

            //- Write the field data for the cloud
            virtual void writeFields() const;

        // References to the carrier gas fields

        //- Return the dissipation rate
        inline const volScalarField& epsilon() const;

        //- Return the turbulence kinetic energy
        inline const volScalarField& k() const;
        
        //- Get end time and dt according to source type
        //- the input shoule be "phase" "mass" "momentum" "heat" "radia" "species"
        
        scalar diffusionDeltaT(word type);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ReactingClusterCloud&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingClusterCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingClusterCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
