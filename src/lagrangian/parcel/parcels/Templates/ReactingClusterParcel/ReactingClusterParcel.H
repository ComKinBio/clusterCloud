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
    Foam::ReactingClusterParcel

Description
    Cluster variant of the reacting parcel class with one/two-way coupling
    with the continuous phase.

SourceFiles
    ReactingClusterParcelI.H
    ReactingClusterParcel.C
    ReactingClusterParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingClusterParcel_H
#define ReactingClusterParcel_H

#include "particle.H"
#include "interpolation.H"
#include "fluidThermo.H"
#include "demandDrivenEntry.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class ReactingClusterParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const ReactingClusterParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                Class ReactingClusterParcelName Declaration
\*---------------------------------------------------------------------------*/

TemplateName(ReactingClusterParcel);


/*---------------------------------------------------------------------------*\
                  Class ReactingClusterParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ReactingClusterParcel
:
    public ParcelType,
    public ReactingClusterParcelName
{
    // Private Data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;


public:

    // IDs of phases in ReactingParcel phase list (Y)

        static const label GAS;
        static const label LIQ;
        static const label SLD;


    //- Class to hold reacting multiphase particle constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {
        // Private Data

            //- Devolatilisation activation temperature [K]
            demandDrivenEntry<scalar> TDevol_;

            //- Latent heat of devolatilisation [J/kg]
            demandDrivenEntry<scalar> LDevol_;

            //- Fraction of enthalpy retained by parcel due to surface
            //  reactions
            demandDrivenEntry<scalar> hRetentionCoeff_;

            //- coefficients to calculate the diffusion coefficient
            demandDrivenEntry<scalar> diffCoeff0_;
            demandDrivenEntry<scalar> diffCoeff1_;
            demandDrivenEntry<scalar> diffCoeff2_;

            //- flag to define if constant heat cpacity should be used
            demandDrivenEntry<bool> constCp_;


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Construct from dictionary
            constantProperties(const dictionary& parentDict);


        // Access

            //- Return const access to the devolatilisation temperature
            inline scalar TDevol() const;

            //- Return const access to the latent heat of devolatilisation
            inline scalar LDevol() const;

            //- Return const access to the fraction of enthalpy retained by
            //  parcel due to surface reactions
            inline scalar hRetentionCoeff() const;

            //- Return const access to the diffusion coefficients
            inline scalar diffCoeff0() const;
            inline scalar diffCoeff1() const;
            inline scalar diffCoeff2() const;

            //- Return a flag, if constant heat capacity should be used for the particle
            inline bool constCp() const;

    };


    //- Use base tracking data
    //typedef typename ParcelType::trackingData trackingData;


    class trackingData
    :
        public ParcelType::trackingData
    {
    private:

        // Private Data

            //- Local copy of carrier turbulent dissipation rate
            //  epsilon not stored on carrier thermo, but returned as tmp<...>
            const volScalarField epsilon_;

            //- Local copy of carrier turbulent kinetic energy field
            //  k not stored on carrier thermo, but returned as tmp<...>
            const volScalarField k_;


            // Interpolators for continuous phase fields

                //- Interpolator for turbulent dissipation rate
                autoPtr<interpolation<scalar>> epsilonInterp_;

                //- Turbulent kinetic energy field interpolator
                autoPtr<interpolation<scalar>> kInterp_;

            // Cached continuous phase properties

                //- turublent dissipation rate 
                scalar epsilonc_;

                //- turbulent kinetic energy
                scalar kc_;



    public:

        // Constructors

            //- Construct from components
            template <class TrackCloudType>
            inline trackingData(const TrackCloudType& cloud);


        // Member Functions

            //- Return access to the locally stored carrier turbulent dissipation rate field
            inline const volScalarField& epsilon() const;

            //- Return access to the locally stored turbulent kinetic energy field
            inline const volScalarField& k() const;

            //- Return const access to the interpolator for continuous phase
            //  turbulent dissipation rate
            inline const interpolation<scalar>& epsilonInterp() const;

            //- Return const access to the interpolator for continuous
            //  turbulent kinetic energy
            inline const interpolation<scalar>& kInterp() const;

            //- Return the continuous phase turbulent dissipation rate
            inline scalar epsilonc() const;

            //- Access the continuous phase turbulent dissipation rate
            inline scalar& epsilonc();

            //- Return the continuous phase turbulent kinetic energy
            inline scalar kc() const;

            //- Access the continuous phase turbulent kinetic energy
            inline scalar& kc();
    };

private:

    // Private Member Functions

        //- Return the mixture effective specific heat capacity
        template<class TrackCloudType>
        scalar CpEff
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Return the mixture effective sensible enthalpy
        template<class TrackCloudType>
        scalar HsEff
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Return the mixture effective latent heat
        template<class TrackCloudType>
        scalar LEff
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Update the mass fractions (Y, YGas, YLiquid, YSolid)
        scalar updateMassFractions
        (
            const scalar mass0,
            const scalarField& dMassGas,
            const scalarField& dMassLiquid,
            const scalarField& dMassSolid
        );


protected:

    // Protected data

        // Parcel properties

            //- Mass fractions of gases []
            scalarField YGas_;

            //- Mass fractions of liquids []
            scalarField YLiquid_;

            //- Mass fractions of solids []
            scalarField YSolid_;

            //- Flag to identify if the particle can devolatilise and combust
            //  Combustion possible only after volatile content falls below
            //  threshold value.  States include:
            //  0 = can devolatilise, cannot combust but can change
            //  1 = can devolatilise, can combust
            // -1 = cannot devolatilise or combust, and cannot change
            label canCombust_;

            //- factor to reduce the mass transport effect in the char conversion
            scalar alphaTilde_;


    // Protected Member Functions

        //- Calculate Devolatilisation
        template<class TrackCloudType>
        void calcDevolatilisation
        (
            TrackCloudType& cloud,
            trackingData& td,
            const scalar dt,           // timestep
            const scalar age,          // age
            const scalar Ts,           // surface temperature
            const scalar d,            // diameter
            const scalar T,            // temperature
            const scalar mass,         // mass
            const scalar mass0,        // mass (initial on injection)
            const scalarField& YGasEff,// gas component mass fractions
            const scalarField& YLiquidEff,// liquid component mass fractions
            const scalarField& YSolidEff,// solid component mass fractions
            label& canCombust,          // 'can combust' flag
            scalarField& dMassDV,      // mass transfer - local to particle
            scalar& Sh,                // explicit particle enthalpy source
            scalar& N,                 // flux of species emitted from particle
            scalar& NCpW,              // sum of N*Cp*W of emission species
            scalarField& Cs            // carrier conc. of emission species
        ) const;

        //- Calculate oxidation reactions
        template<class TrackCloudType>
        void calcClusterOxidationReactions
        (
            TrackCloudType& cloud, 
            trackingData& td,
            const scalar dt,           // timestep
            const scalar d,            // diameter
            const scalar T,            // temperature
            const scalar mass,         // mass
            const label canCombust,     // 'can combust' flag
            const scalar N,            // flux of species emitted from particle
            const scalarField& YMix,   // mixture mass fractions
            const scalarField& YGas,   // gas-phase mass fractions
            const scalarField& YLiquid,// liquid-phase mass fractions
            const scalarField& YSolid, // solid-phase mass fractions
            scalarField& dMassORGas,   // gas-phase mass transfer - local
            scalarField& dMassORLiquid,// liquid-phase mass transfer - local
            scalarField& dMassORSolid, // solid-phase mass transfer - local
            scalarField& dMassORCarrier, // carrier phase mass transfer
            scalar& Sh,                // explicit particle enthalpy source
            scalar& dhsTrans           // sensible enthalpy transfer to carrier
        ) const;

        //- Calculate gasification reactions
        template<class TrackCloudType>
        void calcClusterGasificationReactions
        (
            TrackCloudType& cloud, 
            trackingData& td,
            const scalar dt,           // timestep
            const scalar d,            // diameter
            const scalar T,            // temperature
            const scalar mass,         // mass
            const label canCombust,     // 'can combust' flag
            const scalar N,            // flux of species emitted from particle
            const scalarField& YMix,   // mixture mass fractions
            const scalarField& YGas,   // gas-phase mass fractions
            const scalarField& YLiquid,// liquid-phase mass fractions
            const scalarField& YSolid, // solid-phase mass fractions
            scalarField& dMassORGas,   // gas-phase mass transfer - local
            scalarField& dMassORLiquid,// liquid-phase mass transfer - local
            scalarField& dMassORSolid, // solid-phase mass transfer - local
            scalarField& dMassORCarrier, // carrier phase mass transfer
            scalar& Sh,                // explicit particle enthalpy source
            scalar& dhsTrans           // sensible enthalpy transfer to carrier
        ) const;


        
public:

    // Static Data Members

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " nGas(Y1..YN)"
          + " nLiquid(Y1..YN)"
          + " nSolid(Y1..YN)"
        );


    // Constructors

        //- Construct from mesh, position and topology
        //  Other properties initialised as null
        inline ReactingClusterParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline ReactingClusterParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from Istream
        ReactingClusterParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        ReactingClusterParcel(const ReactingClusterParcel& p);

        //- Construct as a copy
        ReactingClusterParcel
        (
            const ReactingClusterParcel& p,
            const polyMesh& mesh
        );

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new ReactingClusterParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>(new ReactingClusterParcel(*this, mesh));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<ReactingClusterParcel<ParcelType>> operator()
            (
                Istream& is
            ) const
            {
                return autoPtr<ReactingClusterParcel<ParcelType>>
                (
                    new ReactingClusterParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to mass fractions of gases
            inline const scalarField& YGas() const;

            //- Return const access to mass fractions of liquids
            inline const scalarField& YLiquid() const;

            //- Return const access to mass fractions of solids
            inline const scalarField& YSolid() const;

            //- Return const access to the canCombust flag
            inline label canCombust() const;

            //- Return const access to the alphaTilde value
            inline const scalar& alphaTilde() const;

        // Edit

            //- Return access to mass fractions of gases
            inline scalarField& YGas();

            //- Return access to mass fractions of liquids
            inline scalarField& YLiquid();

            //- Return access to mass fractions of solids
            inline scalarField& YSolid();

            //- Return access to the canCombust flag
            inline label& canCombust();

            //- Return const access to the alphaTilde value
            inline scalar& alphaTilde();

        // Main calculation loop

            //- Set cell values
            template<class TrackCloudType>
            void setCellValues(TrackCloudType& cloud, trackingData& td);

            //- Correct cell values using latest transfer information
            template<class TrackCloudType>
            void cellValueSourceCorrection
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Update parcel properties over the time interval
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );


        // I-O

            //- Read
            template<class CloudType, class CompositionType>
            static void readFields
            (
                CloudType& c,
                const CompositionType& compModel
            );

            //- Read - no composition
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write
            template<class CloudType, class CompositionType>
            static void writeFields
            (
                const CloudType& c,
                const CompositionType& compModel
            );

            //- Read - composition supplied
            template<class CloudType>
            static void writeFields(const CloudType& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const ReactingClusterParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingClusterParcelI.H"
#include "ReactingClusterParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingClusterParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
