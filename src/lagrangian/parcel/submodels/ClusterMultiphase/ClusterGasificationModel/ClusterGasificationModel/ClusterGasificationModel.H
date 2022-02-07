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
    Foam::ClusterGasificationModel

Description
    Templated ClusterGasification model class

SourceFiles
    ClusterGasificationModel.C
    ClusterGasificationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef ClusterGasificationModel_H
#define ClusterGasificationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ClusterGasificationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ClusterGasificationModel
:
    public CloudSubModelBase<CloudType>
{
protected:

    // Protected data

        //- Mass of lagrangian phase converted
        scalar dMass_;


public:

    //-Runtime type information
    TypeName("clusterGasificationModel");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        ClusterGasificationModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& cloud
        ),
        (dict, cloud)
    );


    // Constructors

        //- Construct null from owner
        ClusterGasificationModel(CloudType& owner);

        //- Construct from dictionary
        ClusterGasificationModel
        (
            const dictionary& dict,
            CloudType& cloud,
            const word& type
        );

        //- Construct copy
        ClusterGasificationModel(const ClusterGasificationModel<CloudType>& srm);

        //- Construct and return a clone
        virtual autoPtr<ClusterGasificationModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~ClusterGasificationModel();


    //- Selector
    static autoPtr<ClusterGasificationModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& cloud
    );


    // Member Functions

        //- Update surface reactions
        //  Returns the heat of reaction
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
        ) const = 0;

        //- Add to devolatilisation mass
        void addToGasificationMass(const scalar dMass);


        // I-O

            //- Write injection info to stream
            virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeClusterGasificationModel(CloudType)                                \
                                                                               \
    typedef Foam::CloudType::reactingClusterCloudType                          \
        reactingClusterCloudType;                                              \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::ClusterGasificationModel<reactingClusterCloudType>,              \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            ClusterGasificationModel<reactingClusterCloudType>,                \
            dictionary                                                         \
        );                                                                     \
    }


#define makeClusterGasificationModelType(SS, CloudType)                        \
                                                                               \
    typedef Foam::CloudType::reactingClusterCloudType                          \
        reactingClusterCloudType;                                              \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::SS<reactingClusterCloudType>, 0);                               \
                                                                               \
    Foam::ClusterGasificationModel<reactingClusterCloudType>::                 \
        adddictionaryConstructorToTable                                        \
        <Foam::SS<reactingClusterCloudType>>                                   \
        add##SS##CloudType##reactingClusterCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ClusterGasificationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
