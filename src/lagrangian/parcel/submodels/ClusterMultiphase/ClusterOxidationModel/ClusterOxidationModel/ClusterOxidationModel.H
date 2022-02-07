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
    Foam::ClusterOxidationModel

Description
    Templated oxidation model class

SourceFiles
    ClusterOxidationModel.C
    ClusterOxidationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef ClusterOxidationModel_H
#define ClusterOxidationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ClusterOxidationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ClusterOxidationModel
:
    public CloudSubModelBase<CloudType>
{
protected:

    // Protected data

        //- Mass of lagrangian phase converted
        scalar dMass_;


public:

    //-Runtime type information
    TypeName("clusterOxidationModel");


    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        ClusterOxidationModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& cloud
        ),
        (dict, cloud)
    );


    // Constructors

        //- Construct null from owner
        ClusterOxidationModel(CloudType& owner);

        //- Construct from dictionary
        ClusterOxidationModel
        (
            const dictionary& dict,
            CloudType& cloud,
            const word& type
        );

        //- Construct copy
        ClusterOxidationModel(const ClusterOxidationModel<CloudType>& srm);

        //- Construct and return a clone
        virtual autoPtr<ClusterOxidationModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~ClusterOxidationModel();


    //- Selector
    static autoPtr<ClusterOxidationModel<CloudType>> New
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
            scalarField& dMassSRCarrier
        ) const = 0;

        //- Add to devolatilisation mass
        void addToOxidationMass(const scalar dMass);


        // I-O

            //- Write injection info to stream
            virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeClusterOxidationModel(CloudType)                                          \
                                                                               \
    typedef Foam::CloudType::reactingClusterCloudType                           \
        reactingClusterCloudType;                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::ClusterOxidationModel<reactingClusterCloudType>,                         \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            ClusterOxidationModel<reactingClusterCloudType>,                           \
            dictionary                                                         \
        );                                                                     \
    }


#define makeClusterOxidationModelType(SS, CloudType)                                  \
                                                                               \
    typedef Foam::CloudType::reactingClusterCloudType                           \
        reactingClusterCloudType;                                               \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::SS<reactingClusterCloudType>, 0);                                \
                                                                               \
    Foam::ClusterOxidationModel<reactingClusterCloudType>::                            \
        adddictionaryConstructorToTable                                        \
        <Foam::SS<reactingClusterCloudType>>                                    \
        add##SS##CloudType##reactingClusterCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ClusterOxidationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
