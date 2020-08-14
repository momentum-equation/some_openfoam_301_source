/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::MyDevolatilisationModel

Description
    Templated devolatilisation model class

SourceFiles
    MyDevolatilisationModel.C
    MyDevolatilisationModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef MyDevolatilisationModel_H
#define MyDevolatilisationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "SubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class MyDevolatilisationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MyDevolatilisationModel
:
    public SubModelBase<CloudType>
{
protected:

    // Protected data

        //- Mass of lagrangian phase converted
        scalar dMass_;


public:

    //- Runtime type information
    TypeName("mydevolatilisationModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        MyDevolatilisationModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        MyDevolatilisationModel(CloudType& owner);

        //- Construct from dictionary
        MyDevolatilisationModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        MyDevolatilisationModel(const MyDevolatilisationModel<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<MyDevolatilisationModel<CloudType> > clone() const
        {
            return autoPtr<MyDevolatilisationModel<CloudType> >
            (
                new MyDevolatilisationModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~MyDevolatilisationModel();


    //- Selector
    static autoPtr<MyDevolatilisationModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        //- Update model
        virtual void calculate
        (
            const scalar dt,
            const scalar age,
            const scalar mass0,
            const scalar mass,
            const scalar T,
            const scalarField& YGasEff,
            const scalarField& YLiquidEff,
            const scalarField& YSolidEff,
            bool& canCombust,
            scalarField& dMassDV
        ) const;

        //- Add to devolatilisation mass
        void addToDevolatilisationMass(const scalar dMass);


        // I-O

            //- Write injection info to stream
            virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeMyDevolatilisationModel(CloudType)                                  \
                                                                              \
    typedef CloudType::sprayCloudType                            \
        sprayCloudType;                                          \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        MyDevolatilisationModel<sprayCloudType>,                   \
        0                                                                     \
    );                                                                        \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        MyDevolatilisationModel<sprayCloudType>,                   \
        dictionary                                                            \
    );


#define makeMyDevolatilisationModelType(SS, CloudType)                          \
                                                                              \
    typedef CloudType::sprayCloudType                            \
        sprayCloudType;                                          \
    defineNamedTemplateTypeNameAndDebug(SS<sprayCloudType>, 0);  \
                                                                              \
    MyDevolatilisationModel<sprayCloudType>::                      \
        adddictionaryConstructorToTable<SS<sprayCloudType> >     \
        add##SS##CloudType##sprayCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MyDevolatilisationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
