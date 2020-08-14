/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::NewLiquidEvaporation

Description
    Liquid evaporation model
    - uses ideal gas assumption
    - includes boiling model based on:

    \verbatim
        "Studies of Superheated Fuel Spray Structures and Vaporization in
        GDI Engines"

        Zuo, B., Gomes, A. M. and Rutland C. J.

        International Journal of Engine Research, 2000, Vol. 1(4), pp. 321-336
    \endverbatim

\*---------------------------------------------------------------------------*/

#ifndef NewLiquidEvaporation_H
#define NewLiquidEvaporation_H

#include "PhaseChangeModel.H"
#include "liquidMixtureProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                   Class NewLiquidEvaporation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NewLiquidEvaporation
:
    public PhaseChangeModel<CloudType>
{
protected:

    // Protected data

        //- Global liquid properties data
        const liquidMixtureProperties& liquids_;

        //- List of active liquid names
        List<word> activeLiquids_;

        //- Mapping between liquid and carrier species
        List<label> liqToCarrierMap_;

        //- Mapping between local and global liquid species
        List<label> liqToLiqMap_;


    // Protected Member Functions

        //- Sherwood number as a function of Renolds and Schmidt numbers
        scalar Sh(const scalar Re, const scalar Sc, const scalar Gr) const;

        //- Nusselt number as a function of Renolds and Prandtl numbers
        scalar Nusselt(const scalar Re, const scalar Pr, const scalar Gr) const;

        //- Calculate the carrier phase component volume fractions at cellI
        tmp<scalarField> calcXc(const label cellI) const;

        //- Calculate the carrier phase molecular weigth at cellI
        tmp<scalar> calcW(const label cellI) const;


public:

    //- Runtime type information
    TypeName("newLiquidEvaporation");


    // Constructors

        //- Construct from dictionary
        NewLiquidEvaporation(const dictionary& dict, CloudType& cloud);

        //- Construct copy
        NewLiquidEvaporation(const NewLiquidEvaporation<CloudType>& pcm);

        //- Construct and return a clone
        virtual autoPtr<PhaseChangeModel<CloudType> > clone() const
        {
            return autoPtr<PhaseChangeModel<CloudType> >
            (
                new NewLiquidEvaporation<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NewLiquidEvaporation();


    // Member Functions

        //- Update model
        virtual void calculate
        (
            const scalar dt,//
            const label cellI,//
            scalar& Re_No,//
            const scalar Pr,
            const scalar d,//
            const scalar nu,
            const scalar T,//
            const scalar Ts,//
            const scalar pc,//
            const scalar Tc,//
            const scalarField& Yl,//
            scalarField& dMassPC,//
            scalar& bt,//
            const scalar mu
        ) const;

        //- Return the enthalpy per unit mass
        virtual scalar dh
        (
            const label idc,
            const label idl,
            const scalar p,
            const scalar T
        ) const;

        //- Return vapourisation temperature
        virtual scalar Tvap(const scalarField& Y) const;

        //- Return maximum/limiting temperature
        virtual scalar TMax(const scalar p, const scalarField& Y) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NewLiquidEvaporation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
