/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    UreaKineticRate

Description
    Kinetic rate surface reaction model for urea SOLID parcels according to Birkhold parameters.
    Limited to:

        CH4N2O(s) -> NH3(g) + HNCO (g)

\*---------------------------------------------------------------------------*/

#ifndef UreaKineticRate_H
#define UreaKineticRate_H

#include "MySurfaceReactionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
template<class CloudType>
class UreaKineticRate;

/*---------------------------------------------------------------------------*\
            Class UreaKineticRate Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class UreaKineticRate
:
    public MySurfaceReactionModel<CloudType>
{
    // Private data

        // Model constants

            //- Kinetics rate pre-exponential factor
            const scalar A_;

            //- Heat of reaction at standard conditions (298K)
            const scalar Dh_;

            //- Kinetics rate activation energy
            const scalar E_;


        // Addressing

            //- Cs positions in global/local lists
            label CsLocalId_;

            //- NH3 position in global list
            label NH3GlobalId_;

            //- HNCO positions in global list
            label HNCOGlobalId_;


        // Local copies of thermo properties

            //- Molecular weight of CH4N2O [kg/kmol]
            scalar WCH4N2O_;

            //- Molecular weight of NH3 [kg/kmol]
            scalar WNH3_;

            //- Molecular weight of HNCO [kg/kmol]
            scalar WHNCO_;


public:

        //- Phase type
        enum phaseType
        {
            liquid,
            solid
        };

        //- Name representations of enthalpy transfer types
        static const Foam::wordList phaseTypeNames;


protected:

    // Protected data

        //- Phase type enumeration
        phaseType phase_;

        //- Convert word to phase type
            phaseType wordToPhase(const word& etName) const;


public:

    //- Runtime type information
    TypeName("UreaKineticRate");


    // Constructors

        //- Construct from dictionary
        UreaKineticRate
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Construct copy
        UreaKineticRate
        (
            const UreaKineticRate<CloudType>& srm
        );

        //- Construct and return a clone
        virtual autoPtr<MySurfaceReactionModel<CloudType> > clone() const
        {
            return autoPtr<MySurfaceReactionModel<CloudType> >
            (
                new UreaKineticRate<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~UreaKineticRate();


    //- Return the enthalpy transfer type enumeration
        const phaseType& phase() const;


    // Member Functions

        //- Update surface reactions
        virtual scalar calculate
        (
            const scalar dt,
            const label cellI,
            const scalar d,
            const scalar T,
            const scalar Tc,
            const scalar pc,
            const scalar rho,
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
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UreaKineticRate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
