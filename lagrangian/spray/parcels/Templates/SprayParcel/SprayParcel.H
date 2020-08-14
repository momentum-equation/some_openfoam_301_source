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
    Foam::SprayParcel

Description
    Reacing spray parcel, with added functionality for atomization and breakup

\*---------------------------------------------------------------------------*/

#ifndef SprayParcel_H
#define SprayParcel_H

#include "particle.H"
#include "SLGThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class SprayParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const SprayParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                          Class SprayParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class SprayParcel
:
    public ParcelType
{

public:

    // IDs of phases in ReacingParcel phase list (Y)

        static const label GAS;
        static const label LIQ;
        static const label SLD;


    //- Class to hold reacting multiphase particle constant properties
    class constantProperties
    :
        public ParcelType::constantProperties
    {
        // Private data

            //- Latent heat of devolatilisation [J/kg]
            scalar LDevol_;

            //- Fraction of enthalpy retained by parcel due to surface
            //  reactions
            scalar hRetentionCoeff_;


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Constructor from dictionary
            constantProperties
            (
                const dictionary& parentDict,
                const bool readFields = true
            );


        // Access

            //- Return const access to the latent heat of devolatilisation
            inline scalar LDevol() const;

            //- Return const access to the fraction of enthalpy retained by
            //  parcel due to surface reactions
            inline scalar hRetentionCoeff() const;
    };


private:

    // Private Member Functions

        //- Return the mixture effective specific heat capacity
        template<class TrackData>
        scalar CpEff
        (
            TrackData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Return the mixture effective density
        template<class TrackData>
        scalar rhoEff
        (
            TrackData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Return the mixture effective sensible enthalpy
        template<class TrackData>
        scalar HsEff
        (
            TrackData& td,
            const scalar p,
            const scalar T,
            const label idG,
            const label idL,
            const label idS
        ) const;

        //- Return the mixture effective latent heat
        template<class TrackData>
        scalar LEff
        (
            TrackData& td,
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

        // Spray parcel properties

            //- Initial droplet diameter
            scalar d0_;

            //- Injection position
            vector position0_;

            //- Part of liquid core ( >0.5=liquid, <0.5=droplet )
            scalar liquidCore_;

            //- Renolds number
            scalar Renolds_;

            //- Weber number
            scalar Weber_;

            //- Squared oscillation frequency
            scalar om2_;

            //- Index for KH Breakup
            scalar KHindex_;

            //- Spherical deviation
            scalar y_;

            //- Rate of change of spherical deviation
            scalar yDot_;

            //- Characteristic time (used in atomization and/or breakup model)
            scalar tc_;

            //- Stripped parcel mass due to breakup
            scalar ms_;

            //- Injected from injector (needed e.g. for calculating distance from injector)
            scalar injector_;

            //- Momentum relaxation time (needed for calculating parcel acc.)
            scalar tMom_;

            //- Passive scalar (extra variable to be defined by user)
            scalar user_;

            //- Spalding tranfer factor
            scalar Bt_;

            //- Mass fractions of gases []
            scalarField YGas_;

            //- Mass fractions of liquids []
            scalarField YLiquid_;

            //- Mass fractions of solids []
            scalarField YSolid_;

            //- Flag to say that the particle is allowed to combust(Only true after volatile content falls below threshold value)
            bool canCombust_;

	//- Calculate Phase change
        template<class TrackData>
        void newCalcPhaseChange
        (
            TrackData& td,
            const scalar dt,           // timestep
            const label cellI,         // owner cell
            scalar& Re_No,                 // Renolds number
            const scalar Pr,           // Prandtl number
            const scalar Ts,           // Surface temperature
            const scalar nus,          // Surface kinematic viscosity
            const scalar mus,          // Surface dynamic viscosity
            const scalar d,            // diameter
            const scalar T,            // temperature
            const scalar pc,           // pressure
            const scalar mass,         // mass
            const label idPhase,       // id of phase involved in phase change
            const scalar YPhase,       // total mass fraction
            const scalarField& YComponents, // component mass fractions
            scalarField& dMassPC,      // mass transfer - local to particle
            scalar& Sh,                // explicit particle enthalpy source
            scalar& N,                 // flux of species emitted from particle
            scalar& NCpW,              // sum of N*Cp*W of emission species
            scalarField& Cs,           // carrier conc. of emission species
            scalar& bt                  // Spalding mass transfer
        );

	//- Calculate Devolatilisation
        template<class TrackData>
        void newCalcDevolatilisation
        (
            TrackData& td,
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
            bool& canCombust,          // 'can combust' flag
            scalarField& dMassDV,      // mass transfer - local to particle
            scalar& Sh,                // explicit particle enthalpy source
            scalar& N,                 // flux of species emitted from particle
            scalar& NCpW,              // sum of N*Cp*W of emission species
            scalarField& Cs            // carrier conc. of emission species
        ) const;

        //- Calculate surface reactions
        template<class TrackData>
        void newCalcSurfaceReactions
        (
            TrackData& td,
            const scalar dt,           // timestep
            const label cellI,         // owner cell
            const scalar d,            // diameter
            const scalar T,            // temperature
            const scalar mass,         // mass
            const bool canCombust,     // 'can combust' flag
            const scalar N,            // flux of species emitted from particle
            const scalarField& YMix,   // mixture mass fractions
            const scalarField& YGas,   // gas-phase mass fractions
            const scalarField& YLiquid,// liquid-phase mass fractions
            const scalarField& YSolid, // solid-phase mass fractions
            scalarField& dMassSRGas,   // gas-phase mass transfer - local
            scalarField& dMassSRLiquid,// liquid-phase mass transfer - local
            scalarField& dMassSRSolid, // solid-phase mass transfer - local
            scalarField& dMassSRCarrier, // carrier phase mass transfer
            scalar& Sh,                // explicit particle enthalpy source
            scalar& dhsTrans           // sensible enthalpy transfer to carrier
        ) const;

public:

    // Static data members

        //- String representation of properties
        static string propHeader;

        //- Runtime type information
        TypeName("SprayParcel");


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline SprayParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );

        //- Construct from components
        inline SprayParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const scalar dTarget0,
            const vector& U0,
            const vector& f0,
            const vector& angularMomentum0,
            const vector& torque0,
            const scalarField& Y0,
            const scalar liquidCore,
            const scalar Renolds,
            const scalar Weber,
            const scalar om2,
            const scalar KHindex,
            const scalar y,
            const scalar yDot,
            const scalar tc,
            const scalar ms,
            const scalar injector,
            const scalar tMom,
            const scalar user,
            const scalar Bt,
	    const scalarField& YGas0,
	    const scalarField& YLiquid0,
	    const scalarField& YSolid0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        SprayParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        SprayParcel
        (
            const SprayParcel& p,
            const polyMesh& mesh
        );

        //- Construct as a copy
        SprayParcel(const SprayParcel& p);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new SprayParcel<ParcelType>(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>
            (
                new SprayParcel<ParcelType>(*this, mesh)
            );
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

            autoPtr<SprayParcel<ParcelType> > operator()(Istream& is) const
            {
                return autoPtr<SprayParcel<ParcelType> >
                (
                    new SprayParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to initial droplet diameter
            inline scalar d0() const;

            //- Return const access to initial droplet position
            inline const vector& position0() const;

            //- Return const access to liquid core
            inline scalar liquidCore() const;

            //- Return const access to Renolds number
            inline scalar Renolds() const;

            //- Return const access to Weber number
            inline scalar Weber() const;

            //- Return const access to squared oscillation frequency
            inline scalar om2() const;

            //- Return const access to Kelvin-Helmholtz breakup index
            inline scalar KHindex() const;

            //- Return const access to spherical deviation
            inline scalar y() const;

            //- Return const access to rate of change of spherical deviation
            inline scalar yDot() const;

            //- Return const access to atomization characteristic time
            inline scalar tc() const;

            //- Return const access to stripped parcel mass
            inline scalar ms() const;

            //- Return const access to injector id
            inline scalar injector() const;

            //- Return const access to momentum relaxation time
            inline scalar tMom() const;

            //- Return const access to passive user scalar
            inline scalar user() const;

            //- Return const access to Spalding transfer factor
            inline scalar Bt() const;

            //- Return const access to mass fractions of gases
            inline const scalarField& YGas() const;

            //- Return const access to mass fractions of liquids
            inline const scalarField& YLiquid() const;

            //- Return const access to mass fractions of solids
            inline const scalarField& YSolid() const;

            //- Return const access to the canCombust flag
            inline bool canCombust() const;

        // Edit

            //- Return access to initial droplet diameter
            inline scalar& d0();

            //- Return access to initial droplet position
            inline vector& position0();

            //- Return access to liquid core
            inline scalar& liquidCore();

            //- Return const access to Renolds number
            inline scalar& Renolds();

            //- Return const access to Weber number
            inline scalar& Weber();

            //- Return const access to squared oscillation frequency
            inline scalar& om2();

            //- Return access to Kelvin-Helmholtz breakup index
            inline scalar& KHindex();

            //- Return access to spherical deviation
            inline scalar& y();

            //- Return access to rate of change of spherical deviation
            inline scalar& yDot();

            //- Return access to atomization characteristic time
            inline scalar& tc();

            //- Return access to stripped parcel mass
            inline scalar& ms();

            //- Return access to injector id
            inline scalar& injector();

            //- Return access to momentum relaxation time
            inline scalar& tMom();

            //- Return access to passive user scalar
            inline scalar& user();

            //- Return const access to Spalding transfer factor
            inline scalar& Bt();

            //- Return access to mass fractions of gases
            inline scalarField& YGas();

            //- Return access to mass fractions of liquids
            inline scalarField& YLiquid();

            //- Return access to mass fractions of solids
            inline scalarField& YSolid();

            //- Return access to the canCombust flag
            inline bool& canCombust();

        // Main calculation loop

            //- Set cell values
            template<class TrackData>
            void setCellValues
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Correct parcel properties according to atomization model
            template<class TrackData>
            void calcAtomization
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Correct parcel properties according to breakup model
            template<class TrackData>
            void calcBreakup
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Correct cell values using latest transfer information
            template<class TrackData>
            void cellValueSourceCorrection
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Correct surface values due to emitted species
            template<class TrackData>
            void correctSurfaceValues
            (
                TrackData& td,
                const label cellI,
                const scalar T,
                const scalarField& Cs,
                scalar& rhos,
                scalar& mus,
                scalar& Pr,
                scalar& kappa
            );

            //- Update parcel properties over the time interval
            template<class TrackData>
            void calc
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

	    //- Calculate new particle velocity
		template<class TrackData>
		const vector newCalcVelocity
		(
		    TrackData& td,
		    const scalar dt,           // timestep
		    const label cellI,         // owner cell
		    const scalar Re,           // Renolds number
		    const scalar y,	       // droplet distorsion
		    const scalar We,           // Weber number
		    const scalar mu,           // local carrier viscosity
		    const scalar mass,         // mass
		    const scalar YPhase,       // TOTAL mass fraction of liquid phase
		    const vector& Su,          // explicit particle momentum source
		    vector& dUTrans,           // momentum transfer to carrier
		    scalar& Spu,                // linearised drag coefficient
		    const scalar bt      // Spalding evaporation coefficient for drag
		) const;

            //- Calculate the chi-factor for flash-boiling for the atomization model
            template<class TrackData>
            scalar chi
            (
                TrackData& td,
                const scalarField& X
            ) const;

            //- Solve the TAB equation
            template<class TrackData>
            void solveTABEq
            (
                TrackData& td,
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

            //- Write - composition supplied
            template<class CloudType>
            static void writeFields(const CloudType& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const SprayParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SprayParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SprayParcel.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
