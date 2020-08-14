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

\*---------------------------------------------------------------------------*/

#include "UreaKineticRate.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
const Foam::wordList Foam::UreaKineticRate<CloudType>::
phaseTypeNames
(
    IStringStream
    (
        "("
            "liquid "
            "solid"
        ")"
    )()
);


// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class CloudType>
typename Foam::UreaKineticRate<CloudType>::phaseType
Foam::UreaKineticRate<CloudType>::wordToPhase(const word& etName)
const
{
    forAll(phaseTypeNames, i)
    {
        if (etName == phaseTypeNames[i])
        {
            return phaseType(i);
        }
    }

    FatalErrorIn
    (
        "UreaKineticRate<CloudType>::phaseType"
        "UreaKineticRate<CloudType>::wordToPhaseTransfer(const word&) const"
    )   << "Unknown phaseType " << etName << ". Valid selections are:" << nl
        << phaseTypeNames << exit(FatalError);

    return phaseType(0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::UreaKineticRate<CloudType>::UreaKineticRate
(
    const dictionary& dict,
    CloudType& owner
)
:
    MySurfaceReactionModel<CloudType>(dict, owner, typeName),
    A_(readScalar(this->coeffDict().lookup("A"))),
    Dh_(readScalar(this->coeffDict().lookup("Dh"))),
    E_(readScalar(this->coeffDict().lookup("E"))),
    CsLocalId_(-1),
    NH3GlobalId_(owner.composition().globalCarrierId("NH3")),
    HNCOGlobalId_(owner.composition().globalCarrierId("HNCO")),
    WCH4N2O_(0.0),
    WNH3_(0.0),
    WHNCO_(0.0),
    phase_(wordToPhase(this->coeffDict().lookup("phase")))
{
    // Determine specie mass fraction
    label idSolid = owner.composition().idSolid();
    label idLiquid = owner.composition().idLiquid();
    const scalar YSolidTot = owner.composition().YMixture0()[idSolid];
    const scalar YLiquidTot = owner.composition().YMixture0()[idLiquid];

    switch (phase_)
    {
        case (liquid):
        {
	    CsLocalId_ = owner.composition().localId(idLiquid, "CH4N2O"); //LABEL
	    const scalar YUreaLocal = owner.composition().Y0(idLiquid)[CsLocalId_];
	    Info<< "    Urea(L): particle mass fraction = " << YUreaLocal*YLiquidTot << endl;
            break;
        }
        case (solid):
        {
	    CsLocalId_ = owner.composition().localId(idSolid, "CH4N2O"); //LABEL
	    const scalar YUreaLocal = owner.composition().Y0(idSolid)[CsLocalId_];
	    Info<< "    Urea(S): particle mass fraction = " << YUreaLocal*YSolidTot << endl;
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::scalar Foam::UreaKineticRate<CloudType>"
            )   << "Unknown phase type" << abort(FatalError);
        }
    }

    // Set local thermo properties
    WCH4N2O_ = 60.056;
    WNH3_ = owner.thermo().carrier().W(NH3GlobalId_);
    WHNCO_ = owner.thermo().carrier().W(HNCOGlobalId_);
}


template<class CloudType>
Foam::UreaKineticRate<CloudType>::
UreaKineticRate
(
    const UreaKineticRate<CloudType>& srm
)
:
    MySurfaceReactionModel<CloudType>(srm),
    A_(srm.A_),
    Dh_(srm.Dh_),
    E_(srm.E_),
    CsLocalId_(srm.CsLocalId_),
    NH3GlobalId_(srm.NH3GlobalId_),
    HNCOGlobalId_(srm.HNCOGlobalId_),
    WCH4N2O_(srm.WCH4N2O_),
    WNH3_(srm.WNH3_),
    WHNCO_(srm.WHNCO_),
    phase_(srm.phase_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::UreaKineticRate<CloudType>::
~UreaKineticRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::UreaKineticRate<CloudType>::calculate
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
) const
{
    // Mass Fraction of remaining combustible material
    const label idSolid = CloudType::parcelType::SLD;
    const label idLiquid = CloudType::parcelType::LIQ;

    scalar Yurea = 0;
    switch (phase_)
    {
        case (liquid):
        {
	    Yurea = YMixture[idLiquid]*YLiquid[CsLocalId_];
            break;
        }
        case (solid):
        {
	    Yurea = YMixture[idSolid]*YSolid[CsLocalId_];
            break;
        }
    }


    // NO Surface reaction if Urea fraction is small
    if (Yurea < SMALL || T < 425)
    {
        return 0.0;
    }


    // Kinetic rate
    const scalar Rk = A_*exp(-E_/(specie::RR*T));

    // Change in Urea mass [kg]
    scalar dmUrea = constant::mathematical::pi*d*Rk*dt; //Birkhold
//    scalar dmUrea = constant::mathematical::pi*sqr(d)/2*Rk*Yurea*dt; //AbuRamadan surface
//    scalar dmUrea = constant::mathematical::pi*pow(d,3)/6*rho*Rk*Yurea*dt; //Grunwald + AbuRamadan volumetric

    // Limit mass transfer by availability of Urea
    dmUrea = min(mass*Yurea, dmUrea);

    // Molar consumption
    const scalar dOmega = dmUrea/WCH4N2O_;

    // Mass of newly created NH3 [kg]
    const scalar dmNH3 = dOmega*WNH3_;

    // Mass of newly created HNCO [kg]
    const scalar dmHNCO = dOmega*WHNCO_;

    // Update local particle Urea mass
    scalar HsUrea = 0;
    switch (phase_)
    {
        case (liquid):
        {
	    dMassLiquid[CsLocalId_] += dmUrea;
	    HsUrea = dmUrea*(235.8/2*(pow(T,2)-pow(298,2)) - 2.173e4*(T-298) - 0.9121/3*(pow(T,3)-pow(298,3)) + 0.001627/4*(pow(T,4)-pow(298,4)) - 1.08e-6/5*(pow(T,5)-pow(298,5))); //from AdBLue and water Cp (mixture)
            break;
        }
        case (solid):
        {
	    dMassSolid[CsLocalId_] += dmUrea;
	    HsUrea = dmUrea*(38.43*(T-298) + 0.0498/2*(pow(T,2)-pow(298,2)) + 0.000705/3*(pow(T,3)-pow(298,3)) - 0.000000861/4*(pow(T,4)-pow(298,4)))/WCH4N2O_*1000; //http://webbook.nist.gov
            break;
        }
    }

    // Update carrier O2 and CO2 mass
    dMassSRCarrier[NH3GlobalId_] += dmNH3;
    dMassSRCarrier[HNCOGlobalId_] += dmHNCO;

    // carrier sensible enthalpy exchange handled via change in mass
    const scalar HsHNCO = dmHNCO*(7.823681*(T-298) + 15.2521/2000*(pow(T,2)-pow(298,2)) - 9.1925/3000000*(pow(T,3)-pow(298,3)) + 2.325/4000000000*(pow(T,4)-pow(298,4)) + 0.075905*1000000*(1/T-1/298))/WHNCO_*1000*4.184; //http://webbook.nist.gov
    const scalar HsNH3 = dmNH3*(19.99563*(T-298)+49.77119/2000*(pow(T,2)-pow(298,2)) -15.37599/3000000*(pow(T,3)-pow(298,3)) +1.921168/4000000000*(pow(T,4)-pow(298,4)) -0.189174*1000000*(1/T-1/298))/WNH3_*1000; //http://webbook.nist.gov

    
    // Heat of reaction [J]
    return dmUrea*Dh_ + HsHNCO + HsNH3 - HsUrea;
}


// ************************************************************************* //
