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

\*---------------------------------------------------------------------------*/

#include "NewLiquidEvaporation.H"
#include "specie.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::tmp<Foam::scalarField> Foam::NewLiquidEvaporation<CloudType>::calcXc
(
    const label cellI
) const
{
    scalarField Xc(this->owner().thermo().carrier().Y().size());

    forAll(Xc, i)
    {
        Xc[i] = this->owner().thermo().carrier().Y()[i][cellI]/this->owner().thermo().carrier().W(i);
    }

    return Xc/sum(Xc);
}


template<class CloudType>
Foam::scalar Foam::NewLiquidEvaporation<CloudType>::Sh
(
    const scalar Re,
    const scalar Sc,
    const scalar Gr
) const
{
//    return 2.0 + 0.552*Foam::sqrt(Re)*cbrt(Sc); //Shashank,Strom
    return 2.0009 + 0.514*Foam::sqrt(max(Re,sqrt(max(Gr,0))))*cbrt(Sc); //Trujillo,Ebraihimian
}

template<class CloudType>
Foam::scalar Foam::NewLiquidEvaporation<CloudType>::Nusselt
(
    const scalar Re,
    const scalar Pr,
    const scalar Gr
) const
{
//    return 2.0 + 0.552*Foam::sqrt(Re)*cbrt(Pr);
    return 2.0009 + 0.514*Foam::sqrt(max(Re,sqrt(max(Gr,0))))*cbrt(Pr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NewLiquidEvaporation<CloudType>::NewLiquidEvaporation
(
    const dictionary& dict,
    CloudType& owner
)
:
    PhaseChangeModel<CloudType>(dict, owner, typeName),
    liquids_(owner.thermo().liquids()),
    activeLiquids_(this->coeffDict().lookup("activeLiquids")),
    liqToCarrierMap_(activeLiquids_.size(), -1),
    liqToLiqMap_(activeLiquids_.size(), -1)
{
    if (activeLiquids_.size() == 0)
    {
        WarningIn
        (
            "Foam::NewLiquidEvaporation<CloudType>::NewLiquidEvaporation"
            "("
                "const dictionary& dict, "
                "CloudType& owner"
            ")"
        )   << "Evaporation model selected, but no active liquids defined"
            << nl << endl;
    }
    else
    {
        Info<< "Participating liquid species:" << endl;

        // Determine mapping between liquid and carrier phase species
        forAll(activeLiquids_, i)
        {
            Info<< "    " << activeLiquids_[i] << endl;
            liqToCarrierMap_[i] =
                owner.composition().globalCarrierId(activeLiquids_[i]);
        }

        // Determine mapping between model active liquids and global liquids
        const label idLiquid = owner.composition().idLiquid();
        forAll(activeLiquids_, i)
        {
            liqToLiqMap_[i] =
                owner.composition().localId(idLiquid, activeLiquids_[i]);
        }
    }
}


template<class CloudType>
Foam::NewLiquidEvaporation<CloudType>::NewLiquidEvaporation
(
    const NewLiquidEvaporation<CloudType>& pcm
)
:
    PhaseChangeModel<CloudType>(pcm),
    liquids_(pcm.owner().thermo().liquids()),
    activeLiquids_(pcm.activeLiquids_),
    liqToCarrierMap_(pcm.liqToCarrierMap_),
    liqToLiqMap_(pcm.liqToLiqMap_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NewLiquidEvaporation<CloudType>::~NewLiquidEvaporation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::NewLiquidEvaporation<CloudType>::calculate
(
    const scalar dt,
    const label cellI,
    scalar& Re_No,
    const scalar Pr,
    const scalar d,
    const scalar nu,
    const scalar T,
    const scalar Ts,
    const scalar pc,
    const scalar Tc,
    const scalarField& Yl,
    scalarField& dMassPC,
    scalar& bt,
    const scalar mu
) const
{
    // GAS PHASE carrier properties
    const scalarField XcMix(calcXc(cellI));
    scalarField YcMix(this->owner().thermo().carrier().Y().size(), 0.0);
    scalar Wc = 0.0;
    scalar Cpc = 0.0;
    scalar kappac = 0.0;
    scalar muc = 0.0;
    forAll(this->owner().thermo().carrier().Y(), i)
    {
        YcMix[i] = this->owner().thermo().carrier().Y()[i][cellI];
      	scalar Wi = this->owner().thermo().carrier().W(i);
      	scalar kappai = this->owner().thermo().carrier().kappa(i, pc, Ts);
      	scalar mui = this->owner().thermo().carrier().mu(i, pc, Ts);
      	scalar kappafi = 0.0;
      	scalar mufi = 0.0;
        	forAll(this->owner().thermo().carrier().Y(), j)
        	{
        		scalar kappaj = this->owner().thermo().carrier().kappa(j, pc, Ts);
        		scalar muj = this->owner().thermo().carrier().mu(j, pc, Ts);
        		scalar Wj = this->owner().thermo().carrier().W(j);
        		kappafi += XcMix[j]*sqr(1+sqrt(kappai/kappaj)*pow(Wi/Wj,0.25))/sqrt(8*(1+Wi/Wj));
        		mufi += XcMix[j]*sqr(1+sqrt(mui/muj)*pow(Wi/Wj,0.25))/sqrt(8*(1+Wi/Wj));
        	} //Mason&Saxena (Poling)
        kappac += (XcMix[i]*kappai)/kappafi;
        muc += (XcMix[i]*mui)/mufi;
        Wc += XcMix[i]*Wi;
        Cpc += YcMix[i]*this->owner().thermo().carrier().Cp(i, pc, Ts);
    }
    scalar rhoc = pc/(specie::RR*Ts/Wc);
    scalar nuc = muc/rhoc;

    const scalarField X(liquids_.X(Yl));
    if ((liquids_.Tc(X) - T) < SMALL) // Immediately evaporate mass in critical condition (for water is 374°C)
    {
        if (debug)
        {
            WarningIn
            (
                "void Foam::NewLiquidEvaporation<CloudType>::calculate"
                "("
                    "const scalar, "
                    "const label, "
                    "const scalar, "
                    "const scalar, "
                    "const scalar, "
                    "const scalar, "
                    "const scalar, "
                    "const scalar, "
                    "const scalar, "
                    "const scalar, "
                    "const scalarField&, "
                    "scalarField&"
                ") const"
            )   << "Parcel reached critical conditions: "
                << "evaporating all avaliable mass" << endl;
        }
        forAll(activeLiquids_, i)
        {
            const label lid = liqToLiqMap_[i];
            dMassPC[lid] = GREAT;
        }
        return;
    }
    forAll(activeLiquids_, i)
    {
        const label gid = liqToCarrierMap_[i];
        const label lid = liqToLiqMap_[i];
        const scalar TBoil = liquids_.properties()[lid].pvInvert(pc);// boiling temperature at cell pressure for liquid species lid [K]
        const scalar Tfilm = min(Ts, 0.999*TBoil);   // limit surface temperature to boiling temperature (e.g 373K for water at pc=1atm)

	// Concentration of the ACTIVE LIQUID in CARRIER phase
        const scalar Xc = XcMix[gid];
        const scalar Yc = YcMix[gid]; // carrier phase concentration

        // VAPOR properties (ACTIVE liquids)
        const scalar pSat0 = liquids_.properties()[lid].pv(pc, T); //Ebrahimian,Ramadan,Shashank,Sirignano at T
      	const scalar pVap = X[lid]*pSat0;
        const scalar CpVap = liquids_.properties()[lid].Cpg(pVap, Tfilm); // specific heat capacity for lid - Shashank measured at Tfilm
      	const scalar rhoVap = pVap*liquids_.properties()[lid].W()/(specie::RR*Tfilm);
        const scalar muVap = liquids_.properties()[lid].mug(pVap, Tfilm); // vapour dynamic viscosity - Shashank measured at Tfilm
        const scalar kappaVap = liquids_.properties()[lid].Kg(pVap, Tfilm); // vapour thermal conductivity - Shashank measured at Tfilm
        const scalar WVap = liquids_.properties()[lid].W(); // molecular weight
        const scalar DabVap = liquids_.properties()[lid].D(pVap, Tfilm, Wc); // vapour diffusivity in air [m2/s] - Shashank measured at Tfilm

      	const scalar Xs = X[lid]*pSat0/pc; // surface VAPOUR molar fraction - Raoult's Law
      //	const scalar Xs = pc/101325*exp(-liquids_.properties()[lid].hl(pc, Td)/8.3145*liquids_.properties()[lid].W()*(1/Tfilm-1/liquids_.properties()[lid].pvInvert(pc))); //Shashank
      	const scalar Ys = WVap*Xs/(WVap*Xs + Wc*(1 - Xs)); //assuming at FILM only VAPOUR (Xs) & "freestream" gas (1-Xs)

	// REFERENCE (FILM) COMPOSITION
        const scalar Xfilm = (2*Xs + Xc)/3; //Shanshank
        const scalar Yfilm = (2*Ys + Yc)/3;

      	const scalar CpFilm = CpVap*Yfilm + Cpc*(1 - Yfilm); //Shashank,Ramadan,Sirignano
      	const scalar rhoFilm = 1/(Yfilm/rhoVap + (1 - Yfilm)/rhoc); //Sirignano
      	const scalar lfifg = sqr(1+sqrt(kappaVap/kappac)*pow(Wc/WVap,0.25))/sqrt(8*(1+WVap/Wc)); //Shashank,Ramadan (Mason)
      	const scalar lfigf = sqr(1+sqrt(kappac/kappaVap)*pow(WVap/Wc,0.25))/sqrt(8*(1+Wc/WVap)); //Shashank,Ramadan (Mason)
      	const scalar kappaFilm = kappaVap*Xfilm/(Xfilm + (1 - Xfilm)*lfifg) + kappac*(1 - Xfilm)/(1 - Xfilm + Xfilm*lfigf); //Shashank,Ramadan (Mason)

      //	const scalar muFilm = muVap*Xfilm/(Xfilm + (1 - Xfilm)*sqrt(Wc/WVap)) + muc*(1 - Xfilm)/(1 - Xfilm + Xfilm*sqrt(WVap/Wc)); //Strom (Zipperer)
      	const scalar mfifg = sqr(1+sqrt(muVap/muc)*pow(Wc/WVap,0.25))/sqrt(8*(1+WVap/Wc)); //Shashank,Ramadan (Wilke)
      	const scalar mfigf = sqr(1+sqrt(muc/muVap)*pow(WVap/Wc,0.25))/sqrt(8*(1+Wc/WVap)); //Shashank,Ramadan (Wilke)
      	const scalar muFilm = muVap*Xfilm/(Xfilm + (1 - Xfilm)*mfifg) + muc*(1 - Xfilm)/(1 - Xfilm + Xfilm*mfigf); //Shashank,Ramadan (Wilke)

      	const scalar nuFilm = muFilm/rhoFilm;

	// TRANSFER PARAMETERS
        const scalar ScFilm = nuFilm/(DabVap + ROOTVSMALL); // Ebrahimian
//        const scalar ScFilm = nuc/(Dab + ROOTVSMALL); // only Gas
      	const scalar Gr = 9.81*pow(d,3)*(Tc-T)/sqr(nuc)/Tc; //Trujillo,Ebraihimian
      //	const scalar Gr = 0; //altri
      	Re_No = Re_No*mu/muFilm; //Strom,Sirignano
      //        const scalar Prc = muFilm*CpFilm/kappac; //Shashank
        const scalar Prc = muFilm*CpFilm/kappaFilm; //Ramadan
      	const scalar Le = ScFilm/Prc; // Lewis Number
      	const scalar Nu0 = this->Nusselt(Re_No, Prc, Gr); //only gas phase for Nusselt (Ramadan)
        const scalar Sh0 = this->Sh(Re_No, ScFilm, Gr);  // Sherwood number
        const scalar Bm = (Ys - Yc)/max(SMALL, 1.0 - Ys); // mass transfer number
      	scalar Nu = Nu0;
      	scalar Sh = Sh0;
        if (Bm > 0)
        {
            const scalar Fm = pow((1 + Bm), 0.7)/Bm*log(1 + Bm); // Strom - mass transfer correction factor
//                    const scalar Fm = pow((1 + Bm), 0.7); // Ramadan
//                    const scalar Fm = 1; // Shashank
            Sh = 2 + (Sh0-2)/Fm; // Strom - corrected Sherwood
//                    const scalar Sh = Sh0/Fm; // Ramadan,Shashank

//Iteration for Bt:
            scalar Ft = 0.0; // heat transfer correction factor
            Nu = 0.0; // corrected Nusselt
            scalar psiB = 0.0; // Exponent
      	    scalar eps = 1.0; // Relative Error
      	    scalar Btiter = bt; //Bt of previous timestep
    		    for (label iter = 0; iter < 5; iter++)
            {
      		      if (eps > 5e-2)
                {
                    Ft = pow((1 + Btiter), 0.7)/Btiter*log(1 + Btiter); //Strom
      //                        Ft = pow((1 + Btiter), 0.7); //Ramadan
      //                        Ft = 1; //Shashank
              			Nu = 2 + (Nu0-2)/Ft; //Strom
              //			Nu = Nu0/Ft; //Ramadan,Shashank
              			psiB = CpVap/CpFilm*Sh/Nu/Le; //Strom,Shashank
              			eps = std::abs((pow((1 + Bm), psiB)-1)/Btiter-1);
              			Btiter = (pow((1 + Bm), psiB)-1);
                }
            }
    		    bt = Btiter;
		      }

        if (Xc*pc > pSat0)  //specie partial pressure Xc*pc
        {
            		// saturated vapour - no phase change
        }
        else
        {
            if (Xs > 0.999) //pvap=pc
            {		// boiling
                dMassPC[lid] = pi*DabVap*kappac*Nu/CpVap*log(1.0 + bt)*dt; // Birkhold
            }
            else
            {		// evaporation
                dMassPC[lid] = pi*d*Sh*DabVap*rhoFilm*dt*log(1.0 + Bm); // mass transfer [kg] - Strom,Shashank
            }
        }
    }
}


template<class CloudType>
Foam::scalar Foam::NewLiquidEvaporation<CloudType>::dh
(
    const label idc,
    const label idl,
    const scalar p,
    const scalar T
) const
{
    scalar dh = 0;

    scalar TDash = T;
    if (liquids_.properties()[idl].pv(p, T) >= 0.999*p)
    {
        TDash = liquids_.properties()[idl].pvInvert(p);
    }

    typedef PhaseChangeModel<CloudType> parent;
    switch (parent::enthalpyTransfer_)
    {
        case (parent::etLatentHeat):
        {
            dh = liquids_.properties()[idl].hl(p, TDash);
            break;
        }
        case (parent::etEnthalpyDifference):
        {
            scalar hc = this->owner().composition().carrier().Ha(idc, p, TDash);
            scalar hp = liquids_.properties()[idl].h(p, TDash);

            dh = hc - hp;
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::scalar Foam::NewLiquidEvaporation<CloudType>::dh"
                "("
                    "const label, "
                    "const label, "
                    "const scalar, "
                    "const scalar"
                ") const"
            )   << "Unknown enthalpyTransfer type" << abort(FatalError);
        }
    }

    return dh;
}


template<class CloudType>
Foam::scalar Foam::NewLiquidEvaporation<CloudType>::Tvap
(
    const scalarField& Y
) const
{
    const scalarField X(liquids_.X(Y));

    return liquids_.Tpt(X);
}


template<class CloudType>
Foam::scalar Foam::NewLiquidEvaporation<CloudType>::TMax
(
    const scalar p,
    const scalarField& Y
) const
{
    const scalarField X(liquids_.X(Y));

    return liquids_.pvInvert(p, X);
}


// ************************************************************************* //
