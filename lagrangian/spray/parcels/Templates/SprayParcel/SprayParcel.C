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

#include "SprayParcel.H"
#include "CompositionModel.H"
#include "AtomizationModel.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
const Foam::label Foam::SprayParcel<ParcelType>::GAS(0);

template<class ParcelType>
const Foam::label Foam::SprayParcel<ParcelType>::LIQ(1);

template<class ParcelType>
const Foam::label Foam::SprayParcel<ParcelType>::SLD(2);


// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
Foam::scalar Foam::SprayParcel<ParcelType>::CpEff
(
    TrackData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*td.cloud().composition().Cp(idG, YGas_, p, T)
      + this->Y_[LIQ]*td.cloud().composition().Cp(idL, YLiquid_, p, T) //mass fraction (it's not calculated in liquidMixtureProperties)
      + this->Y_[SLD]*td.cloud().composition().Cp(idS, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackData>
Foam::scalar Foam::SprayParcel<ParcelType>::rhoEff
(
    TrackData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{

    scalar rho = 0.0;
    forAll(YGas_, i)
    {
        const phaseProperties& propsGas = td.cloud().composition().phaseProps()[idG];
        label gid = propsGas.globalIds()[i];
        rho += (1/YGas_[i]*td.cloud().composition().thermo().carrier().rho(gid, p, T))/this->Y_[GAS];
    }
    forAll(YLiquid_, i)
    {
        const phaseProperties& propsLiq = td.cloud().composition().phaseProps()[idL];
        label gid = propsLiq.globalIds()[i];
        rho += (1/YLiquid_[i]*td.cloud().composition().thermo().liquids().properties()[gid].rho(p, T))/this->Y_[LIQ];
    }
    forAll(YSolid_, i)
    {
        const phaseProperties& propsSol = td.cloud().composition().phaseProps()[idS];
        label gid = propsSol.globalIds()[i];
        rho += (1/YSolid_[i]*td.cloud().composition().thermo().solids().properties()[gid].rho())/this->Y_[SLD];
    }

    return rho;
}


template<class ParcelType>
template<class TrackData>
Foam::scalar Foam::SprayParcel<ParcelType>::HsEff
(
    TrackData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*td.cloud().composition().Hs(idG, YGas_, p, T)
      + this->Y_[LIQ]*td.cloud().composition().Hs(idL, YLiquid_, p, T)
      + this->Y_[SLD]*td.cloud().composition().Hs(idS, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackData>
Foam::scalar Foam::SprayParcel<ParcelType>::LEff
(
    TrackData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*td.cloud().composition().L(idG, YGas_, p, T)
      + this->Y_[LIQ]*td.cloud().composition().L(idL, YLiquid_, p, T)
      + this->Y_[SLD]*td.cloud().composition().L(idS, YSolid_, p, T);
}


template<class ParcelType>
Foam::scalar Foam::SprayParcel<ParcelType>::updateMassFractions
(
    const scalar mass0,
    const scalarField& dMassGas,
    const scalarField& dMassLiquid,
    const scalarField& dMassSolid
)
{
    scalarField& YMix = this->Y_;

    scalar massGas = this->updateMassFraction(mass0*YMix[GAS], dMassGas, YGas_); //ReactingP
    scalar massLiquid = this->updateMassFraction(mass0*YMix[LIQ], dMassLiquid, YLiquid_);
    scalar massSolid = this->updateMassFraction(mass0*YMix[SLD], dMassSolid, YSolid_);

    scalar massNew = max(massGas + massLiquid + massSolid, ROOTVSMALL);

    YMix[GAS] = massGas/massNew;
    YMix[LIQ] = massLiquid/massNew;
    YMix[SLD] = 1.0 - YMix[GAS] - YMix[LIQ];

    return massNew;
}



// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    ParcelType::setCellValues(td, dt, cellI); //all methods in Template classes
}


template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::cellValueSourceCorrection
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    ParcelType::cellValueSourceCorrection(td, dt, cellI); //only ReactingParcel
}


template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    typedef typename TrackData::cloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition = td.cloud().composition();

    // 1. check if parcel belongs to liquid core
    if (liquidCore() > 0.5)
    {
        td.cloud().forces().setCalcCoupled(false); // these parcels should not experience coupled forces
    }


    // 2. get old properties (check if we have critical or boiling conditions)
    scalarField& YMix = this->Y_; //MASS FRACTION OF EACH PHASE (1x3)
    const scalarField& YMix0 = YMix;
    const label idG = composition.idGas();
    const label idL = composition.idLiquid();
    const label idS = composition.idSolid();
    const scalarField& Ysol0 = this->YSolid(); //MASS FRACTION OF SOLID PHASE (1xNs)
    const scalarField& Yliq0 = this->YLiquid(); //MASS FRACTION OF LIQUID PHASE (1xNl)
    const scalarField& Ygas0 = this->YGas(); ///MASS FRACTION OF GAS PHASE (1xNg)
    const scalar T0 = this->T();
    const scalar pc0 = this->pc_;
    const scalar cp0 = this->Cp();
    const scalar rho0 = this->rho(); //FROM LAST TIMESTEP or SPRAYCLOUD
    const scalar np0 = this->nParticle_;
    const scalar d0 = this->d_;
    const vector& U0 = this->U_;
    const scalar mass0 = this->mass();
    scalar bt = this->Bt_;

    if ((YMix[LIQ] - 1.0) < SMALL)
    {
       scalarField X0(composition.liquids().X(Yliq0));
       scalar TMax = composition.liquids().Tc(X0);
       if (composition.liquids().pv(pc0, T0, X0) >= pc0*0.999)
       {
           TMax = composition.liquids().pvInvert(pc0, X0);  // set TMax to boiling temperature
       }
       td.cloud().constProps().TMax() = TMax;
    }


    // 3. calculate surface properties (!)
    scalar Ts, rhos, mus, Prs, kappas;
    this->calcSurfaceValues(td, cellI, T0, Ts, rhos, mus, Prs, kappas);
    scalar Res = this->Re(U0, d0, this->rhoc_, mus);


    // 4. Initialize Source Terms and Phase change terms
    vector Su = vector::zero; // Explicit momentum source for particle
    scalar Spu = 0.0; // Linearised momentum source coefficient
    vector dUTrans = vector::zero; // Momentum transfer from the particle to the carrier phase
    scalar Sh = 0.0; // Explicit enthalpy source for particle
    scalar Sph = 0.0; // Linearised enthalpy source coefficient
    scalar dhsTrans = 0.0; // Sensible enthalpy transfer from particle to the carrier phase

    scalar Ne = 0.0; // Molar flux of species emitted from the particle (kmol/m^2/s)
    scalar NCpW = 0.0; // Sum Ni*Cpi*Wi of emission species
    scalarField Cs(composition.carrier().species().size(), 0.0); // Surface concentrations of emitted species
    scalarField dMassPC(YLiquid_.size(), 0.0); //only for LIQUIDS
    scalarField dMassDV(YGas_.size(), 0.0);
    scalarField dMassSRGas(YGas_.size(), 0.0);
    scalarField dMassSRLiquid(YLiquid_.size(), 0.0);
    scalarField dMassSRSolid(YSolid_.size(), 0.0);
    scalarField dMassSRCarrier(composition.carrier().species().size(), 0.0);


    // 5. Evaporation
    newCalcPhaseChange
    (
        td,
        dt,
        cellI,
        Res,
        Prs,
        Ts,
        mus/rhos,
        mus,
        d0,
        T0,
	pc0,
        mass0,
        idL, //idPhase
        YMix[LIQ], //TOTAL Liquid fraction of the droplet
        YLiquid_, //VECTOR LIQUID SPECIES
        dMassPC,
        Sh,
        Ne,
        NCpW,
        Cs,
	bt
    );
    this->Bt_ = bt;

    // 6. Devolatilisation
    newCalcDevolatilisation
    (
        td,
        dt,
        this->age_,
        Ts,
        d0,
        T0,
        mass0,
        this->mass0_,
        YMix[GAS]*YGas_,
        YMix[LIQ]*YLiquid_,
        YMix[SLD]*YSolid_,
        canCombust_,
        dMassDV,
        Sh,
        Ne,
        NCpW,
        Cs
    );


    // 7. Surface reactions
    newCalcSurfaceReactions
    (
        td,
        dt,
        cellI,
        d0,
        T0,
        mass0,
        canCombust_,
        Ne,
        YMix,
        YGas_,
        YLiquid_,
        YSolid_,
        dMassSRGas,
        dMassSRLiquid,
        dMassSRSolid,
        dMassSRCarrier,
        Sh,
        dhsTrans
    );


    // 8. Update the parcel mass and mass fractions (YGas_,YLiquid_,YSolid_,YMix) due to phase change
    scalarField dMassGas(dMassDV + dMassSRGas);
    scalarField dMassLiquid(dMassPC + dMassSRLiquid);
    scalarField dMassSolid(dMassSRSolid);

    scalar mass1 = updateMassFractions(mass0, dMassGas, dMassLiquid, dMassSolid);


    // Remove the particle when mass falls below minimum threshold
    if (np0*mass1 < td.cloud().constProps().minParticleMass())
    {
        td.keepParticle = false; // Absorb parcel in carrier
        if (td.cloud().solution().coupled())
        {
            scalar dm = np0*mass0;

            forAll(Ygas0, i) //not YGas_ (already updated!!!!!!!!)
            {
                label gid = composition.localToGlobalCarrierId(GAS, i);
                td.cloud().rhoTrans(gid)[cellI] += dm*YMix0[GAS]*Ygas0[i];
            }
            forAll(Yliq0, i)
            {
                label gid = composition.localToGlobalCarrierId(LIQ, i);
                td.cloud().rhoTrans(gid)[cellI] += dm*YMix0[LIQ]*Yliq0[i];
            }
            /*
            forAll(YSolid_, i)
            {
                label gid = composition.localToGlobalCarrierId(SLD, i);
                td.cloud().rhoTrans(gid)[cellI] += dm*YMix[SLD]*YSolid_[i];
            }*/

            td.cloud().UTrans()[cellI] += dm*U0;
            td.cloud().hsTrans()[cellI] += dm*HsEff(td, pc0, T0, idG, idL, idS);
            td.cloud().phaseChange().addToPhaseChangeMass(np0*mass0);
        }
        return;
    }


    // 9. Update the parcel properties (Cp, rho, d)
    this->Cp_ = CpEff(td, pc0, T0, idG, idL, idS);

    if (td.cloud().constProps().constantVolume())
    {
        this->rho_ = mass1/this->volume();
    }
    else
    {
        /*if ((YMix[LIQ] - 1.0) < SMALL)
        {
	       scalarField X1 = composition.X(idL,YLiquid_);
	       this->rho() = composition.liquids().rho(pc0, T0, X1);
        }
        else
        {
	       scalarField X1 = composition.solids().X(YSolid_);
	       this->rho() = composition.solids().rho(X1);
        }*/ //ALTERNATIVE
        this->rho_ = rhoEff(td, pc0, T0, idG, idL, idS);
        this->d_ = cbrt(mass1/this->rho_*6.0/pi);
    }

    // 10. Correct surface values due to emitted species and calculate T
    ParcelType::correctSurfaceValues(td, cellI, Ts, Cs, rhos, mus, Prs, kappas); //caso aggiornato
    Res = this->Re(U0, this->d_, rhos, mus);
    this->T_ =
        this->calcHeatTransfer
        (
            td,
            dt,
            cellI,
            Res,
            Prs,
            kappas,
            NCpW,
            Sh,
            dhsTrans,
            Sph
        );


    // 11. Calculate new particle velocity
    scalar Wec = 0;
    if ((YMix[LIQ] - 1.0) < SMALL)
    {
	scalarField X0 = composition.liquids().X(Yliq0); //NOT UPDATED
	scalar sigma0 = composition.liquids().sigma(pc0, T0, X0);
	Wec = this->We(U0, this->d_, this->rhoc_, sigma0);
	//scalarField X1(composition().X(idL,YLiquid_)); //UPDATED
	//scalar sigma1 = composition.liquids().sigma(pc0, this->T_, X1);
	//Wec = this->We(U0, this->d_, this->rhoc_, sigma1);
    }
    this->Renolds() = Res;
    this->Weber() = Wec;
    this->U_ = newCalcVelocity(td, dt, cellI, Res, this->y(), Wec, mus, mass1, YMix[LIQ], Su, dUTrans, Spu, bt); //UPDATED
    //this->U_ = newCalcVelocity(td, dt, cellI, Res, this->y(), Wec, mus, mass0, YMix[LIQ], Su, dUTrans, Spu, bt); //NOT UPDATED


    // 12. Update properties due to change in temperature and/or composition
    scalar T1 = this->T();
    this->Cp_ = CpEff(td, this->pc_, this->T(), idG, idL, idS);

    /*scalar rho1 = 0;
    if ((YMix[LIQ] - 1.0) < SMALL)
    {
        scalarField X1 = composition.X(idL,YLiquid_);
        rho1 = composition.liquids().rho(this->pc_, T1, X1);
    }
    else
    {
	scalarField X1 = composition.solids().X(YSolid_);
        rho1 = composition.solids().rho(X1);
    }*/ //ALTERNATIVE
    scalar rho1 = rhoEff(td, this->pc_, T1, idG, idL, idS);

    this->rho() = rho1;
    scalar d1 = cbrt(mass1/rho1*6.0/pi); //mass CHANGES!!!!
    this->d() = d1;


    // 13. Accumulate carrier phase source terms
    if (td.cloud().solution().coupled())
    {
        // Transfer mass lost to carrier mass, momentum and enthalpy sources
        forAll(YGas_, i)
        {
            scalar dm = np0*dMassGas[i];
            label gid = composition.localToGlobalCarrierId(GAS, i);
            scalar hs = composition.carrier().Hs(gid, pc0, T1); //UPDATED T
            td.cloud().rhoTrans(gid)[cellI] += dm;
            td.cloud().UTrans()[cellI] += dm*this->U(); //UPDATED U
            td.cloud().hsTrans()[cellI] += dm*hs;
        }
        forAll(YLiquid_, i)
        {
            scalar dm = np0*dMassLiquid[i];
            label gid = composition.localToGlobalCarrierId(LIQ, i);
            scalar hs = composition.carrier().Hs(gid, pc0, T1); //UPDATED T
            td.cloud().rhoTrans(gid)[cellI] += dm;
            td.cloud().UTrans()[cellI] += dm*this->U(); //UPDATED U
            td.cloud().hsTrans()[cellI] += dm*hs;
        }
	/*// No mapping between solid components and carrier phase
        forAll(YSolid_, i)
        {
            scalar dm = np0*dMassSolid[i];
            label gid = composition.localToGlobalCarrierId(SLD, i);
            scalar hs = composition.carrier().Hs(gid, pc0, T0);
            td.cloud().rhoTrans(gid)[cellI] += dm;
            td.cloud().UTrans()[cellI] += dm*U0;
            td.cloud().hsTrans()[cellI] += dm*hs;
        }*/
        forAll(dMassSRCarrier, i)
        {
            scalar dm = np0*dMassSRCarrier[i];
            scalar hs = composition.carrier().Hs(i, pc0, T1); //ho aggiornato T
            td.cloud().rhoTrans(i)[cellI] += dm;
            td.cloud().UTrans()[cellI] += dm*this->U(); //ho aggiornato U
            td.cloud().hsTrans()[cellI] += dm*hs;
        }
        // Update momentum transfer
        td.cloud().UTrans()[cellI] += np0*dUTrans;
        td.cloud().UCoeff()[cellI] += np0*Spu;
        // Update sensible enthalpy transfer
        td.cloud().hsTrans()[cellI] += np0*dhsTrans;
        td.cloud().hsCoeff()[cellI] += np0*Sph;
        // Update radiation fields
        if (td.cloud().radiation())
        {
            const scalar ap = this->areaP();
            const scalar T4 = pow4(this->T_);
            td.cloud().radAreaP()[cellI] += dt*np0*ap;
            td.cloud().radT4()[cellI] += dt*np0*T4;
            td.cloud().radAreaPT4()[cellI] += dt*np0*ap*T4;
        }
    }


    // 14. Calc breakup
    if ((YMix[LIQ] - 1.0) < SMALL)
    {
      if (td.keepParticle)
      {
        if (liquidCore() > 0.5)
        {
            calcAtomization(td, dt, cellI); //preserve the total mass/volume by increasing the number of particles in parcels due to breakup
            scalar d2 = this->d();
            this->nParticle() *= pow3(d1/d2);
        }
        else
        {
            calcBreakup(td, dt, cellI);
        }
      }
    }


    // 15. restore coupled forces if switched off by liquidCore>0.5
    td.cloud().forces().setCalcCoupled(true);
}


template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::newCalcPhaseChange
(
    TrackData& td,
    const scalar dt,
    const label cellI,
    scalar& Re_No,
    const scalar Pr,
    const scalar Ts,
    const scalar nus,
    const scalar mus,
    const scalar d,
    const scalar T,
    const scalar pc,
    const scalar mass,
    const label idPhase,
    const scalar YPhase,
    const scalarField& YComponents,
    scalarField& dMassPC,
    scalar& Sh,
    scalar& N,
    scalar& NCpW,
    scalarField& Cs,
    scalar& bt
)
{
    typedef typename TrackData::cloudType::reactingCloudType reactingCloudType;
    PhaseChangeModel<reactingCloudType>& phaseChange = td.cloud().phaseChange();

    scalar Tvap = phaseChange.Tvap(YComponents); // = -GREAT, so not useful
    if (!phaseChange.active() || T < Tvap || YPhase < SMALL)
    {
        return;
    }
    const scalar TMax = phaseChange.TMax(pc, YComponents); // = GREAT, so not useful
    const scalar Tdash = min(T, TMax);
    const scalar Tsdash = min(Ts, TMax);

    // Calculate mass transfer due to phase change
    phaseChange.calculate
    (
        dt,
        cellI,
        Re_No,
        Pr,
        d,
        nus,
        Tdash,
        Tsdash,
        pc,
        this->Tc_,
        YComponents,
        dMassPC,
      	bt,
        mus
    );

    dMassPC = min(mass*YPhase*YComponents, dMassPC); // Limit phase change mass by availability of each specie in the PARTICLE
    const scalar dMassTot = sum(dMassPC);

    phaseChange.addToPhaseChangeMass(this->nParticle_*dMassTot); // Add to CUMULATIVE (PARCEL) phase change mass
    const CompositionModel<reactingCloudType>& composition = td.cloud().composition();

    // Calculate energy transfer due to phase change
    forAll(dMassPC, i)
    {
        const label idc = composition.localToGlobalCarrierId(idPhase, i);
        const label idl = composition.globalIds(idPhase)[i];

        const scalar dh = phaseChange.dh(idc, idl, pc, Tdash);
        Sh -= dMassPC[i]*dh/dt;
    }

    // Update molar emissions
    if (td.cloud().heatTransfer().BirdCorrection())
    {
        const scalar Wc = this->rhoc_*specie::RR*this->Tc_/pc; // Average molecular weight of carrier - assumes perfect gas
        forAll(dMassPC, i)
        {
            const label idc = composition.localToGlobalCarrierId(idPhase, i); //ID of the liquid specie in the carrier phase
            const label idl = composition.globalIds(idPhase)[i]; //ID of the phase (liquid)

            const scalar Cp = composition.carrier().Cp(idc, pc, Tsdash); //specie Cp in carrier (NOT vapour)
            const scalar W = composition.carrier().W(idc);
            const scalar Ni = dMassPC[i]/(this->areaS(d)*dt*W);
            const scalar Dab = composition.liquids().properties()[idl].D(pc, Tsdash, Wc); //diff coeff. (DIFFERENT PRESSURE in *LiquidEvaporationModel)

            N += Ni; // Molar flux of species coming from the particle (kmol/m^2/s)
            NCpW += Ni*Cp*W; // Sum of Ni*Cpi*Wi of emission species
            Cs[idc] += Ni*d/(2.0*Dab); // Concentrations of emission species
        }
    }
}


template<class ParcelType>
template<class TrackData>
const Foam::vector Foam::SprayParcel<ParcelType>::newCalcVelocity
(
    TrackData& td,
    const scalar dt,
    const label cellI,
    const scalar Re,
    const scalar y,
    const scalar We,
    const scalar mu,
    const scalar mass,
    const scalar YPhase,
    const vector& Su,
    vector& dUTrans,
    scalar& Spu,
    const scalar bt
) const
{
    typedef typename TrackData::cloudType cloudType;
    typedef typename cloudType::parcelType parcelType;
    typedef typename cloudType::forceType forceType;
    const forceType& forces = td.cloud().forces();

    // Particle forces
    const parcelType& p = static_cast<const parcelType&>(*this);
    const forceSuSp Fcp = forces.calcCoupled(p, dt, mass, Re, mu, y, YPhase, bt, We);
    const forceSuSp Fncp = forces.calcNonCoupled(p, dt, mass, Re, mu);
    const forceSuSp Feff = Fcp + Fncp;
    const scalar massEff = forces.massEff(p, mass);

    // New particle velocity
    const vector abp = (Feff.Sp()*this->Uc() + (Feff.Su() + Su))/massEff;
    const scalar bp = Feff.Sp()/massEff;
    IntegrationScheme<vector>::integrationResult Ures = td.cloud().UIntegrator().integrate(this->U_, dt, abp, bp);
    vector Unew = Ures.value();

    // Single source terms
    Spu = dt*Feff.Sp();
    dUTrans += dt*(Feff.Sp()*(Ures.average() - this->Uc()) - Fcp.Su());

    // Apply correction to velocity and dUTrans for reduced-D cases
    const polyMesh& mesh = td.cloud().pMesh();
    meshTools::constrainDirection(mesh, mesh.solutionD(), Unew);
    meshTools::constrainDirection(mesh, mesh.solutionD(), dUTrans);

    return Unew;

}


template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::calcAtomization
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    typedef typename TrackData::cloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition = td.cloud().composition();
    typedef typename TrackData::cloudType::sprayCloudType sprayCloudType;
    const AtomizationModel<sprayCloudType>& atomization = td.cloud().atomization();
       // cell state info is updated in ReactingParcel calc
    const scalarField& Y = this->YLiquid();
    scalarField X(composition.liquids().X(Y));
    scalar rho = composition.liquids().rho(this->pc(), this->T(), X);
    scalar mu = composition.liquids().mu(this->pc(), this->T(), X);
    scalar sigma = composition.liquids().sigma(this->pc(), this->T(), X);
       // Average molecular weight of carrier mix - assumes perfect gas
    scalar Wc = this->rhoc_*specie::RR*this->Tc()/this->pc();
    scalar R = specie::RR/Wc;
    scalar Tav = atomization.Taverage(this->T(), this->Tc());
       // calculate average gas density based on average temperature
    scalar rhoAv = this->pc()/(R*Tav);
    scalar soi = td.cloud().injectors().timeStart();
    scalar currentTime = td.cloud().db().time().value();
    const vector& pos = this->position();
    const vector& injectionPos = this->position0();
       // disregard the continous phase when calculating the relative velocity (in line with the deactivated coupled assumption)
    scalar Urel = mag(this->U());
    scalar t0 = max(0.0, currentTime - this->age() - soi);
    scalar t1 = min(t0 + dt, td.cloud().injectors().timeEnd() - soi);
       // this should be the vol flow rate from when the parcel was injected
    scalar volFlowRate = td.cloud().injectors().volumeToInject(t0, t1)/dt;
    scalar chi = 0.0;
    if (atomization.calcChi())
    {
        chi = this->chi(td, X);
    }
    atomization.update
    (
        dt,
        this->d(),
        this->liquidCore(),
        this->tc(),
        rho,
        mu,
        sigma,
        volFlowRate,
        rhoAv,
        Urel,
        pos,
        injectionPos,
        td.cloud().pAmbient(),
        chi,
        td.cloud().rndGen()
    );
}


template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::calcBreakup
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    typedef typename TrackData::cloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition = td.cloud().composition();
    typedef typename TrackData::cloudType cloudType;
    typedef typename cloudType::parcelType parcelType;
    typedef typename cloudType::forceType forceType;

    const parcelType& p = static_cast<const parcelType&>(*this);
    const forceType& forces = td.cloud().forces();

    if (td.cloud().breakup().solveOscillationEq())
    {
        solveTABEq(td, dt);
    }

    // cell state info is updated before
    const scalarField& Y = this->YLiquid();
    scalarField X(composition.liquids().X(Y));
    scalar rho = composition.liquids().rho(this->pc(), this->T(), X);
    scalar mu = composition.liquids().mu(this->pc(), this->T(), X);
    scalar sigma = composition.liquids().sigma(this->pc(), this->T(), X);

    // Average molecular weight of carrier mix - assumes perfect gas
    scalar Wc = this->rhoc()*specie::RR*this->Tc()/this->pc();
    scalar R = specie::RR/Wc;
    scalar Tav = td.cloud().atomization().Taverage(this->T(), this->Tc());
    // calculate average gas density based on average temperature
    scalar rhoAv = this->pc()/(R*Tav);
    scalar muAv = this->muc();
    vector Urel = this->U() - this->Uc();
    scalar Urmag = mag(Urel);
    scalar Re = this->Re(this->U(), this->d(), rhoAv, muAv);
    scalar We = this->We(this->U(), this->d_, this->rhoc_, sigma);

    const scalar mass = p.mass();
    const forceSuSp Fcp = forces.calcCoupled(p, dt, mass, Re, mu, this->y(), this->Y_[LIQ], this->Bt(), We);
    const forceSuSp Fncp = forces.calcNonCoupled(p, dt, mass, Re, muAv);
    this->tMom() = mass/(Fcp.Sp() + Fncp.Sp());
    const vector g = td.cloud().g().value();

    scalar massChild = 0.0;
    scalar dChild = 0.0;
    if
    (
        td.cloud().breakup().update
        (
            dt,
            g,
            this->d(),
            this->tc(),
            this->ms(),
            this->nParticle(),
            this->KHindex(),
            this->y(),
            this->yDot(),
            this->d0(),
            rho,
            mu,
            sigma,
            this->U(),
            rhoAv,
            muAv,
            Urel,
            Urmag,
            this->tMom(),
            dChild,
            massChild
        )
    )
    {
        scalar Re = rhoAv*Urmag*dChild/muAv;
        this->mass0() -= massChild; //xkè considerare mass0 e non massa attuale???
           // Add child parcel as copy of parent
        SprayParcel<ParcelType>* child = new SprayParcel<ParcelType>(*this);
        child->mass0() = massChild;
        child->d() = dChild;
        child->nParticle() = massChild/rho*this->volume(dChild);
        const forceSuSp Fcp = forces.calcCoupled(*child, dt, massChild, Re, muAv, td.cloud().breakup().y0(), this->Y_[LIQ], 2.0, We);
        const forceSuSp Fncp = forces.calcNonCoupled(*child, dt, massChild, Re, muAv);
        child->liquidCore() = 0.0;
        child->Renolds() = 0.0;
        child->Weber() = 0.0;
        child->om2() = 0.0;
        child->KHindex() = 1.0;
        child->y() = td.cloud().breakup().y0();
        child->yDot() = td.cloud().breakup().yDot0();
        child->tc() = -GREAT;
        child->ms() = 0.0;
        child->injector() = this->injector();
        child->tMom() = massChild/(Fcp.Sp() + Fncp.Sp());
        child->user() = 1.0;
        child->setCellValues(td, dt, cellI);

        td.cloud().addParticle(child);
    }
}


template<class ParcelType>
template<class TrackData>
Foam::scalar Foam::SprayParcel<ParcelType>::chi
(
    TrackData& td,
    const scalarField& X
) const
{
    // modifications to take account of the flash boiling on primary break-up
    typedef typename TrackData::cloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        td.cloud().composition();
    scalar chi = 0.0;
    scalar T0 = this->T();
    scalar p0 = this->pc();
    scalar pAmb = td.cloud().pAmbient();
    scalar pv = composition.liquids().pv(p0, T0, X);
    forAll(composition.liquids(), i)
    {
        if (pv >= 0.999*pAmb)
        {
            // liquid is boiling - calc boiling temperature
            const liquidProperties& liq = composition.liquids().properties()[i];
            scalar TBoil = liq.pvInvert(p0);
            scalar hl = liq.hl(pAmb, TBoil);
            scalar iTp = liq.h(pAmb, T0) - liq.rho(pAmb, T0);
            scalar iTb = liq.h(pAmb, TBoil) - pAmb/liq.rho(pAmb, TBoil);
            chi += X[i]*(iTp - iTb)/hl;
        }
    }
    chi = min(1.0, max(chi, 0.0));

    return chi;
}


template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::solveTABEq
(
    TrackData& td,
    const scalar dt
)
{
    typedef typename TrackData::cloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition = td.cloud().composition();

    const scalar& TABCmu = td.cloud().breakup().TABCmu();
    const scalar& TABWeCrit = td.cloud().breakup().TABWeCrit();
    const scalar& TABComega = td.cloud().breakup().TABComega();

    scalar r = 0.5*this->d();
    scalar r2 = r*r;
    scalar r3 = r*r2;

    const scalarField& Y = this->YLiquid();
    scalarField X(composition.liquids().X(Y));
    scalar rho = composition.liquids().rho(this->pc(), this->T(), X);
    scalar mu = composition.liquids().mu(this->pc(), this->T(), X);
    scalar sigma = composition.liquids().sigma(this->pc(), this->T(), X);
    scalar rtd = 0.5*TABCmu*mu/(rho*r2); // inverse of characteristic viscous damping time
    scalar omega2 = TABComega*sigma/(rho*r3) - rtd*rtd; // oscillation frequency (squared)
    this->om2() = omega2;

    if (omega2 > 0) // update distortion parameters
    {
        scalar omega = sqrt(omega2);
        scalar rhoc = this->rhoc();
        scalar Wetmp = this->We(this->U(), r, rhoc, sigma)/TABWeCrit;
        scalar y1 = this->y() - Wetmp;
        scalar y2 = this->yDot()/omega;
        scalar c = cos(omega*dt);
        scalar s = sin(omega*dt);
        scalar e = exp(-rtd*dt);
        y2 = (this->yDot() + y1*rtd)/omega;
        this->y() = Wetmp + e*(y1*c + y2*s);
        if (this->y() < 0)
        {
            this->y() = 0.0;
            this->yDot() = 0.0;
        }
        else
        {
            this->yDot() = (Wetmp - this->y())*rtd + e*omega*(y2*c - y1*s);
        }
    }
    else
    {
        // reset distortion parameters
        this->y() = 0;
        this->yDot() = 0;
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::newCalcDevolatilisation
(
    TrackData& td,
    const scalar dt,
    const scalar age,
    const scalar Ts,
    const scalar d,
    const scalar T,
    const scalar mass,
    const scalar mass0,
    const scalarField& YGasEff,
    const scalarField& YLiquidEff,
    const scalarField& YSolidEff,
    bool& canCombust,
    scalarField& dMassDV,
    scalar& Sh,
    scalar& N,
    scalar& NCpW,
    scalarField& Cs
) const
{
    // Check that model is active, and that the parcel temperature is within necessary limits for devolatilisation to occur
    if (!td.cloud().mydevolatilisation().active() || T < td.cloud().constProps().Tvap())
    {
        return;
    }

    typedef typename TrackData::cloudType::reactingCloudType reactingCloudType;
    const CompositionModel<reactingCloudType>& composition =
        td.cloud().composition();

    // Total mass of volatiles evolved
    td.cloud().mydevolatilisation().calculate
    (
        dt,
        age,
        mass0,
        mass,
        T,
        YGasEff,
        YLiquidEff,
        YSolidEff,
        canCombust,
        dMassDV
    );

    scalar dMassTot = sum(dMassDV);

    td.cloud().mydevolatilisation().addToDevolatilisationMass( this->nParticle_*dMassTot );

    Sh -= dMassTot*td.cloud().constProps().LDevol()/dt;

    // Update molar emissions
    if (td.cloud().heatTransfer().BirdCorrection())
    {
        // Molar average molecular weight of carrier mix
        const scalar Wc = max(SMALL, this->rhoc_*specie::RR*this->Tc_/this->pc_);

        // Note: hardcoded gaseous diffusivities for now
        // TODO: add to carrier thermo
        const scalar beta = sqr(cbrt(15.0) + cbrt(15.0));

        forAll(dMassDV, i)
        {
            const label id = composition.localToGlobalCarrierId(GAS, i);
            const scalar Cp = composition.carrier().Cp(id, this->pc_, Ts);
            const scalar W = composition.carrier().W(id);
            const scalar Ni = dMassDV[i]/(this->areaS(d)*dt*W);
            // Dab calc'd using API vapour mass diffusivity function
            const scalar Dab = 3.6059e-3*(pow(1.8*Ts, 1.75))*sqrt(1.0/W + 1.0/Wc)/(this->pc_*beta);
            N += Ni;
            NCpW += Ni*Cp*W;
            Cs[id] += Ni*d/(2.0*Dab);
        }
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::newCalcSurfaceReactions
(
    TrackData& td,
    const scalar dt,
    const label cellI,
    const scalar d,
    const scalar T,
    const scalar mass,
    const bool canCombust,
    const scalar N,
    const scalarField& YMix,
    const scalarField& YGas,
    const scalarField& YLiquid,
    const scalarField& YSolid,
    scalarField& dMassSRGas,
    scalarField& dMassSRLiquid,
    scalarField& dMassSRSolid,
    scalarField& dMassSRCarrier,
    scalar& Sh,
    scalar& dhsTrans
) const
{
    // Check that model is active
    if (!td.cloud().mysurfaceReaction().active()) // || !canCombust)
    {
        return;
    }

    // Calculate ABSOLUTE VALUES of the HEAT and MASS of surface reactions
    const scalar hReaction = td.cloud().mysurfaceReaction().calculate
    (
        dt,
        cellI,
        d,
        T,
        this->Tc_,
        this->pc_,
        this->rho(), //particle density!!!!
        mass,
        YGas,
        YLiquid,
        YSolid,
        YMix,
        N,
        dMassSRGas,
        dMassSRLiquid,
        dMassSRSolid,
        dMassSRCarrier
    );

    td.cloud().mysurfaceReaction().addToSurfaceReactionMass
    (
        this->nParticle_*(sum(dMassSRGas) + sum(dMassSRLiquid) + sum(dMassSRSolid))
    ); //!!! dMass limited in submodel (not as in PhaseChange)

    const scalar coeff = (1.0 - td.cloud().constProps().hRetentionCoeff());

    Sh += coeff*hReaction/dt;

    dhsTrans += (1.0 - coeff)*hReaction;

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::SprayParcel<ParcelType>::SprayParcel(const SprayParcel<ParcelType>& p)
:
    ParcelType(p),
    d0_(p.d0_),
    position0_(p.position0_),
    liquidCore_(p.liquidCore_),
    Renolds_(p.Renolds_),
    Weber_(p.Weber_),
    om2_(p.om2_),
    KHindex_(p.KHindex_),
    y_(p.y_),
    yDot_(p.yDot_),
    tc_(p.tc_),
    ms_(p.ms_),
    injector_(p.injector_),
    tMom_(p.tMom_),
    user_(p.user_),
    Bt_(p.Bt_),
    YGas_(p.YGas_),
    YLiquid_(p.YLiquid_),
    YSolid_(p.YSolid_),
    canCombust_(p.canCombust_)
{}


template<class ParcelType>
Foam::SprayParcel<ParcelType>::SprayParcel
(
    const SprayParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    d0_(p.d0_),
    position0_(p.position0_),
    liquidCore_(p.liquidCore_),
    Renolds_(p.Renolds_),
    Weber_(p.Weber_),
    om2_(p.om2_),
    KHindex_(p.KHindex_),
    y_(p.y_),
    yDot_(p.yDot_),
    tc_(p.tc_),
    ms_(p.ms_),
    injector_(p.injector_),
    tMom_(p.tMom_),
    user_(p.user_),
    Bt_(p.Bt_),
    YGas_(p.YGas_),
    YLiquid_(p.YLiquid_),
    YSolid_(p.YSolid_),
    canCombust_(p.canCombust_)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "SprayParcelIO.C"


// ************************************************************************* //
