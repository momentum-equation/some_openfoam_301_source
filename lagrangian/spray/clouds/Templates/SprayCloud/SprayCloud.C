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

#include "SprayCloud.H"
#include "AtomizationModel.H"
#include "BreakupModel.H"
#include "StochasticCollisionModel.H"
#include "MyDevolatilisationModel.H"
#include "MySurfaceReactionModel.H"
#include <iostream>

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::SprayCloud<CloudType>::setModels()
{
    atomizationModel_.reset
    (
        AtomizationModel<SprayCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    breakupModel_.reset
    (
        BreakupModel<SprayCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    stochasticCollisionModel_.reset
    (
        StochasticCollisionModel<SprayCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    mydevolatilisationModel_.reset
    (
        MyDevolatilisationModel<SprayCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    mysurfaceReactionModel_.reset
    (
        MySurfaceReactionModel<SprayCloud<CloudType> >::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );
}


template<class CloudType>
void Foam::SprayCloud<CloudType>::cloudReset
(
    SprayCloud<CloudType>& c
)
{
    CloudType::cloudReset(c);

    atomizationModel_.reset(c.atomizationModel_.ptr());
    breakupModel_.reset(c.breakupModel_.ptr());
    stochasticCollisionModel_.reset(c.stochasticCollisionModel_.ptr());

    mydevolatilisationModel_.reset(c.mydevolatilisationModel_.ptr());
    mysurfaceReactionModel_.reset(c.mysurfaceReactionModel_.ptr());

    dMassDevolatilisation_ = c.dMassDevolatilisation_;
    dMassSurfaceReaction_ = c.dMassSurfaceReaction_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SprayCloud<CloudType>::SprayCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields
)
:
    CloudType(cloudName, rho, U, g, thermo, false),
    sprayCloud(),
    cloudCopyPtr_(NULL),
    averageParcelMass_(0.0),
    constProps_(this->particleProperties(), this->solution().active()),
    atomizationModel_(NULL),
    breakupModel_(NULL),
    stochasticCollisionModel_(NULL),
    mydevolatilisationModel_(NULL),
    mysurfaceReactionModel_(NULL),
    dMassDevolatilisation_(0.0),
    dMassSurfaceReaction_(0.0)
{
    if (this->solution().active())
    {
        setModels();

        averageParcelMass_ = this->injectors().averageParcelMass();

        if (readFields)
        {
            parcelType::readFields(*this, this->composition());
        }

        Info << "Average parcel mass: " << averageParcelMass_ << endl;
    }

    if (this->solution().resetSourcesOnStartup())
    {
        CloudType::resetSourceTerms();
    }
}


template<class CloudType>
Foam::SprayCloud<CloudType>::SprayCloud
(
    SprayCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    sprayCloud(),
    cloudCopyPtr_(NULL),
    averageParcelMass_(c.averageParcelMass_),
    constProps_(c.constProps_),
    atomizationModel_(c.atomizationModel_->clone()),
    breakupModel_(c.breakupModel_->clone()),
    stochasticCollisionModel_(c.stochasticCollisionModel_->clone()),
    mydevolatilisationModel_(c.mydevolatilisationModel_->clone()),
    mysurfaceReactionModel_(c.mysurfaceReactionModel_->clone()),
    dMassDevolatilisation_(c.dMassDevolatilisation_),
    dMassSurfaceReaction_(c.dMassSurfaceReaction_)
{}


template<class CloudType>
Foam::SprayCloud<CloudType>::SprayCloud
(
    const fvMesh& mesh,
    const word& name,
    const SprayCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    sprayCloud(),
    cloudCopyPtr_(NULL),
    averageParcelMass_(0.0),
    constProps_(),
    atomizationModel_(NULL),
    breakupModel_(NULL),
    stochasticCollisionModel_(NULL),
    mydevolatilisationModel_(NULL),
    mysurfaceReactionModel_(NULL),
    dMassDevolatilisation_(0.0),
    dMassSurfaceReaction_(0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SprayCloud<CloudType>::~SprayCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::SprayCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);

    label idGas = this->composition().idGas();
    label idLiquid = this->composition().idLiquid();
    label idSolid = this->composition().idSolid();

    parcel.YGas() = this->composition().Y0(idGas);
    parcel.YLiquid() = this->composition().Y0(idLiquid);
    parcel.YSolid() = this->composition().Y0(idSolid);

    if (mag(sum(parcel.YLiquid()) - 1.0) < SMALL)
    {
        const liquidMixtureProperties& liqMix = this->composition().liquids();
	scalarField X(liqMix.X(parcel.YLiquid()));
	parcel.Cp() = liqMix.Cp(parcel.pc(), parcel.T(), X); //override rho and Cp from constantProperties
	parcel.rho() = liqMix.rho(parcel.pc(), parcel.T(), X);
    }
}


template<class CloudType>
void Foam::SprayCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);

    // store the injection position and initial drop size
    parcel.position0() = parcel.position();
    parcel.d0() = parcel.d();

    parcel.y() = breakup().y0();
    parcel.yDot() = breakup().yDot0();

    parcel.liquidCore() = atomization().initLiquidCore();

    if (fullyDescribed)
    {
        label idGas = this->composition().idGas();
        label idLiquid = this->composition().idLiquid();
        label idSolid = this->composition().idSolid();

        this->checkSuppliedComposition
        (
            parcel.YGas(),
            this->composition().Y0(idGas),
            "YGas"
        );
        this->checkSuppliedComposition
        (
            parcel.YLiquid(),
            this->composition().Y0(idLiquid),
            "YLiquid"
        );
        this->checkSuppliedComposition
        (
            parcel.YSolid(),
            this->composition().Y0(idSolid),
            "YSolid"
        );
    }
}


template<class CloudType>
void Foam::SprayCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<SprayCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::SprayCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::SprayCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::template
            TrackingData<SprayCloud<CloudType> > td(*this);

        this->solve(td);
    }
}


template<class CloudType>
template<class TrackData>
void Foam::SprayCloud<CloudType>::motion(TrackData& td)
{
    const scalar dt = this->solution().trackTime();

    td.part() = TrackData::tpLinearTrack;
    CloudType::move(td, dt);

    this->updateCellOccupancy();

    scalar collisionCase = 0.0;

    //questa parte andrebbe fatta solo se ho goccia liquida, nel caso solido cambierebbe
    if (mag(sum(this->composition().Y0(this->composition().idLiquid())) - 1.0) < SMALL)
    {
    if (stochasticCollision().active())
    {
        const liquidMixtureProperties& liqMix = this->composition().liquids();

        label i = 0;
        forAllIter(typename SprayCloud<CloudType>, *this, iter)
        {
            label j = 0;
            forAllIter(typename SprayCloud<CloudType>, *this, jter)
            {
                if (j > i)
                {
                    parcelType& p = iter();
                    scalar Vi = this->mesh().V()[p.cell()];
                    scalarField X1(liqMix.X(p.Y()));
                    scalar sigma1 = liqMix.sigma(p.pc(), p.T(), X1);
                    scalar mp = p.mass()*p.nParticle();

                    parcelType& q = jter();
                    scalar Vj = this->mesh().V()[q.cell()];
                    scalarField X2(liqMix.X(q.Y()));
                    scalar sigma2 = liqMix.sigma(q.pc(), q.T(), X2);
                    scalar mq = q.mass()*q.nParticle();

                    bool updateProperties = stochasticCollision().update
                    (
                        dt,
                        this->rndGen(),
                        p.position(),
                        mp,
                        p.d(),
                        p.nParticle(),
                        p.U(),
                        p.rho(),
                        p.T(),
                        p.Y(),
                        sigma1,
                        p.cell(),
                        Vi,
                        q.position(),
                        mq,
                        q.d(),
                        q.nParticle(),
                        q.U(),
                        q.rho(),
                        q.T(),
                        q.Y(),
                        sigma2,
                        q.cell(),
                        Vj
                    );

                    // for coalescence we need to update the density and the diameter cause of the temp/conc/mass-change
                    if (updateProperties)
                    {
			collisionCase += 1;                        
			if (mp > VSMALL)
                        {
                            scalarField Xp(liqMix.X(p.Y()));
                            p.rho() = liqMix.rho(p.pc(), p.T(), Xp);
                            p.Cp() = liqMix.Cp(p.pc(), p.T(), Xp);
                            p.d() =
                                cbrt
                                (
                                    6.0*mp
                                   /(
                                        p.nParticle()
                                       *p.rho()
                                       *constant::mathematical::pi
                                    )
                                );
                        }

                        if (mq > VSMALL)
                        {
                            scalarField Xq(liqMix.X(q.Y()));
                            q.rho() = liqMix.rho(q.pc(), q.T(), Xq);
                            q.Cp() = liqMix.Cp(q.pc(), q.T(), Xq);
                            q.d() =
                                cbrt
                                (
                                    6.0*mq
                                   /(
                                        q.nParticle()
                                       *q.rho()
                                       *constant::mathematical::pi
                                    )
                                );
                        }
                    }
                }
                j++;
            }
            i++;
        }

        // remove coalesced parcels that fall below minimum mass threshold
        forAllIter(typename SprayCloud<CloudType>, *this, iter)
        {
            parcelType& p = iter();
            scalar mass = p.nParticle()*p.mass();

            if (mass < td.cloud().constProps().minParticleMass())
            {
                this->deleteParticle(p);
            }
        }
    }
    Info << "    Coalescence cases = " << collisionCase << nl << endl;
    }
}


template<class CloudType>
void Foam::SprayCloud<CloudType>::info()
{
    CloudType::info();
    scalar d32 = 1.0e+6*this->Dij(3, 2);
    scalar d10 = 1.0e+6*this->Dij(1, 0);
    scalar dMax = 1.0e+6*this->Dmax();
    scalar pen = this->penetration(0.95);

    Info << "    D10, D32, Dmax (mu)             = " << d10 << ", " << d32
         << ", " << dMax << nl
         << "    Liquid penetration 95% mass (m) = " << pen << nl;

    this->mydevolatilisation().info(Info);
    this->mysurfaceReaction().info(Info);
}


// ************************************************************************* //
