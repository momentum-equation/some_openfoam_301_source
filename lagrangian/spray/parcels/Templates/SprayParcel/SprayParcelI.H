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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::SprayParcel<ParcelType>::constantProperties::constantProperties()
:
    ParcelType::constantProperties(),
    LDevol_(0.0),
    hRetentionCoeff_(0.0)
{}


template<class ParcelType>
inline Foam::SprayParcel<ParcelType>::constantProperties::constantProperties
(
    const constantProperties& cp
)
:
    ParcelType::constantProperties(cp),
    LDevol_(cp.LDevol_),
    hRetentionCoeff_(cp.hRetentionCoeff_)
{}


template<class ParcelType>
inline Foam::SprayParcel<ParcelType>::constantProperties::constantProperties
(
    const dictionary& parentDict,
    const bool readFields
)
:
    ParcelType::constantProperties(parentDict, readFields),
    LDevol_(0.0),
    hRetentionCoeff_(0.0)
{
    if (readFields)
    {
        this->dict().lookup("LDevol") >> LDevol_;
        this->dict().lookup("hRetentionCoeff") >> hRetentionCoeff_;

        if ((hRetentionCoeff_ < 0) || (hRetentionCoeff_ > 1))
        {
            FatalErrorIn
            (
                "Spray<ParcelType>::constantProperties::"
                "constantProperties"
            )   << "hRetentionCoeff must be in the range 0 to 1" << nl
                << exit(FatalError) << endl;
        }

        hRetentionCoeff_ = max(1e-06, hRetentionCoeff_);
    }
}


template<class ParcelType>
inline Foam::SprayParcel<ParcelType>::SprayParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const label tetFaceI,
    const label tetPtI
)
:
    ParcelType(mesh, position, cellI, tetFaceI, tetPtI),
    d0_(this->d()),
    position0_(position),
    liquidCore_(0.0),
    Renolds_(0.0),
    Weber_(0.0),
    om2_(0.0),
    KHindex_(0.0),
    y_(0.0),
    yDot_(0.0),
    tc_(0.0),
    ms_(0.0),
    injector_(1.0),
    tMom_(GREAT),
    user_(0.0),
    Bt_(2.0),
    YGas_(0),
    YLiquid_(0),
    YSolid_(0),
    canCombust_(false)
{}


template<class ParcelType>
inline Foam::SprayParcel<ParcelType>::SprayParcel
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
)
:
    ParcelType
    (
        mesh,
        position,
        cellI,
        tetFaceI,
        tetPtI,
        typeId,
        nParticle0,
        d0,
        dTarget0,
        U0,
        f0,
        angularMomentum0,
        torque0,
        Y0,
        constProps
    ),
    d0_(d0),
    position0_(position),
    liquidCore_(liquidCore),
    Renolds_(Renolds),
    Weber_(Weber),
    om2_(om2),
    KHindex_(KHindex),
    y_(y),
    yDot_(yDot),
    tc_(tc),
    ms_(ms),
    injector_(injector),
    tMom_(tMom),
    user_(user),
    Bt_(Bt),
    YGas_(YGas0),
    YLiquid_(YLiquid0),
    YSolid_(YSolid0),
    canCombust_(false)
{}


// * * * * * * * * * constantProperties Member Functions * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::constantProperties::LDevol() const
{
    return LDevol_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::constantProperties::hRetentionCoeff() const
{
    return hRetentionCoeff_;
}


// * * * * * * * * * * SprayParcel Member Functions * * * * * * * * * * * * //

template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::d0() const
{
    return d0_;
}


template<class ParcelType>
inline const Foam::vector& Foam::SprayParcel<ParcelType>::position0() const
{
    return position0_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::liquidCore() const
{
    return liquidCore_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::Renolds() const
{
    return Renolds_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::Weber() const
{
    return Weber_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::om2() const
{
    return om2_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::KHindex() const
{
    return KHindex_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::y() const
{
    return y_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::yDot() const
{
    return yDot_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::tc() const
{
    return tc_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::ms() const
{
    return ms_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::injector() const
{
    return injector_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::tMom() const
{
    return tMom_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::user() const
{
    return user_;
}


template<class ParcelType>
inline Foam::scalar Foam::SprayParcel<ParcelType>::Bt() const
{
    return Bt_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::SprayParcel<ParcelType>::YGas() const
{
    return YGas_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::SprayParcel<ParcelType>::YLiquid() const
{
    return YLiquid_;
}


template<class ParcelType>
inline const Foam::scalarField& Foam::SprayParcel<ParcelType>::YSolid() const
{
    return YSolid_;
}


template<class ParcelType>
inline bool Foam::SprayParcel<ParcelType>::canCombust() const
{
    return canCombust_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::d0()
{
    return d0_;
}


template<class ParcelType>
inline Foam::vector& Foam::SprayParcel<ParcelType>::position0()
{
    return position0_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::liquidCore()
{
    return liquidCore_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::Renolds()
{
    return Renolds_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::Weber()
{
    return Weber_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::om2()
{
    return om2_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::KHindex()
{
    return KHindex_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::y()
{
    return y_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::yDot()
{
    return yDot_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::tc()
{
    return tc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::ms()
{
    return ms_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::injector()
{
    return injector_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::tMom()
{
    return tMom_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::user()
{
    return user_;
}


template<class ParcelType>
inline Foam::scalar& Foam::SprayParcel<ParcelType>::Bt()
{
    return Bt_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::SprayParcel<ParcelType>::YGas()
{
    return YGas_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::SprayParcel<ParcelType>::YLiquid()
{
    return YLiquid_;
}


template<class ParcelType>
inline Foam::scalarField& Foam::SprayParcel<ParcelType>::YSolid()
{
    return YSolid_;
}


template<class ParcelType>
inline bool& Foam::SprayParcel<ParcelType>::canCombust()
{
    return canCombust_;
}


// ************************************************************************* //
