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
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::string Foam::SprayParcel<ParcelType>::propHeader =
    ParcelType::propHeader
  + " d0"
  + " position0"
  + " liquidCore"
  + " Renolds"
  + " Weber"
  + " om2"
  + " KHindex"
  + " y"
  + " yDot"
  + " tc"
  + " ms"
  + " injector"
  + " tMom"
  + " user"
  + " Bt"
  + " YGas"
  + " YLiquid"
  + " YSolid"
  + " canCombust";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::SprayParcel<ParcelType>::SprayParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    ParcelType(mesh, is, readFields),
    d0_(0.0),
    position0_(vector::zero),
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
    Bt_(0.0),
    YGas_(0),
    YLiquid_(0),
    YSolid_(0),
    canCombust_(false)
{
    if (readFields)
    {
        DynamicList<scalar> Yg;
        DynamicList<scalar> Yl;
        DynamicList<scalar> Ys;

        is >> Yg >> Yl >> Ys;

        YGas_.transfer(Yg);
        YLiquid_.transfer(Yl);
        YSolid_.transfer(Ys);

        // scale the mass fractions
        const scalarField& YMix = this->Y_;
        YGas_ /= YMix[GAS] + ROOTVSMALL;
        YLiquid_ /= YMix[LIQ] + ROOTVSMALL;
        YSolid_ /= YMix[SLD] + ROOTVSMALL;

        if (is.format() == IOstream::ASCII)
        {
            d0_ = readScalar(is);
            is >> position0_;
            liquidCore_ = readScalar(is);
            Renolds_ = readScalar(is);
            Weber_ = readScalar(is);
            om2_ = readScalar(is);
            KHindex_ = readScalar(is);
            y_ = readScalar(is);
            yDot_ = readScalar(is);
            tc_ = readScalar(is);
            ms_ = readScalar(is);
            injector_ = readScalar(is);
            tMom_ = readScalar(is);
            user_ = readScalar(is);
            Bt_ = readScalar(is);
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&d0_),
                sizeof(d0_)
              + sizeof(position0_)
              + sizeof(liquidCore_)
              + sizeof(Renolds_)
              + sizeof(Weber_)
              + sizeof(om2_)
              + sizeof(KHindex_)
              + sizeof(y_)
              + sizeof(yDot_)
              + sizeof(tc_)
              + sizeof(ms_)
              + sizeof(injector_)
              + sizeof(tMom_)
              + sizeof(user_)
              + sizeof(Bt_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "SprayParcel<ParcelType>::SprayParcel"
        "("
            "const polyMesh, "
            "Istream&, "
            "bool"
        ")"
    );
}


template<class ParcelType>
template<class CloudType>
void Foam::SprayParcel<ParcelType>::readFields(CloudType& c)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::SprayParcel<ParcelType>::readFields
(
    CloudType& c,
    const CompositionType& compModel
)
{
    if (!c.size())
    {
        return;
    }

    ParcelType::readFields(c, compModel);

    IOField<scalar> d0(c.fieldIOobject("d0", IOobject::MUST_READ));
    c.checkFieldIOobject(c, d0);

    IOField<vector> position0
    (
        c.fieldIOobject("position0", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, position0);

    IOField<scalar> liquidCore(c.fieldIOobject
    (
        "liquidCore", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, liquidCore);

    IOField<scalar> Renolds(c.fieldIOobject("Renolds", IOobject::MUST_READ));
    c.checkFieldIOobject(c, Renolds);

    IOField<scalar> Weber(c.fieldIOobject("Weber", IOobject::MUST_READ));
    c.checkFieldIOobject(c, Weber);

    IOField<scalar> om2(c.fieldIOobject("om2", IOobject::MUST_READ));
    c.checkFieldIOobject(c, om2);

    IOField<scalar> KHindex(c.fieldIOobject("KHindex", IOobject::MUST_READ));
    c.checkFieldIOobject(c, KHindex);

    IOField<scalar> y(c.fieldIOobject("y", IOobject::MUST_READ));
    c.checkFieldIOobject(c, y);

    IOField<scalar> yDot(c.fieldIOobject("yDot", IOobject::MUST_READ));
    c.checkFieldIOobject(c, yDot);

    IOField<scalar> tc(c.fieldIOobject("tc", IOobject::MUST_READ));
    c.checkFieldIOobject(c, tc);

    IOField<scalar> ms(c.fieldIOobject("ms", IOobject::MUST_READ));
    c.checkFieldIOobject(c, ms);

    IOField<scalar> injector(c.fieldIOobject("injector", IOobject::MUST_READ));
    c.checkFieldIOobject(c, injector);

    IOField<scalar> tMom(c.fieldIOobject("tMom", IOobject::MUST_READ));
    c.checkFieldIOobject(c, tMom);

    IOField<scalar> user(c.fieldIOobject("user", IOobject::MUST_READ));
    c.checkFieldIOobject(c, user);

    IOField<scalar> Bt(c.fieldIOobject("Bt", IOobject::MUST_READ));
    c.checkFieldIOobject(c, Bt);

    // Get names and sizes for each Y...
    const label idGas = compModel.idGas();
    const wordList& gasNames = compModel.componentNames(idGas);
    const label idLiquid = compModel.idLiquid();
    const wordList& liquidNames = compModel.componentNames(idLiquid);
    const label idSolid = compModel.idSolid();
    const wordList& solidNames = compModel.componentNames(idSolid);
    const wordList& stateLabels = compModel.stateLabels();
    // Set storage for each Y... for each parcel
    forAllIter(typename Cloud<SprayParcel<ParcelType> >, c, iter)
    {
        SprayParcel<ParcelType>& p = iter();
        p.YGas_.setSize(gasNames.size(), 0.0);
        p.YLiquid_.setSize(liquidNames.size(), 0.0);
        p.YSolid_.setSize(solidNames.size(), 0.0);
    }
    // Populate YGas for each parcel
    forAll(gasNames, j)
    {
        IOField<scalar> YGas
        (
            c.fieldIOobject
            (
                "Y" + gasNames[j] + stateLabels[idGas],
                IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIter
        (
            typename Cloud<SprayParcel<ParcelType> >,
            c,
            iter
        )
        {
            SprayParcel<ParcelType>& p = iter();
            p.YGas_[j] = YGas[i++]/(p.Y()[GAS] + ROOTVSMALL);
        }
    }
    // Populate YLiquid for each parcel
    forAll(liquidNames, j)
    {
        IOField<scalar> YLiquid
        (
            c.fieldIOobject
            (
                "Y" + liquidNames[j] + stateLabels[idLiquid],
                 IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIter
        (
            typename Cloud<SprayParcel<ParcelType> >,
            c,
            iter
        )
        {
            SprayParcel<ParcelType>& p = iter();
            p.YLiquid_[j] = YLiquid[i++]/(p.Y()[LIQ] + ROOTVSMALL);
        }
    }
    // Populate YSolid for each parcel
    forAll(solidNames, j)
    {
        IOField<scalar> YSolid
        (
            c.fieldIOobject
            (
                "Y" + solidNames[j] + stateLabels[idSolid],
                IOobject::MUST_READ
            )
        );

        label i = 0;
        forAllIter
        (
            typename Cloud<SprayParcel<ParcelType> >,
            c,
            iter
        )
        {
            SprayParcel<ParcelType>& p = iter();
            p.YSolid_[j] = YSolid[i++]/(p.Y()[SLD] + ROOTVSMALL);
        }
    }


    label i = 0;
    forAllIter(typename Cloud<SprayParcel<ParcelType> >, c, iter)
    {
        SprayParcel<ParcelType>& p = iter();
        p.d0_ = d0[i];
        p.position0_ = position0[i];
        p.liquidCore_ = liquidCore[i];
        p.Renolds_ = Renolds[i];
        p.Weber_ = Weber[i];
        p.om2_ = om2[i];
        p.KHindex_ = KHindex[i];
        p.y_ = y[i];
        p.yDot_ = yDot[i];
        p.tc_ = tc[i];
        p.ms_ = ms[i];
        p.injector_ = injector[i];
        p.tMom_ = tMom[i];
        p.user_ = user[i];
        p.Bt_ = Bt[i];
        i++;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::SprayParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);
}


template<class ParcelType>
template<class CloudType, class CompositionType>
void Foam::SprayParcel<ParcelType>::writeFields
(
    const CloudType& c,
    const CompositionType& compModel
)
{
    ParcelType::writeFields(c, compModel);

    label np = c.size();


    // Write the composition fractions
    if (np > 0)
    {
        const wordList& stateLabels = compModel.stateLabels();

        const label idGas = compModel.idGas();
        const wordList& gasNames = compModel.componentNames(idGas);
        forAll(gasNames, j)
        {
            IOField<scalar> YGas
            (
                c.fieldIOobject
                (
                    "Y" + gasNames[j] + stateLabels[idGas],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            forAllConstIter
            (
                typename Cloud<SprayParcel<ParcelType> >,
                c,
                iter
            )
            {
                const SprayParcel<ParcelType>& p0 = iter();
                YGas[i++] = p0.YGas()[j]*p0.Y()[GAS];
            }

            YGas.write();
        }

        const label idLiquid = compModel.idLiquid();
        const wordList& liquidNames = compModel.componentNames(idLiquid);
        forAll(liquidNames, j)
        {
            IOField<scalar> YLiquid
            (
                c.fieldIOobject
                (
                    "Y" + liquidNames[j] + stateLabels[idLiquid],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            forAllConstIter
            (
                typename Cloud<SprayParcel<ParcelType> >,
                c,
                iter
            )
            {
                const SprayParcel<ParcelType>& p0 = iter();
                YLiquid[i++] = p0.YLiquid()[j]*p0.Y()[LIQ];
            }

            YLiquid.write();
        }

        const label idSolid = compModel.idSolid();
        const wordList& solidNames = compModel.componentNames(idSolid);
        forAll(solidNames, j)
        {
            IOField<scalar> YSolid
            (
                c.fieldIOobject
                (
                    "Y" + solidNames[j] + stateLabels[idSolid],
                    IOobject::NO_READ
                ),
                np
            );

            label i = 0;
            forAllConstIter
            (
                typename Cloud<SprayParcel<ParcelType> >,
                c,
                iter
            )
            {
                const SprayParcel<ParcelType>& p0 = iter();
                YSolid[i++] = p0.YSolid()[j]*p0.Y()[SLD];
            }

            YSolid.write();
        }
    }


    IOField<scalar> d0(c.fieldIOobject("d0", IOobject::NO_READ), np);
    IOField<vector> position0
    (
        c.fieldIOobject("position0", IOobject::NO_READ),
        np
    );
    IOField<scalar> liquidCore
    (
        c.fieldIOobject("liquidCore", IOobject::NO_READ),
        np
    );
    IOField<scalar> Renolds(c.fieldIOobject("Renolds", IOobject::NO_READ), np);
    IOField<scalar> Weber(c.fieldIOobject("Weber", IOobject::NO_READ), np);
    IOField<scalar> om2(c.fieldIOobject("om2", IOobject::NO_READ), np);
    IOField<scalar> KHindex(c.fieldIOobject("KHindex", IOobject::NO_READ), np);
    IOField<scalar> y(c.fieldIOobject("y", IOobject::NO_READ), np);
    IOField<scalar> yDot(c.fieldIOobject("yDot", IOobject::NO_READ), np);
    IOField<scalar> tc(c.fieldIOobject("tc", IOobject::NO_READ), np);
    IOField<scalar> ms(c.fieldIOobject("ms", IOobject::NO_READ), np);
    IOField<scalar> injector
    (
        c.fieldIOobject("injector", IOobject::NO_READ),
        np
    );
    IOField<scalar> tMom(c.fieldIOobject("tMom", IOobject::NO_READ), np);
    IOField<scalar> user(c.fieldIOobject("user", IOobject::NO_READ), np);
    IOField<scalar> Bt(c.fieldIOobject("Bt", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(typename Cloud<SprayParcel<ParcelType> >, c, iter)
    {
        const SprayParcel<ParcelType>& p = iter();
        d0[i] = p.d0_;
        position0[i] = p.position0_;
        liquidCore[i] = p.liquidCore_;
        Renolds[i] = p.Renolds_;
        Weber[i] = p.Weber_;
        om2[i] = p.om2_;
        KHindex[i] = p.KHindex_;
        y[i] = p.y_;
        yDot[i] = p.yDot_;
        tc[i] = p.tc_;
        ms[i] = p.ms_;
        injector[i] = p.injector_;
        tMom[i] = p.tMom_;
        user[i] = p.user_;
        Bt[i] = p.Bt_;
        i++;
    }

    d0.write();
    position0.write();
    liquidCore.write();
    Renolds.write();
    Weber.write();
    om2.write();
    KHindex.write();
    y.write();
    yDot.write();
    tc.write();
    ms.write();
    injector.write();
    tMom.write();
    user.write();
    Bt.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const SprayParcel<ParcelType>& p
)
{
    scalarField YGasLoc(p.YGas()*p.Y()[0]);
    scalarField YLiquidLoc(p.YLiquid()*p.Y()[1]);
    scalarField YSolidLoc(p.YSolid()*p.Y()[2]);

    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
        << token::SPACE << p.d0()
        << token::SPACE << p.position0()
        << token::SPACE << p.liquidCore()
        << token::SPACE << p.Renolds()
        << token::SPACE << p.Weber()
        << token::SPACE << p.om2()
        << token::SPACE << p.KHindex()
        << token::SPACE << p.y()
        << token::SPACE << p.yDot()
        << token::SPACE << p.tc()
        << token::SPACE << p.ms()
        << token::SPACE << p.injector()
        << token::SPACE << p.tMom()
        << token::SPACE << p.user()
        << token::SPACE << p.Bt()
        << token::SPACE << YGasLoc
        << token::SPACE << YLiquidLoc
        << token::SPACE << YSolidLoc;
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os  << YGasLoc << YLiquidLoc << YSolidLoc;
        os.write
        (
            reinterpret_cast<const char*>(&p.d0_),
            sizeof(p.d0())
          + sizeof(p.position0())
          + sizeof(p.liquidCore())
          + sizeof(p.Renolds())
          + sizeof(p.Weber())
          + sizeof(p.om2())
          + sizeof(p.KHindex())
          + sizeof(p.y())
          + sizeof(p.yDot())
          + sizeof(p.tc())
          + sizeof(p.ms())
          + sizeof(p.injector())
          + sizeof(p.tMom())
          + sizeof(p.user())
          + sizeof(p.Bt())
        );
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const SprayParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
