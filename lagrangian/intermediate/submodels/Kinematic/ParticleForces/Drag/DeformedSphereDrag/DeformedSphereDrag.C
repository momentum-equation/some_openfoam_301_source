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

#include "DeformedSphereDrag.H"
//#include "mathematicalConstants.H"
//#include "dimensionedScalar.H"
//using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::DeformedSphereDrag<CloudType>::CdRe
 (const scalar Re,
  const scalar y,
  const scalar YPhase, 
  const scalar Bt,
  const scalar We) const
{

 scalar Cd = 0;
 scalar bit = Bt;
 if (bit==2) //i.e. no evporation
 {
     bit = 0;
 }

 if ((YPhase - 1.0) < SMALL) 
 {
    if (Re > 1000.0)
    {
        Cd = 0.424*Re; //Kiva
//        Cd = 0.44*Re; //classico
    }
    else if (Re < 20.0)
    {
        Cd = 24.0*(1.0 + 0.13*pow(Re, (0.82-log10(Re))));  //Clift
    }
    else if (20 < Re && Re < 300.0)
    {
//        Cd = 24.0*(1.0 + 1.0/6.0*pow(Re, 2.0/3.0));   //classico Yuen
	  Cd = 24.0/pow((1+bit), 0.2)*(1.0 + 0.2*pow(Re, 0.63))*(1+0.06*pow(Re, -0.12)*pow(We, 1.4));   //Renksizbulut
    }
    else 
    {
        Cd = 24.0*(1.0 + 1.0/6.0*pow(Re, 2.0/3.0));   //Yuen (better, ref.Kralj)
//        Cd = 24.0*(1.0 + 0.15*pow(Re, 0.687));   //Clift,Fritsching
    }
    return Cd;
//    return Cd*(1+2.632*min(1.0,y)); //Reitz,Sirignano (ALTERNATIVE to Renksizbulut!!!!!!!)
 }
 else 
 {
    if (Re > 1000.0)
    {
        Cd = 0.44*Re;
    }
    else 
    {
        Cd = 24.0*(1.0 + 1.0/6.0*pow(Re, 2.0/3.0));   //Yuen (better, ref.Kralj)
//        Cd = 24.0*(1.0 + 0.15*pow(Re, 0.687));   //Clift,Fritsching
    }
    return Cd;
 }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DeformedSphereDrag<CloudType>::DeformedSphereDrag
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, false)
{}


template<class CloudType>
Foam::DeformedSphereDrag<CloudType>::DeformedSphereDrag
(
    const DeformedSphereDrag<CloudType>& df
)
:
    ParticleForce<CloudType>(df)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::DeformedSphereDrag<CloudType>::~DeformedSphereDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::DeformedSphereDrag<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc,
    const scalar y, 
    const scalar YPhase, 
    const scalar Bt, 
    const scalar We
) const
{   

    forceSuSp value(vector::zero, 0.0);

    value.Sp() = mass*0.75*muc*CdRe(Re, y, YPhase, Bt, We)/(p.rho()*sqr(p.d())); 
    return value;
}


// ************************************************************************* //
