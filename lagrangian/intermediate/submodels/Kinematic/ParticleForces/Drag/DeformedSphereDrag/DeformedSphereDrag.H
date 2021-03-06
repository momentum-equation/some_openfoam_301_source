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
    Foam::DeformedSphereDrag

Description
    Drag model based on assumption of solid spheres

\*---------------------------------------------------------------------------*/

#ifndef DeformedSphereDrag_H
#define DeformedSphereDrag_H

#include "ParticleForce.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class DeformedSphereDrag Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DeformedSphereDrag
:
    public ParticleForce<CloudType>
{
    // Private Member Functions

        //- Drag coefficient multiplied by Reynolds number
//        scalar CdRe(const scalar Re, const scalar y, const scalar YPhase, const scalar Bt, const scalar We) const; //temporaneo non uso
        scalar CdRe(const scalar Re, const scalar y, const scalar YPhase, const scalar Bt, const scalar We) const;

public:

    //- Runtime type information
    TypeName("deformedSphereDrag");


    // Constructors

        //- Construct from mesh
        DeformedSphereDrag
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        DeformedSphereDrag(const DeformedSphereDrag<CloudType>& df);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType> > clone() const
        {
            return autoPtr<ParticleForce<CloudType> >
            (
                new DeformedSphereDrag<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~DeformedSphereDrag();


    // Member Functions

        // Evaluation

            //- Calculate the coupled force
            virtual forceSuSp calcCoupled
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
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DeformedSphereDrag.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
