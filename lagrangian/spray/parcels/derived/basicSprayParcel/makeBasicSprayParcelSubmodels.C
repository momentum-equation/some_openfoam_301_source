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

#include "basicSprayCloud.H"

#include "makeParcelCloudFunctionObjects.H"

// Kinematic
#include "makeThermoParcelForces.H" // thermo variant
#include "makeParcelDispersionModels.H"
#include "makeSprayParcelInjectionModels.H" // Spray variant
#include "makeParcelPatchInteractionModels.H"

// Thermodynamic
#include "makeParcelHeatTransferModels.H"

// Reacting
#include "makeReactingParcelPhaseChangeModels.H"
#include "makeReactingParcelSurfaceFilmModels.H"

// Spray
#include "makeSprayParcelAtomizationModels.H"
#include "makeSprayParcelBreakupModels.H"
#include "makeSprayParcelCollisionModels.H"

// New Spray (from ReactingMultiphase)
#include "makeSprayParcelCompositionModels.H"
#include "makeSprayParcelDevolatilisationModels.H"
#include "makeSprayParcelSurfaceReactionModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeParcelCloudFunctionObjects(basicSprayCloud);

    // Kinematic sub-models
    makeThermoParcelForces(basicSprayCloud);
    makeParcelDispersionModels(basicSprayCloud);
    makeSprayParcelInjectionModels(basicSprayCloud);
    makeParcelPatchInteractionModels(basicSprayCloud);

    // Thermo sub-models
    makeParcelHeatTransferModels(basicSprayCloud);

    // Reacting sub-models
    makeReactingParcelPhaseChangeModels(basicSprayCloud);
    makeReactingParcelSurfaceFilmModels(basicSprayCloud);

    // Spray sub-models
    makeSprayParcelAtomizationModels(basicSprayCloud);
    makeSprayParcelBreakupModels(basicSprayCloud);
    makeSprayParcelCollisionModels(basicSprayCloud);

    // New Spray sub-models (from ReactingMultiphase)
    makeSprayParcelCompositionModels(basicSprayCloud);
    makeSprayParcelDevolatilisationModels(basicSprayCloud);
    makeSprayParcelSurfaceReactionModels(basicSprayCloud);
};


// ************************************************************************* //
