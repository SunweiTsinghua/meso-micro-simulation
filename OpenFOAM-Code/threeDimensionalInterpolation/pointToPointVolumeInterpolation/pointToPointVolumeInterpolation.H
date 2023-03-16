/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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
    Foam::pointToPointVolumeInterpolation

Description
    Interpolates between two sets of unstructured points through finding the
    nearest bounding pyramid.

SourceFiles
    pointToPointVolumeInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef pointToPointVolumeInterpolation_H
#define pointToPointVolumeInterpolation_H

#include "FixedList.H"
#include "instantList.H"
#include "vectorField.H"
#include "pointField.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class pointToPointPlanarInterpolation Declaration
\*---------------------------------------------------------------------------*/

class pointToPointVolumeInterpolation
{
    // Private data

        //- Number of source points (for checking)
        label nPoints_;

        //- Current interpolation addressing to the pyramid vertices
        List<FixedList<label, 4> > nearestVertex_;

        //- Current interpolation factors to the pyramid vertices
        List<FixedList<scalar, 4> > nearestVertexWeight_;

    // Private Member Functions

		//- Test if a certain point is in the tetrahedron
		bool inTet
		(
			const FixedList<vector,4>& tetPoints,
			const vector& inPoint
		);
		
        //- Calculate addressing and weights
        void calcWeights
        (
            const pointField& sourcePoints,
            const pointField& destPoints
        );

public:

    // Declare name of the class and its debug switch
    ClassName("pointToPointVolumeInterpolation");


    // Constructors

        //- Construct from 3D locations
        
        pointToPointVolumeInterpolation
        (
            const pointField& sourcePoints,
            const pointField& destPoints
        );


    // Member Functions

        //- Interpolate from field on source points to dest points
        tmp<vectorField> interpolate(const vectorField& sourceFld) const;
        tmp<scalarField> interpolate(const scalarField& sourceFld) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
