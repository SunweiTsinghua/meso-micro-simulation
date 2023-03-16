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

#ifndef trilinearInterpolation_H
#define trilinearInterpolation_H

#include "FixedList.H"
#include "vectorList.H"
#include "labelListIOList.H"
#include "instantList.H"
#include "vectorField.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class pointToPointPlanarInterpolation Declaration
\*---------------------------------------------------------------------------*/

class trilinearInterpolation
{
    // Private data

		//- Wheter use the value at the nearest vertex
		bool Nearest_;
		
        //- Number of source points (for checking)
        label nPoints_;

        //- Current interpolation addressing to the hex vertices
        List<FixedList<label, 8>> boundingVertex_;

        //- Current interpolation factors to the pyramid vertices
        List<FixedList<scalar, 8>> boundingWeight_;

    // Private Member Functions

		//- Find the match cell
		void matchCells
		(
			const pointField& sourcePoint,
			const labelListList& sourceIndices,
			const pointField& destPoints
		);
		
        //- Calculate weights
        void calcWeights
        (
            const pointField& sourcePoints,
            const pointField& destPoints
        );
        
        //- Test if the point is within the hex
        bool inHex
        (
			const FixedList<vector,8>& vs,
			const vector& c0
        );

public:

    // Declare name of the class and its debug switch
    ClassName("trilinearInterpolation");


    // Constructors

        //- Construct from source points and corrresponding indices
        
        trilinearInterpolation
        (
            const pointField& sourcePoints,
            const labelListList& sourceIndices,
            const pointField& destPoints,
            bool nearestSwitch
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
