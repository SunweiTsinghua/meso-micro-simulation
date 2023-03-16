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

\*---------------------------------------------------------------------------*/

#include "pointToPointVolumeInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointToPointVolumeInterpolation, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::pointToPointVolumeInterpolation::inTet
(
	const FixedList<vector,4>& tetPoints,
	const vector& inPoint
)
{
	FixedList<scalarSquareMatrix,5> tetMatrix;
	scalarSquareMatrix tmpMatrix;
	
	for (label mi=0; mi<5; mi++) tetMatrix[mi].setSize(4);
	tmpMatrix.setSize(4);
	
	for (label pi=0; pi<4; pi++)
	{
		tetMatrix[0](pi,0) = tetPoints[pi].x();
		tetMatrix[0](pi,1) = tetPoints[pi].y();
		tetMatrix[0](pi,2) = tetPoints[pi].z();
		tetMatrix[0](pi,3) = 1.0;
	}
	
	tmpMatrix = tetMatrix[0];
		
	if (det(tmpMatrix) < 0) 
	{
		tetMatrix[0](0,0) = tetPoints[1].x();
		tetMatrix[0](0,1) = tetPoints[1].y();
		tetMatrix[0](0,2) = tetPoints[1].z();
		
		tetMatrix[0](1,0) = tetPoints[0].x();
		tetMatrix[0](1,1) = tetPoints[0].y();
		tetMatrix[0](1,2) = tetPoints[0].z();
	}
			
	tetMatrix[1](0,0) = inPoint.x();
	tetMatrix[1](0,1) = inPoint.y();
	tetMatrix[1](0,2) = inPoint.z();
	tetMatrix[1](0,3) = 1.0;
	
	for (label pi=0; pi<4; pi++)
	{
		if (pi==0) continue;
		
		tetMatrix[1](pi,0) = tetPoints[pi].x();
		tetMatrix[1](pi,1) = tetPoints[pi].y();
		tetMatrix[1](pi,2) = tetPoints[pi].z();
		tetMatrix[1](pi,3) = 1.0;
	}
	
	tetMatrix[2](1,0) = inPoint.x();
	tetMatrix[2](1,1) = inPoint.y();
	tetMatrix[2](1,2) = inPoint.z();
	tetMatrix[2](1,3) = 1.0;
	
	for (label pi=0; pi<4; pi++)
	{
		if (pi==1) continue;
		
		tetMatrix[2](pi,0) = tetPoints[pi].x();
		tetMatrix[2](pi,1) = tetPoints[pi].y();
		tetMatrix[2](pi,2) = tetPoints[pi].z();
		tetMatrix[2](pi,3) = 1.0;
	}
	
	tetMatrix[3](2,0) = inPoint.x();
	tetMatrix[3](2,1) = inPoint.y();
	tetMatrix[3](2,2) = inPoint.z();
	tetMatrix[3](2,3) = 1.0;
	
	for (label pi=0; pi<4; pi++)
	{
		if (pi==2) continue;
		
		tetMatrix[3](pi,0) = tetPoints[pi].x();
		tetMatrix[3](pi,1) = tetPoints[pi].y();
		tetMatrix[3](pi,2) = tetPoints[pi].z();
		tetMatrix[3](pi,3) = 1.0;
	}
	
	tetMatrix[4](3,0) = inPoint.x();
	tetMatrix[4](3,1) = inPoint.y();
	tetMatrix[4](3,2) = inPoint.z();
	tetMatrix[4](3,3) = 1.0;
	
	for (label pi=0; pi<4; pi++)
	{
		if (pi==3) continue;
		
		tetMatrix[4](pi,0) = tetPoints[pi].x();
		tetMatrix[4](pi,1) = tetPoints[pi].y();
		tetMatrix[4](pi,2) = tetPoints[pi].z();
		tetMatrix[4](pi,3) = 1.0;
	}	
	
	label volumeIndict = 0;
	
	for (label mi=1; mi<5; mi++) 
	{
		tmpMatrix = tetMatrix[mi];
		volumeIndict += sign(det(tmpMatrix));
	}
	
	if (fabs(volumeIndict) == 4)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void Foam::pointToPointVolumeInterpolation::calcWeights
(
    const pointField& sourcePoints,
    const pointField& destPoints
)
{
	nearestVertex_.setSize(destPoints.size());
    nearestVertexWeight_.setSize(destPoints.size());
    
    // first found three vertices with shortest distances from the target
    
	FixedList<scalar, 3> nearestDist;
	label maxDistIndict;
	scalar maxDist;
	scalar currDist;
	scalar prevDist;
	
	vector v1, v2, v3;
	
	forAll(destPoints,dpi)
	{

		nearestDist[0] = mag(sourcePoints[0]-destPoints[dpi]);
		nearestDist[1] = mag(sourcePoints[1]-destPoints[dpi]);
		
		nearestVertex_[dpi][0] = 0;
		nearestVertex_[dpi][1] = 1;
		
		nearestVertexWeight_[dpi][0] = 1.0/nearestDist[0];
		nearestVertexWeight_[dpi][1] = 1.0/nearestDist[1];
		
		nearestDist[0]>nearestDist[1] ? maxDistIndict=0 : maxDistIndict=1;
		maxDist = nearestDist[maxDistIndict];
		
		for (label spi=2; spi<sourcePoints.size(); spi++)
		{
			currDist = mag(sourcePoints[spi]-destPoints[dpi]);
			
			if (currDist < maxDist)
			{
				nearestDist[maxDistIndict] = currDist;
				nearestVertex_[dpi][maxDistIndict] = spi;
				nearestVertexWeight_[dpi][maxDistIndict] = 1.0/currDist;
				
				nearestDist[0]>nearestDist[1] ? maxDistIndict=0 : maxDistIndict=1;
				maxDist = nearestDist[maxDistIndict];
			}
		}
		
		v1 = sourcePoints[nearestVertex_[dpi][1]] - sourcePoints[nearestVertex_[dpi][0]];
		prevDist = GREAT;
		
		forAll(sourcePoints,spi)
		{
			v2 = sourcePoints[spi] - sourcePoints[nearestVertex_[dpi][0]];
			
			if (mag(v1^v2) < SMALL) continue;
			
			currDist = mag(sourcePoints[spi]-destPoints[dpi]);
			
			if (currDist < prevDist)
			{
				nearestVertex_[dpi][2] = spi;
				nearestVertexWeight_[dpi][2] = 1.0/currDist;
				prevDist = currDist;
			}
			
		}
		
		v2 = sourcePoints[nearestVertex_[dpi][2]] - sourcePoints[nearestVertex_[dpi][0]];
		prevDist = GREAT;
		
		forAll(sourcePoints,spi)
		{
			v3 = sourcePoints[spi] - sourcePoints[nearestVertex_[dpi][0]];
			
			if (mag((v1^v2)&v3) < SMALL) continue;
		
			currDist = mag(sourcePoints[spi]-destPoints[dpi]);
			
			if (currDist < prevDist)
			{
				nearestVertex_[dpi][3] = spi;
				nearestVertexWeight_[dpi][3] = 1.0/currDist;
				prevDist = currDist;
			}
		}
		
		FixedList<vector,4> currTet;
		
		for (label ti=0; ti<4; ti++) currTet[ti] = sourcePoints[nearestVertex_[dpi][ti]];
		
		if (inTet(currTet,destPoints[dpi]))
		{
			scalar sumWeight = nearestVertexWeight_[dpi][0]+nearestVertexWeight_[dpi][1]+nearestVertexWeight_[dpi][2]+nearestVertexWeight_[dpi][3];
			
			nearestVertexWeight_[dpi][0] = nearestVertexWeight_[dpi][0] / sumWeight;
			nearestVertexWeight_[dpi][1] = nearestVertexWeight_[dpi][1] / sumWeight;
			nearestVertexWeight_[dpi][2] = nearestVertexWeight_[dpi][2] / sumWeight;
			nearestVertexWeight_[dpi][3] = nearestVertexWeight_[dpi][3] / sumWeight;
		
		}
		else
		{
			label minDistIndict;
			nearestDist[0]>nearestDist[1] ? minDistIndict=1 : minDistIndict=0;
			nearestVertexWeight_[dpi] = 0.0;
			nearestVertexWeight_[dpi][minDistIndict] = 1.0;
		}
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointToPointVolumeInterpolation::pointToPointVolumeInterpolation
(
    const pointField& sourcePoints,
    const pointField& destPoints
)
:
    nPoints_(sourcePoints.size())
{
    calcWeights(sourcePoints, destPoints);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::pointToPointVolumeInterpolation::interpolate
(
    const vectorField& sourceFld
) const
{
    tmp<vectorField> tfld(new vectorField(nearestVertex_.size()));
    vectorField& fld = tfld.ref();

    forAll(fld, fi)
    {
        const FixedList<label, 4>& verts = nearestVertex_[fi];
        const FixedList<scalar, 4>& w = nearestVertexWeight_[fi];

        fld[fi] =   w[0]*sourceFld[verts[0]]
                  + w[1]*sourceFld[verts[1]]
                  + w[2]*sourceFld[verts[2]]
                  + w[3]*sourceFld[verts[3]];
    }
    return tfld;
}

Foam::tmp<Foam::scalarField> Foam::pointToPointVolumeInterpolation::interpolate
(
    const scalarField& sourceFld
) const
{
    tmp<scalarField> tfld(new scalarField(nearestVertex_.size()));
    scalarField& fld = tfld.ref();

    forAll(fld, fi)
    {
        const FixedList<label, 4>& verts = nearestVertex_[fi];
        const FixedList<scalar, 4>& w = nearestVertexWeight_[fi];

        fld[fi] =   w[0]*sourceFld[verts[0]]
                  + w[1]*sourceFld[verts[1]]
                  + w[2]*sourceFld[verts[2]]
                  + w[3]*sourceFld[verts[3]];
    }
    return tfld;
}

// ************************************************************************* //
