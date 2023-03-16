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

#include "trilinearInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(trilinearInterpolation, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::trilinearInterpolation:: inHex
(
	const FixedList<vector,8>& vs,
	const vector& c0
)
{
	scalar minX = vs[0].x();
	scalar maxX = vs[3].x();
	scalar minY = vs[0].y();
	scalar maxY = vs[1].y();
	scalar minZ = vs[0].z();
	scalar maxZ = vs[5].z();
	
	if (c0.x()<minX || c0.x()>maxX) return false;
	if (c0.y()<minY || c0.y()>maxY) return false;
	if (c0.z()<minZ || c0.z()>maxZ) return false;
	
	return true;
}

void Foam::trilinearInterpolation:: matchCells
(
	const pointField& sourcePoints,
	const labelListList& sourceIndices,
	const pointField& destPoints
)
{
	boundingVertex_.setSize(destPoints.size());
	
	FixedList<label,8> currIndex; 
	FixedList<vector,8> vs;
	
	scalar currDist;
	scalar minDist;
	vector v0;
	
	label si0 = 0;
	bool foundCell;
	
	forAll(destPoints,dpi)
	{	
		minDist = GREAT;
		foundCell = false;
		
		for(label sii=si0; sii<sourceIndices.size(); sii++)
		{		
			currIndex = sourceIndices[sii];
			vs[0] = sourcePoints[currIndex[0]];
			vs[1] = sourcePoints[currIndex[1]];
			vs[2] = sourcePoints[currIndex[2]];
			vs[3] = sourcePoints[currIndex[3]];
			vs[4] = sourcePoints[currIndex[4]];
			vs[5] = sourcePoints[currIndex[5]];
			vs[6] = sourcePoints[currIndex[6]];
			vs[7] = sourcePoints[currIndex[7]];		
			
			v0 = (vs[0]+vs[1]+vs[2]+vs[3]+vs[4]+vs[5]+vs[6]+vs[7]) / 8;	
			
			if (inHex(vs,destPoints[dpi]))
			{
				boundingVertex_[dpi] = sourceIndices[sii];
				foundCell = true;
				si0 = sii;
				break;
			}
			else
			{
				currDist = mag(v0-destPoints[dpi]);
				
				if (currDist < minDist)
				{
					minDist = currDist;
					boundingVertex_[dpi] = sourceIndices[sii];
				}
			}	
		}
		
		if (!foundCell)
		{
			for (label sii=0; sii<si0; sii++)
			{
				currIndex = sourceIndices[sii];
				vs[0] = sourcePoints[currIndex[0]];
				vs[1] = sourcePoints[currIndex[1]];
				vs[2] = sourcePoints[currIndex[2]];
				vs[3] = sourcePoints[currIndex[3]];
				vs[4] = sourcePoints[currIndex[4]];
				vs[5] = sourcePoints[currIndex[5]];
				vs[6] = sourcePoints[currIndex[6]];
				vs[7] = sourcePoints[currIndex[7]];		
				
				v0 = (vs[0]+vs[1]+vs[2]+vs[3]+vs[4]+vs[5]+vs[6]+vs[7]) / 8;	
				
				if (inHex(vs,destPoints[dpi]))
				{
					boundingVertex_[dpi] = sourceIndices[sii];
					foundCell = true;
					si0 = sii;
					break;
				}
				else
				{
					currDist = mag(v0-destPoints[dpi]);
					
					if (currDist < minDist)
					{
						minDist = currDist;
						boundingVertex_[dpi] = sourceIndices[sii];
					}
				}				
			}
		}
	}
}

void Foam::trilinearInterpolation::calcWeights
(
    const pointField& sourcePoints,
    const pointField& destPoints
)
{
	FixedList<label,8> currIndex;
	FixedList<vector,8> vs;
	
	vector c0;
	scalar dom, minDist;
	label minDistIndict;
	
    boundingWeight_.setSize(destPoints.size());
 
	forAll(destPoints,dpi)
	{
		currIndex = boundingVertex_[dpi];
		
		vs[0] = sourcePoints[currIndex[0]];
		vs[1] = sourcePoints[currIndex[1]];
		vs[2] = sourcePoints[currIndex[2]];
		vs[3] = sourcePoints[currIndex[3]];
		vs[4] = sourcePoints[currIndex[4]];
		vs[5] = sourcePoints[currIndex[5]];
		vs[6] = sourcePoints[currIndex[6]];
		vs[7] = sourcePoints[currIndex[7]];
		
		c0 = destPoints[dpi];
		
		if (inHex(vs,c0) && !Nearest_)
		{
			dom = (vs[3]-vs[0]).x()*(vs[1]-vs[0]).y()*(vs[4]-vs[0]).z();	
					
			boundingWeight_[dpi][0] = (vs[3]-c0).x()*(vs[1]-c0).y()*(vs[4]-c0).z() / dom;
			boundingWeight_[dpi][1] = (vs[3]-c0).x()*(c0-vs[0]).y()*(vs[4]-c0).z() / dom;
			boundingWeight_[dpi][2] = (c0-vs[0]).x()*(c0-vs[0]).y()*(vs[4]-c0).z() / dom;
			boundingWeight_[dpi][3] = (c0-vs[0]).x()*(vs[1]-c0).y()*(vs[4]-c0).z() / dom;
			boundingWeight_[dpi][4] = (vs[3]-c0).x()*(vs[1]-c0).y()*(c0-vs[0]).z() / dom;
			boundingWeight_[dpi][5] = (vs[3]-c0).x()*(c0-vs[0]).y()*(c0-vs[0]).z() / dom;
			boundingWeight_[dpi][6] = (c0-vs[0]).x()*(c0-vs[0]).y()*(c0-vs[0]).z() / dom;
			boundingWeight_[dpi][7] = (c0-vs[0]).x()*(vs[1]-c0).y()*(c0-vs[0]).z() / dom;
		}
		else
		{
			minDist = mag(c0-vs[0]);
			minDistIndict = 0;
			
			forAll(vs,vi)
			{
				if (mag(c0-vs[vi]) <= minDist)
				{
					minDistIndict = vi;
					minDist = mag(c0-vs[vi]);
				}
			}
			
			boundingWeight_[dpi] = 0.0;
			boundingWeight_[dpi][minDistIndict] = 1.0;		
		}
		
	}

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::trilinearInterpolation::trilinearInterpolation
(
    const pointField& sourcePoints,
    const labelListList& sourceIndices,
    const pointField& destPoints,
    bool nearestSwitch
)
{
	Nearest_ = nearestSwitch;
	matchCells(sourcePoints,sourceIndices,destPoints);
    calcWeights(sourcePoints, destPoints);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::trilinearInterpolation::interpolate
(
    const vectorField& sourceFld
) const
{
    tmp<vectorField> tfld(new vectorField(boundingVertex_.size()));
    vectorField& fld = tfld.ref();

	forAll(fld, fi)
	{
		const FixedList<label, 8>& verts = boundingVertex_[fi];
		const FixedList<scalar, 8>& w = boundingWeight_[fi];

		fld[fi] =   w[0]*sourceFld[verts[0]]
				  + w[1]*sourceFld[verts[1]]
				  + w[2]*sourceFld[verts[2]]
				  + w[3]*sourceFld[verts[3]]
				  + w[4]*sourceFld[verts[4]]
				  + w[5]*sourceFld[verts[5]]
				  + w[6]*sourceFld[verts[6]]
				  + w[7]*sourceFld[verts[7]];
	}
		
    return tfld;
}

Foam::tmp<Foam::scalarField> Foam::trilinearInterpolation::interpolate
(
    const scalarField& sourceFld
) const
{
    tmp<scalarField> tfld(new scalarField(boundingVertex_.size()));
    scalarField& fld = tfld.ref();

	forAll(fld, fi)
	{
		const FixedList<label, 8>& verts = boundingVertex_[fi];
		const FixedList<scalar, 8>& w = boundingWeight_[fi];

		fld[fi] =   w[0]*sourceFld[verts[0]]
				  + w[1]*sourceFld[verts[1]]
				  + w[2]*sourceFld[verts[2]]
				  + w[3]*sourceFld[verts[3]]
				  + w[4]*sourceFld[verts[4]]
				  + w[5]*sourceFld[verts[5]]
				  + w[6]*sourceFld[verts[6]]
				  + w[7]*sourceFld[verts[7]];
	}
		
    return tfld;
}

// ************************************************************************* //
