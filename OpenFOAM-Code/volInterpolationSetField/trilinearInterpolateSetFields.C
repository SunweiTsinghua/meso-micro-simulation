/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    trilinearInterpolationSetFields

Author
    Sunwei Li, Tsinghua University.  All rights reserved.

Description
    Set the initial fields from the external data from the trilinear interpolation

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "vectorIOField.H"
#include "scalarIOField.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "trilinearInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

	// Reading the configurations from the interpDict and configure the executable
	IOdictionary interpDict
	(
		IOobject
		(
			"interpDict",
			runTime.system(),
			runTime.db(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	bool nearest = interpDict.lookup("nearestSwitch");
	wordList setFields = wordList(interpDict.lookup("setFields"));
	wordList setPatches = wordList(interpDict.lookup("setPatches"));

	// Reading the points where the interpolated data is available
	Info << "Reading points\n" << endl;
	pointIOField inputPoints
	(
		IOobject
		(
			"points",
			runTime.constant(),
			"domainData",
			runTime.db(),
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE,
			false
		)
	);

	Info << "Reading cell labels\n" << endl;
	labelListIOList inputIndices
	(
		IOobject
		(
			"indices",
			runTime.constant(),
			"domainData",
			runTime.db(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE,
			false
		)
	);

	const volVectorField& cellCenters = mesh.C();

	Info << "Calculating the trilinear interpolation\n" << endl;
	trilinearInterpolation volInterp 
	(
		inputPoints,
		inputIndices,		
		cellCenters,
		nearest
	);

	forAll(setFields, fi)
	{
		if(setFields[fi] == "U")
		{
			Info<< "Reading field U\n" << endl;
			volVectorField U
			(
				IOobject
				(
					"U",
					runTime.timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				mesh
			);  

			Info<< "Reading input U\n" << endl;
			vectorIOField inputU
			(
				IOobject
				(
					"inputU",
					runTime.constant(),
					"domainData",
					runTime.db(),
					IOobject::MUST_READ,
					IOobject::NO_WRITE,
					false
				)
			);	

			#include "createPhi.H"
			
			singlePhaseTransportModel laminarTransport(U, phi);	
			autoPtr<incompressible::turbulenceModel> turbulence
			(
				incompressible::turbulenceModel::New(U, phi, laminarTransport)
			);

			turbulence->validate();	

			Info << "Interpolating the internal field for U\n" << endl;
			vectorField interpU = volInterp.interpolate(inputU);
			
			forAll(cellCenters,celli)
			{
				U[celli] = interpU[celli];
			}
			Info << "Interpolating the boundary fields for U\n" << endl;
			typename GeometricField<vector, fvPatchField, volMesh>::Boundary& 
					bdyRefU = U.boundaryFieldRef();

			forAll(setPatches, seti)
			{
				label currPatchID = mesh.boundaryMesh().findPatchID(setPatches[seti]);

				if (currPatchID == -1)
				{
					FatalErrorIn("main()")
					<< "the patch of " << setPatches[seti] << " is not in the mesh, please check\n"
					<< exit(FatalError);	
				}
				else
				{
					const vectorField& faceCenters = mesh.boundary()[currPatchID].Cf();
					trilinearInterpolation patchInterp
					(
						inputPoints,
						inputIndices,
						faceCenters,
						nearest
					);
					vectorField interpPatchU = patchInterp.interpolate(inputU);
					bdyRefU[currPatchID] == interpPatchU;					
				}
			}
			U.correctBoundaryConditions();	
			U.write();				
		}
		else if (setFields[fi] == "p")
		{
			Info<< "Reading field p\n" << endl;
			volScalarField p
			(
				IOobject
				(
					"p",
					runTime.timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				mesh
			); 

			Info<< "Reading input p\n" << endl;
			scalarIOField inputp
			(
				IOobject
				(
					"inputP",
					runTime.constant(),
					"domainData",
					runTime.db(),
					IOobject::MUST_READ,
					IOobject::NO_WRITE,
					false
				)
			);

			Info << "Interpolating the internal field for p\n" << endl;
			scalarField interpp = volInterp.interpolate(inputp);
			
			forAll(cellCenters,celli)
			{
				p[celli] = interpp[celli];
			}

			Info << "Interpolating the boundary fields for p\n" << endl;
			typename GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
					bdyRefp = p.boundaryFieldRef();

			forAll(setPatches, seti)
			{
				label currPatchID = mesh.boundaryMesh().findPatchID(setPatches[seti]);

				if (currPatchID == -1)
				{
					FatalErrorIn("main()")
					<< "the patch of " << setPatches[seti] << " is not in the mesh, please check\n"
					<< exit(FatalError);	
				}
				else
				{
					const vectorField& faceCenters = mesh.boundary()[currPatchID].Cf();
					trilinearInterpolation patchInterp 
					(
						inputPoints,
						inputIndices,
						faceCenters,
						nearest
					);
					scalarField interpPatchp = patchInterp.interpolate(inputp);
					bdyRefp[currPatchID] == interpPatchp;					
				}
			}

			p.write();				
		}		
		else if (setFields[fi] == "k")
		{
			Info<< "Reading field k\n" << endl;
			volScalarField k
			(
				IOobject
				(
					"k",
					runTime.timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				mesh
			); 

			Info<< "Reading input k\n" << endl;
			scalarIOField inputk
			(
				IOobject
				(
					"inputK",
					runTime.constant(),
					"domainData",
					runTime.db(),
					IOobject::MUST_READ,
					IOobject::NO_WRITE,
					false
				)
			);

			Info << "Interpolating the internal field for k\n" << endl;
			scalarField interpk = volInterp.interpolate(inputk);
			
			forAll(cellCenters,celli)
			{
				k[celli] = interpk[celli];
			}

			Info << "Interpolating the boundary fields for k\n" << endl;
			typename GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
					bdyRefk = k.boundaryFieldRef();

			forAll(setPatches, seti)
			{
				label currPatchID = mesh.boundaryMesh().findPatchID(setPatches[seti]);

				if (currPatchID == -1)
				{
					FatalErrorIn("main()")
					<< "the patch of " << setPatches[seti] << " is not in the mesh, please check\n"
					<< exit(FatalError);	
				}
				else
				{
					const vectorField& faceCenters = mesh.boundary()[currPatchID].Cf();
					trilinearInterpolation patchInterp 
					(
						inputPoints,
						inputIndices,
						faceCenters,
						nearest
					);
					scalarField interpPatchk = patchInterp.interpolate(inputk);
					bdyRefk[currPatchID] == interpPatchk;					
				}
			}

			k.write();				
		}
		else if (setFields[fi] == "epsilon")
		{
			Info<< "Reading field epsilon\n" << endl;
			volScalarField epsilon 
			(
				IOobject
				(
					"epsilon",
					runTime.timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				mesh
			);

			Info<< "Reading input epsilon\n" << endl;
			scalarIOField inputEpsilon
			(
				IOobject
				(
					"inputEpsilon",
					runTime.constant(),
					"domainData",
					runTime.db(),
					IOobject::MUST_READ,
					IOobject::NO_WRITE,
					false
				)
			);

			Info << "Interpolating the internal field for epsilon\n" << endl;
			scalarField interpEpsilon = volInterp.interpolate(inputEpsilon);
			
			forAll(cellCenters,celli)
			{
				epsilon[celli] = interpEpsilon[celli];
			}

			Info << "Interpolating the boundary fields for epsilon\n" << endl;
			typename GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
					bdyRefE = epsilon.boundaryFieldRef();

			forAll(setPatches, seti)
			{
				label currPatchID = mesh.boundaryMesh().findPatchID(setPatches[seti]);

				if (currPatchID == -1)
				{
					FatalErrorIn("main()")
					<< "the patch of " << setPatches[seti] << " is not in the mesh, please check\n"
					<< exit(FatalError);	
				}
				else
				{
					const vectorField& faceCenters = mesh.boundary()[currPatchID].Cf();
					trilinearInterpolation patchInterp 
					(
						inputPoints,
						inputIndices,
						faceCenters,
						nearest
					);
					scalarField interpPatchE = patchInterp.interpolate(inputEpsilon);
					bdyRefE[currPatchID] == interpPatchE;					
				}
			}

			epsilon.write();				
		}
		else if (setFields[fi] == "T")
		{
			Info<< "Reading field T\n" << endl;
			volScalarField T
			(
				IOobject
				(
						"T",
						runTime.timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::AUTO_WRITE
				),
				mesh
			);

			Info<< "Reading input T\n" << endl;
			scalarIOField inputT
			(
				IOobject
				(
					"inputT",
					runTime.constant(),
					"domainData",
					runTime.db(),
					IOobject::MUST_READ,
					IOobject::NO_WRITE,
					false
				)
			);

			Info << "Interpolating the internal field for T\n" << endl;
			scalarField interpT = volInterp.interpolate(inputT);
			
			forAll(cellCenters,celli)
			{
				T[celli] = interpT[celli];
			}

			Info << "Interpolating the boundary fields for T\n" << endl;
			typename GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
					bdyRefT = T.boundaryFieldRef();

			forAll(setPatches, seti)
			{
				label currPatchID = mesh.boundaryMesh().findPatchID(setPatches[seti]);

				if (currPatchID == -1)
				{
					FatalErrorIn("main()")
					<< "the patch of " << setPatches[seti] << " is not in the mesh, please check\n"
					<< exit(FatalError);	
				}
				else
				{
					const vectorField& faceCenters = mesh.boundary()[currPatchID].Cf();
					trilinearInterpolation patchInterp 
					(
						inputPoints,
						inputIndices,
						faceCenters,
						nearest
					);
					scalarField interpPatchT = patchInterp.interpolate(inputT);
					bdyRefT[currPatchID] == interpPatchT;					
				}
			}

			T.correctBoundaryConditions();
			T.write();				
		}
		else
		{
			FatalErrorIn("main()")
				<< "the field of " << setFields[fi] << " is currently not support\n"
				<< exit(FatalError);	
		}
	}

	Info << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
	<< nl << endl;

	Info<< "End\n" << endl;

	return(0);
}


// ************************************************************************* //
