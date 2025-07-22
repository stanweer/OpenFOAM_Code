/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Application
  
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    const instantList& timeDirs = timeSelector::select0(runTime, args);
    const fvPatchList& patches = mesh.boundary();

    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
        "nu",
        dimViscosity,
        transportProperties
    );

    volVectorField wallShearStress
    (
        IOobject
        (
            "wallShearStress",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(dimPressure, Zero)
    );

    volVectorField::Boundary& wallShearStressBf = wallShearStress.boundaryFieldRef();

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        Info<< "Reading field U\n" << endl;
        volVectorField U
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        volTensorField gradU(fvc::grad(U));
        volSymmTensorField stress(nu*twoSymm(gradU));

        forAll(patches, patchi)
        {
            if (isA<wallFvPatch>(patches[patchi]))
            {
                const vectorField& n = patches[patchi].nf();
                const symmTensorField& stressPatch = stress.boundaryField()[patchi];
                
                // Calculate wall shear stress: tau = (I - nn) & (n & Reff)
                wallShearStressBf[patchi] = (n & stressPatch) - (n & stressPatch & n)*n;
            }
        }

        wallShearStress.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
