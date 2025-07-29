/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

Application
    baroclinic production of vorticity.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

    dimensionedScalar Ri
    (
        "Ri",
        transportProperties.get<dimensionedScalar>("Ri")
    );

    dimensionedScalar refL
    (
        "refL",
        transportProperties.get<dimensionedScalar>("refL")
    );

    dimensionedVector g
    (
        "g",
        transportProperties.get<dimensionedVector>("g")
    );

    dimensionedScalar Tc
    (
        "Tc",
        transportProperties.get<dimensionedScalar>("Tc")
    );

    dimensionedScalar Th
    (
        "Th",
        transportProperties.get<dimensionedScalar>("Th")
    );

    dimensionedScalar dT("dT", Th - Tc);

    // Select time directories starting with 0
    instantList timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info << "Time = " << runTime.timeName() << endl;

   
        Info << "Reading field T" << endl;
        volScalarField T
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE 
            ),
            mesh
        );

        // Calculate normalized temperature gradient
        volVectorField gradT ("gradT", (refL / dT) * fvc::grad(T));
        

        // Calculate baroclinic vorticity: Ri * (gradT ^ g)
        volVectorField baroclinicVorticity ("baroclinicVorticity",  Ri * (gradT ^ g));
        baroclinicVorticity.write();

        // Clear large fields to free memory
        gradT.clear();
        T.clear();
    }

    Info << "\nEnd\n" << endl;
    return 0;
}

// ************************************************************************* //
