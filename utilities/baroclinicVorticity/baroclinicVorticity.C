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

const scalar Ri = readScalar(transportProperties.lookup("Ri"));

dimensionedVector g
(
    "g",
    transportProperties.lookup("g")
);

dimensionedScalar Tc
(
    "Tc",
    transportProperties.lookup("Tc")
);

dimensionedScalar Th
(
    "Th",
    transportProperties.lookup("Th")
);

dimensionedScalar dT = Th - Tc;
    
//select time directories starting with 0.
const instantList& timeDirs = timeSelector::select0(runTime,args);

    forAll (timeDirs,timeI)
     {   
        //set the time to the value of current time directory.
        runTime.setTime(timeDirs[timeI],timeI);
        Info << "Time" << runTime.timeName() <<endl;


        //Read the temperature field.
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
      
     


         volVectorField gradT("gradT",fvc::grad(T)/dT); 
         volVectorField baroclinicVorticity("baroclinicVorticity",Ri*(gradT^g));    
         baroclinicVorticity.write();
      }

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
