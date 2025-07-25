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


Description
Eigenvalues of strain rate tensor

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
    
//select time directories starting with 0.
const instantList& timeDirs = timeSelector::select0(runTime,args);

     forAll (timeDirs,timeI)
     {   
        //set the time to the value of current time directory.
        runTime.setTime(timeDirs[timeI],timeI);
        Info << "Time" << runTime.timeName() <<endl;

       //Read the U field.
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
         
         
        volSymmTensorField strainRate("strainRate",symm(fvc::grad(U)));
        
        volVectorField EigStrainRate
        (
            IOobject
            (
                 "EigStrainRate",
                 runTime.timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE
             ),
             mesh,
             dimensionedVector("EigStrainRate",strainRate.dimensions(), Zero),
             "zeroGradient"
         );


         
         forAll(strainRate.internalField(),cellI)
         {  
              EigStrainRate[cellI] = eigenValues(strainRate.internalField()[cellI]);
         }
         
        // strainRate.write();
         EigStrainRate.write();

      }
                   

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
