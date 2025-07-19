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
    First extract cells using toposet then access selected cells.
    Calculate fields at selected cells.

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "cellSet.H"
#include "meshTools.H"
#include "treeBoundBox.H"
#include "treeDataCell.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


IOdictionary filterProperties
(
    IOobject
    (
        "filterProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const scalar delta = readScalar(filterProperties.lookup("filterWidth"));
const scalar beta = readScalar(filterProperties.lookup("decayRate")); //decay rate for Gaussian filter
const scalar radius = delta / 2.0;
const scalar radiusSqr = sqr(radius);
const vector rVec(radius, radius, radius);


const pointField& cellCenters = mesh.C();
const scalarField& cellVolumes = mesh.V();
const indexedOctree<treeDataCell>& cellTree = mesh.cellTree();


const cellSet myCellSet(mesh,"filteredCellSet");
const labelList& cellsToFilter = myCellSet.toc();


const instantList& timeDirs = timeSelector::select0(runTime,args);

forAll (timeDirs,timeI)
     {   
        //set the time to the value of current time directory.
        runTime.setTime(timeDirs[timeI],timeI);
        Info << "Time = " << runTime.timeName() <<endl;

       
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
        
        
        volVectorField filteredU
        (
            IOobject
            (
                "filteredU",
                 runTime.timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("filteredU", U.dimensions(), vector::zero)
         );


         
           forAll(cellsToFilter, i)
           {
           	const label cellI = cellsToFilter[i];
                const point& centre = cellCenters[cellI];
                const boundBox box(centre - rVec, centre + rVec);
                const treeBoundBox searchBox(box);
                const labelList nearCells = cellTree.findBox(searchBox);
                
                vector sumPhiVol = vector::zero;
                scalar sumVol = 0.0;
                
                for (label neighborCellI : nearCells)
                {
                    const scalar distSqr = magSqr(cellCenters[neighborCellI] - centre);
                    
                    if (distSqr <= radiusSqr)
                    {                                             
                       const scalar rSqr = distSqr/radiusSqr;
                       const scalar weight = Foam::exp(-beta*rSqr) * cellVolumes[neighborCellI];
                       sumPhiVol += U[neighborCellI] * weight;
                       sumVol += weight;
                    }
                }
         
		if (sumVol > VSMALL)
                {
                    filteredU[cellI] = sumPhiVol / sumVol;
                }
                
                Info << "cell number = " << i << endl;
           }
           
           filteredU.boundaryFieldRef() = U.boundaryField();
           filteredU.correctBoundaryConditions();
           filteredU.write();
	}
        Info<< "\nEnd\n" << endl;
	return 0;
}


// ************************************************************************* //
