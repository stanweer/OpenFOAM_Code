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
    fieldName       U;          // Name of the field to filter
    filterRegion    mySet;     // Name of the cellSet to process
    filterWidth     0.1; 
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

// Read filter properties from constant/filterProperties dictionary
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

const word fieldName(filterProperties.lookup("fieldName"));    // Field to filter
const word filterRegion(filterProperties.lookup("filterRegion"));  // cellSet region to apply filter
const scalar delta = readScalar(filterProperties.lookup("filterWidth"));    // Filter width

 // Precompute constants for Gaussian filter
const scalar deltaSqr = sqr(delta);
const scalar radius = 2.0 * delta / Foam::sqrt(12.0);  //2 sigma (95%)
const scalar radiusSqr = sqr(radius);
const vector radiusVec(radius, radius, radius);
const scalar coeff = -6.0 / deltaSqr;


const pointField& cellCenters = mesh.C();
const scalarField& cellVolumes = mesh.V();
const indexedOctree<treeDataCell>& cellTree = mesh.cellTree();   // Octree for neighbor search

// Get cells where filtering will be applied
const cellSet filterCellSet(mesh,filterRegion);
const labelList& cellsToFilter = filterCellSet.toc();


const instantList& timeDirs = timeSelector::select0(runTime,args);

forAll (timeDirs,timeI)
     {   
     
        runTime.setTime(timeDirs[timeI],timeI);
        Info << "Time = " << runTime.timeName() <<endl;

       
        volVectorField originalField
        (
            IOobject
            (  
                fieldName,
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        
        Info << "Filtering field: " << originalField.name() <<endl;
        
        volVectorField filteredField
        (
            IOobject
            (
                "filteredField",
                 runTime.timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector("filteredField", originalField.dimensions(), vector::zero)
         );



           
           // Loop over all cells in the filter region
           forAll(cellsToFilter, i)
           {
           	const label cellI = cellsToFilter[i];
                const point& centre = cellCenters[cellI];
                
                // Define search box around the current cell
                const boundBox box(centre - radiusVec, centre + radiusVec);
                const treeBoundBox searchBox(box);
                
                // Find neighboring cells inside the search box using octree
                const labelList nearCells = cellTree.findBox(searchBox);
                
                
                vector sumPhiVol = vector::zero;   // Weighted sum of field values
                scalar sumVol = 0.0;              // Sum of weights (normalization)
                
                //for (label neighborCellI : nearCells)
                // Loop over neighboring cells
                forAll(nearCells, n)
                {
                    const label neighborCellI = nearCells[n];
                    const scalar distSqr = magSqr(cellCenters[neighborCellI] - centre);
                    
                    if (distSqr <= radiusSqr)  // Check if inside filter radius
                    {                                             
                       const scalar weight = Foam::exp(coeff*distSqr) * cellVolumes[neighborCellI];
                       sumPhiVol += originalField[neighborCellI] * weight;
                       sumVol += weight;
                    }
                }
         
		if (sumVol > VSMALL)
                {
                    filteredField[cellI] = sumPhiVol / sumVol;
                }
                
                if (i % 10000 == 0) Info << "Processed cell " << i << endl;
           }
           
           filteredField.boundaryFieldRef() = originalField.boundaryField();
           filteredField.correctBoundaryConditions();
           filteredField.write();
	}
        Info<< "\nEnd\n" << endl;
	return 0;
}


// ************************************************************************* //
