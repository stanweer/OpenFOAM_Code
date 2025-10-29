/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "Random.H"
#include "OFstream.H" 
#include "cellSet.H"
#include "basicKinematicCloud.H"
#include "interpolationCellPoint.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//function for particle injection on a patch
void injectParticles(label numParticles, Random& rng, label patchID, const fvMesh& mesh, basicKinematicCloud& kinematicCloud)
{
    const vectorField faceCentres = mesh.boundaryMesh()[patchID].faceCentres();
    const labelList& faceCells = mesh.boundaryMesh()[patchID].faceCells();

    for (label i = 0; i < numParticles; ++i)
    {
        // Select a random face on the patch
        const label randomFaceIndex = rng.position<label>(0, faceCentres.size() - 1);
        const point position = faceCentres[randomFaceIndex];
        const label cellID = faceCells[randomFaceIndex]; // Cell adjacent to the face

        // Create a new particle
        KinematicParcel<particle>* p = new KinematicParcel<particle>(mesh, position, cellID);

        p->U() = vector::zero;  //particle velocity
        p->d() = 1;             //particle diameter  
        p->rho() = 1;           //particle density
        p->nParticle() = 1;     // number of particle in a parcel
        
        kinematicCloud.addParticle(p);
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow,"
        " using the PISO algorithm."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//reading the parameters from particleTrackingDict
const word patchName(particleTrackingDict.lookup("particleInjectionPatch"));
const label numParticles = readLabel(particleTrackingDict.lookup("numParticles"));
const scalar injectionInterval = readScalar(particleTrackingDict.lookup("injectionInterval"));    
const scalar particleWriteInterval = readScalar(particleTrackingDict.lookup("particleWriteInterval"));

label patchID = mesh.boundaryMesh().findPatchID(patchName);
if (patchID == -1)
{
     FatalErrorInFunction
         << "Patch " << patchName << " not found!" << exit(FatalError);
}

scalar injectionTime = runTime.startTime().value();
scalar particleWriteTime = runTime.startTime().value();

Random rng; //for random number
injectParticles(numParticles, rng, patchID, mesh, kinematicCloud);


//-------------------------------------------------------------------------------------------------   
    

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;       
        kinematicCloud.evolve();
        
                
        interpolationCellPoint<vector> UInterp(U);  
               
        forAllIter(typename basicKinematicCloud, kinematicCloud, iter)
        {
            iter().U() = UInterp.interpolate(iter().position(), iter().cell());
        } 
        
        
        
        if (runTime.time().value() >= particleWriteTime)
        {
            //writing particle data in a file
            OStringStream fileNameStream;
            fileNameStream << "particleData_"
                     << setw(4) << setfill('0') << runTime.timeIndex()
                     << ".csv";

            OFstream particleDataFile(runTime.path() / fileNameStream.str());
            particleDataFile << "ID" << "," << "age" << "," << "x_position" << "," << "y_position" << "," << "z_position" << "," << "U_x" << "," << "U_y" << "," << "U_z" << "," << "nut" << "," << "EigStress_x" << "," << "EigStress_y" << "," << "EigStress_z" << "," << "MagEigStress" << "," << "MagTotalStress" << nl;
            
            interpolationCellPoint<scalar> nutInterp(nut);
            interpolationCellPoint<vector> EigStressInterp(EigStress);
            interpolationCellPoint<scalar> MagEigStressInterp(MagEigStress);
            interpolationCellPoint<scalar> MagTotalStressInterp(MagTotalStress);
            


            
            forAllIter(basicKinematicCloud, kinematicCloud, iter)
            {
                                
                scalar nut_particle = nutInterp.interpolate(iter().position(), iter().cell());
                vector EigStress_particle = EigStressInterp.interpolate(iter().position(), iter().cell());
                scalar MagEigStress_particle = MagEigStressInterp.interpolate(iter().position(), iter().cell());
                scalar MagTotalStress_particle = MagTotalStressInterp.interpolate(iter().position(), iter().cell());
                
                particleDataFile 
                    << iter().origId() << ","
                    << iter().age() << ","
                    << iter().position().x() << ","
                    << iter().position().y() << ","
                    << iter().position().z() << ","
                    << iter().U().x() << ","
                    << iter().U().y() << ","
                    << iter().U().z() << ","
                    << nut_particle << ","
                    << EigStress_particle.x() << ","
                    << EigStress_particle.y() << ","
                    << EigStress_particle.z() << ","
                    << MagEigStress_particle << ","
                    << MagTotalStress_particle << nl;
            }
            
            particleWriteTime += particleWriteInterval;
        }
               
                
        
        if (runTime.time().value() >= injectionTime)
        { 
            injectParticles(numParticles, rng, patchID, mesh, kinematicCloud);
            injectionTime += injectionInterval;
        }   

//------------------------------------------------------------------------------------------------------------------------------

        
        runTime.write();
        runTime.printExecutionTime(Info);        
    }

    Info<< "End\n" << endl;

    return 0;
}



// ************************************************************************* //
