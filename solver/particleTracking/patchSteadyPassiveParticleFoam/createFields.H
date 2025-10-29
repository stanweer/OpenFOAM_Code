IOdictionary particleTrackingDict
(
    IOobject
    (
        "particleTrackingDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

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

Info<< "Reading field nut\n" << endl;
volScalarField nut
(
    IOobject
    (
        "nut",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading field EigStress\n" << endl;
volVectorField EigStress
(
    IOobject
    (
        "EigStress",
         runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


Info<< "Reading field MagEigStress\n" << endl;
volScalarField MagEigStress
(
    IOobject
    (
        "MagEigStress",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


Info<< "Reading field MagTotalStress\n" << endl;
volScalarField MagTotalStress
(
    IOobject
    (
        "MagTotalStress",
         runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


#include "createPhi.H"
#include "readGravitationalAcceleration.H"
singlePhaseTransportModel laminarTransport(U, phi);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rho", dimDensity, 1.0)
);

const word kinematicCloudName
(
    args.getOrDefault<word>("cloud", "kinematicCloud")
);

Info<< "Constructing kinematicCloud " << kinematicCloudName << endl;
basicKinematicCloud kinematicCloud
(
    kinematicCloudName,
    rho,
    U,
    laminarTransport.nu(),
    g 
);


