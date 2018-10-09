/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2009-2012 JKU, Linz
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Application
    cfdemSolverPimple

Description
    Transient solver for compressible flow.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver pimpleFoam in OpenFOAM(R) 1.6,
    where additional functionality for CFD-DEM coupling is added.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"
#include "turbulentFluidThermoModel.H"
#include "bound.H"
#include "pimpleControl.H"
#include "pressureControl.H"

    #include "OFversion.H"
    //#if defined(version30)
     //   #include "turbulentTransportModel.H"
    //        #include "pimpleControl.H"
    //#else
    #include "turbulenceModel.H"
    //#endif
    #if defined(versionv1606plus) || defined(version40)
        #include "fvOptions.H"
    #else
        #include "fvIOoptionList.H"
    #endif
    #include "fixedFluxPressureFvPatchScalarField.H"
    #include "cfdemCloud.H"

    #if defined(anisotropicRotation)
        #include "cfdemCloudRotation.H"
    #endif
    #if defined(superquadrics_flag)
        #include "cfdemCloudRotationSuperquadric.H"
    #endif
    #include "implicitCouple.H"
    #include "clockModel.H"
    #include "smoothingModel.H"
    #include "forceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #if defined(version30)
        pimpleControl pimple(mesh);
        #include "createTimeControls.H"
    #endif
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

IOdictionary couplingParameters
(
   IOobject
   (
       "couplingParameters",
       runTime.time().system(),
       runTime,
       IOobject::MUST_READ,
       IOobject::NO_WRITE
   )
);

bool enableCoupling;
(couplingParameters.lookup("enableCoupling")) >> enableCoupling;

// create cfdemCloud
#include "readGravitationalAcceleration.H"
#include "checkImCoupleM.H"
#if defined(anisotropicRotation)
   cfdemCloudRotation particleCloud(mesh);
#elif defined(superquadrics_flag)
    cfdemCloudRotationSuperquadric particleCloud(mesh);
#else
    cfdemCloud particleCloud(mesh);
#endif
#include "checkModelType.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

//version30 not checked
//       #if defined(version30)
//            #include "readTimeControls.H"
//            #include "CourantNo.H"
//            #include "setDeltaT.H"
//        #else
//            #include "readPIMPLEControls.H"
//            #include "CourantNo.H"
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
//        #endif


//      for debugging
//       Info << "Pimple Controls " << endl;
//       Info << "nOuterCorrectors, nCorrectors = " << pimple.nCorrPIMPLE() << " " << pimple.nCorrPISO() << endl;

if (enableCoupling)
{
        // do particle stuff
        particleCloud.clockM().start(1,"Global");
        particleCloud.clockM().start(2,"Coupling");
        bool hasEvolved = particleCloud.evolve(voidfraction,Us,U);

        if(hasEvolved)
        {
            particleCloud.smoothingM().smoothenAbsolutField(particleCloud.forceM(0).impParticleForces());
        }
    
        Ksl = particleCloud.momCoupleM(particleCloud.registryM().getProperty("implicitCouple_index")).impMomSource();
        Ksl.correctBoundaryConditions();
}
surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);

if (enableCoupling)
{
        phi = voidfractionf*phiByVoidfraction;

        //Force Checks
        #include "forceCheckIm.H"

        #include "solverDebugInfo.H"
        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");
}

        volScalarField voidfractionRho = voidfraction * rho;

        if(particleCloud.solveFlow() || !enableCoupling)
        {

            if (pimple.nCorrPIMPLE() <= 1)
            {
                 #include "rhoEqn.H"
                voidfractionRho = voidfraction * rho;
            }

        // --- Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())   
            {
                // Momentum predictor


// WORK AROUND
// - fvm::Sp(fvc::ddt(voidfraction),U) --> cannot figure out how to mult by rho ...
// So create new field that is voidfraction * rho
// Perhaps there is a simpler solution ... but this works.


                tmp<fvVectorMatrix> tUEqn
                (
                    fvm::ddt(voidfraction*rho,U) - fvm::Sp(fvc::ddt(voidfractionRho),U)
                  + fvm::div(phi,U) - fvm::Sp(fvc::div(phi),U)
// TODO: turbulence->divDevReff(U): this is commented out in cfdemSolverPiso.  Why?
                  + turbulence->divDevRhoReff(U)
                  + particleCloud.divVoidfractionTau(U, voidfraction)
                  ==
                  - fvm::Sp(Ksl,U)
                );

                fvVectorMatrix& UEqn = tUEqn.ref();

                UEqn.relax();


                #if defined(version30)
                    if (pimple.momentumPredictor())
                #else
                    if (momentumPredictor)
                #endif
                {
                    if (modelType=="B" || modelType=="Bfull")
                        solve(UEqn == - fvc::grad(p) + Ksl*Us);
                    else
                        solve(UEqn == - voidfraction*fvc::grad(p) + Ksl*Us);
                }

                #include "hEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {

                    #include "pEqn.H"
 
                }
            }// END --- Pressure-velocity PIMPLE corrector loop

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }

        }// end solveFlow
        else
        {
            Info << "skipping flow solution." << endl;
        }

        rho = thermo.rho();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
