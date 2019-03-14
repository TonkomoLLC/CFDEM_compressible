rho = thermo.rho();

volScalarField rAU = 1.0/UEqn.A();
surfaceScalarField rAUf("(1|A(U))", fvc::interpolate(rho*rAU));
volScalarField rAUvoidfraction("(voidfraction2|A(U))",rho*rAU*voidfraction);
surfaceScalarField rAUfvoidfraction("(voidfraction2|A(U)F)", fvc::interpolate(rAUvoidfraction));

U = rAU*UEqn.H();
phi = ( fvc::interpolate(U*rho) & mesh.Sf() )

       + rAUfvoidfraction*fvc::ddtCorr(rho, U, phiByVoidfraction); // FYI,In OF-ext the equivalent is commented out for non-transonic flow


surfaceScalarField phiS(fvc::interpolate(Us) & mesh.Sf());
phi += rAUf*(fvc::interpolate(Ksl) * phiS);

if (modelType=="A")
    rAUvoidfraction = volScalarField("(voidfraction2|A(U))",rho*rAU*voidfraction*voidfraction);
    
// Update the fixedFluxPressure BCs to ensure flux consistency
#include "compressibleFixedFluxPressureHandling.H" 


if (nCorr <= 1)
{
    tUEqn.clear();
}

// transonic part removed from code


    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(voidfractionf*phi)
          - fvm::laplacian(rAUvoidfraction, p)
          ==
            rho*particleCloud.ddtVoidfraction()        
        );

        if
        (
            oCorr == nOuterCorr-1
         && corr == nCorr-1
         && nonOrth == nNonOrthCorr
        )
        {
            pEqn.solve(mesh.solver("pFinal"));
        }
        else
        {
            pEqn.solve();
        }

        if (nonOrth == nNonOrthCorr)
        {
            phiByVoidfraction = phi + pEqn.flux()/voidfractionf;
        }
    }

phi = voidfractionf*phiByVoidfraction;
#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

//if (oCorr != nOuterCorr-1)
{
    // Explicitly relax pressure for momentum corrector
    p.relax();

    rho = thermo.rho();
    rho.relax();
    Info<< "rho max/min : " << max(rho).value()
        << " " << min(rho).value() << endl;
}

if (modelType=="B" || modelType=="Bfull")
    U -= rAU*fvc::grad(p) - Ksl*Us*rAU;
else
    U -= voidfraction*rAU*fvc::grad(p) - Ksl*Us*rAU;
U.correctBoundaryConditions();
fvOptions.correct(U);

DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

//bound(p, pMin);   //ED EDIT

if (pressureControl.limit(p))   // ED EDIT: Think this is the way pressure is limited in OF5
{
    p.correctBoundaryConditions();
    rho = thermo.rho();
}

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
/*
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(psi*p))
        /fvc::domainIntegrate(psi);
}
*/
