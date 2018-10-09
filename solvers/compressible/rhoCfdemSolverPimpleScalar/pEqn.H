rho = thermo.rho();

volScalarField rAU = 1.0/UEqn.A();
surfaceScalarField rAUf("(1|A(U))", fvc::interpolate(rho*rAU));
volScalarField rAUvoidfraction("(voidfraction2|A(U))",rho*rAU*voidfraction);
surfaceScalarField rAUfvoidfraction("(voidfraction2|A(U)F)", fvc::interpolate(rAUvoidfraction));

U = rAU*UEqn.H();
phi = ( fvc::interpolate(U*rho) & mesh.Sf() ) + rAUfvoidfraction*fvc::ddtCorr(rho, U, phiByVoidfraction);
surfaceScalarField phiS(fvc::interpolate(Us) & mesh.Sf());
phi += rAUf*(fvc::interpolate(Ksl) * phiS);

if (modelType=="A")
    rAUvoidfraction = volScalarField("(voidfraction2|A(U))",rho*rAU*voidfraction*voidfraction);
    
// Update the fixedFluxPressure BCs to ensure flux consistency
#include "compressibleFixedFluxPressureHandling.H"    
    
if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

if (pimple.transonic())
// transonic is not functional!
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            (fvc::interpolate(U) & mesh.Sf())
//          + fvc::ddtPhiCorr(rUA, rho, U, phi)
// to use transonic, find the appropriate replacement, for the above commented line, in OpenFoam 5.x
        )
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(voidfractionf*phid, p)
          - fvm::laplacian(rAUvoidfraction, p)
        );

        // Relax the pressure equation to ensure diagonal-dominance
        pEqn.relax();

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi == pEqn.flux();
        }
    }
}
else
{
/*
    phi =
        fvc::interpolate(rho)*
        (
            (fvc::interpolate(U) & mesh.Sf())
        //+ fvc::ddtPhiCorr(rUA, rho, U, phi)
        );
*/


    while (pimple.correctNonOrthogonal())
    {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(voidfractionf*phi)
          - fvm::laplacian(rAUvoidfraction, p)
          + rho*particleCloud.ddtVoidfraction()
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
//            phi += pEqn.flux()/voidfractionf;
            phiByVoidfraction = phi + pEqn.flux()/voidfractionf;
        }
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

DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

//bound(p, pMin);

if (pressureControl.limit(p))
{
    p.correctBoundaryConditions();
    rho = thermo.rho();
}
else if (pimple.SIMPLErho())
{
    rho = thermo.rho();
}

if (thermo.dpdt())
{
    DpDt = fvc::ddt(p);
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
