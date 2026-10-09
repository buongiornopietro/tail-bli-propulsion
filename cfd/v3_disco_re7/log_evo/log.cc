/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  14
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
Build  : 14-7b05503f98a8
Exec   : foamPostProcess -func writeCellCentres -time 4800
Date   : Oct 09 2026
Time   : 05:10:18
Host   : "EVO-X2"
PID    : 39886
I/O    : uncollated
Case   : /home/buongiornopietro/cfd/v3_disco_re7
nProcs : 1
sigFpe : Enabling floating point exception trapping (FOAM_SIGFPE).
fileModificationChecking : Monitoring run-time modified files using timeStampMaster (fileModificationSkew 10)
allowSystemOperations : Allowing user-supplied system call operations

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Create time

Create mesh for time = 4800

--> FOAM Warning: 
    From function virtual void Foam::wedgePolyPatch::calcGeometry(Foam::PstreamBuffers&)
    in file meshes/polyMesh/polyPatches/constraint/wedge/wedgePolyPatch.C at line 231
    Wedge patch 'front' may not be sufficiently planar
At patch face #15070 with centre (0.42947434 0.069884628 -0.0030512287) the normal (1.9250266e-05 -0.04344442 -0.99905585) differs from the average normal (-4.5559499e-09 -0.043619384 -0.99904822) by 3.1041116e-08

--> FOAM Warning: 
    From function virtual void Foam::wedgePolyPatch::calcGeometry(Foam::PstreamBuffers&)
    in file meshes/polyMesh/polyPatches/constraint/wedge/wedgePolyPatch.C at line 231
    Wedge patch 'back' may not be sufficiently planar
At patch face #15070 with centre (0.42947434 0.069884628 0.0030512287) the normal (1.9250266e-05 -0.04344442 0.99905585) differs from the average normal (-4.5559499e-09 -0.043619384 0.99904822) by 3.1041116e-08


Constructing functionObjects from system/functions
    Selecting writeCellCentres with name writeCellCentres
Time = 4800s

Reading fields:

Executing functionObjects
    Writing cell-centre field C to 4800
    Writing the x component field of the cell-centres Ccx to 4800
    Writing the y component field of the cell-centres Ccy to 4800
    Writing the z component field of the cell-centres Ccz to 4800

End

