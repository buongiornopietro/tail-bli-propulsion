/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  14
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
Build  : 14-7b05503f98a8
Exec   : foamPostProcess -func writeCellCentres -time 2400
Date   : Oct 09 2026
Time   : 04:55:09
Host   : "EVO-X2"
PID    : 37602
I/O    : uncollated
Case   : /home/buongiornopietro/cfd/v3_conv_rada
nProcs : 1
sigFpe : Enabling floating point exception trapping (FOAM_SIGFPE).
fileModificationChecking : Monitoring run-time modified files using timeStampMaster (fileModificationSkew 10)
allowSystemOperations : Allowing user-supplied system call operations

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Create time

Create mesh for time = 2400

--> FOAM Warning: 
    From function virtual void Foam::wedgePolyPatch::calcGeometry(Foam::PstreamBuffers&)
    in file meshes/polyMesh/polyPatches/constraint/wedge/wedgePolyPatch.C at line 231
    Wedge patch 'front' may not be sufficiently planar
At patch face #3822 with centre (0.12078275 0.060075749 -0.0026229639) the normal (4.0273997e-06 -0.043639428 -0.99904735) differs from the average normal (1.8780841e-09 -0.043619388 -0.99904822) by 4.1856197e-10

--> FOAM Warning: 
    From function virtual void Foam::wedgePolyPatch::calcGeometry(Foam::PstreamBuffers&)
    in file meshes/polyMesh/polyPatches/constraint/wedge/wedgePolyPatch.C at line 231
    Wedge patch 'back' may not be sufficiently planar
At patch face #3822 with centre (0.12078275 0.060075749 0.0026229639) the normal (4.0273997e-06 -0.043639428 0.99904735) differs from the average normal (1.8780841e-09 -0.043619388 0.99904822) by 4.1856197e-10


Constructing functionObjects from system/functions
    Selecting writeCellCentres with name writeCellCentres
Time = 2400s

Reading fields:

Executing functionObjects
    Writing cell-centre field C to 2400
    Writing the x component field of the cell-centres Ccx to 2400
    Writing the y component field of the cell-centres Ccy to 2400
    Writing the z component field of the cell-centres Ccz to 2400

End

