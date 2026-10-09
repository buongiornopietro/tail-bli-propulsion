/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  14
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
Build  : 14-7b05503f98a8
Exec   : foamPostProcess -func writeCellCentres -time 5600
Date   : Oct 09 2026
Time   : 03:45:50
Host   : "EVO-X2"
PID    : 28526
I/O    : uncollated
Case   : /home/buongiornopietro/cfd/v3_carena_F
nProcs : 1
sigFpe : Enabling floating point exception trapping (FOAM_SIGFPE).
fileModificationChecking : Monitoring run-time modified files using timeStampMaster (fileModificationSkew 10)
allowSystemOperations : Allowing user-supplied system call operations

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Create time

Create mesh for time = 5600

--> FOAM Warning: 
    From function virtual void Foam::wedgePolyPatch::calcGeometry(Foam::PstreamBuffers&)
    in file meshes/polyMesh/polyPatches/constraint/wedge/wedgePolyPatch.C at line 231
    Wedge patch 'front' may not be sufficiently planar
At patch face #6000 with centre (0.064663375 0.045900689 -0.0020040673) the normal (6.7004191e-06 -0.043640264 -0.99904731) differs from the average normal (-4.0911235e-10 -0.043619386 -0.99904822) by 4.8159729e-10

--> FOAM Warning: 
    From function virtual void Foam::wedgePolyPatch::calcGeometry(Foam::PstreamBuffers&)
    in file meshes/polyMesh/polyPatches/constraint/wedge/wedgePolyPatch.C at line 231
    Wedge patch 'back' may not be sufficiently planar
At patch face #6000 with centre (0.064663375 0.045900689 0.0020040673) the normal (6.7004191e-06 -0.043640264 0.99904731) differs from the average normal (-4.0911235e-10 -0.043619386 0.99904822) by 4.8159729e-10


Constructing functionObjects from system/functions
    Selecting writeCellCentres with name writeCellCentres
Time = 5600s

Reading fields:

Executing functionObjects
    Writing cell-centre field C to 5600
    Writing the x component field of the cell-centres Ccx to 5600
    Writing the y component field of the cell-centres Ccy to 5600
    Writing the z component field of the cell-centres Ccz to 5600

End

