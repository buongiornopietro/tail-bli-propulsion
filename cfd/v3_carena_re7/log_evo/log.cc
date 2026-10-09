/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  14
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
Build  : 14-7b05503f98a8
Exec   : foamPostProcess -func writeCellCentres -time 10400
Date   : Oct 09 2026
Time   : 09:47:24
Host   : "EVO-X2"
PID    : 49897
I/O    : uncollated
Case   : /home/buongiornopietro/cfd/v3_carena_re7
nProcs : 1
sigFpe : Enabling floating point exception trapping (FOAM_SIGFPE).
fileModificationChecking : Monitoring run-time modified files using timeStampMaster (fileModificationSkew 10)
allowSystemOperations : Allowing user-supplied system call operations

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Create time

Create mesh for time = 10400

--> FOAM Warning: 
    From function virtual void Foam::wedgePolyPatch::calcGeometry(Foam::PstreamBuffers&)
    in file meshes/polyMesh/polyPatches/constraint/wedge/wedgePolyPatch.C at line 231
    Wedge patch 'front' may not be sufficiently planar
At patch face #28201 with centre (0.61347691 0.040968261 -0.0017887129) the normal (-3.4852354e-05 -0.043810355 -0.99903986) differs from the average normal (-2.0102689e-10 -0.043619394 -0.99904822) by 3.7750644e-08

--> FOAM Warning: 
    From function virtual void Foam::wedgePolyPatch::calcGeometry(Foam::PstreamBuffers&)
    in file meshes/polyMesh/polyPatches/constraint/wedge/wedgePolyPatch.C at line 231
    Wedge patch 'back' may not be sufficiently planar
At patch face #28201 with centre (0.61347692 0.040968261 0.0017887129) the normal (-3.4852354e-05 -0.043810355 0.99903986) differs from the average normal (-2.0102689e-10 -0.043619394 0.99904822) by 3.7750644e-08


Constructing functionObjects from system/functions
    Selecting writeCellCentres with name writeCellCentres
Time = 10400s

Reading fields:

Executing functionObjects
    Writing cell-centre field C to 10400
    Writing the x component field of the cell-centres Ccx to 10400
    Writing the y component field of the cell-centres Ccy to 10400
    Writing the z component field of the cell-centres Ccz to 10400

End

