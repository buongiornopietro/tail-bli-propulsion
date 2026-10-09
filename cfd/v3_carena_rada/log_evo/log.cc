/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  14
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
Build  : 14-7b05503f98a8
Exec   : foamPostProcess -func writeCellCentres -time 800
Date   : Oct 09 2026
Time   : 04:51:26
Host   : "EVO-X2"
PID    : 34803
I/O    : uncollated
Case   : /home/buongiornopietro/cfd/v3_carena_rada
nProcs : 1
sigFpe : Enabling floating point exception trapping (FOAM_SIGFPE).
fileModificationChecking : Monitoring run-time modified files using timeStampMaster (fileModificationSkew 10)
allowSystemOperations : Allowing user-supplied system call operations

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
Create time

Create mesh for time = 800

--> FOAM Warning: 
    From function virtual void Foam::wedgePolyPatch::calcGeometry(Foam::PstreamBuffers&)
    in file meshes/polyMesh/polyPatches/constraint/wedge/wedgePolyPatch.C at line 231
    Wedge patch 'front' may not be sufficiently planar
At patch face #15052 with centre (0.6817054 0.028751736 -0.0012553279) the normal (-4.3539804e-06 -0.043644578 -0.99904712) differs from the average normal (-1.2198791e-09 -0.043619386 -0.99904822) by 6.5480871e-10

--> FOAM Warning: 
    From function virtual void Foam::wedgePolyPatch::calcGeometry(Foam::PstreamBuffers&)
    in file meshes/polyMesh/polyPatches/constraint/wedge/wedgePolyPatch.C at line 231
    Wedge patch 'back' may not be sufficiently planar
At patch face #15052 with centre (0.6817054 0.028751736 0.0012553279) the normal (-4.3539804e-06 -0.043644578 0.99904712) differs from the average normal (-1.2198791e-09 -0.043619386 0.99904822) by 6.5480871e-10


Constructing functionObjects from system/functions
    Selecting writeCellCentres with name writeCellCentres
Time = 800s

Reading fields:

Executing functionObjects
    Writing cell-centre field C to 800
    Writing the x component field of the cell-centres Ccx to 800
    Writing the y component field of the cell-centres Ccy to 800
    Writing the z component field of the cell-centres Ccz to 800

End

