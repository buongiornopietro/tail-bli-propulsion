"""Genera un caso OpenFOAM 14 assialsimmetrico (cuneo) per un corpo di rivoluzione.

Il corpo è scalato a lunghezza L = 1, velocità U = 1, quindi nu = 1/Re_L.
Griglia: gmsh 2D nel piano meridiano (x, r >= 0) con strati prismatici alla
parete, estrusa di 1 cella in z, convertita con gmshToFoam e poi ruotata a
cuneo con extrudeMesh (modello wedge).

uso: python cfd_caso.py <nome_caso> <file_profilo.json> <chiave> <Re_L> [y+ primo strato]
"""
import json, os, sys, textwrap
import numpy as np

HDR = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Version:  14
    \\\\  /    A nd           |
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       {cls};
    location    "{loc}";
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""


def scrivi(caso, rel, cls, testo):
    p = os.path.join(caso, rel)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    loc, obj = os.path.split(rel)
    open(p, "w", newline="\n").write(HDR.format(cls=cls, loc=loc, obj=obj) + textwrap.dedent(testo))


def geo(x, r, Re_L, yplus=1.0, dom=(-4.0, 8.0, 4.0)):
    """Script gmsh: corpo (x, r) con L = 1, dominio [x_in, x_out] x [0, R]."""
    xin, xout, R = dom
    Cf = 0.074 * Re_L ** -0.2
    h0 = yplus / (Re_L * np.sqrt(Cf / 2))           # primo strato per y+ ~ yplus
    delta = 0.37 * Re_L ** -0.2 * 1.3                # spessore strato limite in coda (+30%)
    # punti del profilo (senza estremi, che coincidono con i punti sull'asse)
    L = []
    L.append(f"// corpo: Re_L = {Re_L:.3e}, primo strato h0 = {h0:.3e}, strato limite ~{delta:.3e}")
    L.append("lc_far = 1.0; lc_body = 0.008;")
    L.append(f"Point(1) = {{{xin}, 0, 0, lc_far}};")
    L.append("Point(2) = {0, 0, 0, lc_body};")          # naso
    L.append("Point(3) = {1, 0, 0, lc_body};")          # coda
    L.append(f"Point(4) = {{{xout}, 0, 0, lc_far}};")
    L.append(f"Point(5) = {{{xout}, {R}, 0, lc_far}};")
    L.append(f"Point(6) = {{{xin}, {R}, 0, lc_far}};")
    ids = []
    for k, (xx, rr) in enumerate(zip(x[1:-1], r[1:-1])):
        L.append(f"Point({100 + k}) = {{{xx:.8f}, {rr:.8f}, 0, lc_body}};")
        ids.append(100 + k)
    L.append(f"Spline(2) = {{2, {', '.join(map(str, ids))}, 3}};")   # corpo
    L.append("Line(1) = {1, 2};")                                   # asse a monte
    L.append("Line(3) = {3, 4};")                                   # asse in scia
    L.append("Line(4) = {4, 5};")                                   # uscita
    L.append("Line(5) = {5, 6};")                                   # bordo esterno
    L.append("Line(6) = {6, 1};")                                   # ingresso
    L.append("Curve Loop(1) = {1, 2, 3, 4, 5, 6};")
    L.append("Plane Surface(1) = {1};")
    # strato limite strutturato sul corpo
    L.append("Field[1] = BoundaryLayer;")
    L.append("Field[1].CurvesList = {2};")
    L.append(f"Field[1].Size = {h0:.3e};")
    L.append("Field[1].Ratio = 1.15;")
    L.append(f"Field[1].Thickness = {delta:.4f};")
    L.append("Field[1].Quads = 1;")
    L.append("BoundaryLayer Field = 1;")
    # infittimento della scia e attorno al corpo
    L.append("Field[2] = Box;")
    L.append("Field[2].VIn = 0.04; Field[2].VOut = lc_far;")
    L.append("Field[2].XMin = -0.2; Field[2].XMax = 3.0; Field[2].YMin = 0; Field[2].YMax = 0.35;")
    L.append("Field[2].Thickness = 1.0;")
    L.append("Field[3] = Min; Field[3].FieldsList = {2};")
    L.append("Background Field = 3;")
    L.append("Mesh.MeshSizeExtendFromBoundary = 0;")
    L.append("Mesh.Algorithm = 6;")
    # solo quadrilateri: sull'asse servono quad perché la rotazione a cuneo dia
    # celle prismatiche regolari (i triangoli con un lato sull'asse degenerano)
    L.append("Recombine Surface{1};")
    L.append("Mesh.RecombineAll = 1;")
    L.append("Mesh.RecombinationAlgorithm = 1;")
    L.append("Mesh.SubdivisionAlgorithm = 1;")
    # 1 strato in z -> volume, poi gruppi fisici per gmshToFoam
    L.append("out[] = Extrude {0, 0, 0.01} { Surface{1}; Layers{1}; Recombine; };")
    L.append('Physical Surface("front") = {1};')
    L.append('Physical Surface("back") = {out[0]};')
    L.append('Physical Volume("fluido") = {out[1]};')
    L.append('Physical Surface("asse") = {out[2], out[4]};')
    L.append('Physical Surface("corpo") = {out[3]};')
    L.append('Physical Surface("uscita") = {out[5]};')
    L.append('Physical Surface("esterno") = {out[6]};')
    L.append('Physical Surface("ingresso") = {out[7]};')
    L.append("Mesh.MshFileVersion = 2.2;")
    return "\n".join(L) + "\n", h0


def caso(nome, x, r, Re_L, yplus=1.0, mesher="blockMesh", modello="kOmegaSST", Tu=0.001):
    os.makedirs(nome, exist_ok=True)
    if mesher == "gmsh":
        g, h0 = geo(x, r, Re_L, yplus)
        open(os.path.join(nome, "corpo.geo"), "w", newline="\n").write(g)
    else:
        import cfd_blockmesh
        bm, h0, _ = cfd_blockmesh.blockMeshDict(x, r, Re_L, yplus)
        scrivi(nome, "system/blockMeshDict", "dictionary", bm)
    nu = 1.0 / Re_L
    I = Tu; k0 = 1.5 * I**2; om0 = 10.0              # turbolenza esterna bassa
    if modello == "kOmegaSSTLM":
        om0 = k0 * Re_L                              # nu_t/nu = 1 nel flusso esterno
        tu = 100 * Tu                                # correlazione di Langtry-Menter (Tu in %)
        ReTt = 1173.51 - 589.428 * tu + 0.2196 / tu**2 if tu <= 1.3 else 331.5 * (tu - 0.5658) ** -0.671
    scrivi(nome, "system/controlDict", "dictionary", f"""
    solver          incompressibleFluid;
    startFrom       latestTime;
    startTime       0;
    stopAt          endTime;
    endTime         4000;
    deltaT          1;
    writeControl    timeStep;
    writeInterval   500;
    purgeWrite      2;
    writeFormat     binary;
    writePrecision  8;
    timeFormat      general;
    runTimeModifiable true;

    functions
    {{
        forze
        {{
            type            forceCoeffs;
            libs            ("libforces.so");
            writeControl    timeStep;
            writeInterval   10;
            patches         (corpo);
            rho             rhoInf;
            rhoInf          1;
            CofR            (0 0 0);
            liftDir         (0 1 0);
            dragDir         (1 0 0);
            pitchAxis       (0 0 1);
            magUInf         1;
            lRef            1;
            Aref            1;
        }}
        yPlus {{ type yPlus; libs ("libfieldFunctionObjects.so"); writeControl writeTime; }}
        wallShearStress {{ type wallShearStress; libs ("libfieldFunctionObjects.so"); patches (corpo); writeControl writeTime; }}
    }}
    """)
    scrivi(nome, "system/fvSchemes", "dictionary", """
    ddtSchemes      { default steadyState; }
    gradSchemes     { default Gauss linear; grad(U) cellLimited Gauss linear 1; }
    divSchemes
    {
        default         none;
        div(phi,U)      bounded Gauss linearUpwind grad(U);
        div(phi,k)      bounded Gauss limitedLinear 1;
        div(phi,omega)  bounded Gauss limitedLinear 1;
        div(phi,ReThetat) bounded Gauss limitedLinear 1;
        div(phi,gammaInt) bounded Gauss limitedLinear 1;
        div((nuEff*dev2(T(grad(U))))) Gauss linear;
    }
    laplacianSchemes { default Gauss linear corrected; }
    interpolationSchemes { default linear; }
    snGradSchemes   { default corrected; }
    wallDist        { method meshWave; }
    """)
    scrivi(nome, "system/fvSolution", "dictionary", """
    solvers
    {
        p { solver GAMG; smoother GaussSeidel; tolerance 1e-8; relTol 0.05; }
        "(U|k|omega|ReThetat|gammaInt)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-9; relTol 0.1; }
    }
    SIMPLE
    {
        nNonOrthogonalCorrectors 1;
        consistent      yes;
        residualControl { p 1e-6; U 1e-7; "(k|omega)" 1e-7; }
    }
    relaxationFactors
    {
        equations { U 0.9; ".*" 0.7; }
        fields    { p 0.9; }
    }
    """)
    scrivi(nome, "system/extrudeMeshDict", "dictionary", """
    constructFrom patch;
    sourceCase "$FOAM_CASE";
    sourcePatches (front);
    exposedPatchName back;
    extrudeModel
    {
        type            wedge;
        axisPt          (0 0 0);
        axis            (-1 0 0);
        angle           5;
    }
    flipNormals false;
    mergeFaces false;
    """)
    scrivi(nome, "system/collapseDict", "dictionary", """
    controlMeshQuality off;
    collapseEdgesCoeffs { minimumEdgeLength 1e-10; maximumMergeAngle 180; }
    """)
    scrivi(nome, "constant/physicalProperties", "dictionary", f"""
    viscosityModel  constant;
    nu              [0 2 -1 0 0 0 0] {nu:.6e};
    """)
    scrivi(nome, "constant/momentumTransport", "dictionary", f"""
    simulationType RAS;
    RAS {{ model {modello}; turbulence on; printCoeffs off; }}
    """)
    bc = lambda body, inl, out, ext: f"""
        corpo    {{ {body} }}
        {"ingresso {{ " + inl + " }}" if mesher == "gmsh" else ""}
        uscita   {{ {out} }}
        esterno  {{ {ext} }}
        front    {{ type wedge; }}
        back     {{ type wedge; }}
        asse     {{ type empty; }}
    """
    scrivi(nome, "0/U", "volVectorField", """
    dimensions [0 1 -1 0 0 0 0];
    internalField uniform (1 0 0);
    boundaryField
    {""" + bc("type noSlip;", "type fixedValue; value uniform (1 0 0);",
              "type zeroGradient;", "type freestreamVelocity; freestreamValue uniform (1 0 0);") + "}\n")
    scrivi(nome, "0/p", "volScalarField", """
    dimensions [0 2 -2 0 0 0 0];
    internalField uniform 0;
    boundaryField
    {""" + bc("type zeroGradient;", "type zeroGradient;",
              "type fixedValue; value uniform 0;", "type freestreamPressure; freestreamValue uniform 0;") + "}\n")
    scrivi(nome, "0/k", "volScalarField", f"""
    dimensions [0 2 -2 0 0 0 0];
    internalField uniform {k0:.3e};
    boundaryField
    {{""" + bc("type fixedValue; value uniform 1e-14;", f"type fixedValue; value uniform {k0:.3e};",
               "type zeroGradient;", f"type inletOutlet; inletValue uniform {k0:.3e}; value uniform {k0:.3e};") + "}\n")
    scrivi(nome, "0/omega", "volScalarField", f"""
    dimensions [0 0 -1 0 0 0 0];
    internalField uniform {om0};
    boundaryField
    {{""" + bc("type omegaWallFunction; value uniform 1;", f"type fixedValue; value uniform {om0};",
               "type zeroGradient;", f"type inletOutlet; inletValue uniform {om0}; value uniform {om0};") + "}\n")
    scrivi(nome, "0/nut", "volScalarField", """
    dimensions [0 2 -1 0 0 0 0];
    internalField uniform 0;
    boundaryField
    {""" + bc("type nutLowReWallFunction; value uniform 0;", "type calculated; value uniform 0;",
              "type calculated; value uniform 0;", "type calculated; value uniform 0;") + "}\n")
    if modello == "kOmegaSSTLM":
        scrivi(nome, "0/ReThetat", "volScalarField", f"""
        dimensions [0 0 0 0 0 0 0];
        internalField uniform {ReTt:.2f};
        boundaryField
        {{""" + bc("type zeroGradient;", f"type fixedValue; value uniform {ReTt:.2f};",
                   "type zeroGradient;", f"type inletOutlet; inletValue uniform {ReTt:.2f}; value uniform {ReTt:.2f};") + "}\n")
        scrivi(nome, "0/gammaInt", "volScalarField", """
        dimensions [0 0 0 0 0 0 0];
        internalField uniform 1;
        boundaryField
        {""" + bc("type zeroGradient;", "type fixedValue; value uniform 1;",
                  "type zeroGradient;", "type inletOutlet; inletValue uniform 1; value uniform 1;") + "}\n")
    if mesher != "gmsh":
        open(os.path.join(nome, "Allmesh"), "w", newline="\n").write(textwrap.dedent("""\
            #!/bin/bash
            set -e
            cd "$(dirname "$0")"
            . /opt/openfoam14/etc/bashrc
            blockMesh > log.blockMesh 2>&1
            checkMesh > log.checkMesh 2>&1 || true
            tail -25 log.checkMesh
            """))
        return h0
    open(os.path.join(nome, "Allmesh"), "w", newline="\n").write(textwrap.dedent("""\
        #!/bin/bash
        set -e
        cd "$(dirname "$0")"
        . /opt/openfoam14/etc/bashrc
        gmsh -3 corpo.geo -o corpo.msh -format msh2 > log.gmsh 2>&1
        gmshToFoam corpo.msh > log.gmshToFoam 2>&1
        # tipi di patch dopo la conversione
        foamDictionary constant/polyMesh/boundary -entry entry0/front/type -set patch > /dev/null
        extrudeMesh > log.extrudeMesh 2>&1
        collapseEdges -overwrite > log.collapseEdges 2>&1 || true
        foamDictionary constant/polyMesh/boundary -entry entry0/front/type -set wedge > /dev/null
        foamDictionary constant/polyMesh/boundary -entry entry0/back/type -set wedge > /dev/null
        foamDictionary constant/polyMesh/boundary -entry entry0/corpo/type -set wall > /dev/null
        foamDictionary constant/polyMesh/boundary -entry entry0/asse/type -set empty > /dev/null 2>&1 || true
        checkMesh > log.checkMesh 2>&1 || true
        tail -25 log.checkMesh
        """))
    return h0


if __name__ == "__main__":
    nome, fprof, chiave, Re_L = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
    yplus = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
    modello = sys.argv[6] if len(sys.argv) > 6 else "kOmegaSST"
    Tu = float(sys.argv[7]) if len(sys.argv) > 7 else 0.001
    d = json.load(open(fprof))
    for k in chiave.split("/"):
        d = d[k]
    x, r = np.array(d["x"]), np.array(d["r"])
    L = x[-1]
    h0 = caso(nome, x / L, r / L, Re_L, yplus, modello=modello, Tu=Tu)
    print(f"caso {nome}: L/D = {L / (2 * r.max()):.2f}, Re_L = {Re_L:.3e}, primo strato {h0:.2e} L")
