#!/bin/bash
# Controllo leggero (seriale) delle griglie dei casi v3: copia in /tmp, blockMesh, checkMesh.
# Non lancia foamRun e non tocca ~/cfd. Uso: bash controlla_mesh_v3.sh
. /opt/openfoam14/etc/bashrc
SRC=/mnt/c/Users/PC/Desktop/formaMoto/cfd
T=/tmp/v3check
[ -n "$T" ] && [ "$T" != "/" ] || exit 1
rm -rf "$T"; mkdir -p "$T"
LISTA=${*:-"v3_conv_rada v3_disco_rada v3_carena_rada v3_conv_re7 v3_disco_re7 v3_carena_re7 v3_carena_E v3_carena_F v3_carena_G"}
for c in $LISTA; do
    cp -r "$SRC/$c" "$T/$c"
    cd "$T/$c" || continue
    find . -type f -exec sed -i 's/\r$//' {} +
    if [ -x Allmesh ]; then ./Allmesh > log.mesh 2>&1; else blockMesh > log.mesh 2>&1; fi
    checkMesh > log.checkMesh 2>&1
    celle=$(grep -m1 "cells:" log.checkMesh | awk '{print $2}')
    esito=$(grep -E "Mesh OK|Failed [0-9]+ mesh" log.checkMesh | tail -1)
    quale=$(grep -F '***' log.checkMesh | head -2 | tr -s ' ' | tr '\n' ';')
    nonort=$(grep -m1 "non-orthogonality" log.checkMesh | sed 's/^ *//')
    echo "$c | celle $celle | $esito | $nonort | $quale"
    cd "$T"
done
rm -rf "$T"
