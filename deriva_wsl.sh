#!/bin/bash
# Solo lettura: ultime righe dei forceCoeffs dei casi base e fitti in ~/cfd (deriva iterativa, v3).
for c in conv_ref conv_fine disco_A disco_fine carena_D carena_fine2; do
    for f in $(ls -d ~/cfd/$c/postProcessing/forze/* 2>/dev/null | sort -t/ -k7 -n); do
        grep -v '^#' "$f/forceCoeffs.dat"
    done | awk -v c="$c" '{n++; it[n]=$1; D[n]=0.5*$3}
        END{ if (n < 2) { print c, "nessun dato"; exit }
             k=n; for (j=n; j>0; j--) if (it[n]-it[j] >= 400) { k=j; break }
             printf "%-14s it %d  D %.5e  ultime %d it: %+.3f%%\n", c, it[n], D[n], it[n]-it[k], 100*(D[n]-D[k])/D[n] }'
done
