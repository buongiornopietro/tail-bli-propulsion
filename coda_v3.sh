#!/bin/bash
# Calcoli della versione 3 della nota (prima nota: corpo + propulsore in coda).
# Casi generati da genera_v3.py in formaMoto/cfd/v3_*; ognuno viene copiato in ~/cfd,
# calcolato con Allrun_D (4 processi mpirun) ed esportato con esporta.sh.
# Ordine per utilità: carenature diverse, poi griglia rada, poi secondo Reynolds.
# Da lanciare SOLO quando non gira nessun altro calcolo:
#   wsl -d OpenFOAM14D -- bash /mnt/c/Users/PC/Desktop/formaMoto/coda_v3.sh
# Riprende da dove era rimasta: salta i casi già segnati in ~/cfd/v3_fatti.txt.
W=/mnt/c/Users/PC/Desktop/formaMoto
pgrep -x mpirun > /dev/null && { echo "c'e' gia' un calcolo in corso: esco"; exit 1; }
cd ~/cfd || exit 1
[ -x ./Allrun_D ] || { echo "manca ~/cfd/Allrun_D"; exit 1; }
touch v3_fatti.txt

esegui() {   # $1 caso, poi variabili per Allrun_D
    local caso=$1; shift
    [ -n "$caso" ] || return 1
    grep -qx "$caso" v3_fatti.txt && { echo "$caso gia' fatto"; return 0; }
    echo "$(date '+%d/%m %H:%M') parte $caso"
    env "$@" ./Allrun_D "$caso" "$W/cfd/$caso" > "$caso.log" 2>&1
    bash "$W/esporta.sh" "$caso" "$W" >> "$caso.log" 2>&1
    echo "$caso" >> v3_fatti.txt
    echo "$(date '+%d/%m %H:%M') finito $caso"
}

# 3) carenature diverse (griglia e Re base)
esegui v3_carena_E ITC=3200 IT3=800 NAUTO=3
esegui v3_carena_F ITC=3200 IT3=800 NAUTO=3
esegui v3_carena_G ITC=3200 IT3=800 NAUTO=3
# 1) terza griglia (rada)
esegui v3_conv_rada ITC=2400 NAUTO=0
esegui v3_disco_rada ITC=2400 IT3=800 NAUTO=3
esegui v3_carena_rada ITC=3200 IT3=800 NAUTO=3
# 2) secondo Reynolds (Re_L x10): OMEGA=15 per avere nu_t/nu = 5 all'ingresso come nei casi base
esegui v3_conv_re7 OMEGA=15 ITC=2400 NAUTO=0
esegui v3_disco_re7 OMEGA=15 ITC=2400 IT3=800 NAUTO=3
esegui v3_carena_re7 OMEGA=15 ITC=3200 IT3=800 NAUTO=3
echo "$(date '+%d/%m %H:%M') coda v3 finita"
