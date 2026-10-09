# Nota 1 v3: i 9 casi RANS calcolati sull'EVO-X2 (09/10/2026, Regia)

Controllo di riproducibilita' (conv_fine, stessi input del portatile): Cd identico al portatile a tutte le cifre
stampate a 400/1200/2400 iterazioni (1.62412484e-05, 1.42034946e-05, 1.39450997e-05).
Coda: 02:42-05:10 CEST, due corsie da 4 processi (coda_v3_evo2.sh, stesse variabili di coda_v3.sh).
Risultati copiati sul portatile in formaMoto/cfd/v3_*/ (risultati/, postProcessing/, autopropulsione.txt,
log_evo/ con tutti i log e nut_x.txt). Nessun nan / FOAM FATAL nei log.

| caso | nu_t (nut_x) | autopropulsione (ultima iterazione) | giudizio |
|---|---|---|---|
| v3_carena_E | cresce con x (3e-6 -> 7-8e-5): turbolento | T 7.5164e-6, D 7.5115e-6 (0.06%) | OK |
| v3_carena_F | turbolento | T 7.0666e-6, D 7.0631e-6 (0.05%) | OK |
| v3_carena_G | turbolento | T 6.8972e-6, D 6.8990e-6 (0.03%) | OK |
| v3_conv_rada | turbolento | (NAUTO 0) | OK |
| v3_conv_re7 | turbolento | (NAUTO 0) | OK |
| v3_disco_rada | turbolento | T 7.156026e-6, D 7.156027e-6 | OK |
| v3_disco_re7 | turbolento | T 5.051304e-6, D 5.051312e-6 | OK |
| v3_carena_re7 | turbolento | T 6.0324e-6 -> D 6.6923e-6: oscilla (7.18 -> 6.03 -> 6.69), NON convergente | da rifare |
| v3_carena_rada | irregolare | divergito: Cd 1.66e-5 a 400 it., 6.67e-2 a 800; log.run.cont si ferma a 918/3200; T, D assurdi (-5.7, -8e3, -1e7) | da rifare |

## v3_carena_re7 rifatto con spinta smorzata (09/10 08:48-09:47, Allrun_smorzato NAUTO 6, IT3 800, W 0.5)
| iter. | T | D | (T-D)/D |
|---|---|---|---|
| 1 | 6.362335e-6 | 6.544609e-6 | -2.8% |
| 2 | 6.453472e-6 | 6.442549e-6 | +0.2% |
| 3 | 6.448011e-6 | 6.347704e-6 | +1.6% |
| 4 | 6.397858e-6 | 6.252348e-6 | +2.3% |
| 5 | 6.325103e-6 | 6.165404e-6 | +2.6% |
| 6 | 6.245254e-6 | 6.074662e-6 | +2.8% |
**NON convergente** (criterio < 0.5% sulle ultime due). Non e' piu' un'oscillazione: D scende di ~1.4% ogni
800 iterazioni senza fermarsi, cioe' e' il flusso stesso (a spinta quasi fissa) che non e' assestato.
nu_t (nut_x): 2.8e-6 a x 0.05 ... 3.3e-5 a x 0.75, 7.6e-5 a x 0.85 (carenatura): simile a v3_conv_re7
(2.2e-6 ... 3.5e-5), quindi non sembra un laminare spurio, ma piu' basso della prima corsa (9.97e-5 a x 0.75).
Risultati copiati in formaMoto/cfd/v3_carena_re7 (vecchia autopropulsione in autopropulsione_std.txt).
Da decidere (thread della nota): far assestare prima il flusso a spinta fissa (molte piu' iterazioni, controllando
che D si fermi) e solo dopo l'autopropulsione; oppure dichiarare re7 non convergente.

## I due casi da sistemare (niente rilanciato)
- **v3_carena_rada** (griglia rada): diverge tra 400 e 800 iterazioni ("bounding k" fino a 539). checkMesh:
  11 facce molto non ortogonali (> 70 gradi), skewness 2.2 (le griglie E/F/G non le hanno da verificare).
  Probabile causa: qualita' della griglia rada attorno alla carenatura (vedi trappola "griglia al naso" e
  "impostazioni carenatura"). Va corretta la griglia (genera_v3.py) prima di rilanciare.
- **v3_carena_re7** (Re x10): la stima della spinta di autopropulsione oscilla del 10% fra le 3 iterazioni.
  Servono piu' iterazioni o un aggiornamento della spinta smorzato: decisione di metodo (thread della nota).
