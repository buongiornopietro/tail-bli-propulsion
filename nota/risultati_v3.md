# Versione 3 della nota: cosa dicono i 9 calcoli (2026-10-09)

Fonti: out/v3_risultati_evo.md (Regia, calcoli sull'EVO-X2, nu_t controllato: tutti turbolenti),
misura_v3.py -> out/misura_v3.json, forceCoeffs dei casi v3_conv_*. Nessuna simulazione nuova.

## 1. Carenature diverse (griglia e Reynolds base)

| carenatura | K | velocità al ventilatore | spinta / D_nudo | potenza vs elica ideale | vs elica di pari area |
|---|---|---|---|---|---|
| D (nota v2) | 2.1% | 1.01 U | 1.080 | +9.2% | −2.1% |
| E: corda 0.070, spessore 4% | 1.8% | 0.97 U | 1.080 | **+4.3%** | **−6.2%** |
| F: uscita meno stretta (R_TE 0.0275) | 2.7% | 1.10 U | 1.016 | +12.1% | +1.2% |
| G: uscita quasi come l'ingresso (R_TE 0.0290) | 3.1% | 1.16 U | 0.991 | +14.9% | +4.3% |

- La carenatura **più corta e sottile (E) è la migliore**: meno attrito della carenatura e ventilatore
  su aria più lenta (0.97 U). Guadagna circa 5 punti rispetto a D, ma resta sopra l'elica ideale (+4%).
- F e G vanno nella direzione **sbagliata**, ed è una lezione utile: aprire l'uscita non rallenta l'aria
  al ventilatore, la accelera (l'uscita stretta fa da "rubinetto"; aprendola passa più aria e più veloce).
  Per far vedere al ventilatore aria lenta serve il contrario: uscita più stretta o ventilatore più grande
  a parità di uscita. Questo corregge l'ipotesi con cui avevo progettato G.
- Tutte le carenature hanno K = 2-3%, ben sotto K*: la perdita d'imbocco non è il problema; il problema è
  l'attrito della carenatura e la velocità al ventilatore (conferma la conclusione della v2).

## 2. Terza griglia (rada, fattore 1/1.41)

| grandezza | rada | base | fitta | ordine p | Richardson | GCI (fitta) |
|---|---|---|---|---|---|---|
| resistenza corpo nudo (×1e-6) | 6.9344 | 6.9574 | 6.9725 | 1.2 | 7.0015 | 0.5% |
| disco libero: P vs elica ideale | −7.3% | −8.2% | −8.1% | (non monotono) | — | ~1 punto |
| carenatura D | divergito | +9.2% | +10.0% | — | — | — |

Il corpo nudo converge in modo regolare (errore stimato 0.5%). Il disco libero oscilla entro 1 punto.
La carenatura sulla griglia rada è divergita (griglia di qualità peggiore attorno alla carenatura).

## 3. Reynolds 10 volte più alto (Re_L 4.95e7)

| caso | Re base | Re x10 |
|---|---|---|
| resistenza corpo nudo (×1e-6) | 6.957 | 4.911 (−29%) |
| disco libero: P vs elica ideale | −8.2% | −6.7% |
| disco libero: P vs elica di pari area | −15.3% | −11.9% |
| carenatura D | +9.2% | NON convergente |

A Reynolds più alto il guadagno del disco libero scende (8.2 -> 6.7%): lo strato limite è relativamente
più sottile e porta meno "aria lenta". È la stessa tendenza del modello semplificato (K* scende con Re).

## Due casi da decidere (servono nuovi calcoli: via di Pietro)
- v3_carena_re7: la spinta oscilla (7.2 -> 6.0 -> 6.7e-6), numeri non usabili. Proposta: ripartire
  dall'ultimo campo con aggiornamento smorzato T_new = (T + D)/2 per 4-6 iterazioni (~30 min sull'EVO).
- v3_carena_rada: divergito. Proposta: rinunciare (corpo nudo e disco hanno già tre griglie; la carenatura
  ne ha due, base e fitta, differenza 0.8 punti) e dichiararlo nei limiti. In alternativa correggere la
  griglia rada attorno alla carenatura e rifarlo.

## Controlli chiesti dal revisore (2026-10-09, deriva_v3.py e deriva_wsl.sh)
Variazione della resistenza nelle ultime 400 iterazioni: rada -0.002%, base (conv_ref) -0.001%,
fitta (conv_fine) -0.150% (comparabile con la differenza base->fitta 0.22%: con il valore fitto estrapolato
6.968-6.970e-6 il GCI scende a ~0.2%; in ogni caso errore di griglia < 0.5%).
Re x10: conv_re7 -0.006%, disco_re7 0.000%: convergiti, i numeri -29% e -6.7% non sono provvisori.
Carenature E/F/G: -0.03..-0.05%; carena_D -0.02%, carena_fine2 -0.13%.
carena_re7: nu_t della seconda corsa coerente con conv_re7 (ramo turbolento coerente), 3 volte piu' basso
della prima corsa a x 0.75: dichiarato nella nota insieme alla non convergenza.
