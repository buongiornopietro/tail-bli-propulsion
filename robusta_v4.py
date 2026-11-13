# Regia 11/11/2026 (via di Pietro 16:30 UTC; revisione Fable xhigh punto (d)): forma "garantita" riottimizzata col
# modello corretto (axi.py con e^N di Drela-Giles). Stesse impostazioni di robusta.py (Re_V 1e6 3e6 1e7 3e7, 4 processi).
# Scrive out/robusta_v4.json: il vecchio out/robusta.json non viene sovrascritto. Python di Windows.
import os, sys, time, json
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool
import robusta

if __name__ == "__main__":
    print(time.strftime("%d/%m %H:%M"), "parte robusta v4", flush=True)
    with Pool(4) as pool:
        ris = list(pool.imap_unordered(robusta.caso, [1e6, 3e6, 1e7, 3e7]))
    ris.sort(key=lambda d: d["ReV"])
    json.dump(ris, open("out/robusta_v4.json", "w"), indent=1)
    for d in ris:
        print(f"Re_V {d['ReV']:.0e} ({d['secondi'] / 60:.0f} min)")
        for nome, f in d["forme"].items():
            g = f["guadagno_vel"]
            print(f"  {nome:10s} L/D {f['LsuD']:5.2f} " + " ".join(f"{s} {g[s]:+5.1f}%" for s in robusta.SCENARI)
                  + f" | minimo {min(g.values()):+.1f}%", flush=True)
    print("FATTO", flush=True)
