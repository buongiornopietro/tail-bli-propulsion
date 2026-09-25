"""English versions of the figures (the plotting scripts write Italian labels).

Importing this module patches Figure.savefig: when the environment variable FIG_EN=1 is
set, every text in the figure (titles, axis labels, legends, annotations, tick labels) is
translated with the TRADUZIONI table and the file is saved as <name>_en.png.
With FIG_EN unset the scripts behave exactly as before.
FIG_EN=dump prints every text, to extend the table.
"""
import os
import re
import matplotlib
from matplotlib.figure import Figure
from matplotlib.text import Text

# substrings, applied in order (longer phrases first)
TRADUZIONI = [
    # plot.py (2D lattice Boltzmann)
    ("Scia (vorticità) dietro ogni forma, flusso da sinistra", "Wake (vorticity) behind each shape, flow from the left"),
    ("coefficiente di resistenza (più basso = più veloce)", "drag coefficient (lower = faster)"),
    ("Simulazione Lattice Boltzmann 2D", "2D lattice-Boltzmann simulation"),
    ("resistenza a pari altezza (sezione frontale)", "drag at equal height (frontal area)"),
    ("resistenza a pari area (volume)", "drag at equal area (volume)"),
    ("Resistenza a pari ALTEZZA frontale", "Drag at equal frontal HEIGHT"),
    ("Resistenza a pari AREA (volume in 2D)", "Drag at equal AREA (2D volume)"),
    ("Spirale aurea (ruotata 180°)", "Golden spiral (rotated 180°)"),
    ("Spirale aurea", "Golden spiral"),
    ("Rettangolo aureo", "Golden rectangle"),
    ("Ellisse aurea (1:1,618)", "Golden ellipse (1:1.618)"),
    ("Goccia aurea (L/H 1,618)", "Golden teardrop (L/H 1.618)"),
    ("Goccia L/H 4,5", "Teardrop L/H 4.5"),
    ("Goccia L/H 3", "Teardrop L/H 3"),
    ("Cerchio", "Circle"),
    ("rettangolo aureo", "golden rectangle"),
    ("ellisse aurea", "golden ellipse"),
    ("goccia aurea", "golden teardrop"),
    ("spirale aurea 180°", "golden spiral 180°"),
    ("spirale aurea", "golden spiral"),
    ("cerchio", "circle"),
    # plot_axi.py
    ("Migliore resistenza ottenibile a ogni snellezza (★ = ottimo globale)",
     "Best drag achievable at each fineness ratio (★ = global optimum)"),
    ("resistenza in più rispetto all'ottimo [%]", "drag above the optimum [%]"),
    ("Resistenza in più rispetto all'ottimo\n(* flusso staccato: valore stimato)",
     "Drag above the optimum\n(* separated flow: estimated value)"),
    ("Forme ottime a PARI VOLUME (sezione del solido di rivoluzione, flusso da sinistra)",
     "Optimal shapes at EQUAL VOLUME (meridian section, flow from the left)"),
    ("Pari SEZIONE FRONTALE (es. veicolo di altezza data)", "Equal FRONTAL AREA (e.g. a vehicle of given height)"),
    ("Corpi di rivoluzione 3D a Reynolds reali — pannelli + strato limite integrale",
     "Axisymmetric bodies at realistic Reynolds numbers: panels + integral boundary layer"),
    ("ellissoide aureo (stesso volume)", "golden ellipsoid (same volume)"),
    ("ellissoide aureo", "golden ellipsoid"),
    ("spessore max a", "max thickness at"),
    ("snellezza L/D", "fineness ratio L/D"),
    ("Volume fisso, strato limite turbolento", "Fixed volume, turbulent boundary layer"),
    ("Volume fisso, transizione naturale", "Fixed volume, natural transition"),
    ("Sezione frontale fissa, strato limite turbolento", "Fixed frontal area, turbulent boundary layer"),
    ("Sezione frontale fissa", "Fixed frontal area"),
    ("ottimo D", "optimum D"),
    ("scen.", "case"),
    ("goccia\naurea", "golden\nteardrop"),
    ("ellissoide\naureo", "golden\nellipsoid"),
    ("goccia\nL/D", "teardrop\nL/D"),
    ("sfera", "sphere"),
    ("φ = 1,618", "φ = 1.618"),
    # plot_mappa_finale.py
    ("Soglia di convenienza del propulsore in coda\n", "Break-even duct loss of the tail propulsor\n"),
    ("sotto la curva conviene il propulsore in coda, sopra l'elica convenzionale",
     "below the curve the tail propulsor pays, above it the conventional propeller"),
    ("turbolento (superficie ruvida o aria agitata)", "turbulent (rough surface or turbulent free stream)"),
    ("transizione naturale: fascia di incertezza del criterio", "natural transition: spread between criteria"),
    ("transizione naturale, criterio e^N (stima preferita)", "natural transition, e^N criterion"),
    ("transizione naturale, criterio di Michel (ottimistico)", "natural transition, Michel criterion (optimistic)"),
    ("numero di Reynolds volumetrico  Re_V", "volumetric Reynolds number  Re_V"),
    ("perdita massima ammissibile nel condotto  K*  [% di ½ρU²]", "maximum tolerable duct loss  K*  [% of ½ρU²]"),
    # plot_cfd.py
    ("Velocità e linee di corrente nel piano meridiano (asse in basso)",
     "Velocity and streamlines in the meridian plane (axis at the bottom)"),
    ("Zoom: strato limite in coda e scia", "Close-up: tail boundary layer and wake"),
    ("Viscosità turbolenta: strato limite e scia", "Eddy viscosity: boundary layer and wake"),
    ("pannelli (flusso potenziale)", "panels (potential flow)"),
    ("Pressione sulla superficie", "Surface pressure"),
    ("y⁺ primo strato", "first-cell y⁺"),
    ("Risoluzione di parete", "Wall resolution"),
    # plot_carena.py
    ("Disco attuatore libero (nessuna carenatura)", "Free actuator disc (no cowl)"),
    ("Propulsore intubato: carenatura con labbro arrotondato", "Ducted propulsor: cowl with rounded lip"),
    ("velocità assiale u / U", "axial velocity u / U"),
    ("ventilatore", "fan"),
    ("disco", "disc"),
    # plot_robusta.py
    ("Forma ottima per la transizione naturale: nominale vs garantita (stessa potenza)",
     "Optimal shape with natural transition: nominal vs guaranteed (same power)"),
    ("velocità vs goccia [%]", "speed vs teardrop [%]"),
    ("turbolenta (", "turbulent optimum ("),
    ("quieto", "low Tu"),
    ("mosso", "medium Tu"),
    ("agitato", "high Tu"),
    ("rugoso", "rough"),
    ("michel", "Michel"),
    ("nominale", "nominal"),
    ("garantita", "guaranteed"),
    ("robusta", "robust"),
    # generic words last, after the phrases that contain them
    ("turbolento", "turbulent"),
    ("ottimo", "optimum"),
    ("goccia", "teardrop"),
]


# tick labels that are whole words (robusta.py transition scenarios)
TICK = {"quieto": "Tu 0.05%", "nominale": "Tu 0.1%", "mosso": "Tu 0.3%", "agitato": "Tu 1%",
        "michel": "Michel", "rugoso": "rough"}


def traduci(s):
    for it, en in TRADUZIONI:
        s = s.replace(it, en)
    return re.sub(r"(\d),(\d)", r"\1.\2", s)          # decimal comma -> point


def _etichette_fisse(fig):
    """Tick labels given as fixed strings and table cells are not plain Text children."""
    fig.canvas.draw()                                   # makes the tick labels exist
    for ax in fig.axes:
        for asse in (ax.xaxis, ax.yaxis):
            for minor in (False, True):
                testi = [t.get_text() for t in asse.get_ticklabels(minor=minor)]
                nuovi = [TICK.get(s, traduci(s)) for s in testi]
                if nuovi != testi:
                    asse.set_ticks(asse.get_ticklocs(minor=minor), nuovi, minor=minor)
        for tab in ax.tables:
            for cella in tab.get_celld().values():
                t = cella.get_text()
                t.set_text(traduci(t.get_text()))


_originale = Figure.savefig


def _savefig(self, fname, *a, **k):
    modo = os.environ.get("FIG_EN")
    if modo:
        _etichette_fisse(self)
        testi = [t for t in self.findobj(Text) if t.get_text().strip()]
        if modo == "dump":
            for t in testi:
                print(repr(t.get_text()))
        for t in testi:
            t.set_text(traduci(t.get_text()))
        base, ext = os.path.splitext(str(fname))
        fname = base + "_en" + ext
    return _originale(self, fname, *a, **k)


Figure.savefig = _savefig
