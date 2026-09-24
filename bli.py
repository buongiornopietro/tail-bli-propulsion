"""Propulsore in coda / corpo di Goldschmied: ingestione dello strato limite (BLI).

Modello v3 (bilancio di energia, cfr. Drela 2009 "Power balance", Smith 1993):
  * fan intubato: presa ad anello in x_s larga Q/(2π r U) (ingresso a ~U) e
    scarico subito a valle, dimensionato per un getto a ~U; nel flusso
    potenziale la funzione di corrente sul corpo sale di Q/2π lungo la presa e
    torna a 0 lungo lo scarico, con rampe lisce (niente singolarità ai bordi);
    pannelli infittiti attorno alla presa;
  * la pompa restituisce all'aria ingerita la pressione totale persa per attrito:
        P_ripristino = φ · ½ρ Ue³ θ* · 2π r   (θ* spessore di energia all'imbocco)
    con φ = min(1, Q/Q_BL) frazione di strato limite ingerita; Q_BL è misurata
    a monte, fuori dall'accelerazione locale verso la presa;
  * la coda a valle del fan è lambita dal getto e resta attaccata: resistenza di
    attrito turbolento di lastra piana alla velocità del getto (fino alla stessa
    stazione di Young del corpo convenzionale) + scia dello strato limite non
    ingerito, compensate da un getto più veloce:
        V_j = U + D_coda / (ρQ),   P_extra = D_coda (V_j + U) / 2
  * perdite: pompa con rendimento η, condotto con perdita K·½ρU² su tutta Q.

Riferimento convenzionale: corpo ottimo con elica in aria libera, P = D·U / η.
Coefficiente di potenza: C_P = P / (½ρU³ V^(2/3)).
"""
import numpy as np
from scipy.interpolate import CubicSpline
import axi

PERDITA = "U"   # "U": K·½ρU²·Q   |   "Ue": K·½ρUe²·Q (velocità locale a monte della fessura)


# ------------------------------------------------ flusso potenziale con aspirazione
def surface_velocity_q(x, r, xs_frac=None, Q=0.0, slot_w=0.012, exit_w=0.02):
    """Come axi.surface_velocity, con un fan intubato: presa larga slot_w·L
    centrata in xs_frac·L (portata Q) e scarico largo exit_w·L subito a valle."""
    xc, rc = (x[1:] + x[:-1]) / 2, (r[1:] + r[:-1]) / 2
    L = x[-1]
    target = np.zeros_like(xc)
    if xs_frac is not None and Q > 0:
        q = Q / (2 * np.pi)
        a, b = (xs_frac - slot_w / 2) * L, (xs_frac + slot_w / 2) * L
        # rampe lisce (smoothstep): aspirazione e getto partono da zero senza spigoli,
        # evitando singolarità logaritmiche della velocità ai bordi della presa
        sm = lambda t: (lambda u: u * u * (3 - 2 * u))(np.clip(t, 0, 1))
        c_ = b + exit_w * L                 # scarico del fan subito a valle della presa
        target = np.where(xc < b, q * sm((xc - a) / (b - a)), q * (1 - sm((xc - b) / (c_ - b))))
    A = influence(x, r)
    g = np.linalg.solve(A, target - 0.5 * rc**2)
    if not (xs_frac is not None and Q > 0):
        return xc, rc, np.abs(g)
    # con portata interna il salto γ non è più la velocità esterna:
    # u_est = (1/r) ∂ψ/∂n appena fuori dalla superficie (differenze 2° ordine)
    ds = np.hypot(np.diff(x), np.diff(r))
    nx, nr = -np.diff(r) / ds, np.diff(x) / ds
    eps = 0.25 * ds
    p1 = psi_at(x, r, g, xc + eps * nx, rc + eps * nr)
    p2 = psi_at(x, r, g, xc + 2 * eps * nx, rc + 2 * eps * nr)
    dpsi = (-3 * target + 4 * p1 - p2) / (2 * eps)
    return xc, rc, np.abs(dpsi) / np.maximum(rc, 1e-9)


def psi_at(x, r, g, px, pr):
    """Funzione di corrente totale nei punti (px, pr)."""
    ds = np.hypot(np.diff(x), np.diff(r))
    t, wg = axi.GL16
    tt = (t + 1) / 2
    xq = x[:-1, None] + tt[None] * np.diff(x)[:, None]       # (N, 16)
    rq = r[:-1, None] + tt[None] * np.diff(r)[:, None]
    G = axi.ring_psi(px[:, None, None], pr[:, None, None], xq[None], rq[None])
    return (G * wg[None, None] / 2).sum(2) @ (g * ds) + 0.5 * pr**2


_cache = {}
def influence(x, r):
    key = (x.tobytes(), r.tobytes())
    if key not in _cache:
        _cache.clear()
        xc, rc = (x[1:] + x[:-1]) / 2, (r[1:] + r[:-1]) / 2
        ds = np.hypot(np.diff(x), np.diff(r))
        N = len(xc); A = np.zeros((N, N))
        t, wg = axi.GL8; t2, w2 = axi.GL16
        for j in range(N):
            tt = (t + 1) / 2
            xq = x[j] + tt * (x[j + 1] - x[j]); rq = r[j] + tt * (r[j + 1] - r[j])
            A[:, j] = (axi.ring_psi(xc[:, None], rc[:, None], xq[None], rq[None]) * wg / 2).sum(1) * ds[j]
            acc = 0.0
            for a0, a1 in ((0.0, 0.5), (0.5, 1.0)):
                tt = a0 + (t2 + 1) / 2 * (a1 - a0)
                xq = x[j] + tt * (x[j + 1] - x[j]); rq = r[j] + tt * (r[j + 1] - r[j])
                acc += (axi.ring_psi(xc[j], rc[j], xq, rq) * w2).sum() * (a1 - a0) / 2
            A[j, j] = acc * ds[j]
        _cache[key] = A
    return _cache[key]


# ------------------------------------------------------------- strato limite
def Hstar(H):
    """Spessore di energia / quantità di moto, da profilo a legge di potenza."""
    n = 2 / np.maximum(H - 1, 1e-3)
    return 2 * (n + 2) / (n + 3)


def march(s, Ue, dUe, R, dR, nu, transition, theta0=None, Rex_max=None):
    """Strato limite da s[0] a s[-1]. Se theta0 è None parte dal punto di
    ristagno (laminare, poi transizione); altrimenti parte turbolento con θ0.
    Restituisce θ, H alla fine, indice di distacco (o None), indice di transizione."""
    Rex_max = axi.REX_MAX if Rex_max is None else Rex_max
    if theta0 is None:
        integ = np.concatenate([[0], np.cumsum(0.5 * (Ue[1:]**5 * R[1:]**2 + Ue[:-1]**5 * R[:-1]**2) * np.diff(s))])
        th = np.sqrt(0.45 * nu * integ / (Ue**6 * R**2) + 1e-30)
        lam = th**2 * dUe / nu
        Reth = Ue * th / nu; Rex = Ue * (s - s[0]) / nu + 1e-9
        if transition == "forzata":
            itr = min(int(0.03 * len(s)), len(s) - 2)
        else:
            itr = axi.indice_transizione(s, Ue, th, lam, nu, Rex_max)
        if itr >= len(s) - 1:                           # tutto laminare
            Hl = float(axi.thwaites_H(lam[-1]))
            return th[-1], Hl, None, itr
        theta = th[itr]
    else:
        itr = 0; theta = theta0
    H = 1.4; H1 = float(axi.H1_of_H(H))
    for i in range(itr, len(s) - 1):
        h = s[i + 1] - s[i]
        def rhs(i_, th_, H1_):
            H_ = float(axi.H_of_H1(H1_))
            Rt = max(Ue[i_] * th_ / nu, 10.0)
            Cf = 0.246 * 10 ** (-0.678 * H_) * Rt ** -0.268
            dth = Cf / 2 - th_ * ((H_ + 2) * dUe[i_] / Ue[i_] + dR[i_] / R[i_])
            F = 0.0306 * max(H1_ - 3.0, 1e-6) ** -0.6169
            dH1 = F / th_ - H1_ * (dUe[i_] / Ue[i_] + dth / th_ + dR[i_] / R[i_])
            return dth, dH1
        a1, b1 = rhs(i, theta, H1)
        a2, b2 = rhs(i + 1, theta + h * a1, H1 + h * b1)
        theta += h * (a1 + a2) / 2; H1 += h * (b1 + b2) / 2
        H = float(axi.H_of_H1(H1))
        if H > 2.4 or not np.isfinite(theta):
            return theta, min(H, 2.4), i + 1, itr
    return theta, H, None, itr


def _splines(x, r, xc, ue, s_lo, s_hi, npts):
    s_nodes = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(r)))])
    sc = (s_nodes[1:] + s_nodes[:-1]) / 2
    return s_nodes, sc


# ------------------------------------------------------------- potenza
EXIT = 0.985        # inizio dell'uscita del getto (frazione di L)
U_PRESA = 1.0       # velocità media d'ingresso nella presa (unità di U)


def power(x, r, ReV, transition, xs_frac, Q, x_end=None, detail=False, npts=2000,
          penal=True, eta_pompa=1.0, K_condotto=0.0):
    """Coefficiente di potenza C_P del corpo autopropulso con BLI in xs_frac.

    La presa è dimensionata perché l'aria entri a velocità U_PRESA·U:
        larghezza = Q / (2π r_s U_PRESA)   (minimo 1,2% di L)
    eta_pompa: rendimento della pompa; K_condotto: perdita di pressione totale
    nel condotto interno, in unità di ½ρU² (o ½ρUe² se PERDITA = "Ue"), su tutta Q.
    Lo scarico del fan è subito a valle della presa, dimensionato per un getto ~U."""
    nu = 1 / ReV
    L = x[-1]
    xs = xs_frac * L
    r_s = float(np.interp(xs, x, r))
    w = max(0.012, Q / (2 * np.pi * max(r_s, 1e-6) * U_PRESA) / L)
    # scarico: getto a velocità ~U attraverso l'anello a valle della presa
    r_b = float(np.interp(xs + w * L / 2, x, r))
    we = max(0.012, Q / (2 * np.pi * max(r_b, 1e-6) * U_PRESA) / L)
    pen = 0.0
    fine = xs_frac + w / 2 + we
    if fine > EXIT:                            # presa + scarico non stanno sul corpo
        pen += 1.0 * (fine - EXIT)
    xc, rc, ue = surface_velocity_q(x, r, xs_frac, Q, slot_w=w, exit_w=we)
    s_nodes = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(r)))])
    sc = (s_nodes[1:] + s_nodes[:-1]) / 2
    Rs = CubicSpline(s_nodes, r); Xs = CubicSpline(s_nodes, x)
    a = (xs_frac - w / 2) * L; b = (xs_frac + w / 2) * L
    # --- avancorpo: dal naso all'inizio della presa
    kf = xc < a
    Sf = CubicSpline(np.concatenate([[0], sc[kf]]), np.concatenate([[0], ue[kf]]))
    s_a = np.interp(a, xc, sc)
    s = np.linspace(0, s_a, npts)[1:]
    Ue = np.maximum(Sf(s), 1e-6); dUe = Sf(s, 1)
    R = np.maximum(Rs(s), 1e-9); dR = Rs(s, 1)
    th, H, sep_f, itr = march(s, Ue, dUe, R, dR, nu, transition)
    out = dict(sep_avancorpo=None if sep_f is None else float(Xs(s[sep_f]) / L),
               x_trans=float(Xs(s[min(itr, len(s) - 1)]) / L), larghezza_presa=float(w))
    if sep_f is not None:           # distacco prima della presa: forma non valida
        pen += 0.2 * np.pi * R[sep_f]**2 * (1 + (xs - Xs(s[sep_f])) / L * 10)
        iend = sep_f
    else:
        iend = len(s) - 1
    Ue_s, r_a = Ue[iend], R[iend]
    # energia persa dallo strato limite fino all'imbocco della presa
    E_tot = 0.5 * Ue_s**3 * Hstar(H) * th * 2 * np.pi * r_a     # ρ = 1
    # portata dello strato limite misurata a monte, fuori dall'accelerazione
    # locale verso la presa (è l'aria che arriverà alla presa)
    x_q = max(0.05, xs_frac - w / 2 - max(0.05, 2 * w)) * L
    k0 = int(np.clip(np.searchsorted(s, np.interp(x_q, xc, sc)), 10, len(s) - 1))
    thq, Hq, _, _ = march(s[:k0], Ue[:k0], dUe[:k0], R[:k0], dR[:k0], nu, transition)
    nq = 2 / max(Hq - 1, 1e-3)
    Q_BL = 2 * np.pi * R[k0 - 1] * Ue[k0 - 1] * thq * ((nq + 1) * (nq + 2) / nq - Hq)
    # ingestione parziale: se Q < Q_BL entra solo la frazione phi dello strato
    # limite; il resto prosegue sulla coda (niente penalità a gradino)
    phi = min(1.0, Q / max(Q_BL, 1e-12))
    E_ing = phi * E_tot
    th_res = (1 - phi) * th                 # deficit residuo, per unità di perimetro
    # --- coda a valle del fan: immersa nel getto, che la mantiene attaccata.
    #     Resistenza = attrito turbolento di lastra piana alla velocità del getto
    #     sulla superficie bagnata fino alla stazione di Young (come il corpo
    #     convenzionale), + scia dello strato limite non ingerito (Young).
    x_end = axi.X_END if x_end is None else x_end
    xa0, xa1 = b + we * L, x_end * L
    D_res = 0.5 * 4 * np.pi * r_a * th_res * Ue_s ** ((H + 5) / 2)
    sep_c = None
    if xa1 > xa0:
        xx = np.linspace(xa0, xa1, 200)
        S_aft = np.trapezoid(2 * np.pi * Rs(np.interp(xx, x, s_nodes)) *
                             np.sqrt(1 + np.gradient(np.interp(xx, x, r), xx) ** 2), xx)
        ell = xa1 - xa0
    else:
        S_aft, ell = 0.0, 0.0
    Vj = 1.0
    for _ in range(4):                       # V_j dipende dalla resistenza della coda
        Cf = 0.074 * max(Vj * ell / nu, 1e3) ** -0.2 if ell > 0 else 0.0
        D_coda = 0.5 * Vj**2 * Cf * S_aft + D_res
        Vj = 1 + D_coda / max(Q, 1e-9)
    q_din = 0.5 * (Ue_s ** 2 if PERDITA == "Ue" else 1.0)
    P = (E_ing + D_coda * (Vj + 1) / 2 + K_condotto * q_din * Q) / eta_pompa
    CP = 2 * P + (pen if penal else 0.0)
    if not detail:
        return CP
    out.update(CP=float(CP), CP_ingestione=float(2 * E_ing), CP_coda=float(2 * D_coda * (Vj + 1) / 2),
               Q=float(Q), Q_BL=float(Q_BL), phi=float(phi), Ue_fessura=float(Ue_s), sep_coda=sep_c,
               penalita=float(pen), Vj=float(Vj), xc=xc, ue=ue)
    return out


def conventional(x, r, ReV, transition, A_disk=None):
    """Potenza del corpo liscio con elica in aria libera.
    A_disk None -> elica ideale (P = D U); altrimenti disco attuatore."""
    CDV = axi.drag(x, r, ReV, transition)
    if A_disk is None:
        return CDV
    D = 0.5 * CDV                    # ρ = U = 1, V = 1
    CT = D / (0.5 * A_disk)
    return CDV * (1 + np.sqrt(1 + CT)) / 2
