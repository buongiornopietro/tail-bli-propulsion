"""Resistenza di corpi di rivoluzione (3D assialsimmetrici) a Reynolds reali.

Metodo classico (Young 1939, Myring 1976, Lutz & Wagner 1998):
  1. flusso potenziale: pannelli ad anelli vorticosi sulla superficie
     (condizione di Dirichlet: funzione di corrente ψ = 0 sul corpo)
  2. strato limite integrale lungo il meridiano:
       laminare   - Thwaites assialsimmetrico (Rott-Crabtree)
       transizione - criterio di Michel (o transizione forzata al naso)
       turbolento - metodo di entrainment di Head + Ludwieg-Tillmann
  3. resistenza totale (attrito + forma) dalla formula di Young applicata
     allo spessore di quantità di moto vicino alla coda.

Unità: U = 1, volume del corpo V = 1.  Re_V = U V^(1/3) / ν.
C_DV = D / (½ρU² V^(2/3)).
"""
import numpy as np
from math import comb
from scipy.special import ellipk, ellipe
from scipy.interpolate import CubicSpline

REX_MAX = 5e6        # Re_x oltre cui la transizione avviene comunque (disturbi reali)
CRITERIO = "eN"      # "eN" (Drela-Giles, tiene conto del gradiente) oppure "michel"
TU = 0.001           # turbolenza esterna (frazione): fissa la soglia N del metodo e^N
X_END = 0.95         # stazione della formula di Young (frazione di L)
N_PANNELLI = 140

# ------------------------------------------------------------- geometria
def nodes(n, x_fitto=None, extra=0, sigma=0.04):
    """Nodi in ξ ∈ [0,1]: distribuzione a coseno (fitta a naso e coda) più,
    se x_fitto è dato, 'extra' nodi concentrati con una gaussiana attorno a x_fitto.
    Costruita invertendo una densità continua: la dimensione dei pannelli
    varia gradualmente (niente salti bruschi)."""
    if x_fitto is None or extra <= 0:
        return (1 - np.cos(np.linspace(0, np.pi, n + 1))) / 2
    g = np.linspace(0, 1, 20001)
    dens = n / (np.pi * np.sqrt(np.clip(g * (1 - g), 1e-6, None)))
    dens += extra * np.exp(-0.5 * ((g - x_fitto) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
    cdf = np.concatenate([[0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(g))])
    cdf /= cdf[-1]
    return np.interp(np.linspace(0, 1, n + extra + 1), cdf, g)


def body(params, n=None, xi=None):
    """Profilo CST: r(ξ) = ξ^0.5 (1-ξ) Σ w_k B_k(ξ), scalato a L/D = f e V = 1.

    params = (f, w_1/w_0, w_2/w_0, ...)  (valori diretti, non logaritmi)
    xi: nodi in [0,1] (di default distribuzione a coseno con N_PANNELLI pannelli)"""
    n = N_PANNELLI if n is None else n
    f = params[0]
    w = np.concatenate([[1.0], np.asarray(params[1:], float)])
    xi = nodes(n) if xi is None else np.asarray(xi, float)
    nn = len(w) - 1
    B = sum(w[k] * comb(nn, k) * xi**k * (1 - xi)**(nn - k) for k in range(nn + 1))
    r = np.sqrt(xi) * (1 - xi) * B
    r = r / r.max() / (2 * f)          # L = 1, D = 1/f
    return scale_to_unit_volume(xi, r)


def from_profile(x, r):
    """Porta un profilo qualunque (x da 0 a L, r>=0) a V = 1."""
    return scale_to_unit_volume(np.asarray(x, float) / x[-1], np.asarray(r, float) / x[-1])


def scale_to_unit_volume(x, r):
    V = np.pi * np.trapezoid(r**2, x)
    s = V ** (-1 / 3)
    return x * s, r * s


def ellipsoid(f, n=140):
    t = np.linspace(0, np.pi, n + 1)
    return from_profile((1 - np.cos(t)) / 2, np.sin(t) / (2 * f))


# ------------------------------------------------------ flusso potenziale
GL8 = np.polynomial.legendre.leggauss(8)
GL16 = np.polynomial.legendre.leggauss(16)


def ring_psi(x, r, x0, a):
    """ψ nel punto (x, r) di un anello vorticoso di raggio a in x0, Γ = 1."""
    d2 = (x - x0) ** 2 + (r + a) ** 2
    m = np.clip(4 * r * a / d2, 0, 1 - 1e-15)
    k = np.sqrt(m)
    k = np.where(k < 1e-12, 1e-12, k)
    return np.sqrt(r * a) / (2 * np.pi) * ((2 / k - k) * ellipk(m) - 2 / k * ellipe(m))


def surface_velocity(x, r):
    """Velocità tangenziale sul corpo (U∞ = 1) ai punti medi dei pannelli."""
    xc, rc = (x[1:] + x[:-1]) / 2, (r[1:] + r[:-1]) / 2
    ds = np.hypot(np.diff(x), np.diff(r))
    N = len(xc)
    A = np.zeros((N, N))
    t, wg = GL8
    for j in range(N):
        tt = (t + 1) / 2
        xq = x[j] + tt * (x[j + 1] - x[j]); rq = r[j] + tt * (r[j + 1] - r[j])
        A[:, j] = (ring_psi(xc[:, None], rc[:, None], xq[None], rq[None]) * wg / 2).sum(1) * ds[j]
        # pannello su se stesso: singolarità logaritmica -> due metà, 16 punti ciascuna
        t2, w2 = GL16
        acc = 0.0
        for a0, a1 in ((0.0, 0.5), (0.5, 1.0)):
            tt = a0 + (t2 + 1) / 2 * (a1 - a0)
            xq = x[j] + tt * (x[j + 1] - x[j]); rq = r[j] + tt * (r[j + 1] - r[j])
            acc += (ring_psi(xc[j], rc[j], xq, rq) * w2).sum() * (a1 - a0) / 2
        A[j, j] = acc * ds[j]
    g = np.linalg.solve(A, -0.5 * rc**2)
    return xc, rc, np.abs(g)


# ----------------------------------------------------------- strato limite
def H1_of_H(H):
    return np.where(H <= 1.6, 3.3 + 0.8234 * np.maximum(H - 1.1, 1e-6) ** -1.287,
                    3.3 + 1.5501 * np.maximum(H - 0.6778, 1e-6) ** -3.064)


def H_of_H1(H1):
    H1 = np.maximum(H1, 3.3 + 1e-6)
    return np.where(H1 >= 5.3, 1.1 + ((H1 - 3.3) / 0.8234) ** (-1 / 1.287),
                    0.6778 + ((H1 - 3.3) / 1.5501) ** (-1 / 3.064))


def thwaites_H(lam):
    lam = np.clip(lam, -0.09, 0.25)
    return np.where(lam >= 0, 2.61 - 3.75 * lam + 5.24 * lam**2,
                    2.088 + 0.0731 / (lam + 0.14))


def n_crit(Tu=None):
    """Soglia del fattore di amplificazione (Mack): N = -8.43 - 2.4 ln(Tu)."""
    Tu = TU if Tu is None else Tu
    return -8.43 - 2.4 * np.log(max(Tu, 1e-6))


def dn_dRetheta(H):
    """Pendenza dell'inviluppo di amplificazione (Drela & Giles 1987)."""
    return 0.01 * np.sqrt((2.4 * H - 3.7 + 2.5 * np.tanh(1.5 * H - 4.65)) ** 2 + 0.25)


def Retheta_crit(H):
    """Re_theta critico: inizio dell'amplificazione (Drela & Giles 1987)."""
    h = np.maximum(H - 1, 0.05)
    lg = (1.415 / h - 0.489) * np.tanh(20 / h - 12.9) + 3.295 / h + 0.44
    return 10 ** lg


def indice_transizione(s, Ue, th, lam, nu, Rex_max=None, criterio=None, Tu=None):
    """Indice di transizione lungo la marcia laminare.

    criterio "eN": si integra dn/ds = (dn/dRe_theta)(H) dRe_theta/ds dove
    Re_theta supera il valore critico; transizione quando n raggiunge N(Tu).
    Il distacco laminare (Thwaites lambda < -0.09) e il tetto Re_x restano
    come innesco aggiuntivo.
    """
    criterio = CRITERIO if criterio is None else criterio
    Rex_max = REX_MAX if Rex_max is None else Rex_max
    Reth = Ue * th / nu
    Rex = Ue * (s - s[0]) / nu + 1e-9
    s0 = s[0] + 0.02 * (s[-1] - s[0])
    if criterio == "michel":
        trig = Reth > 1.174 * (1 + 22400 / Rex) * Rex**0.46
    else:
        H = thwaites_H(lam)
        amp = (Reth > Retheta_crit(H)) & (np.gradient(Reth, s) > 0)
        dn = np.where(amp, dn_dRetheta(H) * np.gradient(Reth, s), 0.0)
        n = np.concatenate([[0], np.cumsum(0.5 * (dn[1:] + dn[:-1]) * np.diff(s))])
        trig = n >= n_crit(Tu)
    cand = np.flatnonzero((trig | (lam < -0.09) | (Rex > Rex_max)) & (s > s0))
    return int(cand[0]) if len(cand) else len(s) - 1


def drag(x, r, ReV, transition="naturale", x_end=None, detail=False, npts=2500,
         Rex_max=None):
    """Coefficiente di resistenza volumetrico C_DV del corpo (V = 1)."""
    Rex_max = REX_MAX if Rex_max is None else Rex_max
    x_end = X_END if x_end is None else x_end
    nu = 1.0 / ReV
    xc, rc, ue = surface_velocity(x, r)
    L = x[-1]
    # ascissa curvilinea lungo il meridiano
    s_nodes = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(r)))])
    sc = (s_nodes[1:] + s_nodes[:-1]) / 2
    keep = xc <= x_end * L + 1e-12
    S = CubicSpline(np.concatenate([[0], sc[keep]]), np.concatenate([[0], ue[keep]]))
    Rs = CubicSpline(s_nodes, r)
    Xs = CubicSpline(s_nodes, x)
    s_end = np.interp(x_end * L, xc, sc)
    s = np.linspace(0, s_end, npts)[1:]
    Ue = np.maximum(S(s), 1e-6); dUe = S(s, 1)
    R = np.maximum(Rs(s), 1e-9); dR = Rs(s, 1)

    # --- laminare: Thwaites/Rott-Crabtree
    integ = np.concatenate([[0], np.cumsum(0.5 * (Ue[1:]**5 * R[1:]**2 + Ue[:-1]**5 * R[:-1]**2) * np.diff(s))])
    th = np.sqrt(0.45 * nu * integ / (Ue**6 * R**2) + 1e-30)
    lam = th**2 * dUe / nu
    Reth = Ue * th / nu; Rex = Ue * s / nu
    if transition == "forzata":
        itr = min(int(0.03 * len(s)), len(s) - 2)           # inciampo vicino al naso
    else:
        itr = indice_transizione(s, Ue, th, lam, nu, Rex_max)
    laminar_sep = transition != "forzata" and lam[itr] < -0.09

    # --- turbolento: Head
    theta = th[itr]; H = 1.4; H1 = float(H1_of_H(H))
    sep_at = None
    for i in range(itr, len(s) - 1):
        h = s[i + 1] - s[i]
        def rhs(i_, theta_, H1_):
            H_ = float(H_of_H1(H1_))
            Rt = max(Ue[i_] * theta_ / nu, 10.0)
            Cf = 0.246 * 10 ** (-0.678 * H_) * Rt ** -0.268
            dth = Cf / 2 - theta_ * ((H_ + 2) * dUe[i_] / Ue[i_] + dR[i_] / R[i_])
            F = 0.0306 * max(H1_ - 3.0, 1e-6) ** -0.6169
            dH1 = F / theta_ - H1_ * (dUe[i_] / Ue[i_] + dth / theta_ + dR[i_] / R[i_])
            return dth, dH1
        a1, b1 = rhs(i, theta, H1)
        a2, b2 = rhs(i + 1, theta + h * a1, H1 + h * b1)
        theta += h * (a1 + a2) / 2; H1 += h * (b1 + b2) / 2
        H = float(H_of_H1(H1))
        if H > 2.4 or not np.isfinite(theta):
            sep_at = i + 1
            break
    iend = sep_at if sep_at is not None else len(s) - 1
    ue_e = Ue[iend]; r_e = R[iend]
    theta = theta if np.isfinite(theta) else th[iend]
    CDV = 4 * np.pi * r_e * theta * ue_e ** ((min(H, 2.4) + 5) / 2)
    if sep_at is not None:
        # flusso staccato: Young non vede la scia -> resistenza di base
        # con Cp_base ≈ -0.2 sulla sezione di distacco (stima prudente)
        CDV += 0.2 * np.pi * r_e**2
    if not detail:
        return CDV
    return dict(CDV=float(CDV),
                x_trans=float(Xs(s[itr]) / L), sep_lam=bool(laminar_sep),
                x_sep_turb=None if sep_at is None else float(Xs(s[sep_at]) / L),
                H_end=float(H), theta_end=float(theta), L=float(L), D=float(2 * r.max()),
                f=float(L / (2 * r.max())), xc=xc, ue=ue)


def CD_frontal(CDV, x, r):
    """Converte C_DV in coefficiente riferito alla sezione frontale."""
    return CDV / (np.pi * r.max() ** 2)


def hoerner_frontal(f, ReL):
    """Hoerner, Fluid-Dynamic Drag (1965) cap. 6: corpi di rivoluzione affusolati,
    strato limite turbolento. C_D (frontale) = Cf [3 f + 4.5 f^-0.5 + 21 f^-2]."""
    Cf = 0.455 / np.log10(ReL) ** 2.58          # Prandtl-Schlichting, lastra piana
    return Cf * (3 * f + 4.5 / np.sqrt(f) + 21 / f**2)
