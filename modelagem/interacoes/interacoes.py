"""M8 — Interações intermoleculares cafeína–epicatequina (estudo Edson).

Pergunta: no extrato etanólico de guaraná, cafeína e epicatequina (unidade dos taninos condensados) coexistem.
Elas se associam em água? Por empilhamento π (cafeína sobre o anel A ou B) ou por ligação de hidrogênio?
E o complexo deixa livres os OH do catecol (anel B), que são os doadores de H• da atividade antioxidante (M2)?

Protocolo
  1. Monômeros: geometrias do M2 (GFN2-xTB/ALPB metanol) reotimizadas em GFN2-xTB/ALPB(água).
  2. Busca de complexos (partidas determinísticas; semente fixa nas aleatórias):
       empilhamento: cafeína paralela ao anel A ou B, nas duas faces, 6 rotações × 2 orientações, a 3,4 Å
                     (3,9 ou 4,5 Å se houver choque);
       ligação H:    um aceptor da cafeína (O2, O6, N9) a 1,9 Å de cada H de OH da epicatequina, na direção O–H;
       aleatórias:   24 orientações a ≥ 5 Å entre centroides (afastadas até não haver choque).
     Cada partida é otimizada em GFN2-xTB/ALPB(água) (geomeTRIC). Duplicatas removidas por energia + RMSD.
  3. Energia de ligação: ΔE = E(complexo) − E(cafeína) − E(epicatequina), monômeros relaxados (inclui deformação).
     Para os 3 complexos mais estáveis de cada tipo: Hessiana xTB → ΔH e ΔG (RRHO com entropia quase-harmônica
     de Grimme, ω0 = 100 cm⁻¹), estado-padrão 1 M (−1,89 kcal/mol para a associação bimolecular, 298,15 K).
  4. NCI (índice de interações não covalentes; JOHNSON et al., 2010) sobre densidade promolecular
     (soma de densidades atômicas esféricas, SAD/def2-SVP no PySCF): s(ρ) = |∇ρ| / [2(3π²)^(1/3) ρ^(4/3)],
     cor por sign(λ2)ρ. Filtro intermolecular: só pontos em que nenhum fragmento contribui com > 95 % de ρ.
     Gera volumes para o visualizador 3D e o gráfico s × sign(λ2)ρ.

Limitações: GFN2-xTB superestima um pouco a dispersão em empilhamentos; o solvente é contínuo (ALPB), sem
águas explícitas competindo pelas ligações H; ΔG por RRHO é aproximado (ordem de grandeza e ranking).

Uso
  "<venv>/bin/python" interacoes.py              tudo (retomável: cada partida vira resultados/partidas/<id>.json)
  "<venv>/bin/python" interacoes.py --analise    só reanalisa (energias, classificação, NCI, termo)
"""
from __future__ import annotations

import argparse
import base64
import json
import math
import os
import time
from pathlib import Path

for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(v, "1")

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

from quimiolib import qm

AQUI = Path(__file__).resolve().parent
EDSON = AQUI.parents[1]
GEO_M3 = EDSON / "modelagem" / "tddft_uv" / "geometrias"
RES = AQUI / "resultados"
PART = RES / "partidas"
SOLV = "water"
KCAL = qm.HARTREE_KCAL
PAR = ("cafeina", "epicatequina")
NOMES = {"cafeina": "Cafeína", "epicatequina": "(−)-Epicatequina"}
CONC_1M = -1.89   # kcal/mol: 1 atm → 1 M para A + B → AB a 298,15 K (RT ln 24,46)


# ------------------------------------------------------------------ geometria e topologia
def ler(mid):
    d = json.loads((GEO_M3 / f"{mid}.json").read_text())
    return list(d["simbolos"]), np.asarray(d["xyz"], float)


def rdmol(simb, xyz):
    blk = f"{len(simb)}\n\n" + qm.xyz_str(simb, xyz)
    m = Chem.MolFromXYZBlock(blk)
    rdDetermineBonds.DetermineConnectivity(m)
    Chem.FastFindRings(m)
    return m


def plano(xyz):
    c = xyz.mean(0)
    _, _, vt = np.linalg.svd(xyz - c)
    return c, vt[2]


def aneis_epi(m, xyz):
    """Anel A (fundido ao C, resorcinol) e anel B (catecol, pendente) da epicatequina."""
    ri = m.GetRingInfo()
    seis = [list(r) for r in ri.AtomRings() if len(r) == 6]
    arom = [r for r in seis if all(m.GetAtomWithIdx(i).GetSymbol() == "C" for i in r)
            and abs(np.linalg.svd(xyz[r] - xyz[r].mean(0))[1][2]) < 0.15]
    todos = [list(r) for r in ri.AtomRings()]
    fundido = lambda r: any(set(r) & set(o) for o in todos if set(o) != set(r))
    A = next(r for r in arom if fundido(r))
    B = next(r for r in arom if not fundido(r))
    return {"A": A, "B": B}


def anel_caf(m):
    ri = m.GetRingInfo()
    at = sorted(set(i for r in ri.AtomRings() for i in r))
    return at   # 9 átomos do sistema purínico


def oh_epi(m):
    """[(idx_O, idx_H, rótulo)] dos OH da epicatequina; rótulo pelo anel a que o O está ligado."""
    return [(a.GetIdx(), n.GetIdx()) for a in m.GetAtoms() if a.GetSymbol() == "O"
            for n in a.GetNeighbors() if n.GetSymbol() == "H"]


def aceptores_caf(m):
    """[(idx_aceptor, idx_vizinho_para_direção, rótulo)]: carbonilas O2/O6 e o N9 (imidazol, sem H/CH3)."""
    out = []
    for a in m.GetAtoms():
        if a.GetSymbol() == "O":
            out.append((a.GetIdx(), a.GetNeighbors()[0].GetIdx()))
        if a.GetSymbol() == "N" and a.GetDegree() == 2:
            out.append((a.GetIdx(), [n.GetIdx() for n in a.GetNeighbors()]))
    return out


def rot_eixo(u, ang):
    u = u / np.linalg.norm(u)
    K = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])
    return np.eye(3) + math.sin(ang) * K + (1 - math.cos(ang)) * K @ K


def alinhar(a, b):
    """Rotação que leva o vetor unitário a em b."""
    a, b = a / np.linalg.norm(a), b / np.linalg.norm(b)
    v, c = np.cross(a, b), float(a @ b)
    if np.linalg.norm(v) < 1e-8:
        if c > 0:
            return np.eye(3)
        p = np.array([1, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1, 0])
        return rot_eixo(np.cross(a, p), math.pi)
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + K + K @ K * (1 / (1 + c))


def _dmin(a, b):
    return float(np.linalg.norm(a[:, None] - b[None], axis=2).min())


def partidas(caf, epi):
    (sc, xc), (se, xe) = caf, epi
    mc, me = rdmol(sc, xc), rdmol(se, xe)
    cc, nc = plano(xc[anel_caf(mc)])
    xc0 = xc - cc
    out = []
    # empilhamento
    for nome, r in aneis_epi(me, xe).items():
        ce, ne = plano(xe[r])
        for face in (1, -1):
            for flip in (0, 1):
                R0 = alinhar(nc, ne) @ (rot_eixo(np.cross(nc, [1, 0, 0]) if abs(nc[0]) < .9 else
                                                 np.cross(nc, [0, 1, 0]), math.pi) if flip else np.eye(3))
                for k in range(6):
                    R = rot_eixo(ne, k * math.pi / 3) @ R0
                    for alt in (3.4, 3.9, 4.5):   # sobe até não haver choque com o resto da epicatequina
                        pos = xc0 @ R.T + ce + face * alt * ne
                        if _dmin(pos, xe) > 1.6:
                            break
                    out.append({"id": f"pi_{nome}{'+' if face > 0 else '-'}_{flip}{k}", "tipo0": f"π {nome}",
                                "xyz": pos})
    # ligação H
    for io, ih in oh_epi(me):
        u = xe[ih] - xe[io]; u /= np.linalg.norm(u)
        for ia, viz in aceptores_caf(mc):
            d = xc[ia] - (xc[viz].mean(0) if isinstance(viz, list) else xc[viz])
            for k in range(2):
                R = rot_eixo(u, k * math.pi) @ alinhar(d, -u)
                pos = (xc - xc[ia]) @ R.T + xe[ih] + 1.9 * u
                out.append({"id": f"hb_O{io}_A{ia}_{k}", "tipo0": "ligação H", "xyz": pos})
    # aleatórias (semente fixa)
    rng = np.random.default_rng(7)
    ce = xe.mean(0)
    for k in range(24):
        q = rng.normal(size=4); q /= np.linalg.norm(q)
        a, b, c, d = q
        R = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)], [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                      [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])
        u = rng.normal(size=3); u /= np.linalg.norm(u)
        dist = rng.uniform(5, 6)
        while _dmin(xc0 @ R.T + ce + dist * u, xe) < 2.2:   # afasta até sair do choque
            dist += 0.3
        out.append({"id": f"rnd_{k:02d}", "tipo0": "aleatória", "xyz": xc0 @ R.T + ce + dist * u})
    # descarta choques graves
    ok = []
    for p in out:
        if _dmin(p["xyz"], xe) > 1.5:
            ok.append(p)
    return ok


# ------------------------------------------------------------------ otimização (um processo por partida)
def otimizar(id, simb, xyz):
    t = time.time()
    x, e, conv = qm.xtb_otimizar(simb, np.asarray(xyz), solvente=SOLV, maxiter=500)
    return {"id": id, "xyz": np.asarray(x).round(5).tolist(), "e_h": e, "convergiu": conv,
            "t_s": round(time.time() - t, 1)}


# ------------------------------------------------------------------ termoquímica (qRRHO)
def termo(simb, xyz, T=298.15):
    """Hessiana numérica xTB → H e G (RRHO; entropia vibracional quase-harmônica de Grimme). kcal/mol."""
    pt = Chem.GetPeriodicTable()
    nums = [pt.GetAtomicNumber(s) for s in simb]
    massas = np.array([qm.MASSAS.get(z, pt.GetMostCommonIsotopeMass(z)) for z in nums])
    x0 = np.asarray(xyz, float).ravel() / qm.BOHR
    n3, h = len(x0), 0.005
    H = np.zeros((n3, n3))
    for i in range(n3):
        xp, xm = x0.copy(), x0.copy(); xp[i] += h; xm[i] -= h
        H[i] = (qm.xtb_eg(nums, xp.reshape(-1, 3), solvente=SOLV)[1].ravel()
                - qm.xtb_eg(nums, xm.reshape(-1, 3), solvente=SOLV)[1].ravel()) / (2 * h)
    H = 0.5 * (H + H.T)
    sm = np.sqrt(np.repeat(massas, 3))
    P, ntr = qm._proj_tr_rot(x0.reshape(-1, 3), massas)
    w = np.linalg.eigvalsh(P @ (H / np.outer(sm, sm)) @ P)
    w = np.sort(w[np.argsort(np.abs(w))[ntr:]])
    fr = np.sign(w) * np.sqrt(np.abs(w)) * qm.FREQ_CONV
    n_imag = int((fr < -50).sum())
    nu = np.abs(fr[np.abs(fr) > 1.0]); nu = np.where(nu < 20, 20, nu)   # pequenos/imag. → 20 cm⁻¹ (tratados como rotores)
    R, kB, hP, c = 1.98720425e-3, 1.380649e-23, 6.62607015e-34, 2.99792458e10
    kT = qm.KB_KCAL * T
    e = nu * qm.CM1_KCAL
    zpe, evib = 0.5 * e.sum(), (e / np.expm1(e / kT)).sum()
    x = e / kT
    s_ho = R * (x / np.expm1(x) - np.log1p(-np.exp(-x)))
    # rotor livre (Grimme 2012): Bav = 1e-44 kg m², μ = h/(8π²ν)
    mu = hP / (8 * math.pi ** 2 * nu * c)
    mu_ef = mu * 1e-44 / (mu + 1e-44)
    s_r = R * (0.5 + np.log(np.sqrt(8 * math.pi ** 3 * mu_ef * kB * T / hP ** 2)))
    wgt = 1 / (1 + (100 / nu) ** 4)
    s_vib = (wgt * s_ho + (1 - wgt) * s_r).sum()
    # translação (1 atm) e rotação (rotor rígido assimétrico, σ = 1)
    M = massas.sum() * 1.66053907e-27
    V = kB * T / 101325
    s_tr = R * (math.log((2 * math.pi * M * kB * T / hP ** 2) ** 1.5 * V) + 2.5)
    xyzm = np.asarray(xyz) * 1e-10
    com = (xyzm * massas[:, None]).sum(0) / massas.sum()
    r = xyzm - com
    I = sum(m * (np.dot(v, v) * np.eye(3) - np.outer(v, v)) for m, v in zip(massas * 1.66053907e-27, r))
    Ia, Ib, Ic = np.linalg.eigvalsh(I)
    s_rot = R * (math.log(math.sqrt(math.pi) * (8 * math.pi ** 2 * kB * T / hP ** 2) ** 1.5
                          * math.sqrt(Ia * Ib * Ic)) + 1.5)
    h_corr = zpe + evib + 3 * kT + kT
    s_tot = s_vib + s_tr + s_rot
    return {"h_corr": float(h_corr), "TS": float(T * s_tot), "g_corr": float(h_corr - T * s_tot),
            "n_imag": n_imag, "freq_min": float(fr.min())}


# ------------------------------------------------------------------ análise geométrica
def classificar(simb, xyz, nc_caf, mc, me, aneis):
    xc, xe = xyz[:nc_caf], xyz[nc_caf:]
    cc, ncv = plano(xc[anel_caf(mc)])
    pi = []
    for nome, r in aneis.items():
        ce, ne = plano(xe[r])
        d = cc - ce
        ang = math.degrees(math.acos(min(1, abs(float(ncv @ ne)))))
        h = abs(float(d @ ne))
        desl = float(np.linalg.norm(d - (d @ ne) * ne))
        pi.append({"anel": nome, "dist": round(float(np.linalg.norm(d)), 2), "altura": round(h, 2),
                   "deslocamento": round(desl, 2), "angulo": round(ang, 1),
                   "empilhado": bool(np.linalg.norm(d) < 4.6 and ang < 30 and desl < 2.5)})
    hb = []
    se = simb[nc_caf:]
    for io, ih in oh_epi(me):
        for j, s in enumerate(simb[:nc_caf]):
            if s in "ON":
                dHA = float(np.linalg.norm(xe[ih] - xc[j]))
                v1, v2 = xe[io] - xe[ih], xc[j] - xe[ih]
                ang = math.degrees(math.acos(float(v1 @ v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))))
                if dHA < 2.4 and ang > 120:
                    hb.append({"doador": f"O{io}", "aceptor": f"{s}{j}", "H_A": round(dHA, 2), "ang": round(ang, 1),
                               "sitio": rotulo_oh(io, me, aneis)})
    emp = [p for p in pi if p["empilhado"]]
    tipo = (f"π-empilhamento (anel {emp[0]['anel']})" if emp else "ligação H" if hb else "contato vdW")
    if emp and hb:
        tipo += " + ligação H"
    return tipo, pi, hb


def rotulo_oh(io, me, aneis):
    c = me.GetAtomWithIdx(io).GetNeighbors()
    cs = [n.GetIdx() for n in c if n.GetSymbol() == "C"]
    if not cs:
        return "?"
    for nome, r in aneis.items():
        if cs[0] in r:
            return f"OH anel {nome}"
    return "3-OH (anel C)"


def rmsd(a, b):
    """RMSD após sobreposição ótima (Kabsch), mesma ordem de átomos."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a - a.mean(0), b - b.mean(0)
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    R = u @ np.diag([1, 1, d]) @ vt
    return float(np.sqrt(((a @ R - b) ** 2).sum(1).mean()))


# ------------------------------------------------------------------ NCI promolecular
def nci(simb, xyz, nc_caf, passo=0.12, margem=3.0):
    from pyscf import gto, scf
    B = 1 / qm.BOHR

    def frag(s, x):
        m = gto.M(atom=[(a, tuple(p)) for a, p in zip(s, x)], basis="def2-svp", verbose=0)
        return m, scf.hf.init_guess_by_atom(m)

    fa, fb = frag(simb[:nc_caf], xyz[:nc_caf]), frag(simb[nc_caf:], xyz[nc_caf:])
    # caixa: só a região entre as moléculas (átomos a < 5 Å do outro fragmento)
    xa, xb = xyz[:nc_caf], xyz[nc_caf:]
    D = np.linalg.norm(xa[:, None] - xb[None], axis=2)
    perto = np.vstack([xa[D.min(1) < 5], xb[D.min(0) < 5]])
    lo, hi = perto.min(0) - margem, perto.max(0) + margem
    n = np.ceil((hi - lo) / passo).astype(int) + 1
    g = [lo[i] + passo * np.arange(n[i]) for i in range(3)]
    pts = np.stack(np.meshgrid(*g, indexing="ij"), -1).reshape(-1, 3) * B   # bohr

    def rho(fr):
        m, dm = fr
        out = np.empty(len(pts))
        for i in range(0, len(pts), 20000):
            ao = m.eval_gto("GTOval", pts[i:i + 20000])
            out[i:i + 20000] = np.einsum("pi,ij,pj->p", ao, dm, ao, optimize=True)
        return out.reshape(n)

    ra, rb = rho(fa), rho(fb)
    r = ra + rb
    hb = passo * B
    gr = np.gradient(r, hb)
    Hs = np.empty(tuple(n) + (3, 3))
    for i in range(3):
        gi = np.gradient(gr[i], hb)
        for j in range(3):
            Hs[..., i, j] = gi[j]
    Hs = 0.5 * (Hs + np.swapaxes(Hs, -1, -2))
    lam2 = np.linalg.eigvalsh(Hs)[..., 1]
    gnorm = np.sqrt(sum(x ** 2 for x in gr))
    rs = np.maximum(r, 1e-12)
    s = gnorm / (2 * (3 * math.pi ** 2) ** (1 / 3) * rs ** (4 / 3))
    slr = np.sign(lam2) * r
    inter = np.maximum(ra, rb) / rs < 0.95
    valido = inter & (r < 0.05) & (r > 1e-5)
    s_v = np.where(valido, s, 2.0)
    # gráfico s × sign(λ2)ρ (pontos com s < 1)
    sel = valido & (s < 1.0)
    pts2 = np.stack([slr[sel], s[sel]], 1)
    if len(pts2) > 4000:
        pts2 = pts2[np.random.default_rng(0).choice(len(pts2), 4000, replace=False)]
    # integrais sobre a região s < 0,5 (descritor por tipo)
    reg = valido & (s < 0.5)
    dv = hb ** 3
    integ = {"atrativa_forte": float(r[reg & (slr < -0.01)].sum() * dv),
             "vdw": float(r[reg & (np.abs(slr) <= 0.01)].sum() * dv),
             "repulsiva": float(r[reg & (slr > 0.01)].sum() * dv),
             "volume_A3": float(reg.sum() * passo ** 3)}
    # recorta a caixa à interface (s < 0,8) + 0,8 Å, para o volume do painel ficar leve
    idx = np.argwhere(s_v < 0.8)
    pad = int(round(0.8 / passo))
    a0 = np.maximum(idx.min(0) - pad, 0); a1 = np.minimum(idx.max(0) + pad + 1, n)
    cut = tuple(slice(a0[i], a1[i]) for i in range(3))
    s_v, slr = s_v[cut], slr[cut]
    lo = lo + a0 * passo; n = a1 - a0
    q_s = np.clip(np.round(s_v / 2.0 * 255), 0, 255).astype(np.uint8)
    q_c = np.clip(np.round(slr / 0.05 * 127), -127, 127).astype(np.int8)
    vol = {"n": n.tolist(), "origem_bohr": (lo * B).tolist(), "passo_bohr": hb,
           "s_escala": 2.0 / 255, "c_escala": 0.05 / 127,
           "s": base64.b64encode(q_s.tobytes()).decode(), "c": base64.b64encode(q_c.tobytes()).decode()}
    return vol, pts2.round(4).tolist(), integ


# ------------------------------------------------------------------ pipeline
def main(so_analise=False, workers=4):
    RES.mkdir(exist_ok=True); PART.mkdir(exist_ok=True)
    caf, epi = ler(PAR[0]), ler(PAR[1])
    mono = {}
    fm = RES / "monomeros.json"
    if fm.exists():
        mono = json.loads(fm.read_text())
    for mid, (s, x) in zip(PAR, (caf, epi)):
        if mid not in mono:
            print(f"monômero {mid} (ALPB água) …", flush=True)
            xo, e, ok = qm.xtb_otimizar(s, x, solvente=SOLV, maxiter=500)
            mono[mid] = {"simbolos": s, "xyz": np.asarray(xo).round(5).tolist(), "e_h": e, "convergiu": ok}
            fm.write_text(json.dumps(mono, indent=1))
    caf = (mono["cafeina"]["simbolos"], np.asarray(mono["cafeina"]["xyz"]))
    epi = (mono["epicatequina"]["simbolos"], np.asarray(mono["epicatequina"]["xyz"]))
    simb = caf[0] + epi[0]
    nc = len(caf[0])
    ps = partidas(caf, epi)
    print(f"{len(ps)} partidas", flush=True)
    if not so_analise:
        tarefas = [{"id": p["id"], "simb": simb, "xyz": np.vstack([p["xyz"], epi[1]]).tolist()} for p in ps]
        qm.rodar_paralelo(otimizar, tarefas, workers=workers, pasta=PART, chave="id")

    # ---- análise
    mc, me = rdmol(*caf), rdmol(*epi)
    aneis = aneis_epi(me, epi[1])
    e0 = mono["cafeina"]["e_h"] + mono["epicatequina"]["e_h"]
    t0 = {p["id"]: p["tipo0"] for p in ps}
    lst = []
    for f in sorted(PART.glob("*.json")):
        r = json.loads(f.read_text())
        x = np.asarray(r["xyz"])
        # complexo desfeito (moléculas afastadas) não conta
        dmin = float(np.linalg.norm(x[:nc, None] - x[None, nc:], axis=2).min())
        if dmin > 3.5:
            continue
        # conectividade preservada (nada reagiu)
        if rdmol(simb, x).GetNumBonds() != mc.GetNumBonds() + me.GetNumBonds():
            continue
        tipo, pi, hb = classificar(simb, x, nc, mc, me, aneis)
        lst.append({"id": r["id"], "partida": t0.get(r["id"], "?"), "tipo": tipo, "dE": (r["e_h"] - e0) * KCAL,
                    "convergiu": r["convergiu"], "pi": pi, "hb": hb, "xyz": x})
    lst.sort(key=lambda d: d["dE"])
    unicos = []
    for d in lst:
        if any(abs(d["dE"] - u["dE"]) < 0.15 and rmsd(d["xyz"], u["xyz"]) < 0.6 for u in unicos):
            continue
        unicos.append(d)
    print(f"{len(lst)} complexos válidos → {len(unicos)} únicos", flush=True)

    # ---- termoquímica: 3 melhores de cada família
    fam = lambda t: "π" if t.startswith("π") else ("H" if t.startswith("ligação") else "vdW")
    escolhidos, cont = [], {}
    for u in unicos:
        k = fam(u["tipo"])
        if cont.get(k, 0) < 3:
            escolhidos.append(u); cont[k] = cont.get(k, 0) + 1
    ft = RES / "termo.json"
    tm = json.loads(ft.read_text()) if ft.exists() else {}
    for mid, (s, x) in zip(PAR, (caf, epi)):
        if mid not in tm:
            print(f"termo {mid} …", flush=True); tm[mid] = termo(s, x); ft.write_text(json.dumps(tm, indent=1))
    for u in escolhidos:
        if u["id"] not in tm:
            print(f"termo {u['id']} …", flush=True); tm[u["id"]] = termo(simb, u["xyz"])
            ft.write_text(json.dumps(tm, indent=1))
    for u in escolhidos:
        a, b, c = tm[u["id"]], tm["cafeina"], tm["epicatequina"]
        u["dH"] = u["dE"] + a["h_corr"] - b["h_corr"] - c["h_corr"]
        u["TdS"] = -(a["TS"] - b["TS"] - c["TS"])
        u["dG"] = u["dE"] + a["g_corr"] - b["g_corr"] - c["g_corr"] + CONC_1M
        u["n_imag"] = a["n_imag"]

    # ---- NCI dos 3 melhores no geral + melhor de cada família
    alvos = unicos[:3] + [next((u for u in unicos if fam(u["tipo"]) == k), None) for k in ("π", "H")]
    alvos = [a for i, a in enumerate(alvos) if a and a not in alvos[:i]]
    vols = {}
    for u in alvos:
        fn = RES / "nci" / f"{u['id']}.json"
        fn.parent.mkdir(exist_ok=True)
        if not fn.exists():
            print(f"NCI {u['id']} …", flush=True)
            v, sc, ig = nci(simb, u["xyz"], nc)
            fn.write_text(json.dumps({"vol": v, "scatter": sc, "integrais": ig}))
        vols[u["id"]] = json.loads(fn.read_text())

    # ---- resumo
    sitios_livres = {}
    for u in unicos:
        usados = {h["sitio"] for h in u["hb"]}
        u["catecol_livre"] = not any("anel B" in s for s in usados)
    resumo = {
        "nivel": "GFN2-xTB/ALPB(água); NCI promolecular (SAD/def2-SVP)",
        "monomeros": {k: {"e_h": v["e_h"], "xyz": v["xyz"], "simbolos": v["simbolos"]} for k, v in mono.items()},
        "n_partidas": len(ps), "n_validos": len(lst), "n_caf": nc, "simbolos": simb,
        "complexos": [{k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in u.items()
                       if k != "xyz"} | {"xyz": np.asarray(u["xyz"]).round(4).tolist()}
                      for u in unicos],
        "nci": vols,
        "boltzmann_dE": boltz([u["dE"] for u in unicos], [fam(u["tipo"]) for u in unicos]),
    }
    (RES / "resumo.json").write_text(json.dumps(resumo, ensure_ascii=False))
    for u in unicos[:12]:
        print(f"  {u['id']:<18} {u['tipo']:<32} ΔE={u['dE']:6.2f}"
              + (f"  ΔG={u['dG']:6.2f}" if "dG" in u else "") + f"  HB={len(u['hb'])}", flush=True)


def boltz(es, fams, T=298.15):
    kT = qm.KB_KCAL * T
    e = np.array(es); w = np.exp(-(e - e.min()) / kT); w /= w.sum()
    out = {}
    for f, x in zip(fams, w):
        out[f] = out.get(f, 0) + float(x)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--analise", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    main(a.analise, a.workers)
