"""M8b — Refino DFT dos melhores complexos cafeína–epicatequina (estudo Edson).

Pergunta: o empate entre empilhamento sobre o anel A e sobre o anel B (catecol), visto com GFN2-xTB, se mantém com
um funcional com dispersão de alta qualidade?

Protocolo (DFT//xTB, como no M2)
  Geometrias: as do M8 (GFN2-xTB/ALPB água), sem reotimizar.
  Complexos:  o melhor de cada tipo de contato π distinto do M8 (ex.: anel A + ligação H, anel B, anel A puro).
  Método:     ωB97X-V/def2-SVP (VV10 não local, nativo do PySCF; dispensa D3), ajuste de densidade,
              malha XC nível 3, malha VV10 nível 1.
  Solvente:   IEF-PCM água (ε = 78,3553) para complexo e monômeros.
  BSSE:       contrapeso de Boys–Bernardi em fase gasosa (monômeros com átomos fantasma do parceiro), somado à
              energia em PCM: a correção é um efeito da base, praticamente independente do solvente.
  Energias:   ΔE_lig = E(AB) − E(A relaxado) − E(B relaxado) + δ_BSSE        (inclui a deformação)
              ΔE_int = E(AB) − E(A na geom. do complexo) − E(B na geom. do complexo) + δ_BSSE
              ΔG(1 M) = ΔE_lig + [G_corr(AB) − G_corr(A) − G_corr(B)]_xTB/qRRHO − 1,89

Uso
  "<venv>/bin/python" refino_dft.py [--workers 2]     retomável: cada cálculo vira resultados/dft/<chave>.json
  "<venv>/bin/python" refino_dft.py --analise         só monta resultados/refino_dft.json
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(v, "1")

import numpy as np

from quimiolib import qm

AQUI = Path(__file__).resolve().parent
RES = AQUI / "resultados"
DFT = RES / "dft"
XC, BASE, EPS = "wb97x-v", "def2-svp", 78.3553
KCAL = qm.HARTREE_KCAL
CONC_1M = -1.89


def sp(chave, atomos, pcm, max_memory=1800):
    """Single point. atomos: [(símbolo, (x, y, z))]; símbolos 'ghost-X' são átomos fantasma (só funções de base)."""
    from pyscf import dft, gto
    t = time.time()
    m = gto.M(atom=atomos, basis=BASE, verbose=0, max_memory=max_memory)
    mf = dft.RKS(m, xc=XC).density_fit()
    mf.grids.level = 3
    mf.nlcgrids.level = 1
    mf.conv_tol = 1e-8
    mf.max_cycle = 200
    if pcm:
        mf = mf.PCM()
        mf.with_solvent.method = "IEF-PCM"
        mf.with_solvent.eps = EPS
    e = mf.kernel()
    if not mf.converged:
        mf.level_shift = 0.2
        e = mf.kernel(mf.make_rdm1())
    return {"chave": chave, "e_h": float(e), "convergiu": bool(mf.converged), "nao": int(m.nao),
            "pcm": pcm, "t_s": round(time.time() - t, 1)}


def escolher(r):
    """Melhor complexo de cada tipo π distinto (a lista do M8 já vem ordenada por ΔE)."""
    vistos, out = set(), []
    for c in r["complexos"]:
        if c["tipo"].startswith("π") and c["tipo"] not in vistos:
            vistos.add(c["tipo"]); out.append(c)
    return out[:3]


def tarefas():
    r = json.loads((RES / "resumo.json").read_text())
    mono = json.loads((RES / "monomeros.json").read_text())
    simb, nc = r["simbolos"], r["n_caf"]
    at = lambda s, x: [(a, tuple(p)) for a, p in zip(s, x)]
    fant = lambda s, x: [("ghost-" + a, tuple(p)) for a, p in zip(s, x)]
    T = [{"chave": "relax__cafeina", "atomos": at(mono["cafeina"]["simbolos"], mono["cafeina"]["xyz"]), "pcm": True},
         {"chave": "relax__epicatequina", "atomos": at(mono["epicatequina"]["simbolos"], mono["epicatequina"]["xyz"]),
          "pcm": True}]
    for c in escolher(r):
        x = c["xyz"]
        sa, xa, sb, xb = simb[:nc], x[:nc], simb[nc:], x[nc:]
        k = c["id"]
        T += [{"chave": f"{k}__AB", "atomos": at(simb, x), "pcm": True},
              {"chave": f"{k}__A", "atomos": at(sa, xa), "pcm": True},
              {"chave": f"{k}__B", "atomos": at(sb, xb), "pcm": True},
              {"chave": f"{k}__A_gas", "atomos": at(sa, xa), "pcm": False},
              {"chave": f"{k}__A_fant", "atomos": at(sa, xa) + fant(sb, xb), "pcm": False},
              {"chave": f"{k}__B_gas", "atomos": at(sb, xb), "pcm": False},
              {"chave": f"{k}__B_fant", "atomos": at(sb, xb) + fant(sa, xa), "pcm": False}]
    # os maiores primeiro (equilibra os processos)
    T.sort(key=lambda t: -len(t["atomos"]))
    return T, r


def analise(r):
    E = {f.stem: json.loads(f.read_text()) for f in DFT.glob("*.json")}
    tm = json.loads((RES / "termo.json").read_text())
    out = []
    for c in escolher(r):
        k = c["id"]
        need = [f"{k}__{s}" for s in ("AB", "A", "B", "A_gas", "A_fant", "B_gas", "B_fant")] + \
               ["relax__cafeina", "relax__epicatequina"]
        if not all(n in E for n in need):
            out.append({"id": k, "tipo": c["tipo"], "dE_xtb": c["dE"], "pronto": False,
                        "feitos": sum(n in E for n in need), "total": len(need)})
            continue
        e = lambda n: E[n]["e_h"]
        bsse = (e(f"{k}__A_gas") - e(f"{k}__A_fant")) + (e(f"{k}__B_gas") - e(f"{k}__B_fant"))
        dE_int = (e(f"{k}__AB") - e(f"{k}__A") - e(f"{k}__B") + bsse) * KCAL
        dE_lig = (e(f"{k}__AB") - e("relax__cafeina") - e("relax__epicatequina") + bsse) * KCAL
        d = {"id": k, "tipo": c["tipo"], "pronto": True, "dE_xtb": c["dE"], "dG_xtb": c.get("dG"),
             "dE_int": dE_int, "dE_lig": dE_lig, "bsse": bsse * KCAL, "deformacao": dE_lig - dE_int,
             "convergiu": all(E[n]["convergiu"] for n in need)}
        if k in tm:
            g = tm[k]["g_corr"] - tm["cafeina"]["g_corr"] - tm["epicatequina"]["g_corr"]
            d["dG"] = dE_lig + g + CONC_1M
        out.append(d)
    res = {"nivel": f"{XC}/{BASE} (VV10), IEF-PCM água; BSSE contrapeso em fase gasosa; geometrias GFN2-xTB/ALPB",
           "complexos": out}
    (RES / "refino_dft.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))
    for d in out:
        if d["pronto"]:
            print(f"  {d['id']:<16} {d['tipo']:<36} xTB {d['dE_xtb']:6.2f} │ DFT ΔE_lig {d['dE_lig']:6.2f}"
                  f"  ΔE_int {d['dE_int']:6.2f}  BSSE {d['bsse']:5.2f}" + (f"  ΔG {d['dG']:5.2f}" if "dG" in d else ""))
        else:
            print(f"  {d['id']:<16} pendente ({d['feitos']}/{d['total']})")
    return res


def main(so_analise=False, workers=2):
    DFT.mkdir(parents=True, exist_ok=True)
    T, r = tarefas()
    # Hessiana xTB para complexos escolhidos que ainda não têm termoquímica (rápido)
    tm_f = RES / "termo.json"
    tm = json.loads(tm_f.read_text())
    import interacoes as I
    for c in escolher(r):
        if c["id"] not in tm:
            print(f"termo xTB {c['id']} …", flush=True)
            tm[c["id"]] = I.termo(r["simbolos"], np.asarray(c["xyz"]))
            tm_f.write_text(json.dumps(tm, indent=1))
    if not so_analise:
        print(f"{len(T)} cálculos ωB97X-V/def2-SVP · {workers} processos", flush=True)
        qm.rodar_paralelo(sp, T, workers=workers, pasta=DFT, chave="chave")
    analise(r)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--analise", action="store_true")
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    main(a.analise, a.workers)
